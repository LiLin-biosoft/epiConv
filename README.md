# epiConv tutorial
Li Lin<br>

EpiConv is a novel algorithm to cluster scATAC-seq data and detect differentially accessible (DE) peaks. We will demonstrate how to analyze scATAC-seq data by epiConv in this tutorial.

## Installation
EpiConv is developed under the following environments:
1. CentOS release 6.9 with perl 5.10.1
2. bedtools v2.27.1
3. MACS2
4. R 3.5.1
    1. Matrix package
    2. bigmemory package
    3. biganalytics package
    4. umap or uwot packages

All source codes of epiConv can be found in the folder `source/`. Here assume we put them into the folder `~/epiConv/`. Two source files need to be compiled.<br>
```
g++ ~/epiConv/matrix_sum.c -o ~/epiConv/matrix_sum
g++ ~/epiConv/matrix_sampl.c -o ~/epiConv/matrix_sampl
```
## Usage
We will use one dataset of [PBMCs](https://support.10xgenomics.com/single-cell-atac/datasets/1.2.0/atac_pbmc_5k_v1) from 10X Genomics as example. Two files are required: `Peak by cell matrix (filtered)` and `Fragments (TSV)`. Here we put them into folder `pbmc5k/`, extract the first file and rename the second file as `pbmck5k_frag.bed.gz`.

### epiConv-simp
There are two versions of epiConv: epiConv-full and epiConv-simp. EpiConv-full calculates the similarities between cells from raw Tn5 insertion profiles and epiConv-simp calculates the similarities from binary matrix. We first show the analysis pipeline for epiConv-simp. It is an implemention in R. First we read the source file of epiConv and the data:
```
source("~/epiConv/epiConv_functions.R")
mat<-readMM(file="pbmc5k/matrix.mtx")
barcode<-read.table(file="pbmc5k/barcodes.tsv",colClass="character")[,1]
peak<-read.table(file="pbmc5k/peaks.bed")
colnames(peak)<-c("seqnames","start","end")
rownames(mat)<-paste(peak$seqnames,":",
                     peak$start,"-",
                     peak$end,sep="")
colnames(mat)<-barcode
lib.estimate<-function(mat){
  mat@x[mat@x>1]<-1
  return(Matrix::colSums(mat))}
lib_size<-lib.estimate(mat)
freq.estimate<-function(mat){
  mat@x[mat@x>1]<-1
  return(Matrix::rowSums(mat)/ncol(mat))
}
freq<-freq.estimate(mat[,lib_size>1000])

res_epiConv<-create.epiconv(meta.features=data.frame(barcode=barcode[lib_size>1000],
                                                     lib_size=lib_size[lib_size>1000]))
res_epiConv@mat[["peak"]]<-mat[freq!=0,lib_size>1000]
```
We create an object `res_epiConv` containing the raw data and results. Low quality cells (library size <1000) are removed. 
Function `create.epiConv`: create an epiConv object.<br>
+ `meta.features`: the data.frame that contains meta features of single cells (e.g. barcodes and library size).<br>
The meta features of cells can be obtained from the object using the form such as `res_epiConv$barcode` or `res_epiConv$lib_size`.<br>
Next, we normalize the matrix by TF-IDF transformation:
```
mat<-tfidf.norm(mat=res_epiConv@mat[["peak"]],lib_size=res_epiConv$lib_size)
infv<-inf.estimate(mat[,sample(1:ncol(mat),size=500)],
                   sample_size=0.125,nsim=30)
```
Function `tfidf.norm`: perform the TF-IDF transformation.<br>
+ `mat`: the sparse matrix in MatrixMarket format.
+ `lib_size`: library size used in normalization. In the scripts above, we used the total number of accessible regions in each cell as library size.<br>
Function `inf.estimate`: learn a small value to replace infinite value in the analysis below.<br>
+ `mat`: the matrix after TF-IDF transformation. We randomly sample a small fraction of cells from the full matrix to save the running time.
+ `sample_size`: the fraction of peaks used in each bootstrap.<br>
+ `nsim`: the number of bootstraps.<br>
Generally the settings above is suitable for most data.

The similarities between single cells is based on a bootstrap approach. In each bootstrap we randomly sample some peaks and calculate the similarites between single cells, the final similarities are calculated by averging the results from bootstraps:
```
sample_size<-floor(nrow(mat)/8)
nsim<-30
sample_matrix<-lapply(1:nsim,function(x) sample(1:nrow(mat),size=sample_size))

Smat<-matrix(0,res_epiConv@ncell,res_epiConv@ncell)
for(i in sample_matrix){
  Smat<-Smat+epiConv.matrix(mat=mat[i,],inf_replace=infv)
}
Smat<-Smat/nsim
res_epiConv<-add.similarity(obj=res_epiConv,x=Smat,name="sampl")
```
In the script above, `sample_matrix` contains the peaks for each bootstrap. Function `epiConv.matrix` is used to calculate the similarites. `mat` specifies the matrix used to calculate the similarities. `inf_replaces` should be specified by the value calculated above or an empirical value (e.g. -8). Function `add.simlarity` adds results to `res_epiConv@similarity[["sampl"]]`. `obj` specifies the epiConv object, `x` specifies the similarity matrix and `name` specifies the list name used to store the data. We can simply use the form `res_epiConv[["sampl"]]` to obtain the similarity matrix from the object.

Next we blur the similarities between single cells to denoise the data.
```
Smat<-sim.blur(Smat=res_epiConv[["sampl"]],
               weight_scale=log10(res_epiConv$lib_size),
               neighbor_frac=0.25,
               knn=20)
res_epiConv<-add.similarity(res_epiConv,x=Smat,name="samplBlurred")
```
`sim.blur` is used to blur the similarity matrix. `Smat` specifies the similarity matrix, `weight_scale` specifies the weight for each cell. Generally we think cells with high library size are more reliable and use log10 library size as weights. `neighbor_frac` specifies the fraction of information used from the neighbors of each cell. It should be within 0~1. Higher value means strong denoising while lower value means weak denoising. The default value 0.25 is suitable for most datasets. If you find the downstream result is poor, you can try higher values (e.g. 0.5). `knn` specifies the number of neighbors for each cell. There is no need to change it unless your data set is very small (e.g. <200 cells).

Finally we use umap to learn the low-dimensional embedding of the data:
```
umap_settings<-umap::umap.defaults
umap_settings$input<-"dist"
umap_settings$n_components<-2
umap_res<-umap::umap(max(Smat)-Smat,config=umap_settings)$layout
res_epiConv<-add.embedding(obj=res_epiConv,x=umap_res,name="samplBlurred")
plot(res_epiConv@embedding[["samplBlurred"]],pch="+")
```
The distance is calculated by `max(Smat)-Smat`. `add.embedding` is used to add the embeddiing to the epiConv object. `obj` specifies the epiConv object, `x` specifies the embedding matrix and `name` specifies the list name used to store the embedding.
To prepare the input for epiConv-full, we save the list of high-quality barcodes to file:
```
temp<-data.frame(res_epiConv@meta.features$barcode,1)
write.table(temp,file="pbmc5k/pbcm5k_ident.tsv",row.names=F,col.names=F,quote=F,sep="\t")
saveRDS(res_epiConv,file="res_epiConv.rds")
```

### epiConv-full
In order to accelerate the running speed, epiConv-full runs in bash shell but some common steps shared by epiConv-simp is performed in R. Its input is a compressed bed file named `<prefix>_frag.bed.gz` (can be read by `zcat`; e.g `pbmc5k/pbmc5k_frag.bed.gz`) with the following format: first column, chromsome; second column, starting site of the fragment, third column, ending stie of the fragment; fourth column, cell barcodes. Other columns are ignored. Generally this file will be provided by low level processing tools (e.g. cellranger from 10X Genomics). <br><br>
  Also you need to prepare a file containing valid barocdes named `<prefiex>_ident.tsv` (e.g. `pbmc5k/pbmc5k_ident.tsv` generated above). The first column of the file should be valid cell barcodes and the second column can be its identies (e.g. cell type, batch, experimental condition, or simply use 1 for all cells if there are not any information on the identities of cells). Only valid barcodes will be processed.<br><br>
  First we use `peak_calling.sh` to call high density regions of Tn5 insertions:
```
peak_calling.sh <prefix> <extsize> <fraction of data retained>
```
`<prefix>`: the prefix of data.<br>
`<extsize>`: the same parameter in `MACS2 pileup`. For example, set `<extsize>` to 100 will make MACS2 extend each insertion from both directions by 100bp before pileup insertions.<br>
`<fraction of data retained>`: specify the fraction of data you want to use in downstream analysis. Based on our preliminary analysis, this parameter need not to be accurately specified but should be close to the fraction of fragments from nucleosome-free regions. For example, you can check the histogram of insertion length from 10X Genomics reporting summary to learn this parameter (the fraction of fragments in first peak in the histogram).<br>
For the PBMC dataset, we use the following command:
```
~/epiConv/peak_calling.sh pbmc5k/pbmc5k 100 0.7
```
When it is finished (~4 hours), there will be a new file `pbmc5k/pbmc5k_peak.bed` in bed format containing the peaks. In order to run epiConv in parallel, we need to split the job using the linux command `split`:
```
split -a2 -d -l127000 pbmc5k/pbmc5k_peak.bed pbmc5k/pbmc5k_peak.run
```
Here we split the peak file into 10 jobs, each containing 127000 peaks. Like epiConv-simp, we need to randomly sample some peaks to perform bootstraps:
```
~/epiConv/peak_sampl.sh <peak file> <number of bootstraps> <fraction of peaks in each bootstrap> <random seed>
```
`<peak file>`: the peak file generated in previous step.<br>
`<number of bootstraps>`: number of bootstraps.<br>
`<fraction of peaks in each bootstrap>`: fraction of peaks in each bootstrap.<br>
`<random seed>`: random seed.<br>
`peak_sampl.sh` will directly print to the standard output. For the PBMC dataset, we use the following command:
```
~/epiConv/peak_sampl.sh pbmc5k/pbmc5k_peak.bed 30 0.125 12345 >pbmc5k/pbmc5k_sampl.mtx
```


Then we use `convolution.sh` to calculate the similarites between single cells:
```

```
