# epiConv tutorial
Li Lin<br>

EpiConv is a novel algorithm to cluster scATAC-seq data and detect differentially accessible (DE) peaks. We will demonstrate how to analyze scATAC-seq data by epiConv in this tutorial.

## Installation
EpiConv is developed under the following environments:
1. CentOS release 6.9 with perl 5.10.1
2. bedtools v2.27.1
3. MACS2
4. R 3.5.1 with:
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
res_epiConv<-add.mat(obj=res_epiConv,x=mat[freq!=0,lib_size>1000],name="peak")
```
We create an object `res_epiConv` containing the raw data and results. Low quality cells (library size <1000) are removed. 
+ `create.epiConv`: create an epiConv object.
  - `meta.features`: the data.frame that contains meta features of single cells (e.g. barcodes and library size).
  - The meta features of cells can be obtained from the object using the form such as `res_epiConv$barcode` or `res_epiConv$lib_size`.
+ `add.mat`: add matrix to the epiConv object.
  - `obj`: the epiConv object.
  - `x`: the matrix.
  - `name`: the list name used to store the matrix.

Next, we normalize the matrix by TF-IDF transformation:
```
mat<-tfidf.norm(mat=res_epiConv@mat[["peak"]],lib_size=res_epiConv$lib_size)
infv<-inf.estimate(mat[,sample(1:ncol(mat),size=500)],
                   sample_size=0.125,nsim=30)
```
+ `tfidf.norm`: perform the TF-IDF transformation.
  - `mat`: the matrix. It must be a Matrix object (created by readMM() or Matrix()).
  - `lib_size`: library size used in normalization. In the scripts above, we used the total number of accessible regions in each cell as library size.
+ `inf.estimate`: learn a small value to replace infinite value in the analysis below.
  - `mat`: the matrix after TF-IDF transformation. We randomly sample a small fraction of cells from the full matrix to save the running time.
  - `sample_size`: the fraction of peaks used in each bootstrap.
  - `nsim`: the number of bootstraps.
  - The settings above is suitable for most data.

The similarities between single cells is calculated based on bootstrap approach. In each replicate we randomly sample some peaks and calculate the similarites between single cells, the final similarities are calculated by averging the results from bootstraps:
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
In the script above, `sample_matrix` contains the peaks for each bootstrap. `Smat` contains the similarities between single cells.
+ `epiConv.matrix`: the function that calculates the similarites.
  - `mat`: the matrix used to calculate the similarities.
  - `inf_replace`: the value used to replace infinity. It can be calculated above or used as an empirical value (e.g. -8). 
+ `add.simlarity`: add similarity matrix to the epiConv object.
  - `obj`: the epiConv object.
  - `x`: the similarity matrix.
  - `name`: the list name used to store the similarity matrix.
  - The similarity matrix can be obtained from the object using the form such as `res_epiConv[["sampl"]]`.

Next we blur the similarities between single cells to denoise the data.
```
Smat<-sim.blur(Smat=res_epiConv[["sampl"]],
               weight_scale=log10(res_epiConv$lib_size),
               neighbor_frac=0.25,
               knn=20)
res_epiConv<-add.similarity(res_epiConv,x=Smat,name="samplBlurred")
```
+ `sim.blur`: the function used to blur the similarity matrix.
  - `Smat`: the similarity matrix for denoising.
  - `weight_scale`: the weight for each cell. Generally we think cells with high library size are more reliable and use log10 library size as weights.
  - `neighbor_frac`: the fraction of information used from the neighbors of each cell. It should be within 0~1. Higher value means strong denoising while lower value means weak denoising. The default value 0.25 is suitable for most datasets. If you find the downstream result is poor, you can try higher values (e.g. 0.5).
  - `knn`: the number of neighbors for each cell. There is no need to change it unless your data set is very small (e.g. <200 cells).

Finally we use umap to learn the low-dimensional embedding of the data:
```
umap_settings<-umap::umap.defaults
umap_settings$input<-"dist"
umap_settings$n_components<-2
umap_res<-umap::umap(max(Smat)-Smat,config=umap_settings)$layout
res_epiConv<-add.embedding(obj=res_epiConv,x=umap_res,name="samplBlurred")
plot(res_epiConv@embedding[["samplBlurred"]],pch="+")
```
The distance is calculated by `max(Smat)-Smat`. 
+ `add.embedding`: add the embeddiing to the epiConv object.
  - `obj`: the epiConv object.
  - `x`: the embedding matrix.
  - `name`: the list name used to store the embedding.

To prepare the input for epiConv-full, we save the list of high-quality barcodes and the epiConv object:
```
temp<-data.frame(res_epiConv@meta.features$barcode,1)
write.table(temp,file="pbmc5k/pbcm5k_ident.tsv",row.names=F,col.names=F,quote=F,sep="\t")
saveRDS(res_epiConv,file="res_epiConv_simp.rds")
```
#### Tips for large datasets
As R does not support long vectors, error will occur when the dataset is large (e.g. >80,000 cells). The function `epiConv.matrix` and `sim.blur` have the `bin` parameter with default value of 10,000. When the dataset contains more than 10,000 cells, the function will split the matrix into several parts, each containing 10,000 cells. Generally there is no need to change `bin`, but if some memory errors occured, you can tried small values such as 5,000 (e.g. `epiConv.matrix(mat=mat,inf_replace=infv,bin=5000)`). Based on our tests, epiConv requires 520GB memory for dataset with 81,173 cells and 436,206 peaks.

### epiConv-full
In order to accelerate the running speed, epiConv-full runs in bash shell but some common steps shared by epiConv-simp is performed in R. Its input is a compressed bed file named `<prefix>_frag.bed.gz` (can be read by `zcat`; e.g `pbmc5k/pbmc5k_frag.bed.gz`):<br>
| | |
|-|-|
|1st column|chromsome|
|2nd column|starting site of the fragment|
|3rd column|ending site of the fragment|
|4th column|cell barcode|

An example file:
```
chr1    10073   10198   CCGCATTGTGTTTCTT-1
chr1    10085   10198   GAGACTTCAGCAACAG-1
chr1    10085   10302   GGGCCATAGAGATTAC-1
chr1    10126   10309   CCACGTTCAGTAGTCT-1
...
```
Other columns will be ignored. Generally this file will be provided by low level processing tools (e.g. cellranger from 10X Genomics). <br>
Also you need to prepare a Tab-delimited file containing valid barocdes named `<prefiex>_ident.tsv` (e.g. `pbmc5k/pbmc5k_ident.tsv` generated above). The first column of the file should be valid cell barcodes and the second column can be its identies (e.g. cell type, batch, experimental condition, or simply use 1 for all cells if there are not any information on the identities of cells). 
| | |
|-|-|
|1st column|cell barcode|
|2nd column|custom information (e.g batch, experimental condition, or 1 for all cells if there are not any information. But you can NOT use "NA" as epiConv will ignore all barcodes with "NA" identy.)|

An example file:
```
AAACGAAAGACACTTC-1      1
AAACGAAAGCATACCT-1      1
AAACGAAAGCGCGTTC-1      1
AAACGAAAGGAAGACA-1      1
AAACGAACAGGCATCC-1      1
...
```
Only valid barcodes will be processed.<br><br>
First we use `peak_calling.sh` to call high density regions of Tn5 insertions:
```
peak_calling.sh <prefix> <extsize> <fraction of data retained>
```
- `<prefix>`: the prefix of data.<br>
- `<extsize>`: the same parameter in `MACS2 pileup`. For example, set `<extsize>` to 100 will make MACS2 extend each insertion from both directions by 100bp before pileup insertions.<br>
- `<fraction of data retained>`:  the fraction of data you want to use in downstream analysis. Based on our preliminary analysis, this parameter need not to be accurately specified but should be close to the fraction of fragments from nucleosome-free regions. For example, you can check the histogram of insertion length from 10X Genomics reporting summary to learn this parameter (the fraction of fragments in first peak in the histogram).

For the PBMC dataset, we use the following command:
```
~/epiConv/peak_calling.sh pbmc5k/pbmc5k 100 0.7
```
When it is finished (~4 hours), there will be a new file `pbmc5k/pbmc5k_peak.bed` in bed format containing the peaks. In order to run epiConv in parallel, we need to split the job using the linux command `split`:
```
split -a2 -d -l127000 pbmc5k/pbmc5k_peak.bed pbmc5k/pbmc5k_peak.run
```
Here we split the peak file into 10 jobs, each containing 127000 peaks. Like epiConv-simp, we need to randomly sample some peaks to perform bootstraps. We generate the sampling file using `peak_sampl.sh`:
```
~/epiConv/peak_sampl.sh <peak file> <number of bootstraps> <fraction of peaks in each bootstrap> \
                        <random seed> > <prefix>_sampl.mtx
```
- `<peak file>`: the peak file generated in previous step.<br>
- `<number of bootstraps>`: number of bootstraps.<br>
- `<fraction of peaks in each bootstrap>`: fraction of peaks in each bootstrap.<br>
- `<random seed>`: random seed.<br>
- `peak_sampl.sh` will directly print to the standard output. 

#### Tips for large datasets
In order to reduce the required memory for each thread, we can adjust settings of bootstrap. For example, we can set `<number of bootstraps>` to 10 and increase `<fraction of peaks in each bootstrap>` to 0.25 to make sure that each peak is still be sampled > 2x times. Actually the bootstrap step is aimed to reduce the noise for low-sequencing data (e.g. < 2,000 fragment for most cells). Based on our analysis, the results remain similar even without any bootstrap for most datasets (number of bootstrap=1; fraction of peaks=1). So reducing the number of bootstraps won't affect the results but can reduce the required memory.

For the PBMC dataset, we use the following command:
```
~/epiConv/peak_sampl.sh pbmc5k/pbmc5k_peak.bed 30 0.125 12345 >pbmc5k/pbmc5k_sampl.mtx
```
Then we use `convolution.sh` to calculate the similarites between single cells:
```
convolution.sh convolution.sh <prefix> <suffix> <standard deviation of normal distribution> 
```
- `<prefix>`: the prefix of data.
- `<suffix>`: the suffix of the data, depends on the `split` command called above (e.g. run00, run01,run02...)
- `<standard deviation of normal distribution>`: decide the interactions between insertions. We think insertions from two cells with distance < 4σ suggest an active regulatory element shared by these two cells.

For the PBMC dataset,we run `convolution.sh` as follows:
```
~/epiConv/convolution.sh data/pbmc5k run00 100
~/epiConv/convolution.sh data/pbmc5k run01 100
~/epiConv/convolution.sh data/pbmc5k run02 100
...
```
This step can be run in parallel to save the running time. `convolution.sh` will automatically read the inputs. So all input files should be properly named as described above. Each thread requires approximately (n cells)^2/2*(n bootstraps)*4/2^30 GB RAM (e.g. 1.4 GB RAM for 5,000 cells).<br>
After running (each job requires ~4 hours), the script will generate two files: `<prefix>_cmat.<suffix>` and `<prefix>_sampled.<suffix>`. `<prefix>_cmat.<suffix>` is a binary files contains the pairwise similarities bewtween cells for each peak. `<prefix>_sampled.<suffix>` is a Tab-delimited file contains (ncells)*(ncells-1)/2 rows and nbootstraps columns, with each element containing the similarity between two cells in one bootstrap. An example for `pbmc5k_sampled.run00`:
```
10.5113	9.8516	18.8105	6.5908	7.2311 ......
3.9813	5.4448	1.8991	3.7418	4.0857 ......
6.8450	12.2343	1.5356	6.6680	2.7696 ......
0.8218	4.7227	11.9622	1.1415	8.0251 ......
27.5526	15.5030	29.7605	15.5724	20.6794 ......
......
```
In order to acquire the summed similarites between cells, we need to sum the corressponding elements for each job using `paste` and `gawk`:
```
paste -d " " pbmc5k/pbmc5k_sampled.run?? |\
	gawk -f ~/epiCOnv/run_merge.gawk ncol=30 >pbmc5k/pbmc5k_sampled.mat
```
In the script above, `ncol` should be equal to number of columns in `pbmc5k/pbmc5k_sampled.run??`. If you split the running into many small jobs (e.g. 100), this step can also run in parallel as follows (assuming we split it into 100 jobs with suffix from run00 to run99):
```
paste -d " " pbmc5k/pbmc5k_sampled.run0? |\
	gawk -f ~/epiCOnv/run_merge.gawk ncol=30 >pbmc5k/pbmc5k_sampled0.mat
paste -d " " pbmc5k/pbmc5k_sampled.run1? |\
	gawk -f ~/epiCOnv/run_merge.gawk ncol=30 >pbmc5k/pbmc5k_sampled1.mat
.......
paste -d " " pbmc5k/pbmc5k_sampled?.mat |\
	gawk -f ~/epiCOnv/run_merge.gawk ncol=30 >pbmc5k/pbmc5k_sampled.mat
```
After all these is done, we summarize the results from bootstraps and transform the data into a ncellXncell square matrix:
```
gawk -f ~/epiConv/rep_merge.gawk  pbmc5k/pbmc5k_sampled.mat \
	>pbmc5k/pbmc5k_smat.txt
```
An example for `pbmc5k/pbmc5k_smat.txt`:
```

1.85311
1.64461,2.03432
1.79543,2.22283,1.82333
1.78524,2.26929,1.97568,2.18579
.......
```
Note that the first line is blank. Only lower triangle elements are stored. The matrix is un-normalized by library size. In order to normalize the matrix, we also need to obtain the number of insertions falling into the peaks for each single cell:
```
paste pbmc5k/pbmc5k_lib.run* |\
	gawk '{sum=0;for(j=1;j<=NF;j++) sum+=$j;print sum}' |\
	paste pbmc5k/pbmc5k_ident.tsv - > pbmc5k/pbmc5k.info
```
In the script above, we calculated the library sizes of single cells and combined them with the identity file. `pbmc5k/pbmc5k.info` shall be like this:
```
AAACGAAAGACACTTC-1	1	13154
AAACGAAAGCATACCT-1	1	22474
AAACGAAAGCGCGTTC-1	1	11189
AAACGAAAGGAAGACA-1	1	35276
AAACGAACAGGCATCC-1	1	17233
```
The 3rd column is the library size for single cells. After that, all steps in bash shell are finished. The following steps are performed in R.
First we still need to read the scource file and data:
```
source("~/epiConv/epiConv_functions.R")
cell_info<-read.table(file="pbmc5k/pbmc5k.info")
mfeature<-data.frame(barcode=as.character(cell_info[,1]),
                     lib_size=cell_info[,3],
                     ident=factor(cell_info[,2]))
res_epiConv<-create.epiconv(meta.features=mfeature)
Smat<-read.csv(file="pbmc5k/pbmc5k_smat.txt",col.names=res_epiConv$barcode,
               header=F,fill=T,blank.lines.skip=F)
Smat<-as.matrix(Smat)
Smat[is.na(Smat)]<-0
Smat<-Smat+t(Smat)
Smat<-t(t(Smat)-log10(res_epiConv$lib_size))-log10(res_epiConv$lib_size)
res_epiConv<-add.similarity(res_epiConv,x=Smat,name="sampl")
```
In the script above, we read the un-normalized similarity matrix and normalize it by library size of single cells. Then we blur the similarity matrix and perform low-dimensional embedding:
```
Smat<-sim.blur(Smat=res_epiConv[["sampl"]],
               weight_scale=log10(res_epiConv$lib_size),
               neighbor_frac=0.25,
               knn=20)
res_epiConv<-add.similarity(res_epiConv,x=Smat,name="samplBlurred")
umap_settings<-umap::umap.defaults
umap_settings$input<-"dist"
umap_settings$n_components<-2
umap_res<-umap::umap(max(Smat)-Smat,config=umap_settings)$layout
res_epiConv<-add.embedding(res_epiConv,x=umap_res,name="samplBlurred")
plot(res_epiConv@embedding[["samplBlurred"]],pch="+")

saveRDS(res_epiConv,file="res_epiConv_full.rds")
```




