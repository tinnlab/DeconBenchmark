# Deconvolution Benchmark

This package provides a common interface to run 46 deconvolution methods.

## Installation

This package requires `docker` to be installed and `R` version 4.0 or later.

### Install docker

Follow instructions at [https://docs.docker.com/engine/install/](https://docs.docker.com/engine/install/) to install `docker`.
To check docker is installed and R can communicate with it, run:

```R
if (!requireNamespace("babelwhale", quietly = TRUE)) {
   install.packages("babelwhale")
}
babelwhale::test_docker_installation(detailed = TRUE)

```

### Install DeconBenchmark
```R
if (!requireNamespace("devtools", quietly = TRUE)) {
   install.packages("devtools")
}
devtools::install_github("tinnlab/DeconBenchmark")

library(DeconBenchmark)
```
## Supported methods

List of supported methods, their required inputs beside bulk data and their original package.

| Method            | Inputs                                                                 | Publication                                    | Original package                                     |
|-------------------|------------------------------------------------------------------------|------------------------------------------------|------------------------------------------------------|
| AdRoit            | single cell expression, single cell labels                             | https://doi.org/10.1038/s42003-021-02739-1     | https://doi.org/10.5281/zenodo.5495546               |
| ARIC              | cell type expression, isMethylation                                    | https://doi.org/10.1093/bib/bbab362            | https://xwanglabthu.github.io/ARIC                   |
| AutoGeneS         | cell type expression                                                   | https://doi.org/10.1016/j.cels.2021.05.006     | https://github.com/theislab/AutoGeneS                |
| BayCount          | number of cell types                                                   | https://doi.org/10.1214/17-AOAS1123            | https://fangzheng-xie.github.io/publication/         |
| BayesCCE          | cell type expression                                                   | https://doi.org/10.1186/s13059-018-1513-2      | https://github.com/cozygene/BayesCCE                 |
| BayICE            | single cell expression, single cell labels                             | https://doi.org/10.1214/20-AOAS1376            | https://github.com/AshTai/BayICE/                    |
| BisqueMarker      | markers                                                                | https://doi.org/10.1038/s41467-020-15816-6     | https://cran.r-project.org/package=BisqueRNA         |
| CellDistinguisher | number of cell types                                                   | https://doi.org/10.1371/journal.pone.0193067   | https://github.com/GeneralElectric/CellDistinguisher |
| CIBERSORT         | signature                                                              | https://doi.org/10.1038/s41587-019-0114-2      | https://cibersortx.stanford.edu/                     |
| CPM               | single cell expression, single cell labels                             | https://doi.org/10.1038/s41592-019-0355-5      | https://cran.r-project.org/package=scBio             |
| DAISM             | single cell expression, single cell labels                             | https://doi.org/10.1016/j.patter.2022.100440   | https://github.com/xmuyulab/DAISM-XMBD               |
| debCAM            | number of cell types                                                   | https://doi.org/10.1093/bioinformatics/btaa205 | http://bioconductor.org/packages/debCAM              |
| Deblender         | number of cell types, markers                                          | https://doi.org/10.1186/s12859-018-2442-5      | https://github.com/kondim1983/Deblender/             |
| DeCompress        | cell type expression, significant genes,                               | https://doi.org/10.1093/nar/gkab031            | https://github.com/bhattacharya-a-bt/DeCompress      |
| deconf            | number of cell types                                                   | https://doi.org/10.1186%2F1471-2105-11-27      | http://web.cbio.uct.ac.za/~renaud/CRAN/              |
| DeconICA          | number of cell types                                                   |                                                | https://github.com/UrszulaCzerwinska/DeconICA        |
| DeconPeaker       | single cell expression, single cell labels                             | https://doi.org/10.3389/fgene.2020.00392       | https://github.com/lihuamei/DeconPeaker              |
| DeconRNASeq       | signature                                                              | https://doi.org/10.1093/bioinformatics/btt090  | https://bioconductor.org/packages/DeconRNASeq/       |
| deconvSeq         | single cell expression, single cell labels, significant genes          | https://doi.org/10.1093/bioinformatics/btz444  | https://github.com/rosedu1/deconvSeq                 |
| DecOT             | single cell expression, single cell labels, single cell subject labels | https://doi.org/10.3389/fgene.2022.825896      | https://github.com/lg-ustb/DecOT                     |
| DeMixT            | single cell expression, single cell labels                             | https://doi.org/10.1016/j.isci.2018.10.028     | http://bioinformatics.mdanderson.org/main/DeMixT     |
| DESeq2            | cell type expression                                                   | https://doi.org/10.1186/s13059-014-0550-8      | https://bioconductor.org/packages/DESeq2/            |
| digitalDLSorter   | single cell expression, single cell labels,                            | https://doi.org/10.3389/fgene.2019.00978       | https://github.com/cartof/digitalDLSorter            |
| DSA               | markers                                                                | https://doi.org/10.1186/1471-2105-14-89        | https://github.com/zhandong/DSA                      |
| dtangle           | cell type expression                                                   | https://doi.org/10.1093/bioinformatics/bty926  | https://cran.r-project.org/package=dtangle           |
| DWLS              | single cell expression, single cell labels                             | https://doi.org/10.1038/s41467-019-10802-z     | https://github.com/dtsoucas/DWLS                     |
| EMeth             | cell type expression                                                   | https://doi.org/10.1038/s41598-021-84864-9     | https://github.com/Hanyuz1996/EMeth                  |
| EPIC              | cell type expression, significant genes                                | https://doi.org/10.7554/elife.26476            | http://epic.gfellerlab.org/                          |
| FARDEEP           | signature                                                              | https://doi.org/10.1371/journal.pcbi.1006976   | https://github.com/YuningHao/FARDEEP                 |
| ImmuCellAI        | signature, markers                                                     | https://doi.org/10.1002/advs.201902880         | http://bioinfo.life.hust.edu.cn/ImmuCellAI           |
| LinDeconSeq       | signature                                                              | https://doi.org/10.1186/s12864-020-06888-1     | https://github.com/lihuamei/LinDeconSeq              |
| Linseed           | number of cell types                                                   | https://doi.org/10.1038/s41467-019-09990-5     | https://github.com/ctlab/linseed                     |
| MCPcounter        | markers                                                                | https://doi.org/10.1186/s13059-016-1070-5      | https://doi.org/10.5281/zenodo.61372                 |
| MethylResolver    | signature                                                              | https://doi.org/10.1038/s42003-020-01146-2     | https://github.com/darneson/MethylResolver           |
| MIXTURE           | signature                                                              | https://doi.org/10.1093/bib/bbaa317            | https://github.com/elmerfer/MIXTURE.App              |
| MOMF              | single cell expression, single cell labels                             | https://doi.org/10.3390/cells8101161           | https://github.com/biostat0903/MOMF/                 |
| MuSic             | single cell expression, single cell labels                             | https://doi.org/10.1038/s41467-018-08023-x     | https://github.com/xuranw/MuSiC                      |
| NITUMID           | signature                                                              | https://doi.org/10.1093/bioinformatics/btz748  | https://github.com/tdw1221/NITUMID                   |
| PREDE             | cell type expression                                                   | https://doi.org/10.1371/journal.pcbi.1008452   | https://xiaoqizheng.github.io/PREDE                  |
| quanTIseq         | signature                                                              | https://doi.org/10.1186/s13073-019-0638-6      | http://icbi.at/quantiseq                             |
| ReFACTor          | number of cell types                                                   | https://doi.org/10.1038/nmeth.3809             | https://github.com/cozygene/refactor                 |
| RNA-Sieve         | single cell expression, single cell labels                             | https://doi.org/10.1101/gr.272344.120          | https://github.com/songlab-cal/rna-sieve             |
| scaden            | single cell expression, single cell labels                             | https://doi.org/10.1126/sciadv.aba2619         | https://github.com/KevinMenden/scaden                |
| SCDC              | single cell expression, single cell labels                             | https://doi.org/10.1093/bib/bbz166             | https://github.com/meichendong/SCDC                  |
| TOAST             | number of cell types, markers                                          | https://doi.org/10.1186/s13059-019-1778-0      | https://bioconductor.org/packages/TOAST              |

Please cite the original publications if you use the methods in your work.

## Examples

```R
# load example data
data(BloodExample)
print(names(BloodExample)) # c("bulk", "singleCellExpr", "singleCellLabels")
bulk <- BloodExample$bulk
singleCellExpr <- BloodExample$singleCellExpr
singleCellLabels <- BloodExample$singleCellLabels
```

### Run an individual method
First, check the required inputs for the method.
```R
methodsToRun <- "AdRoit"

# get required inputs
inputTypes <- getMethodsInputs(methodsToRun)
print(inputTypes) # list (AdRoif = c("bulk", "singleCellExpr", "singleCellLabels"))
```
We already had all inputs needed for the method to run the deconvolution.
```R
deconvolutionResult <- deconvolution(methods = methodsToRun, bulk = bulk, singleCellExpr = singleCellExpr, singleCellLabels = singleCellLabels)
proportion <- deconvolutionResult$AdRoit$P
print(head(proportion))
```


### Run multiple methods
```R
methodsToRun <- c("AdRoit", "ARIC")

# get required inputs
inputTypes <- getMethodsInputs(methodsToRun)
print(inputTypes) # list (AdRoif = c("bulk", "singleCellExpr", "singleCellLabels"), ARIC = c("bulk", "cellTypeExpr", "isMethylation"))
```
The `ARIC` method requires the extra inputs `cellTypeExpr` and `isMethylation`.
Since our data is RNA data, so `isMethylation` is `FALSE`.
Now we generate the cell type expression `cellTypeExpr`

```R
reference <- generateReference(singleCellExpr, singleCellLabels, type="cellTypeExpr")
```
Run the deconvolution using two methods.
```R
deconvolutionResults <- deconvolution(methodsToRun, bulk = bulk, singleCellExpr = singleCellExpr, singleCellLabels = singleCellLabels, cellTypeExpr = reference$cellTypeExpr, isMethylation = F) # Run deconvolution
proportionAdRoit <- deconvolutionResults$AdRoit$P
proportionARIC <- deconvolutionResults$AdRoit$P
print(proportionAdRoit)
print(proportionARIC)
```

