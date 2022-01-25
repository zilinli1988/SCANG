# SCANG (SCAN the Genome)
This is an R package for performing SCANG procedure in whole genome sequencing studies.
## Description
SCANG is an R package for performing a flexible and computationally efficient scan statistic procedure (SCANG) that uses the p-value of a variant set-based test as a scan statistic of each moving window, to detect rare variant association regions for both continuous and dichotomous traits.
## Dependencies
SCANG links to R packages <a href="https://cran.r-project.org/web/packages/Rcpp/index.html">Rcpp</a> and <a href="https://cran.r-project.org/web/packages/RcppArmadillo/index.html">RcppArmadillo</a>, and also imports R packages <a href="https://cran.r-project.org/web/packages/Rcpp/index.html">Rcpp</a>, <a href="https://cran.r-project.org/web/packages/GMMAT/index.html">GMMAT</a>, <a href="https://bioconductor.org/packages/release/bioc/html/GENESIS.html">GENESIS</a>, <a href="https://cran.r-project.org/web/packages/Matrix/index.html">Matrix</a>, <a href="https://cran.r-project.org/web/packages/kinship2/index.html">kinship2</a>. These dependencies should be installed before installing SCANG.
## Installation
```
library(devtools)
devtools::install_github("zilinli1988/SCANG")
```
## Docker Image
A [docker image for SCANG](https://hub.docker.com/repository/docker/zilinli/staarpipeline), including R (version 3.6.1) built with Intel MKL and all STAAR-related packages (STAAR, SCANG, STAARpipeline, STAARpipelineSummary) pre-installed, is located in the Docker Hub. The docker image can be pulled using
```
docker pull zilinli/staarpipeline
```
## Usage
Please see the <a href="doc/SCANG-manual-v1.0.3.pdf">**SCANG** user manual</a> for detailed usage of SCANG package. Please see the <a href="https://htmlpreview.github.io/?https://github.com/zilinli1988/SCANG/blob/master/doc/SCANG_Example_v1.0.3.html">**SCANG** tutorial</a> for an example of analyzing sequencing data using SCANG procedure.

For optimal computational performance, it is recommended to use an R version configured with the Intel Math Kernel Library (or other fast BLAS/LAPACK libraries). See the <a href="https://software.intel.com/en-us/articles/using-intel-mkl-with-r">instructions</a> on building R with Intel MKL.
## Citation
If you use **SCANG** for your work, please cite:

Zilin Li, Xihao Li, Yaowu Liu, Jincheng Shen, Han Chen, Hufeng Zhou, Alanna C. Morrison, Eric Boerwinkle, and Xihong Lin (2019) "Dynamic Scan Procedure for Detecting Rare-Variant Association Regions in Whole-Genome Sequencing Studies". _The American Journal of Human Genetics_, 104(5), 802-814. PMID: <a href="https://www.ncbi.nlm.nih.gov/pubmed/30982610">30982610</a>. PMCID: <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6507043/">PMC6507043</a>. DOI: <a href="https://doi.org/10.1016/j.ajhg.2019.03.002">10.1016/j.ajhg.2019.03.002</a>.
## Version
v1.0.2: Allow incorporating multiple functional annotations in SCANG procedure through <a href="https://github.com/xihaoli/STAAR">STAAR</a>.  
v1.0.2.1: Add one more parameter in SCANG function.<br>
v1.0.3: Clean up the C++ functions to avoid duplication with STAAR package.<br>
v1.0.3.1: Define ARMA_64BIT_WORD 1 in C++ code.
## License
This software is licensed under GPLv3.
