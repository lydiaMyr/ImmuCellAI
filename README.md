# ImmuCellAI
ImmuCellAI (Immune Cell Abundance Identifier, http://bioinfo.life.hust.edu.cn/ImmuCellAI/) is a tool to estimate the abundance of 24 immune cells from gene expression dataset including RNA-Seq and microarray data

## Installation
``` bash
install.packages("devtools")
library(devtools)
install_github("lydiaMyr/ImmuCellAI@main")
#if the "/bin/gtar: not found" error occured, please run the following command "export TAR="/bin/tar" before installation.
```
## Geting started
### new version
``` bash
#The prediction process of ImmuCellAI is updated, which simulated the flow cytometry process to predict cell type abundance by hierarchical strategy. All 24 cell types were divided into two layers, layer1:DC, B cell, Monocyte, Macrophage, NK, Neutrophil, CD4 T, CD8 T, NKT, Tgd; layer2:CD4 naive, CD8 naive, Tc, Tex, Tr1, nTreg, iTreg, Th1, Th2, Th17, Tfh, Tcm, Tem, MAIT.
#Parameters
#sample_expression: Sample expression profile in FPKM, TPM format by RNA-seq or log2-transformed signal by microarray.
#datatype: One of "rnaseq" and "microarray"
#group_tag: One of 0 and 1, if there is the need to perform the comparision between different groups. If the value is 1, users need to add a group tag row in the input epxression matrix to explain the group of each sample.
#response_tag: One of 0 and 1, if there is the need to predict the ICB response of each sample.

ImmuCellAI_new(sample_expression,data_type,group_tag,response_tag)
#output
The output of the function is a list with three variables including sample immune cell abundance, group comparison result and ICB response prediction result.

``` 

### old version
``` bash
#Parameters
#sample_expression: Sample expression profile in FPKM, TPM format by RNA-seq or log2-transformed signal by microarray.
#datatype: One of "rnaseq" and "microarray"
#group_tag: One of 0 and 1, if there is the need to perform the comparision between different groups. If the value is 1, users need to add a group tag row in the input epxression matrix to explain the group of each sample.
#response_tag: One of 0 and 1, if there is the need to predict the ICB response of each sample.
#customer: One of 0 and 1, if there is the need to upload the self-build reference file. if the value = 1, users need to provide the gene signature (by list format in R) and reference expression matrix with rownames is gene, colnames is cell type and separated by tab. 
ImmuCellAI(sample_expression,data_type,group_tag,response_tag,customer)

#output
The output of the function is a list with three variables including sample immune cell abundance, group comparison result and ICB response prediction result.
```
## Citations
``` bash
Miao, Y.-R., Zhang, Q., Lei, Q., Luo, M., Xie, G.-Y., Wang, H., Guo, A.-Y., ImmuCellAI: A Unique Method for Comprehensive T-Cell Subsets Abundance Prediction and its Application in Cancer Immunotherapy. Adv. Sci. 2020, 7, 1902880. https://doi.org/10.1002/advs.201902880
Ya-Ru Miao, Mengxuan Xia, Mei Luo, Tao Luo, Mei Yang, An-Yuan Guo, ImmuCellAI-mouse: a tool for comprehensive prediction of mouse immune cell abundance and immune microenvironment depiction, Bioinformatics, 2021;, btab711. 
```
