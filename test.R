# devtools::create("~/UNSW/VafaeeLab/HostPathogenInteraction/dRNASb")

# BiocManager::install("Mfuzz")
# install.packages("ClueR")

setwd("~/UNSW/VafaeeLab/HostPathogenInteraction/dRNASb")

devtools::document()
devtools::install()



# result <- dRNAR::seq_analyze("Alldat.csv", "pheno.csv", "GO_Function.csv")
#
# hou2_result <- result[[1]]
dRNASb::get_package_info("Abhishek")



dat <- read.csv("Data/Alldat.csv",row.names = 1)
pheno <- read.csv("Data/pheno.csv", row.names = 1)
fun<-read.csv("Data/GO_Function.csv")
ppi<-read.csv("Data/PPI.csv")
ann<-read.csv("Data/Annotation.csv")


data_file_path = "Data/Alldat.csv"
phenotype_file_path = "Data/pheno.csv"
go_function_file_path = "Data/GO_Function.csv"
ppi_file_path = "Data/PPI.csv"
annotation_file_path = "Data/Annotation.csv"
