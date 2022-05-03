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

