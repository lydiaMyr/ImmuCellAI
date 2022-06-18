## code to prepare `DATASET` dataset goes here

#usethis::use_data("DATASET")

#load("data/ImmuCellAI_marker_7_23.Rdata")
#load("data/ImmuCellAI_compensation_mat_7_16.Rdata")
#load("data/marker_exp_T.Rdata")
#load("data/immune_infiltrate_marker.Rdata")
##load("data/train_data_new.Rdata")
#load("data/train_tag_new.Rdata")
#marker_exp=read.table("data/ImmuCellAI_marker_exp.txt",header=T,row.names = 1,sep="\t")

usethis::use_data(paper_marker,overwrite = TRUE)
usethis::use_data(compensation_matrix,overwrite = TRUE)
usethis::use_data(immune_infiltate_marker,overwrite = TRUE)
usethis::use_data(marker_exp,overwrite = TRUE)
usethis::use_data(marker_exp_T,overwrite = TRUE)
usethis::use_data(train_data,overwrite = TRUE)
usethis::use_data(train_tag,overwrite = TRUE)
