filename <- "testing_file.rds"
a<- data.frame(Vals = rep(1,10))

saveRDS(res, file =  paste0(parent,"/Outputs/epidemia_fits/", filename))