filename <- "testing_file.rds"
a<- data.frame(Vals = rep(1,10))

saveRDS(a, file =  paste0(parent,"/Outputs/epidemia_fits/", filename))