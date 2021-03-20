library(epidemia)
library(dplyr)
library(here)
library(readr)
library(optparse)
#options(mc.cores = parallel::detectCores())

job.index <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
job.id <- as.numeric(Sys.getenv("PBS_JOBID"))

short_state <- c("AL", "AK","AZ", "AR", "CA","CO","CT","DE","FL","GA","ID","IL","IN","IA","KS","KY",
                 "LA","ME","MD","MA","MI","MN","MS","MO","MT","NE","NV","NH","NJ","NM","NY","NC","ND",
                 "OH","OK","OR","PA","RI","SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY") 

long_state <- c("alabama", "alaska", "arizona", "arkansas", "california", "colorado", "connecticut",
                "delaware", "florida", "georgia", "idaho", "illinois","indiana","iowa", "kansas",
                "kentucky", "louisiana", "maine", "maryland", "massachussets", "michigan",
                "minnesota", "mississipi", "missouri","montana", "nebraska", "nevada", "new_hampshire",
                "new_jersey", "new_mexico", "new_york", "north_carolina", "north_dakota", "ohio",
                "oklahoma", "oregon", "pennsylvania", "rhode_island", "south_carolina",
                "south_dakota", "tennessee", "texas", "utah", "vermont", "virginia", "wahsington",
                "west_virginia","wisconsin", "wyoming")

#job.id <- 38
this.short <- short_state[job.index]
this.long <- long_state[job.index]

group.file.name<- paste0(this.long, "_groups.rds")

groups <- readRDS(here::here("data", "grouped_counties", group.file.name))
data <- readRDS(here::here("data", "covid", "prepped_data_Jan21.rds"))

state_data <- filter(data, state == this.short)

N <- nrow(groups)


for (i in (seq (1:N))){
  
  #Obtaining the data
  fips <- unlist(groups$county_fips[i])
  
  if (length(fips)==1){
    data_sub <- filter(state_data, countyFIPS == fips)
  }else{
    data_sub <- filter(state_data, countyFIPS %in% fips)
    tmp <- aggregate(cbind(cases, deaths) ~ date, data=data_sub, FUN=sum, na.action = NULL)
    data_sub <- data_sub[1:376,]
    data_sub <- subset(data_sub, select = -c(countyFIPS, date,cases, deaths, state))
    data_sub <- cbind(data_sub, tmp)
    n <- nrow (data_sub)
    data_sub$county <- rep(groups$county[i], n)
  }
  #Doing the models
  cases_no_NA = data_sub$cases
  cases_no_NA[is.na(cases_no_NA)]=0
  idx <- which(cumsum(cases_no_NA) >= 10)[1]
  
  if (length(idx) == 0) {
    stop(paste0("Fewer than 10 cumulative cases in entire epidemic. Not modeling."))
  }
  
  start_date <- data_sub$date[idx] - 7
  data_sub <- filter(data_sub, date >= start_date)
  
  data_sub <- mutate(data_sub, week = as.integer(format(date, "%V")))
  
  pops <- data.frame("county" = groups$county[i], "population"= groups$population[i])
  
  args <- list(data= data_sub, si = EuropeCovid$si)
  
  deaths <- epiobs(
    formula = deaths(county, date) ~ 1,
    family = "neg_binom", # overdispersion for daily counts
    i2o = EuropeCovid$obs$deaths$i2o * 0.01,
    prior_intercept = rstanarm::normal(1, 0.5),
    link = "identity"
  )
  
  cases <- epiobs(
    formula = cases(county, date) ~ 1,
    i2o = c(rep(0, 3), rep(1/7, 7)) * 0.2,
    family = "neg_binom",
    prior_intercept = rstanarm::normal(1, 0.5),
    link = "identity"
  )
  
  args$obs <- list(deaths = deaths)
  
  args$rt <- epirt(
    formula = R(county, date) ~ rw(time=week)
  )
  
  args$pops <- pops
  
  args$algorithm <- "sampling"
  args$pop_adjust <- FALSE
  args$init_run <- TRUE
  args$sampling_args <- list(iter = 1e3, seed=12345, chains = 1)
  
  filename <- paste0(gsub("\\s", "_", groups$county[i]), "_", "PA_",job.id,".rds")
  
  fit <- do.call("epim", args)
  
  
  res <- list(
    fit = fit,
    county = groups$county[i],
    state = this.short
  )
  
  wd <- getwd()
  setwd("..")
  parent <- getwd()
  setwd(wd)
  
  saveRDS(res, file =  paste0(parent,"/Outputs/epidemia_fits", filename))

  
  rt <- posterior_rt(res$fit)
  draws <- rt$draws
  rt_medians <- apply(draws,2, median)
  rt_df <- data.frame("Date"= rt$time, "Rt_medians"= rt_medians)
  
  
  filename2 <- paste0(gsub("\\s", "_", groups$county[i]), "_", "PA_",job.id,"_medians.rds")
  write.csv(rt_medians, file =  paste0(parent,"/Outputs/rt_medians", filename2))
}




