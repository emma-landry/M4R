library(epidemia)
library(dplyr)
library(here)
library(readr)

options(mc.cores = parallel::detectCores())
options(mc.cores = 4)

job.index <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
print(job.index)
job.id <- Sys.getenv("PBS_JOBID")


short_state <- c("AL", "AK","AZ", "AR", "CA","CO","CT","DE","FL","GA","ID","IL","IN","IA","KS","KY",
                 "LA","ME","MD","MA","MI","MN","MS","MO","MT","NE","NV","NH","NJ","NM","NY","NC","ND",
                 "OH","OK","OR","PA","RI","SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY") 

long_state <- c("alabama", "alaska", "arizona", "arkansas", "california", "colorado", "connecticut",
                "delaware", "florida", "georgia", "idaho", "illinois","indiana","iowa", "kansas",
                "kentucky", "louisiana", "maine", "maryland", "massachusetts", "michigan",
                "minnesota", "mississippi", "missouri","montana", "nebraska", "nevada", "new_hampshire",
                "new_jersey", "new_mexico", "new_york", "north_carolina", "north_dakota", "ohio",
                "oklahoma", "oregon", "pennsylvania", "rhode_island", "south_carolina",
                "south_dakota", "tennessee", "texas", "utah", "vermont", "virginia", "washington",
                "west_virginia","wisconsin", "wyoming")


this.short <- short_state[job.index]
this.long <- long_state[job.index]

data <- readRDS(here::here("data", "covid", "prepped_data_state.rds"))
populations <- readRDS(here::here("data", "covid", "state_populations.rds"))

state_data <- filter(data, State == this.short)

#Doing the models
cases_no_NA = state_data$cases
cases_no_NA[is.na(cases_no_NA)]=0
idx <- which(cumsum(cases_no_NA) >= 10)[1]

if (length(idx) == 0) {
  stop(paste0("Fewer than 10 cumulative cases in entire epidemic. Not modeling."))
}

start_date <- state_data$date[idx] - 7
state_data <- filter(state_data, date >= start_date)

state_data <- mutate(state_data, week = as.integer(format(date, "%V")))
new_year <- min(which(state_data$week ==1))
state_data[new_year:nrow(state_data),]$week <- state_data[new_year:nrow(state_data),]$week +53

pops <- subset(populations, populations$State== this.long)$Population
state_data$pop<- pops

args <- list(data=state_data)

inf <- epiinf( gen = EuropeCovid$si,
               pop_adjust = FALSE,
               susceptibles = pop) #check prior_tau

deaths <- epiobs(
  #formula = deaths(county, date) ~ 1,
  formula = deaths ~ 1,
  family = "neg_binom", # overdispersion for daily counts
  i2o = EuropeCovid$inf2death,
  prior_intercept = rstanarm::normal(1, 0.5),
  link = "identity"
)


args$obs <- deaths

args$rt <- epirt(
  formula = R(State, date) ~ rw(time=week)
)

args$algorithm <- "sampling"
args$init_run <- TRUE
args$inf <- inf
args$iter <- 2.5e3
args$chains <- 4
args$seed <- 12345

filename <- paste0(this.long, "_", job.id,".rds")

fit <- do.call("epim", args)

res <- list(
  fit = fit,
  county = groups$county[i],
  state = this.short
)

wd <- getwd()
setwd("..")
setwd("..")
parent <- getwd()
setwd(wd)

saveRDS(res, file =  paste0(parent,"/Outputs/epidemia_fits/state1/", filename))


rt <- posterior_rt(res$fit)
draws <- rt$draws
rt_medians <- apply(draws,2, median)
rt_df <- data.frame("Date"= rt$time, "Rt_medians"= rt_medians)


filename2 <- paste0(this.long,"_",job.id,"_medians.rds")
write.csv(rt_medians, file =  paste0(parent,"/Outputs/rt_medians/state1/", filename2))