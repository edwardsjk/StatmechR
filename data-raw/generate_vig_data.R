## CREATE VIGNETTE DATA ##

####################
#### Vignette 1 ####
####################

n <- 30
id <- c(1:n)
x <- ifelse(id %% 2 == 0, 0, 1)

k <- 1000
spread <- 10
peaktime <- 30 + id * 5
tau <- 225
daily_count <- list()
dcnts <- list()
for(i in 1:n){
  dcnts[[i]] <- round(StatmechR::normmdl(c(k, peaktime[i], log(spread)), tau))
  daily_count[[i]] <- cbind(rep(id[i], tau), c(1:tau), x[i], dcnts[[i]])
}

counts <- do.call(rbind, daily_count)
counts <- as.data.frame(counts)
names(counts) <- c("location", "time", "variable_1", "cases")

vig1_data <- counts
save(vig1_data, file = "data/vig1_data.RData")
