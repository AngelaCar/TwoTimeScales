###########################################################
# Create the dataset reccolontts from the raw data colon  #
# in package survival                                     #
###########################################################

# ---- Part 1: load and prepare the data ----
## Load the library and the data
set.seed(2803)
library(data.table)
library(usethis)

## Prepare the data: select only individuals who had a recurrence
##                   and create the two time scales variables
recur <- subset(survival::colon, status == 1 & etype == 1)
reccolon2ts <- subset(survival::colon, id %in% recur$id)

## We can sort by id and etype using data.table

reccolon2ts <- as.data.table(reccolon2ts)
setorder(reccolon2ts, id, etype)

## now we need to shift the variable "time" from the first
## row of each individual to the second
## timer : time since randomization AT recurrence

reccolon2ts[, timer := shift(time, type = "lag"), by = id]

## we can rename the variable "time" to "timedc" and get
## rid of the rows where "etype = 1" (recurrence)

setnames(reccolon2ts, "time", "timedc")
reccolon2ts <- reccolon2ts[etype == 2]

## The first time scale of interest is time since
## randomization, which is our "timedc"
## the second time scale is time since recurrence
## We call this new variable timesr = timedc - timer

## Entry times (on the time since recurrence scale (s))
## everyone starts at 0
reccolon2ts$entrys <- 0

## Exit times (on the time since recurrence scale)
reccolon2ts[, timesr := timedc - timer]

## Introduce (semi-)random left truncation for 40 individuals (a bit less than 10%)
temp <- reccolon2ts[timer <= 500]
tempsampid <- sample(temp$id, 40)
reccolon2ts[id %in% tempsampid,
            entrys := round(runif(n=1, min=1, max=timesr),0), by = id]
reccolon2ts[, entryt := entrys + timer]

## The time variables are measured in days
## we can transform them in years

reccolon2ts[, timedc_y := timedc / 365.25]
reccolon2ts[, timesr_y := timesr / 365.25]
reccolon2ts[, entrys_y := entrys / 365.25]
reccolon2ts[, entryt_y := entryt / 365.25]
reccolon2ts[, timer_y := timer / 365.25]

# covariates as factors
reccolon2ts[, sex := factor(sex, labels = c("female", "male"))]
reccolon2ts[, differ := factor(differ, labels = c("well", "moderate", "poor"))]
reccolon2ts[, extent := factor(extent, labels = c("submucosa", "muscle", "serosa", "contstruct"))]
reccolon2ts[, surg := factor(surg, labels = c("short", "long"))]
reccolon2ts <- as.data.frame(reccolon2ts)
## Save the data
#save(reccolon2ts, file = "../data/reccolon2ts")
usethis::use_data(reccolon2ts, overwrite = T)
