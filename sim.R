################################################################################
# Updated version of the R code for an example of the analysis in:
#
#   "Investigating uncertainty in the minimum temperature mortality:
#     methods and application to 52 Spanish cities"
#   Armstrong B, Tobias A, Gasparrini A
#   http://www.ag-myresearch.com/2017_tobias_epidem.html
#
# Update: 15 January 2017
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/2017_tobias_Epidem_Rcodedata
################################################################################

################################################################################
# SIMULATION STUDY ON THE FUNCTION findmin()
################################################################################

library(dlnm) ; library(splines)
source("findmin.R")

################################################################################
# DEFINE A KNOWN EXPOSURE-RESPONSE ASSOCIATION

# LOAD A SUBSET OF THE DATA: LONDON JULY-DEC 2005
# (CONVENIENCE SAMPLE, LARGISH ASIMMETRIC CI)
data <- subset(read.csv("london.csv"), year==2005 & month>6)
head(data)

# A UNIDIMENSIONAL 4DF SPLINE 
b <- onebasis(data$tmean,df=4)

# SIMPLE MODEL, WITH PREDICTED OUTCOME
m <- glm(death~b,family=poisson,data)
deathpred <- predict(m,type="response")

# REAL MINIMUM
(min <- findmin(b,m))

################################################################################
# RUN SIMULATION

# NUMBER OF SIMULATIONS
nsim <- 100

# GRID
summary(data$tmean)
at <- seq(0,24.3,by=0.1)

# CREATE THE OBJECT TO STORE THE INFO
res <- matrix(NA,nsim,6,dimnames=list(seq(nsim),
  c("est_min","bias","covered","true_at_CI_edge","boundary","est_SE")))

# RUN THE LOOP
set.seed(13041975)
for(i in seq(nsim)) {
  
  cat(i,"")
  
  # GENERATE SIMULATED DATA
  deathsim <- rpois(length(deathpred),deathpred)
  
  # FIT THE MODEL
  msim <- glm(deathsim~b,family=poisson,data)
  
  # FIND MIN, SAMPLE OF BOOTSTRAP ESTIMATES, AND CI AND SE FROM THOSE
  minsim   <- findmin(b,msim,at=at)
  minsimbs <- findmin(b,msim,,at=at,sim=T)
  mincisim <- quantile(minsimbs,c(2.5,97.5)/100)
  minsesim <- sd(minsimbs)
  
  # STORE THE DATA
  res[i,1] <- minsim
  res[i,2] <- min-minsim
  res[i,3] <- mincisim[1]<=min & mincisim[2]>=min
  res[i,4] <- mincisim[1]==min | mincisim[2]==min
  res[i,5] <- any(mincisim == range(at))
  res[i,6] <- minsesim
}

################################################################################
# ASSESS RESULTS

# BIAS
mean(res[,2])

# COVERAGE
mean(res[,3])*100

# % OF CIS WITH TRUE MIN AT EDGE (DISCRETE SPACE ISSUE)
mean(res[,4])*100

# PERCENTAGE OF CI AT BOUNDARY OF X SPACE
mean(res[,5])*100

# TRUE SAMPLING SD OF ESTIMATED MINIMA (EMPIRICAL SE)
# AND DISTRIBUTION OF ESTIMATED SES 
sd(res[,1])
summary(res[,6])
mean(res[,6])



#
