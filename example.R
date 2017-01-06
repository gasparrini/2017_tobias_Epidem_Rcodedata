################################################################################
# Updated version of the R code for an example of the analysis in:
#
#   "Investigating uncertainty in the minimum temperature mortality:
#     methods and application to 52 Spanish cities"
#   Armstrong B, Tobias A, Gasparrini A
#
# Update: 18 October 2016
# For any problem with this code, please contact antonio.gasparrini@lshtm.ac.uk
# Please refer to the original code for any copyright issue
#
#  See www.ag-myresearch.com for future updates
################################################################################

################################################################################
# EXAMPLES OF ESTIMATION OF THE MMT
################################################################################

library(dlnm) ; library(splines)
source("findmin.R")

################################################################################
# FIRST EXAMPLE: SUBSET OF DATA, SIMPLE MODEL, UNIDIMENSIONAL EXPOSURE-RESPONSE

# LOAD A SUBSET OF THE DATA: LONDON JULY-DEC 2005
# (CONVENIENCE SAMPLE, LARGISH ASIMMETRIC CI)
data <- subset(read.csv("london.csv"), year==2005 & month>6)
head(data)

# A UNIDIMENSIONAL 4DF SPLINE 
b1 <- onebasis(data$tmean,df=4)

# SIMPLE MODEL WITH NO CONTROL FOR CONFOUNDING
m1 <- glm(death~b1,family=quasipoisson,data)

# ESTIMATE MMT, WITH CI AND STANDARD ERROR
(min1 <- findmin(b1,m1))
(min1ci <- quantile(findmin(b1,m1,sim=T),c(2.5,97.5)/100))
(min1se <- sd(findmin(b1,m1,sim=T)))

# ESTIMATE THE MINIMUM WITHIN A SPECIFIED RANGE (15-16: MEANINGLESS ILLUSTRATION)
(min1b <- findmin(b1,m1,from=15,to=16))

# PLOT
plot(crosspred(b1,m1),ylab="RR",xlab="Temperature",xlim=c(0,25),
  ylim=c(0.9,1.3),lwd=1.5)
abline(v=min1)
abline(v=min1ci,lty=2)

################################################################################
# SECOND EXAMPLE: WHOLE DATA, FULL MODEL, BI-MENSIONAL EXPOSURE-LAG-RESPONSE

# LOAD DATA: LONDON 1993-2006
data <- read.csv("london.csv")
head(data)

# BI-DIMENSIONAL EXPOSURE-LAG-RESPONSE SPLINE
vk <- equalknots(data$tmean,fun="bs",df=4,degree=2)
lk <- logknots(25,3)
cb2 <- crossbasis(data$tmean, lag=25, argvar=list(fun="bs",degree=2,knots=vk),
  arglag=list(knots=lk))

# FULL MODEL WITH CONTROL FOR CONFOUNDING
m2 <- glm(death~cb2+ns(time,10*14)+dow,family=quasipoisson(),data)

# ESTIMATE MMT, WITH CI AND STANDARD ERROR
(min2 <- findmin(cb2,m2))
(min2ci <- quantile(findmin(cb2,m2,sim=T),c(2.5,97.5)/100))
(min2se <- sd(findmin(cb2,m2,sim=T)))

# IN PERCENTILE SCALE
sum(data$tmean<min2)/nrow(data)*100

# PLOT 
cb2 <- crossbasis(data$tmean, lag=25, argvar=list(fun="bs",degree=2,knots=vk),
  arglag=list(knots=lk))
plot(crosspred(cb2,m2,cen=min2),"overall",ylab="RR",xlab="Temperature",
  xlim=c(-5,35),ylim=c(0.5,3.5),lwd=1.5)
abline(v=min2)
abline(v=min2ci,lty=2)

#
