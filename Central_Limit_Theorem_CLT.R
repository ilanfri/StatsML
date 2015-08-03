# THIS CODE PROVIDES A DEMONSTRATION OF THE CENTRAL LIMIT THEOREM IN ACTION:
#   IT GENERATES SYNTHETIC DATA USING A MIXTURE OF DISTRIBUTIONS AND THEN COMPUTES ITS MEAN
#   AND STANDARD DEVIATION TO SHOW THAT THE MEAN TENDS TO A GAUSSIAN DISTRIBUTION WITH 
#   WIDTH SIGMA/SQRT(N), REGARDLESS OF WHAT THE MIXTURE DISTRIBUTIONS WERE

#install.packages(c( "foreach", "doParallel","doMC") ) 
# To avoid needing these libraries comment 'foreach' line below and uncomment the for loop below it
library(foreach)
library(doMC) # These two lines fix 'no parallel backend registered' error at the foreach
registerDoMC(cores=4)

# Sample points from a mixture of different distributions (with randomly set parameters).
points <- c(runif(10000), rbeta(10000,5*runif(1),5*runif(1)), rweibull(10000,runif(1),runif(1)),rchisq(10000,5*runif(1)))

# Plot a histogram of the points to see what it looks like, it will likely be very far from Normal
hist(points)

## http://www.exegetic.biz/blog/2013/08/the-wonders-of-foreach/
## %do% = serial, %dopar% = parallel

# Repeatedly take samples of a given size and compute the mean of the sample
samplesize <- 100
numsamples <- length(points)/samplesize
#numsamples <- 1000
meanvect <- numeric()
meanvect <- foreach(n = 1:numsamples, .combine = c) %dopar% mean(sample(points,samplesize,replace=TRUE))

#This uses samples of increasing size
# meanvect <- foreach(n = 10:length(points), .combine = c) %dopar% mean(sample(points,n,replace=TRUE))

# The loop below works without the foreach package but is not parallel and therefore slower
#for(i in 10:length(points)){
#  meanvect <- c(meanvect, mean(sample(points,i,replace=FALSE)))
#}

summary(points)

# Plot a histogram of these computed means
h <- hist(meanvect,xlim=c(mean(points)-0.5*sd(points), mean(points)+0.5*sd(points)),main="Histogram of sample means",xlab="Sample means")
xgauss=seq(min(meanvect),max(meanvect),length=length(meanvect))
# (The normalisation of the Gaussian set by multiplying by bin width and number of total entries)
ygauss=dnorm(xgauss,mean=mean(points),sd=sd(meanvect)/sqrt(length(samplesize)))*diff(h$mids[1:2])*length(meanvect)
lines(xgauss,ygauss,col="darkblue",lwd=2)
abline(v=mean(points),col="red",lwd=2)

# Despite the fact that the points were sampled from non-Normal distributions, and even a mix of them
#    the means of these samples perfectly follows a Normal distribution about the 'population' mean,
#    and with a sd of sd/sqrt(n)... just as the Central Limit Theorem says it should.

# Note that the tails will take more points to settle into a Gaussian than the central parts

# Re-run this code as many times as you like and watch the CLT in action every time