import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt
from scipy.stats import norm, poisson, binom, uniform

# Computes the percentage confidence intervals corresponding to one standard deviation from the mean for different distributions, to see how they differ
#      e.g. verify the well-known 1 SD = ~65% CL for a Gaussian distribution
#           and compare it to other distributions

# Choose the distribution of interest

npoints = 1E5

# Gaussian
setmean = 3
setsd = 1
x = npr.normal(setmean, setsd, npoints)
## Expected [mean-1SD,mean+1SD] CDF: The full integral (1) minus twice the tail integral between (mean + 1 SD) and infinity
expected = 1 - 2 * norm.cdf(x=setmean - setsd, loc=setmean, scale=setsd)

# Binomial
#n = 100
#p = 0.5
#x = npr.binomial(n, p, npoints)
#setmean = n*p
#setsd = np.sqrt(n*p*(1-p))
##   For one-tailed (e.g. Poisson distribution) expected [mean-1SD,mean+1SD] CDF: the tail integral from leftmost bound (setmean+setsd) minus tail integral from rightmost bound (setmean-setsd)
#expected = binom.cdf(setmean+setsd, n=n, p=p) - binom.cdf(setmean-setsd, n=n, p=p)

# Uniform
#a = 0
#b = 1
#x = npr.uniform(a,b,npoints)
#setmean = 0.5*(a+b)
#setsd = np.sqrt(1./12. * (b-a)**2)
#expected = uniform.cdf(setmean+setsd, loc=a, scale=(b-a)) - uniform.cdf(setmean-setsd, loc=a, scale=(b-a))

# Poisson
#mu = 10
#x = npr.poisson(mu,npoints)
#setmean = mu
#setsd = np.sqrt(mu)
#expected = poisson.cdf(setmean+setsd,setmean) - poisson.cdf(setmean-setsd,setmean)

mean = np.mean(x)
se = np.std(x)

inCI = 0.
for i in x:
    if i >= (mean - se) and i <= (mean + se):
        inCI += 1.

print "\nMean: {0:.3g}".format(setmean)
print "Standard error: {0:.3g}".format(setsd)
print "Percentage of points {0:.3g} sampled points within 1 SD of the mean ( [{1:.3g}, {2:.3g}] ):\n{3:.3g}%\n".format(
    npoints, mean - se, mean + se, 100 * inCI / float(len(x)))
print "Exact CDF within +/- 1 SD of mean:\n{0:.3g}%".format(100 * expected)
