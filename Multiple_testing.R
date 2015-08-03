# Inspired by:
# http://www.r-bloggers.com/using-and-abusing-data-visualization-anscombes-quartet-and-cheating-bonferroni/

library(dplyr)
library(ggplot2)

set.seed(10)
nfeats = 20
n = 10

#factor(rep(1:nfeats, each = n))

# Generate n observations each with nfeats features and with label vector y.
# Generate both x and y randomly from Normal distributions (this is the null hypothesis: Gaussianity)
data = data_frame(
  # Create a vector defining sets: each set is a list of all the observations for a given feature
  set = factor(rep(1:nfeats, each = n)),
  # Generate the x data values
  x = rnorm(n * nfeats),
  # Generate the y labels (one for each observation x)
  y = rep(rnorm(n), nfeats))
data

# Create a matrix of scatter plots plotting each feature (from all the observations) against y
#  (say we are hunting for the features which might correlate with y)
ggplot(data, aes(x, y)) + geom_point() + geom_smooth(method = lm) + facet_wrap(~set)

# Now let's (wrongly!!) cherrypick the ones which look like they have some (anti-)correlation
# with set.seed(10), nfeats = 20 and n = 10 these are for example features 12 and 17
pickedset1 = data[data$set==12,]
pickedset2 = data[data$set==17,]

mod1 = lm(y ~ x, data=pickedset1)
mod2 = lm(y ~ x, data=pickedset2)

print(anova(mod1))
summary(mod1)

print(anova(mod2))
summary(mod2)

# We see that y anti-correlates with feature 12 at level p=0.0555
# and y correlates with feature 17 at level p=0.0119 (!!)

# HOWEVER WE KNOW THIS MUST BE WRONG SINCE WE CREATED THE DATA SET BY SAMPLING
# FROM THE NULL HYPOTHESIS DISTRIBUTION! (IN THIS CASE ARBITRARILY CHOSEN TO BE GAUSSIANS)

# The correct way to analyse this data set is to correct for multiple testing
# This will rescale the p values to compensate for the fact that a correlation could be found
#     purely randomly in any of the features we have. The more features we have, this is more likely.

# Possibilities are: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"

# Use dplyr's pipes to do the analysis in 'one' line
?"%>%"
results = data %>%
  group_by(set) %>%
  do(mod = lm(y ~ x, data = .)) %>%
  summarize(set = set, p = anova(mod)$"Pr(>F)"[1]) %>%
  mutate(bon = p.adjust(p, method = "bonferroni")) %>%
  mutate(holm = p.adjust(p, method = "holm")) %>%
  mutate(fdr = p.adjust(p, method = "fdr"))
results

# After correcting for multiple comparisons we see that correlations with previously
# significant (small) p values for (anti-)correlation with y are now much larger and close to 1
# and no (anti-)correlation is found at any remotely statistically significant level
# (i.e. all are completely consistent with the null hypothesis, as expected)

# In particular, for the two features we had singled out:
results %>% filter(set %in% c(12, 17))

# To see that all features are indeed consistent with the null hypothesis we make a QQ of p-values
# If they are uniformly distributed then the null hypothesis (a Gaussian distribution) holds. This
#   is due to the Probability Integral Transform:
#   p-values are the CDF of a Gaussian, therefore if the distribution of the input values (data) is 
#   Gaussian as well, then the p-values should be uniformly distributed.
#   F_X(X) = U   where F_X is the CDF of the function f, according to which X are distributed.

# p-values should be distributed uniformly if the null hypothesis holds
# Alternative proof:
# p-value = int_x^infty 1/sqrt(2pi) e^(-x^2)
# cdf of p-vales is integral of this wrt p
# Integrand on RHS has no p-dependence => its integral is linear in p
# The only distribution for which its cdf is linear in its argument is the uniform distribution
# Therefore p-values are uniformly distributed

# http://stats.stackexchange.com/questions/52192/how-to-interpret-false-positive-from-qq-plot-in-genome-wide-association-studies

# If the observed p-values are not distributed uniformly this means our null hypothesis does not hold(?)
# If the left-most (in the log10 plots) are uniform, but the rightmost (small p) deviate significantly
#    then the features with these p-values may contain a (anti-)correlation(?)



# # We now construct the qq  by hand
# n=length(results$p)
# # We need the quantiles of our data, and those of the Uniform distribution
# # The quantiles of our data are given by dividing by the number of entries,
# #   but quantiles do not include 0 and 1. Therefore, the sorted values of x
# #   are considered the 0.5/n, 1.5/n, 2.5/n, ..., (n - 0.5)/n quantiles of the sample distribution.
# # We want the quantiles of the Uniform distribution associated with the following probabilities:
# ps = ((1:n) - 0.5)/n
# # Use qunif to get the corresponding quantiles for the Uniform distribution:
# qunif(ps)
# # Plot our QQ plot comparing the Uniform quantiles vs the data ones (i.e. the sorted data)
# # Since it's the smallest p values we are interested in take a log to make them the largest points
# plot(-log10(qunif(ps)),-log10(sort(results$p)))
# abline(a=0, b=1, col='red')




# Another way to do the QQ plot, this time with confidence intervals included as well:
#  http://www.gettinggeneticsdone.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
obs <- -log10(results$p)
N <- length(obs) ## number of p-values

## create the null distribution
## (-log10 of the uniform)
null <- -log10(ppoints(N))
#null <- -log10(1:N/N)  # As given on website, incorrect since quantiles should have no 0 or 1
#null <- -log10(((1:N) - 0.5)/N)
#?runif
#null <- -log10(runif(N, min(obs), max(obs))) # We want exact unif dist, don't use this
MAX <- max(c(obs,null))

## create the confidence intervals
# creates vectors of size N, filled with zeros
c95 <- rep(0,N)
c05 <- rep(0,N)

## the jth order statistic from a
## uniform(0,1) sample
## has a beta(j,n-j+1) distribution
## (Casella & Berger, 2002,
## 2nd edition, pg 230, Duxbury, Example 5.4.5)

for(i in 1:N){
  c95[i] <- qbeta(0.95,i,N-i+1)
  c05[i] <- qbeta(0.05,i,N-i+1)
}

## plot the two confidence lines
plot(null, -log10(c95), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
par(new=T)
plot(null, -log10(c05), ylim=c(0,MAX), xlim=c(0,MAX), type="l",
     axes=FALSE, xlab="", ylab="")
par(new=T)

# add shading in the confidence intervals
# First arg: x values of boundaries, forwards and reverse for bottom and top lines
polygon(c(rev(null), null), c(rev(-log10(c95)), -log10(c05)), col = 'grey80', border = NA)
par(new=T) # Make sure these go into the next plot, not plotted now by themselves

## add the diagonal
abline(a=0,b=1,col="red")
par(new=T)

#if(any(results$fdr<=0.05))
#sort(results$fdr)[1:5]

# Take -log10 of p-values, order them ascending, this is their order in plot, set column is labels
labels = results %>%  do(transform(. , p = -log10(p))) %>% arrange(p) %>% .[,1] %>% sapply(., as.character)
#labels = sapply(results$set, as.character)
results %>%  do(transform(. , p = -log10(p))) %>% arrange(p)

text(null, sort(obs,decreasing=TRUE), rev(labels), pos=1, cex=0.4, col="red") 
par(new=T)

## add the qqplot
qqplot(null,obs, ylim=c(0,MAX),xlim=c(0,MAX),
       main="QQ plot of p-values against Uniform Distribution",
       xlab=expression("-log"[10]*"(uniform) quantiles"), ylab=expression("-log"[10]*"(p-values), NOT quantiles"))

# WE SEE THAT ALL FEATURES LIE WITHIN THE CONFIDENCE INTERVAL, EXACTLY AS EXPECTED HERE
#    SINCE WE KNOW THE NULL HYPOTHESIS TO BE TRUE