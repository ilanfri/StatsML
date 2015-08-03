#!/usr/bin/python

# IMPLEMENTS THE CLASSICAL BOOTSTRAP METHOD AND USES IT TO ESTIMATE
# A SECOND-ORDER ERROR: THE STANDARD ERROR ON THE STANDARD DEVIATION
# OF THE MEAN

import numpy as np
import numpy.random as npr
#import pylab
import matplotlib.pyplot as plt
import scipy
import scikits.bootstrap as btstrp

#from pandas.tools.plotting import bootstrap_plot


def bootstrap(data, num_samples, statistic, alpha):
    """Returns bootstrap estimate of 100.0*(1-alpha) CI for statistic."""
    # Generate the indices for the required number of permutations/(resamplings with replacement) required
    idx = npr.randint(0, len(data), (num_samples, len(data)))
    # Generate the multiple resampled data set from the single original one
    samples = data[idx]
    # Apply the 'statistic' function given to each of the data sets produced by the resampling and order the resulting statistic by decreasing size
    stat = np.sort(statistic(samples, 1))
    # Return the value of the computed statistic at the upper and lower percentiles specified by the alpha parameter given. These are, by definition, the boundaries of the Confidence Interval for that value of alpha. E.g. alpha=0.05 --> CI 95%
    return (stat[int((alpha / 2.0) * num_samples)],
            stat[int((1 - alpha / 2.0) * num_samples)])


if __name__ == '__main__':
    # data of interest is bimodal and obviously not normal
    x = np.concatenate([npr.normal(3, 1, 100), npr.normal(6, 2, 200)])

    print "Original sample's mean: ", np.mean(x)
    print "Original sample's SD: ", np.std(x)

    nbootstraps = 10000
    alpha = 0.05

    # find mean 95% CI and 100,000 bootstrap samples
    lowmean, highmean = bootstrap(x, nbootstraps, np.mean, alpha)
    lowsd, highsd = bootstrap(x, nbootstraps, np.std, alpha)

    print "Bootstrapped 95% confidence intervals on the mean\nLow:", lowmean, "\nHigh:", highmean
    print "Bootstrapped 95% confidence intervals on the SD\nLow:", lowsd, "\nHigh:", highsd

    # Scipy implementation of the Bootstrap, for validation purposes
    # compute 95% confidence intervals around the mean 
    # CIsmean = btstrp.ci(data=x, statfunction=np.mean, alpha=alpha, n_samples=nbootstraps)
    # print "Scipy Bootstrapped 95% confidence intervals on the mean\nLow:", CIsmean[0], "\nHigh:", CIsmean[1]

    # CIsSD = btstrp.ci(data=x, statfunction=np.std, alpha=alpha, n_samples=nbootstraps)
    # print "Scipy Bootstrapped 95% confidence intervals on the SD\nLow:", CIsSD[0], "\nHigh:", CIsSD[1]

    # make plots
    plt.figure(figsize=(8, 4))
    plt.subplot(121)
    plt.hist(x, 50, histtype='step')
    plt.title('Histogram of data')
    plt.subplot(122)
    plt.plot([-0.03, 0.03], [np.mean(x), np.mean(x)], 'r', linewidth=2)
    plt.scatter(0.1 * (npr.random(len(x)) - 0.5), x)
    plt.plot([0.19, 0.21], [lowmean, lowmean], 'r', linewidth=2)
    plt.plot([0.19, 0.21], [highmean, highmean], 'r', linewidth=2)
    plt.plot([0.2, 0.2], [lowmean, highmean], 'r',
             linewidth=2,
             label="Mean CI")

    # The SD interval
    plt.plot([0.25, 0.27], [np.mean(x) - np.std(x),
                            np.mean(x) - np.std(x)], 'g',
             linewidth=2)
    plt.plot([0.25, 0.27], [np.mean(x) + np.std(x),
                            np.mean(x) + np.std(x)], 'g',
             linewidth=2)
    plt.plot([0.26, 0.26], [np.mean(x) - np.std(x),
                            np.mean(x) + np.std(x)], 'g',
             linewidth=2)

    # CI bars on each of the downwards and upwards 1 SD bands about the mean
    plt.plot([0.29, 0.31], [np.mean(x) - np.std(x) - (np.std(x) - lowsd),
                            np.mean(x) - np.std(x) - (np.std(x) - lowsd)], 'g',
             linewidth=2)
    plt.plot([0.29,
              0.31], [np.mean(x) - np.std(x) + (highsd - np.std(x)),
                      np.mean(x) - np.std(x) + (highsd - np.std(x))], 'g',
             linewidth=2)
    plt.plot([0.3, 0.3], [np.mean(x) - np.std(x) - (np.std(x) - lowsd),
                          np.mean(x) - np.std(x) + (highsd - np.std(x))], 'g',
             linewidth=2)
    plt.plot([0.29, 0.31], [np.mean(x) + np.std(x) - (np.std(x) - lowsd),
                            np.mean(x) + np.std(x) - (np.std(x) - lowsd)], 'g',
             linewidth=2)
    plt.plot([0.29,
              0.31], [np.mean(x) + np.std(x) + (highsd - np.std(x)),
                      np.mean(x) + np.std(x) + (highsd - np.std(x))], 'g',
             linewidth=2)
    plt.plot([0.3, 0.3], [np.mean(x) + np.std(x) - (np.std(x) - lowsd),
                          np.mean(x) + np.std(x) + (highsd - np.std(x))], 'g',
             linewidth=2,
             label="SD CI")

    plt.xlim([-0.2, 0.35])
    plt.title('Bootstrap 95% CIs')
    plt.legend()
    plt.show()
    #plt.savefig('examples/boostrap.png')
