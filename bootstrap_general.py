#!/usr/bin/python

# IMPLEMENTS THE CLASSICAL BOOTSTRAP METHOD FOR NON-PARAMETRIC ESTIMATION
#    TESTS THE IMPLEMENTATION ON A SMALL DATA SAMPLE

import numpy as np
import numpy.random as npr
#import pylab
import matplotlib.pyplot as plt
import scipy
import scikits.bootstrap as btstrp
import pandas as pan
import os, shutil
#from pandas.tools.plotting import bootstrap_plot
from scipy import stats

#import pickle_csv as pcsv


def bootstrap(data, num_samples, statistic, alpha):
    """Returns the results from num_samples bootstrap samples for an input test statistic, its standard deviation, and its 100*(1-alpha)% confidence level interval."""
    # Generate the indices for the required number of permutations/(resamplings with replacement) required
    idx = npr.randint(0, len(data), (num_samples, len(data)))
    # Generate the multiple resampled data set from the original one
    samples = data[idx]
    # Apply the 'statistic' function given to each of the data sets produced by the resampling and order the resulting statistic by decreasing size.
    stats = np.sort(statistic(samples, 1))
    stat = stats.mean()
    # Return the value of the computed statistic at the upper and lower percentiles specified by the alpha parameter given. These are, by definition, the boundaries of the Confidence Interval for that value of alpha. E.g. alpha=0.05 ==> CI 95%
    low_ci = stats[int((alpha / 2.0) * num_samples)]
    high_ci = stats[int((1 - alpha / 2.0) * num_samples)]

    #sd = np.std(stat)
    # To include Bessel's correction for unbiased standard deviation:
    sd = np.sqrt(len(data) / (len(data) - 1)) * np.std(stats)

    return stat, sd, low_ci, high_ci


if __name__ == '__main__':

    filename = '2008_test.csv'
    na_symbol = 'NA'

    # Read input data file in, create serialised copy if none exists
    #pklfilename=(str(filename).replace('.csv',''))+'.pkl'
    pklfilename = (str(filename).replace('.csv', '')) + '.h5'
    rootdir = os.getcwd()
    filepath = os.path.join(rootdir, pklfilename)
    if (os.path.exists(filepath) == True):
        #data=pan.read_pickle(pklfilename)
        data = pan.read_hdf(pklfilename, 'data')
    else:
        data = pan.read_csv(filename, na_values=na_symbol)
        #data.to_pickle(pklfilename)
        data.to_hdf(pklfilename, 'data', mode='w')

    # Choose two arbitrary (mutually exclusive) subsets of data for statistical comparison
    testdist1 = data["ArrDelay"].dropna()[data["Distance"] < 1000]
    testdist2 = data["ArrDelay"].dropna()[(data["Distance"] >= 1000) &
                                          (data["Distance"] < 2000)]
    #testdist2=data[data["Distance"]>=1000]["ArrDelay"]
    z, pmww = stats.ranksums(testdist1, testdist2)
    print "\n\nThe MWW RankSum p-value for flights less than 1000 miles vs [1000,2000) miles:\n{0:.3g}\n".format(
        pmww)

    #testdist3=data[data["Distance"]>2000]["ArrDelay"]
    testdist3 = data["ArrDelay"].dropna()[data["Distance"] >= 2000]
    f, panova = stats.f_oneway(testdist1, testdist2, testdist3)
    print "One-way ANOVA p-value for <1000, [1000,2000), and >= 2000 miles flight data sets:\n{0:.3g}".format(
        panova)

    nbootstraps = 10000
    alpha = 0.05

    # Observed sample mean and standard deviation
    print "\nOriginal-sample mean:\n{0:.3g}".format(testdist3.mean())
    print "Original-sample standard deviation:\n{0:.3g}".format(
        stats.sem(testdist3))

    # Compute bootstrap mean, SD of mean, and CI of mean
    meanbtstrp, meanbtstrpsd, lowmean, highmean = bootstrap(
        np.array(testdist3), nbootstraps, np.mean, alpha)

    print "\nBootstrap mean:\n{0:.3g}".format(meanbtstrp)
    print "Bootrstrap SD of the mean:\n{0:.3g}".format(meanbtstrpsd)
    print "\nBootstrapped 95% confidence intervals on the mean:\n[{0:.3g}, {1:.3g}]".format(
        lowmean, highmean)

    # Use Scipy's implementation of bootstrapping to do the same
    CIs = btstrp.ci(data=testdist3,
                    statfunction=np.mean,
                    n_samples=nbootstraps,
                    alpha=alpha)
    print "Scipy Bootstrapped 95% confidence intervals on the mean:\n[{0:.3g}, {1:.3g}]\n\n".format(
        CIs[0], CIs[1])

    # Options to deal with NA entries:
    # Drop all of the entries which contain a NA: data.dropna()
    # Drop NA entries in a particular column: data["Distance"].dropna()
    # Replace NA with a given value: data.fillna(0.0)["Distance"]

    # Fields are:
    # Year,Month,DayofMonth,DayOfWeek,DepTime,CRSDepTime,ArrTime,CRSArrTime,UniqueCarrier,FlightNum,TailNum,ActualElapsedTime,CRSElapsedTime,AirTime,ArrDelay,DepDelay,Origin,Dest,Distance,TaxiIn,TaxiOut,Cancelled,CancellationCode,Diverted,CarrierDelay,WeatherDelay,NASDelay,SecurityDelay,LateAircraftDelay

    # # make plots
    # plt.figure(figsize=(8,4))
    # plt.subplot(121)
    # plt.hist(x, 50, histtype='step')
    # plt.title('Histogram of data')
    # plt.subplot(122)
    # plt.plot([-0.03,0.03], [np.mean(x), np.mean(x)], 'r', linewidth=2)
    # plt.scatter(0.1*(npr.random(len(x))-0.5), x)
    # plt.plot([0.19,0.21], [lowmean, lowmean], 'r', linewidth=2)
    # plt.plot([0.19,0.21], [highmean, highmean], 'r', linewidth=2)
    # plt.plot([0.2,0.2], [lowmean, highmean], 'r', linewidth=2, label="Mean CI")

    # # The SD interval
    # plt.plot([0.25,0.27], [np.mean(x)-np.std(x), np.mean(x)-np.std(x)], 'g', linewidth=2)
    # plt.plot([0.25,0.27], [np.mean(x)+np.std(x), np.mean(x)+np.std(x)], 'g', linewidth=2)
    # plt.plot([0.26,0.26], [np.mean(x)-np.std(x), np.mean(x)+np.std(x)], 'g', linewidth=2)

    # # CI bars on each of the downwards and upwards 1 SD bands about the mean
    # plt.plot([0.29,0.31], [np.mean(x)-np.std(x)-(np.std(x)-lowsd), np.mean(x)-np.std(x)-(np.std(x)-lowsd)], 'g', linewidth=2)
    # plt.plot([0.29,0.31], [np.mean(x)-np.std(x)+(highsd-np.std(x)), np.mean(x)-np.std(x)+(highsd-np.std(x))], 'g', linewidth=2)
    # plt.plot([0.3,0.3], [np.mean(x)-np.std(x)-(np.std(x)-lowsd), np.mean(x)-np.std(x)+(highsd-np.std(x))], 'g', linewidth=2)
    # plt.plot([0.29,0.31], [np.mean(x)+np.std(x)-(np.std(x)-lowsd), np.mean(x)+np.std(x)-(np.std(x)-lowsd)], 'g', linewidth=2)
    # plt.plot([0.29,0.31], [np.mean(x)+np.std(x)+(highsd-np.std(x)), np.mean(x)+np.std(x)+(highsd-np.std(x))], 'g', linewidth=2)
    # plt.plot([0.3,0.3], [np.mean(x)+np.std(x)-(np.std(x)-lowsd), np.mean(x)+np.std(x)+(highsd-np.std(x))], 'g', linewidth=2, label="SD CI")

    # plt.xlim([-0.2, 0.35])
    # plt.title('Bootstrap 95% CIs')
    # plt.legend()
    # plt.show()
    # #plt.savefig('examples/boostrap.png')
