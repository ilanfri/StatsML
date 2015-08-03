__author__ = "Ilan Fridman Rojas"

# Please report bugs/suggestions/comments and any other feedback to:
#        ilanfri@mac.com

from numpy import log, sqrt, inf, linspace
from scipy.stats import poisson, norm, chi2
from scipy.special import erfcinv, betainc, erfc, gammainc
from scipy.optimize import minimize
from scipy.integrate import quad
import argparse
import matplotlib.pyplot as plt

# BASIC USAGE:
# python significance.py -n 10 -b 3 --berror 1.3

parser = argparse.ArgumentParser(
    description=
    'Compute the statistical significance (p-value and number of sigmas) for an observation of n events, when a known background of b events (with error berror) was expected.\n\nAlso compute the corresponding upper bound on number of signal events allowed at the 100*(1-alpha)% confidence level (default: alpha=0.05, 95% CL) using both the standard confidence levels (CL_(s+b)), and the CL_s prescription.\n\nIf a number of signal events s is specified, the s+b model itself instead becomes the null hypothesis and the statistical significance computed is that with which the background-only hypothesis can be rejected.')
parser.add_argument('-n', '--observed',
                    help='Observed number of events',
                    required=True,
                    type=int)
parser.add_argument('-b', '--background',
                    help='Expected number of background events',
                    required=True,
                    type=float)
parser.add_argument('-s', '--signal',
                    help='Expected number of signal events',
                    required=False,
                    type=float,
                    default=0)
parser.add_argument('--berror',
                    help='Error on expected number of background events',
                    required=False,
                    default=0.,
                    type=float)
parser.add_argument(
    '--alpha',
    help=
    'Specify the significance (type-I error rate) desired for upper bounds on the signal',
    required=False,
    default=0.05,
    type=float)

#parser.add_argument('--trialsfactor', help = 'The trials factor/Bonferroni correction for multiple comparisons testing', required = False, default=1.,type=float)

# Based on:
# [1] arxiv:physics/0702156
# [2] arxiv:1007.1727

# References to equations are from paper [1] unless otherwise stated.

# All three methods (Z_(Bi), Z_N, Z_(PL) ) of computing statistical significance validated against table Table 1 of [1]. Very good agreement found.
#        (n_on = our n, hat mu_b = our b, sigma_b = our eb).

# Further discussion of three methods used here to compute statistical significance of Poisson process with Gaussian error on the background:
# arXiv:physics/0312059
# arXiv:physics/0511028

# Brief description of profile Likelihood method: section 33.3.2.3 of PDG's statistics section:
#     http://pdg.lbl.gov/2011/reviews/rpp2011-rev-statistics.pdf

# Of the methods used here ATLAS and CMS use the profile Likelihood (PL) method, but with the Gaussian replaced by complicated priors dependent on many nuisance parameters used to introduce systematic uncertainties. For details see:
#     https://cds.cern.ch/record/1375842/files/ATL-PHYS-PUB-2011-011.pdf 
#     arXiv:1108.2288

# Sources on CLs method:
# www.pp.rhul.ac.uk/~cowan/stat/cls/CLsInfo.ps
# Brief summary: Section 33.3.2.2 from PDG section on statistics
#          http://pdg.lbl.gov/2011/reviews/rpp2011-rev-statistics.pdf
# http://iopscience.iop.org/0954-3899/28/10/313/


def TestStatisticPL(non, s, estmub, eb, test):
    def L(mus, mub, non, estmub, eb):
        return poisson.pmf(non,
                           mu=mub + mus) * norm.pdf(x=estmub,
                                                    loc=mub,
                                                    scale=eb)
    # Define function to be minimised over mu_b and mu_s
    def max_sb_logL(c):
        return -2 * log(poisson.pmf(non,
                                    mu=c[0] + c[1]) * norm.pdf(x=estmub,
                                                               loc=c[0],
                                                               scale=eb))
    # Define function to be minimised over mu_b only
    def max_b_logL(c):
        return -2 * log(poisson.pmf(non,
                                    mu=c + s) * norm.pdf(x=estmub,
                                                         loc=c,
                                                         scale=eb))
    # Minimise over mub and mus to get arguments for numerator of likelihood ratio (eqn 22)
    result = minimize(max_sb_logL, [estmub, s])
    maxmub = result.x[0]
    maxmus = result.x[1]

    # Minimise over mub only to get arguments for denominator of likelihood ratio (eqn 22)
    result1 = minimize(max_b_logL, estmub)
    maxmubonly = result1.x[0]

    # # Plot the Likelihood as a function of mus
    # from numpy import linspace
    # plt.figure()
    # x = linspace(0,b,100)
    # plt.plot(x, L(x,b,non,b,eb) )
    # plt.show()
    # # Plot the Likelihood as a function of mub
    # plt.figure()
    # x = linspace(b-eb,b+eb,100)
    # plt.plot(x, L(0,x,non,b,eb) )
    # plt.show()

    if (test == 'q0'):  # eqn 12 of [2]
        if (maxmus >= 0):
            result = -2. * log(L(0, maxmubonly, non, estmub, eb) /
                               L(maxmus, maxmub, non, estmub, eb))
        elif (maxmus < 0):
            result = 0.
    if (test == 'qmu'):  # eqn 14 of [2]
        if (maxmus <= s):
            result = -2. * log(L(s, maxmubonly, non, estmub, eb) /
                               L(maxmus, maxmub, non, estmub, eb))
        elif (maxmus > s):
            result = 0.

    return result


def main():

    args = parser.parse_args()
    n = args.observed
    b = args.background
    s = args.signal
    eb = args.berror
    alpha = args.alpha
    # TO DO: implement trials factors for local p-value --> global p-value (arXiv:1005.1891)
    #    nt = args.trialsfactor

    if (eb ==
        0):  # For speed if error is zero compute Poisson p-value directly
        if (b > 0 and s == 0):
            pval = poisson.sf(n - 1, mu=b)
            print "\nBackground-only null hypothesis p-value = {0:.4g}\t({1:.3g} sigmas)\n".format(
                pval, sqrt(2.) * erfcinv(2. * pval))
        elif (b > 0 and s > 0):
            pval = poisson.sf(n - 1, mu=b + s)
            print "\nBackground+signal null hypothesis \tp-value = {0:.4g}\t({1:.3g} sigmas)\n".format(
                pval, sqrt(2.) * erfcinv(2. * pval))
            print "(Significance according to ATLAS prescription = {0:.3g} sigmas)\n".format(
                sqrt(2 * ((s + b) * log(1 + s / b) - s))
            )  # eqn 97 of [2]

    if (b >= 0 and eb > 0 and s == 0.):

        def computeZBi(n, b, eb):
            #pBi
            tau = b / (eb ** 2)  #eqn 8
            noff = tau * b  #eqn 6
            rho = 1 / (1 + tau)  #above eqn 13

            pBi = betainc(
                n, 1 + noff, rho
            )  # eqn 14 (scipy's beta incomplete already includes the factor of 1/beta)
            ZBi = sqrt(2.) * erfcinv(2. * pBi)
            return pBi, ZBi

        pBi, ZBi = computeZBi(n, b, eb)
        print "\np_(Bi) = {0:.4g},\tp_(Bi) sigmas = {1:.3g}".format(pBi, ZBi)

        # The Bayesian-frequentist hybrid approach takes a weighted average of the p-values: Poisson likelihood times Gaussian prior. Evaluate the Poisson p-value using background values sampled from a Gaussian, weight each one with the Gaussian probability for that background value, and sum (see eqns 16 & 17)
        def computeZN(n, b, eb):
            # Eqns 16 and 17:
            pdensity = lambda x: gammainc(n, x) * norm.pdf(x, loc=b, scale=eb)
            pN = quad(pdensity, 0, inf)[0]
            #print quad(pdensity, b, inf)[1] # The error on the integration
            ZN = sqrt(2.) * erfcinv(2. * pN)
            if (ZN >= b / eb):
                print "(Warning: Z*eb > b : this value of Z_N is unreliable)"
            return pN, ZN

        pN, ZN = computeZN(n, b, eb)
        print "\np_N    = {0:.4g},\tp_N sigmas    = {1:.3g}".format(pN, ZN)

        # The profile Likelihood method if the background were Poissonian and not using the Gaussian approximation, we use the result in eqn 25
        #ntot = non + noff
        #ZPL = sqrt( 2*(non*log((non*(1+tau)/ntot))  + noff*log((noff*(1+tau))/ntot*tau) )   )
        #print ZPL

        ZPL = sqrt(TestStatisticPL(n, s, b, eb,
                                   test='q0'))  # eqn 24 of [1], eqn 53 of [2]
        pPL = 0.5 * erfc(ZPL / sqrt(2))
        print "\np_(PL) = {0:.4g},\tp_(PL) sigmas = {1:.3g}\n".format(pPL, ZPL)

        # Work out the upper limit on the signal
        clsignalbound = 0
        for i in xrange(0, int(10 * n)):
            pval = 1. - norm.cdf(sqrt(TestStatisticPL(n, i, b, eb,
                                                      test='qmu'))
                                 )  # eqn 59 of [2]
            if (pval >= alpha):
                clsignalbound = i
            else:
                break

        clssignalbound = 0
        bpval = 1. - norm.cdf(sqrt(TestStatisticPL(n, 0, b, eb, test='qmu')))
        #bpval = pPL
        for i in xrange(0, int(10 * n)):
            sbpval = 1. - norm.cdf(sqrt(TestStatisticPL(n, i, b, eb,
                                                        test='qmu'))
                                   )  # eqn 59 of [2]

            cls = sbpval / (bpval)
            if (cls >= alpha):
                clssignalbound = i
            else:
                break

        print "Up to {0:.4g} signal events allowed at {1:.0f}% CL\nUp to {2:.4g} signal events allowed at {1:.0f}% CL (with the CLs prescription)\n".format(
            clsignalbound, 100. * (1 - alpha), clssignalbound)

    if (b > 0 and eb > 0 and s > 0.):

        ZPL = sqrt(TestStatisticPL(n, s, b, eb,
                                   test='q0'))  # eqn 24 of [1], eqn 53 of [2]
        pPL = 0.5 * erfc(ZPL / sqrt(2))

        print "\nAssuming {0:.3g} background events and {1:.3g} signal events the (profile likelihood) statiscal significance of observing {2:.3g} events under the background-only null hypothesis is:\n\np-value = {3:.4g}\t sigmas = {4:.3g}\n".format(
            b, s, n, pPL, ZPL)

        # Signal event counts allowed at 100*(1-alpha) CLs level are those for which CLs < alpha
        # Work out the upper limit on the signal
        clsignalbound = 0
        for i in xrange(0, int(10 * n)):
            pval = 1. - norm.cdf(sqrt(TestStatisticPL(n, i, b, eb,
                                                      test='qmu'))
                                 )  # eqn 59 of [2]
            if (pval >= alpha):
                clsignalbound = i
            else:
                break

        clssignalbound = 0
        bpval = 1. - norm.cdf(sqrt(TestStatisticPL(n, 0, b, eb, test='qmu')))
        for i in xrange(0, int(10 * n)):
            sbpval = 1. - norm.cdf(sqrt(TestStatisticPL(n, i, b, eb,
                                                        test='qmu'))
                                   )  # eqn 59 of [2]
            cls = sbpval / (bpval)
            if (cls >= alpha):
                clssignalbound = i
            else:
                break

        print "Up to {0:.4g} signal events allowed at {1:.0f}% CL\nUp to {2:.4g} signal events allowed at {1:.0f}% CL (with the CLs prescription)\n".format(
            clsignalbound, 100. * (1 - alpha), clssignalbound)


if __name__ == "__main__":
    main()

# Extra validation points (all in good agreement with present implementation):
# obs=80, ept=50, pval=5.66x10-5, 3.86 sigmas    https://indico.cern.ch/event/277650/session/15/contribution/275/material/slides/1.pdf

# obs=7, ept=3.2+-1.1 pval=0.072    http://www.desy.de/~blobel/blobel_errors.pdf

# obs=16, ept=6.4+-2.2, 2.03 sigmas  http://www.science20.com/a_quantum_diaries_survivor/spring_flukes_new_3sigma_signals_from_lhcb_and_atlas-154210

# obs=13, ept=4.2+-1.6, 2.26 sigmas  http://www.science20.com/a_quantum_diaries_survivor/spring_flukes_new_3sigma_signals_from_lhcb_and_atlas-154210

# For discussions of general issues surrounding these calculations of statistical significance, etc see:
# 1st CERN workshop on confidence intervals:
#         http://cds.cern.ch/record/411537/files/p29.pdf
#         http://cds.cern.ch/record/411537/files/CERN-2000-005.pdf?version=2
#  Phystat 2007 proceedings:
#         http://phystat-lhc.web.cern.ch/phystat-lhc/2008-001.pdf
#  Phystat 2011 proceedings:
#         http://www.hep.ph.ic.ac.uk/~dauncey/proceedings.pdf
