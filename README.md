# StatsML
Python and R scripts related mostly to Statistics and Machine Learning

* Bayesian_bootstrap.ipynb: A Python implementation of Rasmus Baath's (@rasmusab) R Bayesian bootstrap code.

* lead_generation.ipynb: A prototype for a Business-to-Business lead generation recommendation engine. Uses a database of companies (the Companies House listing) to search for companies by name (using PyGoogle), identify their 'About Us' webpage, extract the body of the text, preprocess it, and perform topic modelling on it, such that a query can be made and LDA topic projections can be used to identify potential client companies from the companies analysed.

* portfolios.R: Given a fixed total investment and fractions into a variety of equities (all denominated in a single currency and with tickers available from Google Finance or Yahoo Finance) for a given timeframe it creates plots of the overall and single-stock performances over that period, as well as evaluating the CAPM parameters.

* statistical_significance.py: Computes the statistical significance of a Poisson process in the presence of a Gaussian uncertainty in its measurement by three different methods, including profile likelihood. For a given number of measured events, with a background process with a known Gaussian error, it can also provide an estimated upper limits on the number of potential signal events which would be excluded at the 95% C.L.

* Breit-Wigner_mapping.py: A demonstration of a variable transformation to efficiently sample from a Breit-Wigner distribution.

* KDD_1999_network_intrusion_ML_analysis.ipynb: A simple machine learning model selection exercise using the KDD 1999 network intrusion data set. Performs cross-validation with care to preserve class prevalence in the folds, and produces a full set of performance metrics and confusion matrix for a set of classifiers.

* Central_Limit_Theorem_CLT.R: A demonstration of the Central Limit Theorem in action, takes the mean of samples from a mixture of distributions and compares it to a the Gaussian it should asymptote to under the Central Limit Theorem in the limit of infinite sample size.

* Multiple_testing.R: Demonstrates the importance of applying multiple-comparison correction factors (e.g. Bonferroni factor a.k.a. trials factor) in the case where correlations are being examined in a multi-variable data set.

* bootstrap_general.py: A simple Python implementation of Bootstrap sampling for non-parametric estimation.

* bootstrap_of_SD.py: An application of Bootstrap sampling for second-order error estimation (estimating the error on the standard deviation of the mean).

* confidence_interval_percentages.py: Extremely simple simulation to compute the confidence interval percentages corresponding to 1 standard deviation for various distributions (e.g. the well-known ~68% CI corresponding to +/- 1 SD for the Gaussian).

* circle through point: A simple script to compute and plot a semi-circle through a given point in the (+,+) quadrant of a plane (used originally for a phase diagram plot).
