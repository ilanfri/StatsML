# Author: Ilan Fridman Rojas

# Based on (i.e. useful code snippets/functions were taken from):
#  http://www.r-bloggers.com/monitoring-an-etf-portfolio-in-r/
#  http://blog.revolutionanalytics.com/2014/01/quantitative-finance-applications-in-r-plotting-xts-time-series.html

library(quantmod, warn.conflicts = FALSE, quietly = TRUE)
library(PerformanceAnalytics, warn.conflicts = FALSE, quietly = TRUE)
library(knitr, warn.conflicts = FALSE, quietly = TRUE)
#library(TFX)

getData<-function(tickers,datasrc){
  for (i in 1:length(tickers)){
    cat(tickers[i],i,"\n")
    getSymbols(tickers[i],src=datasrc,
               auto.assign=getOption("getSymbols.auto.assign",TRUE),
               env=parent.frame())
  }
}

makeIndex<-function(x,inv,ret){
  # Takes an xts object x and returns an index starting at 100 and evolving as the log returns of x.
  # The inv flag tells whether or not to invert the series before calculating returns.
  # The ret flag tells whether or not we have been passed a series of returns already.
  init.val<-100
  dts<-index(x,0)
  if (inv==TRUE) data<-1/x else data<-x
  if (ret==TRUE){ # we have a series of returns...
    ret.series<-x
  } else {
    ret.series<-periodReturn(data,period="daily",subset=NULL,type="arithmetic")
    dts<-index(ret.series,0)
  }
  n<-length(ret.series)
  new.series<-ret.series
  new.series[1]<-init.val
  
  for (i in 2:n){
    new.series[i]<-(1+ret.series[i-1])*new.series[i-1]
  }
  names(new.series)<-c("index")
  return(new.series)
} # My custom index funtion for converting indices to 100 based at inception.


calcWeights<-function(prices,numshares,initial){
  ret<-NULL
  for (i in 1:length(numshares)){
    sh<-numshares[i]
    ret<-cbind(ret,sh*prices[,i]/initial)
  }
  return(ret)
}


getOHLC<-function(assets,OHLC){
  # Takes a list of assets and returns either the Open, High, Low, or Close depending
  # on the passed value of HLOC. Return value is of type xts/zoo.
  ret<-NULL
  for (i in 1:length(assets)){
    if (OHLC=="O" || OHLC=="Open"){
      ret<-cbind(ret,assets[[i]][,1])
    } else {
      if (OHLC=="H" || OHLC=="High"){
        ret<-cbind(ret,assets[[i]][,2])
      } else {
        if (OHLC=="L" || OHLC=="Low"){
          ret<-cbind(ret,assets[[i]][,3])
        } else {
          if (OHLC=="C" || OHLC=="Close"){
            ret<-cbind(ret,assets[[i]][,4])
          }
        }
      }
    }
  }
  return(ret)
}


#dev.off()

# Uncomment this to redirect output plots into a PDF file
#pdf("Portfolios_performance.pdf")


# Set the initial date and investment
#first.date <- c("2015-02-01/") # The "/" indicates all dates up to current date
fdate = as.Date("2015-03-01")
ldate = as.Date("2015-04-01")
first.date = seq(fdate, ldate, by = "1 day")

# Set the initial total investment amount (for each portfolio)
init.invest <- 1000000

# TO DO: IMPLEMENT EXCHANGE RATES FOR NON-USD-DENOMINATED STOCKS
# GetExchangeRates <- function(from, to, dt=Sys.Date()) {
#   require(quantmod)
#   obj.names <- getFX(paste0(from, "/", to), from=dt, to=dt)
#   result <- numeric(length(obj.names))
#   names(result) <- obj.names
#   for (obj.name in obj.names) {
#     result[obj.name] <- as.numeric(get(obj.name))[1]
#     # Clean up    
#     rm(obj.name)
#   }
#   return(result)
# }
# 
# getQuote(paste0(from="GBP", to="USD", "=X"),from="2014-02-05",to="2014-02-08")
# getFX(paste0(from="GBP", "/", to="USD", "=X"),from="2014-02-05",to="2014-02-08")
# getSymbols("USD/GBP",src="oanda")


# Set up the data source and the tickers for all the ETFs...
#data.source = c("yahoo")
data.source = c("google")

# Define the stock to be used as a benchmark for comparison, e.g. the S&P500
tickers.bench = c("SPY")

# Set the tickers in each portfolio
denis_tickers.etf = c("TSLA","EPAM","MSFT","CSIQ","TSL","TWTR","SFUN","BABA")
denis_tickers.human = c("Tesla","EPAM","Microsoft","Canadian Solar","Trina Solar","Twitter","SouFun","Alibaba")

ilan_tickers.etf = c("FXF", "NVS", "AAPL","FB","USO","TSLA","PALL","GOOGL","TM")
ilan_tickers.human = c("SwissFrancTrust", "Novartis", "Apple","Facebook","USoilETF","Tesla","PalladiumETF","Google","Toyota")

pauly_tickers.etf = c("OTCMKTS:ETFM","SWKS","ACHN","OTCMKTS:JOES","FUEL","NBG","FLWS","MMM")
pauly_tickers.human = c("2050motors","Skyworks","Achillion","EatAtJoes","RocketFuel","NatBankGreece","1-800Flowers.com","3M")

max_tickers.etf = c("AAPL","TSLA","YHOO","SCTY")
max_tickers.human = c("Apple","Tesla","Yahoo","SolarCity")

# Get the data from Google Finance
suppressWarnings(getData(tickers.bench, data.source))

suppressWarnings(getData(denis_tickers.etf, data.source))
suppressWarnings(getData(ilan_tickers.etf, data.source))
suppressWarnings(getData(pauly_tickers.etf, data.source))
suppressWarnings(getData(max_tickers.etf, data.source))

denis_etfs <- list(TSLA,EPAM,MSFT,CSIQ,TSL,TWTR,SFUN,BABA)
ilan_etfs <- list(FXF, NVS, AAPL,FB,USO,TSLA,PALL,GOOGL,TM)
pauly_etfs <- list(`OTCMKTS:ETFM`,SWKS,ACHN,`OTCMKTS:JOES`,FUEL,NBG,FLWS,MMM)
max_etfs <- list(AAPL,TSLA,YHOO,SCTY)


# Get all the closing prices into a vector starting from the inception date...
denis_etfs_close <- getOHLC(denis_etfs, "C")
denis_etfs_close <- denis_etfs_close[first.date]
denis_etfs_close <- na.omit(denis_etfs_close)

ilan_etfs_close <- getOHLC(ilan_etfs, "C")
ilan_etfs_close <- ilan_etfs_close[first.date]
ilan_etfs_close <- na.omit(ilan_etfs_close)

pauly_etfs_close <- getOHLC(pauly_etfs, "C")
pauly_etfs_close <- pauly_etfs_close[first.date]
pauly_etfs_close <- na.omit(pauly_etfs_close)

max_etfs_close <- getOHLC(max_etfs, "C")
max_etfs_close <- max_etfs_close[first.date]
max_etfs_close <- na.omit(max_etfs_close)

# Print subsections if necessary to identify missing data, etc.
#max_etfs_close$AAPL["2014-06"]

denis_port.weights <- c(0.15,0.15,0.1,0.15,0.15,0.1,0.1,0.1)
ilan_port.weights <- c(0.1,0.1,0.15,0.1,0.1,0.1,0.15,0.1,0.1)
pauly_port.weights <- c(0.1,0.15,0.1,0.1,0.1,0.15,0.1,0.2)
max_port.weights <- c(0.25,0.25,0.25,0.25)

# Make sure that the weights in each portfolio sum to 1
#sum(denis_port.weights)
#sum(ilan_port.weights)
#sum(pauly_port.weights)
#sum(max_port.weights)

denis_invest.amounts <- init.invest*denis_port.weights
denis_shares <- denis_invest.amounts/first(denis_etfs_close)
denis_shares <- as.numeric(coredata(denis_shares))

ilan_invest.amounts <- init.invest*ilan_port.weights
ilan_shares <- ilan_invest.amounts/first(ilan_etfs_close)
ilan_shares <- as.numeric(coredata(ilan_shares))

pauly_invest.amounts <- init.invest*pauly_port.weights
pauly_shares <- pauly_invest.amounts/first(pauly_etfs_close)
pauly_shares <- as.numeric(coredata(pauly_shares))

max_invest.amounts <- init.invest*max_port.weights
max_shares <- max_invest.amounts/first(max_etfs_close)
max_shares <- as.numeric(coredata(max_shares))


# The total value of the portfolios at the end date
denis_mtm_last_vector <- last(denis_etfs_close) * denis_shares
denis_mtm_last <- sum(denis_mtm_last_vector)

ilan_mtm_last_vector <- last(ilan_etfs_close) * ilan_shares
ilan_mtm_last <- sum(ilan_mtm_last_vector)

pauly_mtm_last_vector <- last(pauly_etfs_close) * pauly_shares
pauly_mtm_last <- sum(pauly_mtm_last_vector)

max_mtm_last_vector <- last(max_etfs_close) * max_shares
max_mtm_last <- sum(max_mtm_last_vector)


# Total value of portfolios at the start date
denis_mtm_first_vector <- first(denis_etfs_close) * denis_shares
denis_mtm_first<- sum(denis_mtm_first_vector)

ilan_mtm_first_vector <- first(ilan_etfs_close) * ilan_shares
ilan_mtm_first<- sum(ilan_mtm_first_vector)

pauly_mtm_first_vector <- first(pauly_etfs_close) * pauly_shares
pauly_mtm_first<- sum(pauly_mtm_first_vector)

max_mtm_first_vector <- first(max_etfs_close) * max_shares
max_mtm_first<- sum(max_mtm_first_vector)

# Print individual entries for checks
# max_last_vector$AAPL.Close["2015-02-05"]


# Line charts of the portfolio values in USD
denis_port_value <- as.xts(apply(denis_etfs_close, 1, FUN = function(x) sum(x * denis_shares)))
plot(denis_port_value, mar=c(4.1,5.1,1.1,1.1), yaxis.right=FALSE, main="1's portfolio value") 

ilan_port_value <- as.xts(apply(ilan_etfs_close, 1, FUN = function(x) sum(x * ilan_shares)))
plot.xts(ilan_port_value, mar=c(4.1,5.1,1.1,1.1), yaxis.right=FALSE,main="2's portfolio value") 

pauly_port_value <- as.xts(apply(pauly_etfs_close, 1, FUN = function(x) sum(x * pauly_shares)))
plot.xts(pauly_port_value, mar=c(4.1,5.1,1.1,1.1), yaxis.right=FALSE, main="3's portfolio value") 

max_port_value <- as.xts(apply(max_etfs_close, 1, FUN = function(x) sum(x * max_shares)))
plot.xts(max_port_value, mar=c(4.1,5.1,1.1,1.1), yaxis.right=FALSE, main="4's portfolio value") 


#install.packages("xtsExtra", repos='http://r-forge.r-project.org')
suppressWarnings(library(xtsExtra))

# Calculate the returns and index the portfolios to 100 at the start date
bench_ret <- periodReturn(SPY[first.date][, 4], period = "daily", subset = NULL, type = "arithmetic")
bench_index <- makeIndex(bench_ret, inv = FALSE, ret = TRUE)

denis_port_ret <- periodReturn(denis_port_value, period = "daily", subset = NULL, type = "arithmetic") 
denis_port_index <- makeIndex(denis_port_ret, inv = FALSE, ret = TRUE) 
ilan_port_ret <- periodReturn(ilan_port_value, period = "daily", subset = NULL, type = "arithmetic") 
ilan_port_index <- makeIndex(ilan_port_ret, inv = FALSE, ret = TRUE) 
pauly_port_ret <- periodReturn(pauly_port_value, period = "daily", subset = NULL, type = "arithmetic")
pauly_port_index <- makeIndex(pauly_port_ret, inv = FALSE, ret = TRUE) 
max_port_ret <- periodReturn(max_port_value, period = "daily", subset = NULL, type = "arithmetic") 
max_port_index <- makeIndex(max_port_ret, inv = FALSE, ret = TRUE) 


#  Put multiple plots in a single graphic
#par(mfrow=c(2,1))


# DENIS'S PORTFOLIO
Ra <- denis_port_ret
Rb <- bench_ret
dates <- index(Rb[index(Ra,0)]) # Find the subset of dates both data sets contain, otherwise plotting is problematic
Ra <- denis_port_ret[dates]
Rb <- bench_ret[dates]
Rf <- as.numeric(last(SPY)/100/252) # Divide by 252, the annual number of trading days

# Performance relative to the benchmark:
chart.RelativePerformance(Ra, as.vector(Rb), main = "1's Relative Performance vs. S&P 500", xaxis = TRUE)

# Put the returns for all the stocks in the portfolio into a data frame 
denis_etfs.ret <- NULL 
for (i in 1:length(denis_tickers.etf)) {
  temp <- as.xts(periodReturn(denis_etfs_close[, i], period = "daily", type = "arithmetic"))
  denis_etfs.ret <- cbind(denis_etfs.ret, temp)
}

# Input the company names instead of the tickers (for more easily readable plots)
names(denis_etfs.ret) <- denis_tickers.human

#par(mfrow = c(1, 1))
# Plot the relative performance of the portfolio vs the benchmark, per asset
denis_etfs.ret <- denis_etfs.ret[dates]
chart.RelativePerformance(denis_etfs.ret, as.vector(Rb), main = "1's Relative Performance (per asset) vs. S&P 500", 
                          colorset = c(1:8), xaxis = TRUE, ylog = FALSE, legend.loc = "topleft", 
                          cex.legend = 0.8)

# Generate a table of CAPM parameters
Rf <- 0
table.CAPM(Ra, Rb, scale = 252, Rf = Rf, digits = 4)

"1's portfolio gave a total return of:"
denis_mtm_last-denis_mtm_first




# ILAN'S PORTFOLIO
Ra <- ilan_port_ret
Rb <- bench_ret
dates <- index(Rb[index(Ra,0)]) 
Ra <- ilan_port_ret[dates]
Rb <- bench_ret[dates]
Rf <- as.numeric(last(SPY)/100/252)

chart.RelativePerformance(Ra, as.vector(Rb), main = "2's Relative Performance vs. S&P 500", xaxis = TRUE)

ilan_etfs.ret <- NULL  # A data frame that holds all of the etf return streams...
for (i in 1:length(ilan_tickers.etf)) {
  temp <- as.xts(periodReturn(ilan_etfs_close[, i], period = "daily", type = "arithmetic"))
  ilan_etfs.ret <- cbind(ilan_etfs.ret, temp)
}
names(ilan_etfs.ret) <- ilan_tickers.human

#par(mfrow = c(1, 1))
ilan_etfs.ret <- ilan_etfs.ret[dates]
chart.RelativePerformance(ilan_etfs.ret, as.vector(Rb), main = "2's Relative Performance (per asset) vs. S&P 500", 
                          colorset = c(1:8), xaxis = TRUE, ylog = FALSE, legend.loc = "bottomleft", 
                          cex.legend = 0.8)

Rf <- 0
table.CAPM(Ra, Rb, scale = 252, Rf = Rf, digits = 4)

"2's Portfolio gave a total return of:"
ilan_mtm_last-ilan_mtm_first




# PAUL'S PORTFOLIO
Ra <- pauly_port_ret
Rb <- bench_ret
dates <- index(Rb[index(Ra,0)])
Ra <- pauly_port_ret[dates]
Rb <- bench_ret[dates]
Rf <- as.numeric(last(SPY)/100/252)

chart.RelativePerformance(Ra, as.vector(Rb), main = "3's: Relative Performance vs. S&P 500", xaxis = TRUE)

pauly_etfs.ret <- NULL
for (i in 1:length(pauly_tickers.etf)) {
  temp <- as.xts(periodReturn(pauly_etfs_close[, i], period = "daily", type = "arithmetic"))
  pauly_etfs.ret <- cbind(pauly_etfs.ret, temp)
}
names(pauly_etfs.ret) <- pauly_tickers.human

#par(mfrow = c(1, 1))
pauly_etfs.ret <- pauly_etfs.ret[dates]
chart.RelativePerformance(pauly_etfs.ret, as.vector(Rb), main = "3's Relative Performance (per asset) vs. S&P 500", 
                          colorset = c(1:8), xaxis = TRUE, ylog = FALSE, legend.loc = "topleft", 
                          cex.legend = 0.8)

Rf <- 0
table.CAPM(Ra, Rb, scale = 252, Rf = Rf, digits = 4)

"3's Portfolio gave a total return of:"
pauly_mtm_last-pauly_mtm_first





# MAX'S PORTFOLIO
Ra <- max_port_ret
Rb <- bench_ret
dates <- index(Rb[index(Ra,0)]) 
Ra <- max_port_ret[dates]
Rb <- bench_ret[dates]
Rf <- as.numeric(last(SPY)/100/252)

act.premium <- ActivePremium(Ra, Rb, scale = 252)
act.premium

chart.RelativePerformance(Ra, as.vector(Rb), main = "4's Relative Performance vs. S&P 500", xaxis = TRUE)

max_etfs.ret <- NULL 
for (i in 1:length(max_tickers.etf)) {
  temp <- as.xts(periodReturn(max_etfs_close[, i], period = "daily", type = "arithmetic"))
  max_etfs.ret <- cbind(max_etfs.ret, temp)
}
names(max_etfs.ret) <- max_tickers.human

#par(mfrow = c(1, 1))
max_etfs.ret <- max_etfs.ret[dates]
chart.RelativePerformance(max_etfs.ret, as.vector(Rb), main = "4's Relative Performance (per asset) vs. S&P 500", 
                          colorset = c(1:8), xaxis = TRUE, ylog = FALSE, legend.loc = "topleft", 
                          cex.legend = 0.8)

Rf <- 0
table.CAPM(Ra, Rb, scale = 252, Rf = Rf, digits = 4)

"4's portfolio gave a total return of:"
max_mtm_last-max_mtm_first





#par(mfrow=c(1,1))

#port_comp_names <- c("Denis", "Ilan", "Paul Y", "Max")
port_comp_names <- c("1", "2", "3", "4")
port_comp <- cbind(denis_port_value, ilan_port_value, pauly_port_value, max_port_value)#/init.invest

# This works without the xtsExtra library but produces much uglier plots
#plot(as.zoo(port_comp), screens=1)
#plot(as.zoo(port_comp), plot.type='single')

names(port_comp) <- port_comp_names
plot(port_comp, screens=1, main="Comparison of portfolios", mar=c(4.1,5.1,1.1,1.1), yaxis.right=FALSE, legend.loc="topleft")

# This ends the pdf() command at the beginning of the code which produces the pdf
#dev.off()
#dev.cur()