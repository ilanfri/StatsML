setwd("/home/ilanfri/Desktop/Hoy_no_circula")

# http://www.aire.df.gob.mx/default.php?opc=%27aKBhnmI=%27&opcion=aw==
#download.file("http://148.243.232.112:8080/opendata/IndiceCalidadAire_gz/indice_2015.csv.gz", "indice_2015.csv.gz")
#download.file("http://148.243.232.112:8080/opendata/IndiceCalidadAire_gz/indice_2016.csv.gz", "indice_2016.csv.gz")
# Manually unzip as these don't work:
#unzip("indice_2015.csv.gz")
#ind15 = read.table(gzfile("indice_2015.csv.gz"))

# Manually delete non-tabular data at head of files

ind15 = read.csv("indice_2015.csv")
ind16 = read.csv("indice_2016.csv")
#colnames(ind15)
#colnames(ind16)

ind15$Fecha <- as.Date(as.character(ind15$Fecha),format="%d/%m/%Y")
ind16$Fecha <- as.Date(as.character(ind16$Fecha),format="%d/%m/%Y")

# Check percentage of rows with missing data
sum(is.na(ind15))/nrow(ind15)
sum(is.na(ind16))/nrow(ind16)

# See examples of missing data rows
#ind15[!complete.cases(ind15), ]
#ind16[!complete.cases(ind16), ]

#install.packages("data.table")
library(data.table)

# Make a table showing percentages of missing values per column
nmisvaltabrows = length(colnames(ind15))
colmisvals = data.table(colnum=integer(nmisvaltabrows), pcmissing=double(nmisvaltabrows))
for(i in seq(length(colnames(ind15)))) {
  colmisvals[i, "colnum"] <- i
  colmisvals[i, "pcmissing"] <- sum(is.na(ind15[,i]))/nrow(ind15)
  }
colmisvals[order(-pcmissing)]

sum(is.na(ind15[,"Centro.mon.xido.de.carbono"]))/nrow(ind15)
length(ind15[,"Centro.mon.xido.de.carbono"])
length(ind15[!is.na(ind15[,"Centro.mon.xido.de.carbono"]),"Centro.mon.xido.de.carbono"])

# Should in principle do the above missing value checks for ind16 as well, but can't be asked

# Check that column names are the same in both tables so they can be merged
#all.equal(colnames(ind15), colnames(ind16))
alldata = merge(ind15,ind16, all=TRUE)
alldata = alldata[!duplicated(alldata),]
alldata = alldata[order(alldata[,"Fecha"], alldata[,"Hora"]), ]
# Check that ordering worked as expected:
#alldata[1:30, c("Fecha", "Hora")]

data = alldata[!is.na(alldata[,"Centro.mon.xido.de.carbono"]), c("Fecha","Hora","Centro.mon.xido.de.carbono")]

library(dplyr)

# CausalImpact doesn't like multiple entries with the same timestamp, summarise by taking mean over hours to get single value per day
data_mean = data %>% group_by(Fecha) %>% summarise_each(funs(mean(., na.rm=TRUE)), -Hora)

# To use full data set from start of 2015:
#before_period = c(head(ozonedata_mean,1)$Fecha, as.Date("2015-07-08"))
#after_period = c(as.Date("2015-07-09"), tail(ozonedata_mean,1)$Fecha)

# To use symmetric time period before and after intervention, using all data up to latest
intervention_date = as.Date("2015-07-09")
intrvtn_to_latest_dt = tail(data_mean,1)$Fecha - intervention_date
intrvtn_minus_dt = intervention_date - intrvtn_to_latest_dt
print(sprintf("Considering symmetric period of 2*%s days between %s and %s (intervention: %s)", intrvtn_to_latest_dt, intrvtn_minus_dt, tail(data_mean,1)$Fecha, intervention_date))

before_period = c(seq(tail(data_mean,1)$Fecha, length=2, by="-1 years")[2], as.Date("2015-07-08"))
after_period = c(as.Date("2015-07-09"), tail(data_mean,1)$Fecha)

# Convect data to time series
library(zoo)
tsdata = read.zoo(data.table(data_mean))
#tsdata

library(xts)
plot(as.xts(tsdata), xlab="Date", ylab="CO concentration (unknown units)", major.format = "%Y-%m")
#title(main = "", adj = 0, outer = TRUE, line = -1.5)

library(CausalImpact)
impact = CausalImpact(tsdata, before_period, after_period)

summary(impact)
summary(impact, "report")
