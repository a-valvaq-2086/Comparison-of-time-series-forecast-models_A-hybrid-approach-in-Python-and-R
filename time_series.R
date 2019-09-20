library(tidyverse)
library(ggfortify)
library(lubridate)
library(xts)
library(forecast)
library(zoo)
library(tseries)
library(urca)

# Importing the data, and parsing to data
pjme = read.csv("data/PJME_hourly.csv")
pjme[1] <- lapply(pjme[1], as.character)
pjme[,1]  <-  parse_date_time(pjme[,1], "ymd HMS")
pjme <- distinct(pjme,Datetime, .keep_all = TRUE)
pjme <- arrange(pjme,Datetime)
pjme <- pjme[order(pjme$Datetime),]

pjme_mini <- pjme %>% 
  filter(Datetime > ymd_hms("2016-01-01 00:00:00"))

train <- pjme %>% 
  filter(Datetime < ymd_hms("2017-12-01 00:00:00")) %>% 
  filter(Datetime > ymd_hms("2015-12-31 23:00:00"))
test  <- pjme %>% 
  filter(Datetime > ymd_hms("2017-11-30 23:00:00"))

firstHour    <- 24*(as.Date("2017-11-30 23:00:00")-as.Date("2017-01-01 00:00:00"))
pjme_ts      <- msts(pjme$PJME_MW, start = c(2002,0), ts.frequency = 24*365.25, seasonal.periods = c(24,24*365.25))
train_ts     <- msts(train$PJME_MW, start = c(2016,0), ts.frequency = 24*365.25, seasonal.periods = c(24,24*365.25))
test_ts      <- msts(test$PJME_MW, start = c(2017,firstHour), ts.frequency = 24*365.25, seasonal.periods = c(24,24*365.25))

autoplot(train_ts) +
  autolayer(train_ts, series = "Train", color="Black") +
  autolayer(test_ts, series  = "Test") + 
  ggtitle("Test") + xlab("Date time") + ylab("Demand in MW")

ggseasonplot(pjme_mini_ts, year.labels = TRUE, year.labels.left = TRUE) +
  ylab("Demand in MW") +
  ggtitle("Seasonal plot: Power Demand")

# ============================================
# STL decomposition for multiple seasonality
# ============================================
pjme_ts %>% mstl() %>% 
  autoplot() + xlab("Year")
# Note the seasonal component is half of the magnitud of the daily and annual components

# ====================
# Stationarity check
# ====================
train_ts %>% ur.kpss() %>% summary() # Non-differenced data

train_ts %>% diff() %>% ur.kpss() %>% summary()

ndiffs(train_ts) # We should take a first difference for the data

train_ts %>% diff(lag=24) %>% ndiffs() # No seasonal difference.

train_ts %>% log() %>% nsdiffs()

# =========================================
# Method 2.1: STL with two seasonal periods - STL
# =========================================
train_ts %>% mstl() %>% 
  autoplot(include=945*24, PI=FALSE) +
  autolayer(test_ts, series="Test", alpha=0.8) +
  ylab("Demand in MW") + xlab("Date")

# =========================================
# Method 2: STL with two seasonal periods - STL + ETS
# =========================================
train_ts %>% stlf(h = 24*245, level=80) %>% +
  autoplot(include=945*24, PI=FALSE) + xlab("Year") +
  autolayer(test_ts, series="Test")

# ===============================================================
# Method 3: Dynamic Harmonic regression. With daily and yearly seasonalities
# ===============================================================
fit <- auto.arima(train_ts, seasonal=TRUE, D=0, lambda = 0, xreg = fourier(train_ts, K=c(12,2)))

fit %>% 
  forecast(xreg=fourier(train_ts, K=c(12,2), h=24*245), level = 80) %>% 
  autoplot(include=945*24, PI=FALSE) +
  autolayer(test_ts, series="Test", alpha=0.8) +
  ylab("Demand in MW") + xlab("Date")

checkresiduals(fit)

# =====================
# TBATS forecasting
# =====================
train_ts %>% 
  tbats(use.arma.errors = TRUE,
        use.parallel = TRUE) -> fit2

fc2 <- forecast(fit2, h=24*245, level = 80)
autoplot(fc2, include=945*24, PI=FALSE) +
  autolayer(test_ts, series="Test") +
  ylab("Demand in MW") + xlab("Date")

checkresiduals(fit2)

# ====================
# Predicted vs. test
# ====================
