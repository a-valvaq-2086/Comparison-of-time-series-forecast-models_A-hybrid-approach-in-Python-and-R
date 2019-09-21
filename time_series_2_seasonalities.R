setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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



# =========================================
# Methods 1: Simple methods - Mean, Naive, Seasonal Naive
# =========================================
pjme_mean    <- meanf(train_ts, h=24*245)
pjme_naive   <- naive(train_ts, h=24*245)
pjme_drift   <- rwf(train_ts, h=24*245, drift=TRUE)
pjme_s_naive <- snaive(train_ts, h=24*245)

# Plot: Mean method  
autoplot(train_ts) +
  autolayer(test_ts, series="Test") +
  autolayer(pjme_mean, series="Mean", PI=FALSE) +
  ggtitle("Power Demand (MW) over time [2016 - 2018]") +
  xlab("Year") + ylab("Demand in MW") +
  guides(colour=guide_legend(title="Forecast"))

# Plot: Naive method
autoplot(train_ts) +
  autolayer(test_ts, series="Test") +
  autolayer(pjme_naive, series="Naïve", PI=FALSE) +
  ggtitle("Power Demand (MW) over time [2016 - 2018]") +
  xlab("Year") + ylab("Demand in MW") +
  guides(colour=guide_legend(title="Forecast"))

# Plot: Naive with drift method
autoplot(train_ts) +
  autolayer(test_ts, series="Test") +
  autolayer(pjme_drift, series="Drift", PI=FALSE) +
  ggtitle("Power Demand (MW) over time [2016 - 2018]") +
  xlab("Year") + ylab("Demand in MW") +
  guides(colour=guide_legend(title="Forecast"))

# Plot: Naive seasonal method
autoplot(train_ts) +
  autolayer(test_ts, series="Test") +
  autolayer(pjme_s_naive, series="S. naïve", PI=FALSE, alpha=0.7) +
  ggtitle("Power Demand (MW) over time [2016 - 2018]") +
  xlab("Year") + ylab("Demand in MW") +
  guides(colour=guide_legend(title="Forecast"))

# Plot: Naive seasonal zoomed in
autoplot(train_ts) +
  autolayer(test_ts, series="Test") +
  autolayer(pjme_s_naive, series="S. naïve", PI=FALSE, alpha=0.7) +
  ggtitle("Power Demand (MW) zoomed [2017 - 2018]") +
  xlab("Year") + ylab("Demand in MW") +
  guides(colour=guide_legend(title="Forecast")) +
  coord_cartesian(xlim = c(2017.94, 2018.6))

# =========================================
# Method 2.1: STL with two seasonal periods - STL - NOT WORKING on short T.S.
# =========================================
train_ts %>% mstl() %>% 
  autoplot(include=945*24, PI=FALSE) +
  autolayer(test_ts, series="Test", alpha=0.8) +
  ylab("Demand in MW") + xlab("Date")

# =========================================
# Method 2: STL with two seasonal periods - STL + ETS - NOT WORKING on short T.S.
# =========================================
train_ts %>% stlf(h = 24*245, level=80) %>% +
  autoplot(include=945*24, PI=FALSE) + xlab("Year") +
  autolayer(test_ts, series="Test")

# ===============================================================
# Method 3: Dynamic Harmonic regression. With daily and yearly seasonalities
# ===============================================================
fit <- auto.arima(train_ts, seasonal=TRUE, lambda = 0, xreg = fourier(train_ts, K=c(12,2)))

autoplot(train_ts) + 
  autolayer(test_ts, series="Test") +
  fit %>% 
  forecast(xreg=fourier(train_ts, K=c(12,2,4), h=24*245), level = 99) %>% 
  autolayer(series="Dyn. Reg. + Fourier", PI=FALSE, alpha=0.7) +
  ggtitle("Power Demand (MW) over time [2016 - 2018]") +
  xlab("Year") + ylab("Demand in MW") +
  guides(colour=guide_legend(title="Forecast"))

# Zoomed plot
  autoplot(train_ts) + 
  autolayer(test_ts, series="Test") +
  fit %>% 
  forecast(xreg=fourier(train_ts, K=c(12,2,4), h=24*245), level = 99) %>% 
  autolayer(series="Dyn. Reg. + Fourier", PI=FALSE, alpha=0.7) +
    ggtitle("Power Demand (MW) zoomed [2017 - 2018]") +
    xlab("Year") + ylab("Demand in MW") +
    guides(colour=guide_legend(title="Forecast")) +
    coord_cartesian(xlim = c(2017.94, 2018.6))

checkresiduals(fit)

# ===========================
# Method 4: TBATS forecasting
# ===========================
train_ts %>% 
  tbats(use.arma.errors = TRUE,
        use.parallel = TRUE) -> fit2

fc2 <- forecast(fit2, h=24*245, level = 80)
autoplot(train_ts) +
  autolayer(test_ts, series="Test") +
  autolayer(fc2, series="TBATS", PI=FALSE, alpha=0.7) +
  ggtitle("Power Demand (MW) over time [2016 - 2018]") +
  xlab("Year") + ylab("Demand in MW") +
  guides(colour=guide_legend(title="Forecast"))

# Zoomed plot
fc2 <- forecast(fit2, h=24*245, level = 80)
autoplot(train_ts) +
  autolayer(test_ts, series="Test") +
  autolayer(fc2, series="TBATS", PI=FALSE, alpha=0.7) +
  ggtitle("Power Demand (MW) zoomed [2017 - 2018]") +
  xlab("Year") + ylab("Demand in MW") +
  guides(colour=guide_legend(title="Forecast")) +
  coord_cartesian(xlim = c(2017.94, 2018.6))

checkresiduals(fit2)

# =============================
# Accuracy Predicted vs. Test
# =============================
acc1 <- round(accuracy(pjme_mean, test_ts),2)[2,5]
acc2 <- round(accuracy(pjme_naive, test_ts),2)[2,5]
acc3 <- round(accuracy(pjme_drift, test_ts),2)[2,5]
acc4 <- round(accuracy(pjme_s_naive, test_ts),2)[2,5]

sprintf("Mean method MAPE: %s %%", acc1)
sprintf("Naive method MAPE: %s %%", acc2)
sprintf("Naive method with drift MAPE: %s %%", acc3)
sprintf("Seasonal Naive method MAPE: %s %%", acc4)

#For ARMA, ARIMA and TBATS no second argument needed
accuracy(fit)
accuracy(fit2)

sprintf("Dynamic Reg. MAPE: %s %%", round(accuracy(fit)[5],2))
sprintf("TBATS MAPE: %s %%", round(accuracy(fit2)[5],2))
