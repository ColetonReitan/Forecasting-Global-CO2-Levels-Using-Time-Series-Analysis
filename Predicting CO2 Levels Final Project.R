library(MASS)
library(tseries)
library(forecast)
library(vrtest)
library(stats)
library(rugarch)
library(FinTS)
library(ggplot2)
library(spectral)
library(zoo)
library(lubridate)
library(tidyverse)
library(car)
library(scales)
library(patchwork)
library(kableExtra)
library(matlab)
library(Metrics)
library(Latex)
library(tinytex)


### For methods, do naive, sarima, spectral. 
#Spectral for the 

# ML Original Data Plot
df <- read.csv("C:/Users/colet/Documents/Wesleyan QAC/QAC 320 Time Series Analysis/Mauna Loa data.csv")
mldf<- df$average
plot(mldf)
#ML TS Plot
mlts <- ts(mldf, start = c(1960, 1), end = c(2021, 8), frequency = 12)
tsdisplay(mlts)
plot(mlts)+
  title("Time series of Atmospheric CO2 Above Mauna Loa") + 
summary(mlts)
#Moving Average is seen within data
trend_mlts <- ma(mlts, order = 12, centre = T)
ma_line <- lines(trend_mlts, col="red")

#Decomposition of ML TS
mlts_decomp <- decompose(mlts)
plot(mlts_decomp) 
#Clearly, there is a seasonal additive trend in the data, with what appears to be constant variance. This tells me that holt-winters
#forecasting may be beneficial to start with, along with a seasonal naive model. Even as these two models seem to be the best simple
#forecasting models to start with, I am going to start by examining a drift model, to understand where about the models should be 
#forecasting. 

acf((mlts))
#Acf shows that the data is non-stationary and has non-constant mean. This makes sense because the original time series plot
#showed a continual additive upward trend. 
adf.test((mlts))
#The adf test tells me that the data is non-stationary. 

#Because the data is non-stationary, differencing the data may show itself useful. 
#Since there appears to be a strong monthly seasonal trend, I will use 12th differencing for the data. 
mlts_diff <- diff(mlts, differences = 1)
plot.ts(mlts_diff, main = "ML one diff")

mlts_diff112 <- diff(mlts_diff, differences = 1)
plot.ts(mlts_diff112, main ="1st and 12th diff")

mlts_diff113 <- diff(mlts_diff112, differences = 1)
plot.ts(mlts_diff112, main ="1st and 12th diff")

mlts_diff114 <- diff(mlts_diff113, differences = 1)
plot.ts(mlts_diff112, main ="1st and 12th diff")

mlts_diff115 <- diff(mlts_diff114, differences = 1)
plot.ts(mlts_diff112, main ="1st and 12th diff")



mlts_diff12 <- diff(mlts, differences = 12)
plot.ts(mlts_diff12, main = "ML Seasonaly Differenced")


mlts_diff112 <- diff(mlts_diff, differences = 4)
plot.ts(mlts_diff112, main ="1st and 12th diff")

#Now run adf test to check for stationary data
adf.test(mlts_diff115)
kpss.test(mlts_diff)
#p-value for the adf test is .01, data is stationary with constant mean. 
#p-value for kpss test is .1, data is stationary with constant mean.

#Check the acf to help set the parameters for sarima model. 
acf(mlts_diff115, lag = 48)
#Acf shows significant correlation until the third lag, so an MA(3) model may be appropriate. 
pacf(mlts_diff115, lag = 48)
#Pacf shows significance until the second lag, so an AR(2) model should be appropriate. 

#Now I will start to create a SARIMA model.
arima213 <- arima(mlts, order=c(2,1,3), seasonal=c(2,1,3))
summary(arima213) #AIC 394.41

arima615 <- arima(mlts, order=c(6,1,5), seasonal=c(2,1,3))
summary(arima615) #AIC 392.12

sarima615 <- arima(mlts, order=c(1,1,1), seasonal=c(6,1,5))
summary(sarima615) #AIC 391.87

arima212 <- arima(mlts, order=c(2,1,2), seasonal=c(2,1,3))
summary(arima212) #AIC 391.44

arima211 <- arima(mlts, order=c(2,1,1), seasonal=c(2,1,3))
summary(arima211) #AIC 389.83

arima113 <- arima(mlts, order=c(1,1,3), seasonal=c(1,1,2))
summary(arima113) #AIC 389.02
arima113_forecast <- forecast(arima113, h=28)
autoplot(arima113_forecast)
summary(arima113_forecast)
checkresiduals(arima113_forecast)
checkresiduals(arima113)


auot_mlts111 <- auto.arima(mlts, trace = T, approximation = F, stepwise = F) 
summary(auot_mlts111) #AIC 383.08
#auto arima gives ARIMA(1,1,1)(0,1,1) as best model
arima111_forecast <- forecast(auot_mlts111, h=28)
autoplot(arima111_forecast)
summary(arima111_forecast)
checkresiduals(auot_mlts111)




checkresiduals(arima211)
checkresiduals(auot_mlts111)
#acf plot of residuals from both arima models shown above shows that all autocorrelations are within the threshold limit
#Indicating residual behavior is similar to white noise. 

#Plotting the arima model onto the data
arima111_forecast <- forecast(auot_mlts111, h=28)
autoplot(arima111_forecast)
summary(arima111_forecast)

#Now lets look into other simple methods to see how they are forecasted. 
#Drift Method
drift_mlts <- rwf(mlts, h=100, drift=T)
plot(drift_mlts)
autoplot(drift_mlts, h=40)
ma_line <- lines(trend_mlts, col="red")
summary(drift_mlts)


#Holt-Winters Forecasting
holt_ml <- HoltWinters(mlts, seasonal = "additive")
plot(holt_ml)
plot(fitted(holt_ml))
holt_pred<- forecast(holt_ml, h=100)
plot(holt_pred)
ma_line <- lines(trend_mlts, col="red")

#holt_window <- window(holt_pred, start = 2000, end = 2025)
#holt_window <- window()

#The moving average line clearly shows a path that would go straight through the center of the forecasted values from the Holt-Winters
#Method, so the HW method seems like a reliable method. 
summary(holt_pred)


#seasonal naive method
mlts2 <- window(mlts, start=1960, end=c(2021,8))
naive_mlts <- snaive(mlts2, h=100)
autoplot(mlts2) +
  autolayer(snaive(mlts2, h=100),
            series="Seasonal Naive", PI=FALSE) +
  ggtitle("Naive forecast for monthly atmospheric CO2 levels") +
  xlab("Year") + ylab("PPM") +
  guides(colour=guide_legend(title="Forecast"))
summary(naive_mlts)

rmse(df,holt_pred)


#SARIMAX Modeling

#Import exogenous data
#First variable is income
df2 <- read.csv("C:/Users/colet/Desktop/QAC 320/income data.csv")
incomedf<- df2$PI
plot(incomedf)
#Income TS Plot
income_ts <- ts(incomedf, start = c(1960, 1), end = c(2021, 8), frequency = 12)
tsdisplay(income_ts)
plot(income_ts)+
  title("Time series of Average Income")
adf.test(income_ts) #The income ts is not stationary, must do some transformation to it

#Second Variable is Unemployment
df3 <- read.csv("C:/Users/colet/Desktop/QAC 320/Unemployment stats.csv")
unemdf <- df3$UNRATE
plot(unemdf)
#Unemployment TS Plot
unem_ts <- ts(unemdf, start = c(1976, 1), end = c(2021, 8), frequency = 12)
plot(unem_ts)+
  title("Time Series of Unemployment Rate")
adf.test(unem_ts) #Unemployment data is also not stationary. 

#Combine Variables into single ts
exog_ts <- cbind(mlts, income_ts, unem_ts)
exog_ts <- na.remove(exog_ts)
plot(exog_ts)

#Fit model with income
#exog_fit1 <- auto.arima(exog_ts[, "mlts"], xreg = exog_ts[,"income_ts"])
#checkresiduals(exog_fit1)
#summary(exog_fit1)
#Residuals seemed to be decent for this model, but the AIC stood out at 280.67

#Fit model with unemployment 
#exog_fit2 <- auto.arima(exog_ts[,"mlts"], xreg = exog_ts[,"unem_ts"])
#checkresiduals(exog_fit2)
#summary(exog_fit2)
#Residuals seemed good as well, AIC stood at 281.05, slightly worse than first exog model.

#Now to create model and plot 
# split the data into train and test sets 
exog_train <- window(exog_ts, end=2015)
exog_test <- window(exog_ts, start=2016)

# retrain model only on train data 
exog_arimaxfit1 <- auto.arima(exog_train[,"mlts"], 
                              xreg = exog_train[,"income_ts"], 
                              trace = FALSE,  
                              seasonal= TRUE, 
                              stepwise= FALSE,
                              approximation=FALSE)

exog_arimaxfit2 <- auto.arima(exog_train[,"mlts"], 
                              xreg = exog_train[,"unem_ts"], 
                              trace = FALSE,  
                              seasonal= TRUE, 
                              stepwise= FALSE,
                              approximation=FALSE)
## produce forecasts with arimax
exog_forecast1 <- forecast::forecast(exog_arimaxfit1, xreg=rep(exog_test[,"income_ts"],28))
exog_forecast2 <- forecast::forecast(exog_arimaxfit2, xreg=rep(exog_test[,"unem_ts"],28))
## plot the forecasts
autoplot(exog_forecast1) + autolayer(exog_ts[,"mlts"]) + xlim(1976, 2023) + ylim(325, 425) 
autoplot(exog_forecast2) + autolayer(exog_ts[,"mlts"]) + xlim(1976, 2023) + ylim(325, 425)

summary(exog_arimaxfit1)
#AIC for exog. var. income is 233.49
summary(exog_arimaxfit2)
#AIC for exog. var. unemployment is 232

checkresiduals(exog_forecast1)
checkresiduals(exog_arimaxfit2)
#for both residual tests, the p-value for the ljung box test is below .05. 







##############################################
#spectral analysis
#I am very unsure how to code for spectral analysis, so I will put together an ARIMAX model
spectrum(df)

mspect <- spectrum(df$average, log="no", spans=c(2,2), plot=FALSE)
delta <- 1/12
specx <- mspect$freq/delta
specy <- 2*mspect$spec
plot(specx, specy, xlab="Period (years)", ylab="Spectral Density", type="l")
#zero and one years show clear peaks


#Two models for differenced data
ssp=spectrum(mlts)
per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]

chocSM <- lm(mlts ~ sin(2*pi/per*index(mlts))+cos(2*pi/per*index(mlts))+sin(4*pi/per*index(mlts))+cos(4*pi/per*index(mlts)))
plot(chocSM)

######################################################
###Trying to plot differently
#Import data
data <- 
  read.delim('ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt', comment.char = '#', header = F, sep = '', col.names = c('Year','Month','Time','Co2_Concentration','Interpolated','Trend','Days','X'))
#Mutate the data
mauna_cc <- data %>% 
  mutate(
    Co2_Con = case_when(
      Co2_Concentration == -99.99 ~ Interpolated,
      TRUE ~ Co2_Concentration
    )
  )
#Transform data into ymd data
mauna_cc$Date <- ymd(paste0(data$Year, " ", data$Month, " ", "15"))
#Create data with only necessary columns
mauna_cc_sel <- mauna_cc %>% 
  select(Year, Month, Date, Co2_Con )
mauna_cc_sel<-subset(mauna_cc_sel, Year!="1958" & Year!="1959")
#Set training/testing dates. using 2015 on as test, everything before is train.
mauna_cc_sel_test <- mauna_cc_sel %>% 
  filter(Year > 2015)
mauna_cc_sel_train <- mauna_cc_sel %>% 
  filter(Year <= 2015)
#Now to plot the data. 
ggplot(mauna_cc_sel,aes(Date, Co2_Con)) +
  geom_line(color='red') +
  xlab("Year, Month") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "5 year") +
  theme(axis.text.x = element_text(face = "bold", color = "#993333", 
                                   size = 12, angle = 45, hjust = 1)) +
  ylab("CO2 Concentration (ppm)") +
  #scale_x_continuous(breaks = trans_breaks(identity, identity, n = 10))
  scale_y_continuous() +
  theme(axis.text.y = element_text(face = "bold", color = "#993333", 
                                   size = 10, hjust = 1),axis.title.y = element_text(size = 10))












#Look at ACF PACF of the training set
mauna_co2_train <- ts(mauna_cc_sel_train$Co2_Con, start = c(1960,1), frequency = 12)
mauna_co2_train %>% ggtsdisplay()
#Take the differnece to correct for non-stationarity
mauna_co2_train %>% diff(lag=12) %>% diff() %>% ggtsdisplay()
#Now with the differencing done, we can start setting parameters for the arima model
#Since we took a difference, d=1, D=1.
#There is only one significant spike in the seasonal lags, Q=1, while the non-differenced data (before 12 in ACF) showed 3 spikes, q=3
#First model that will be used is ARIMA(0,1,3)(3,1,1)[12]
#Run multiple arima models, see which works, and report the AICs.
mauna_aicsvalue <- function(p,q,P,Q) {
  fit <- Arima(mauna_co2_train, order=c(p,1,q), seasonal = list(order=c(P,1,Q),period=12),
               lambda = "auto"
  )
  return(fit$aicc)
}
mauna_eva <- data.frame(Model_name=c("ARIMA(0,1,3)(3,1,1)[12]","ARIMA(0,1,1)(3,1,1)[12]","ARIMA(1,1,0)(1,1,0)[12]",
                                     "ARIMA(1,1,2)(1,1,0)[12]","ARIMA(1,1,3)(0,1,1)[12]","ARIMA(1,1,1)(1,1,0)[12]",
                                     "ARIMA(1,1,1)(1,1,0)[12]","ARIMA(1,1,0)(1,1,1)[12]","ARIMA(1,1,1)(0,1,1)[12]"),
                        AIC=c(mauna_aicsvalue(0,3,3,1),mauna_aicsvalue(0,1,3,1),mauna_aicsvalue(1,0,1,0),
                               mauna_aicsvalue(1,2,1,0),mauna_aicsvalue(1,3,0,1),mauna_aicsvalue(1,1,1,0), 
                               mauna_aicsvalue(1,1,1,0),mauna_aicsvalue(1,0,1,1),mauna_aicsvalue(1,1,0,1)))
mauna_eva
#so, the arima(1,1,1)(0,1,1)[12] model had the lowest aic at -1101.2607
#my chosen parameters arima model had the third lowest aic at -1095.1323

#check the residuals
(mauna_fit_minaicc <- Arima(mauna_co2_train, order=c(1,1,1),seasonal=list(order=c(0,1,1),period=12),
                      lambda = "auto"
))
checkresiduals(mauna_fit_minaicc, lag=36)
mauna_fit_minaicc$aicc
#Test test data
mauna_co2_test <- ts(mauna_cc_sel_test$Co2_Con, start = c(2015,1), frequency = 12)
mauna_acc <- accuracy(forecast(mauna_fit_minaicc,h=35)$mean, mauna_co2_test )


#Now let's plot the forecasted model!
mauna_co2_train %>%
  Arima(order=c(1,1,1),seasonal=list(order=c(0,1,1),period=12),
        lambda = "auto"
  ) %>%
  forecast(h=400) %>%
  autoplot() +
  ylab("CO2 Concentration (PPM)") + xlab("Year")
#  autolayer(mauna_co2_test)
#Zoom in to see visual performance. What's shown is the testing vs. predicted values. 
prediction <- forecast(mauna_fit_minaicc,h=71) 
mauna_cc_sel_test$prediction <- prediction$mean
mauna_test_pre_tidy <- gather(mauna_cc_sel_test, "type", "Co2", -Year,-Month,-Date)
ggplot(mauna_test_pre_tidy,aes(Date, Co2,color=type)) +
  geom_line() +
  xlab("Year, Month") +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "1 year") +
  theme(axis.text.x = element_text(face = "bold", color = "#993333", 
                                   size = 12, angle = 45, hjust = 1)) +
  ylab("CO2 Concentration (ppm)") +
  #scale_x_continuous(breaks = trans_breaks(identity, identity, n = 10))
  scale_y_continuous() +
  theme(axis.text.y = element_text(face = "bold", color = "#993333", 
                                   size = 10, hjust = 1), axis.title.y = element_text(size = 8))




rmse_eva <- function(p,d,q,P,D,Q) {
  fit <- Arima(mauna_co2_train, order=c(p,d,q),seasonal=list(order=c(P,D,Q),period=12),
               lambda = "auto", include.drift = T
  )
  mm <- accuracy(forecast(fit,h=35)$mean, mauna_co2_test)
  return(mm[2])
}
rmse_eva <- data.frame(Model_name=c(
  "ARIMA(0,1,3)(3,1,1)[12]",
  "ARIMA(1,1,3)(0,1,1)[12]",
  "ARIMA(1,1,1)(0,1,1)[12]"
), RMSE=c(                        
  rmse_eva(0,1,3,3,1,1),
  rmse_eva(1,1,3,0,1,1),
  rmse_eva(1,1,1,0,1,1)))








###############################################################################################
#Remove trend from ts, seasonality exposed. 
detrend_mlts = mlts - trend_mlts
plot(as.ts(detrend_mlts))

#Averaging seasonality
m_mlts = t(matrix(data = detrend_mlts, nrow = 12))
seasonal_mlts = colMeans(m_mlts, na.rm = T)
plot(as.ts(rep(seasonal_mlts,16)))

#Examine Randomness
random_mlts = mlts - trend_mlts - seasonal_mlts
plot(as.ts(random_mlts))

#Reconstruct Original
recomposed_mlts = trend_mlts+seasonal_mlts+random_mlts
plot(as.ts(recomposed_mlts))

#Histogram Plot, not very necessary
hist(mlts, main = "Histogram of GSPC", freq = FALSE, col = "yellow") 

#decomposition plot
decomp <- decompose(mlts)
plot(decomp)


#Differencing Plot to check variance
diff1 <- diff(mlts, lag = 12, ylab = "change in Atmospheric CO2") #see variance not constant
plot(diff1)
title("First difference of Atmospheric CO2")


#step 1. test for normality
shapiro.test(mlts[0:5000]) #using first 5000 points to get rid of error
# p <0.05 therefore not normal, difference using log
df1 <- log(mlts)
plot(df1)
shapiro.test(diff1[0:5000])

#stationarity, mean, and variance
#check for variance
Auto.VR(df1)
#variance is not constant given how but $stat value is 
#check for stationarity
adf.test(df1) #data is not stationary
df2 <- diff(mlts)
df2 <- na.remove(df2)
#do another test
adf.test(df2) #it is significant for stationary

Auto.VR(df2) #negative value so constant

acf(mlts)
acf(df2)
acf(df2^2) 
# Seasonal pattern, sarima model would likely be preferred.


#AR and MA terms
#manual
acf.df2 <- acf(df2, main = "ACF of Atmoshperic CO2", lag.max = 50) #3 lags for MA -> q=3
pacf.df2 <- pacf(df2, main = "PACF of Atmospheric CO2", lag.max = 50) #1 lag for AR -> p=1
arima103 <- arima(df2, order = c(1,0,3), seasonal = c(0,2,1)) #ARIMA(1,0,3)(0,2,1) 
arima103 #AIC is 810.88

#automatic
auto.arima(df2, trace = T, approximation = F, stepwise = F) 
auto.arima(df2) #ARIMA(1,0,2)(0,1,1) 376.97

ar1<-arima(df2,c(1,0,2), seasonal = c(0,1,1))
ar1 # AIC is 388.56

#Let's go with the auto-ARIMA (1,0,2)(0,1,1) in df2 (after differencing, residuals were not great before differencing). AIC was significantly better (for autoarima model)

#(S)ARIMA model 
arima102 <- arima(df2, order = c(1,0,2), seasonal = c(0,1,1))
summary(arima102)
checkresiduals(arima102)
#Residuals show some upside; only one line on ACF that is significant, residual graph looks smooth.

arimar <- arima102$residuals
ggtsdisplay(arimar) #note, ACF and PACF only have one significant spike, both at the 3rd spot. 
