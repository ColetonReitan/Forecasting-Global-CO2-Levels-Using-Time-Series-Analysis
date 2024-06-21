# Forecasting Global CO2 Levels - A time series analysis

This was a semester long project conducted in the fall of 2021. This analysis includes data ETL, EDA & supervised machine learning modeling for time series analysis.    
Completed fully within R

## Project Overview
To forecast the rate at which atmospheric CO2 levels are increasing globally, I performed a
time series analysis on the atmospheric concentration of CO2 based on data from an observatory located
in Mauna Loa, Hawaii. The use of simple forecasting models such as Holt-Winters seasonal smoothing and
moderately advanced forecasting models such as an autoregressive integrated moving average (ARIMA)
model was used to determine the forecasted values. As a result of this analysis, there is an
understanding of not only the rate of increase of which CO2 levels are rising globally, but also where these
levels could potentially rise to in coming years.

## Methodology

### EDA
Numerous techniques of exploratory data anlysis was used to understand the data more deeply than what is seen on the surface. Of these techniques, the most important was the data decomposition technique, which broke the data down into four different categories: Noise, Trend, Seasonality, & observed. This decomposition gives the base understanding as to where the data has patterns and what techniques should be used going forward. 

### Modeling
Several models were used to gain a better understanding of the data. 
Simple forecasting models, such as a Holt-Winters, Naive, and Drift were used as benchmarking models.
More advanced forecasting models such as SARIMA and SARIMAX were the final models. 

## Data
The data is from the Global Monitoring Lab located in Mauna Loa, Hawaii, taken from the FRED. The dating of the data goes back as far as 1960 through 2021. There are roughly 750 observations across 5 features. 

### Variables
Year: The year in which the data is from  
Month: The month in which the data is from (1-12)  
Decimal.Date: The date in which the data is from in numeric form  
Average: Average Atmospheric Level of CO2 (PPM)  
Interpolated: Secondary reading of average atmospheric level of CO2 (PPM)  
