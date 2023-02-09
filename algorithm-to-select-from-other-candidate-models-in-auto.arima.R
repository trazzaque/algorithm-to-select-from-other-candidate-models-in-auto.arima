library(hts)
library(forecast)
library(miscTools)
# Read Data
DistrictGenericNameData=read.csv('G:/Porashuna/MS Thesis/District Generic data for time series final.csv')

# Convert to time series object
DistrictGenericNameDataTimeSeries=ts(DistrictGenericNameData,start = c(2013,12), end = c(2015,6),frequency = 12)

# Convert to hierarchical series
DistrictGenericNameDataHts=hts(DistrictGenericNameDataTimeSeries,nodes=list(10,rep(64,10)))

# Change names of node 1
DistrictGenericNameDataHts$`labels`$`Level 1`[1]="AZITHROMYCIN"
DistrictGenericNameDataHts$`labels`$`Level 1`[2]="CALCIUM.CARBONATE.VIT.D"
DistrictGenericNameDataHts$`labels`$`Level 1`[3]="CEFIXIME"
DistrictGenericNameDataHts$`labels`$`Level 1`[4]="DOMPERIDONE"
DistrictGenericNameDataHts$`labels`$`Level 1`[5]="ESOMEPRAZOLE"
DistrictGenericNameDataHts$`labels`$`Level 1`[6]="MELITRACEN.FLUPENTHIXOL"
DistrictGenericNameDataHts$`labels`$`Level 1`[7]="OMEPRAZOLE"
DistrictGenericNameDataHts$`labels`$`Level 1`[8]="PANTOPRAZOLE"
DistrictGenericNameDataHts$`labels`$`Level 1`[9]="PARACETAMOL"
DistrictGenericNameDataHts$`labels`$`Level 1`[10]="VITAMIN.B1.B6.B12"

# Seperate all of the series from the three nodes
allSeries=aggts(DistrictGenericNameDataHts)


##DistrictGenericNameDataTraining = window(DistrictGenericNameDataHts, start = c(2013,12), end = c(2015,5))
##DistrictGenericNameDataTest = window(DistrictGenericNameDataHts, start = c(2015,6), end = c(2015,6))
##DistrictGenericNameDataForecast <- forecast(DistrictGenericNameDataTraining, h=1, method="bu")
##DistrictGenericNameDataAccuracy=accuracy.gts(DistrictGenericNameDataForecast, DistrictGenericNameDataTest)

#chittagongCefixime = allSeries[,'CHITTAGONG.CEFIXIME']
#chittagongCefiximeTraining = window(chittagongCefixime,start=c(2013,12),end=c(2015,5))
#chittagongCefiximeTest = window(chittagongCefixime,start=c(2015,6),end=c(2015,6))
#chittagongCefiximeArima = auto.arima(chittagongCefiximeTraining) 
#chittagongCefiximeForecast = forecast(chittagongCefiximeArima,h=1)
#chittagongCefiximeAccuracy = accuracy(chittagongCefiximeForecast,chittagongCefiximeTest)


# Model Selection Algorithm

allSeriesNames=colnames(as.data.frame(allSeries))                          # Store series names in a variable
ChosenModels <- list()                                                     # Declare a blank list for storing all models
for(j in 1:ncol(allSeries)){                                                
  tempModel <- auto.arima(allSeries[,j])                                   # apply auto.arima on a series
if (sum(arimaorder(tempModel))==0){                                        # check if it returns any parameter or not
  trace <- capture.output({                                                # If doesn't return any parameters, store all other possible models
    model <- auto.arima(allSeries[,j], trace = TRUE)
  })
  con    <- textConnection(trace)
  models <- read.table(con, sep=":")                                      # Cointains all other possible models
  close(con)
  models$V1=as.character(models$V1)                                       # Contains model details in string format
  
  modelsRmse=data.frame(P=as.numeric(),D=as.numeric(),Q=as.numeric(),RMSE=as.numeric())           # Declare a blank data frame to hold parameter and MSE values
  for (i in 1:(length(models$V1)-1)){
    param=as.numeric(unlist(regmatches(models$V1[i], gregexpr("[[:digit:]]+", models$V1[i]))))    # Extract model details from the string format
    forecastFunction <- function(x,h){forecast(Arima(x, order=param,method = 'ML') ,h=h)}         # Declare a funtion to forecast ARIMA model
    CrossValidate=tsCV(allSeries[,j],forecastFunction,h=1)                                        # Cross validate for the ARIMA model and its forecast
    tempModeldetails=c(param,sqrt(mean(CrossValidate^2, na.rm=TRUE)))                             # Store model details and RMSE in a variable
    modelsRmse[i,]=tempModeldetails                                                               # Insert a row with model details and RMSE to that data frame declared before
  }
  bestModelOrder=as.numeric(unique(modelsRmse[modelsRmse$RMSE==min(modelsRmse$RMSE),1:3]))        # Select the best ARIMA model with lowest RMSE after cross validation
  expForecastFunction=function(x, h){forecast(ets(allSeries[,j]), h=h)}                           # Declare a funtion to forecast Exponential model
  expCrossValidate=tsCV(allSeries[,j],expForecastFunction,h=1)                                    # Cross validate for the Exponential model and its forecast
  expMse=sqrt(mean(expCrossValidate^2, na.rm=TRUE))                                               # Store RMSE of the Exponential Model
  if (expMse>min(modelsRmse$RMSE)){                                                               # If RMSE of Exponential is bigger than the ARIMA RMSE,
  ChosenModels[[j]]=Arima(allSeries[,j],order = bestModelOrder,method='ML')                       # Run the selected ARIMA model
  names(ChosenModels)[j]=paste(allSeriesNames[j])}                                                # Store it in the declared list
  else{
    ChosenModels[[j]]=ets(allSeries[,j])                                                          # If Exponential model has the least RMSE,
    names(ChosenModels)[j]=paste(allSeriesNames[j])}                                              # Run Exponential model and store in the list
}
  else{
    ChosenModels[[j]]=tempModel                                                                   # If auto.arima return at least one parameter, store it
    names(ChosenModels)[j]=paste(allSeriesNames[j])}                                              # in the list declare above 
  print(j)
  }


ChosenModelsForecastFuntion = function(x){forecast(x,h=1)}
ChosenModelsForecastBoxTestFunction=function(x){Box.test(x$residuals,type = "Ljung-Box")}
ChosenModelsForecast = lapply(ChosenModels,ChosenModelsForecastFuntion)
ChosenModelsForecastBoxTest=lapply(ChosenModelsForecast,ChosenModelsForecastBoxTestFunction)

BoxTestPValues=set.seed(651)
for (i in 12:length(ChosenModelsForecastBoxTest)){
  BoxTestPValues[i]=ChosenModelsForecastBoxTest[[i]]$p.value}
write.csv(BoxTestPValues,'G:/Porashuna/MS Thesis/BoxTestPValues.csv')


ChosenModelsForecastAcfDetect=list()
for (i in 1:length(ChosenModelsForecast)){
  AcfValues=Acf(ChosenModelsForecast[[i]]$residuals,plot=F)
  UnlistedAcfValues=unlist(AcfValues[[1]])
  ChosenModelsForecastAcfDetect[[i]]=which(!(UnlistedAcfValues<=(1.96/sqrt(19)) & UnlistedAcfValues>=-(1.96/sqrt(19))))
  names(ChosenModelsForecastAcfDetect)[i]=paste(allSeriesNames[i])
}

ChosenModelsForecastAcfDetectOrders=NULL
for(i in 12:length(ChosenModelsForecastAcfDetect)){
  tempModelsForecastAcfDetectOrders=ChosenModelsForecastAcfDetect[[i]]
  if(length(tempModelsForecastAcfDetectOrders)==1){
    ChosenModelsForecastAcfDetectOrders=append(ChosenModelsForecastAcfDetectOrders,0)
  }
  else{
  ChosenModelsForecastAcfDetectOrders=append(ChosenModelsForecastAcfDetectOrders,tempModelsForecastAcfDetectOrders[2])}
}

ChosenModelsForecastPacfDetect=list()
for (i in 1:length(ChosenModelsForecast)){
  PacfValues=Pacf(ChosenModelsForecast[[i]]$residuals,plot=F)
  UnlistedPacfValues=unlist(PacfValues[[1]])
  ChosenModelsForecastPacfDetect[[i]]=which(!(UnlistedPacfValues<=(1.96/sqrt(19)) & UnlistedPacfValues>=-(1.96/sqrt(19))))
  names(ChosenModelsForecastPacfDetect)[i]=paste(allSeriesNames[i])
}

ChosenModelsForecastPacfDetectOrders=NULL
for(i in 12:length(ChosenModelsForecastPacfDetect)){
  tempModelsForecastPacfDetectOrders=ChosenModelsForecastPacfDetect[[i]]
  if(length(tempModelsForecastPacfDetectOrders)==0)
  {
    ChosenModelsForecastPacfDetectOrders=append(ChosenModelsForecastPacfDetectOrders,0)  
  }
  else{
  ChosenModelsForecastPacfDetectOrders=append(ChosenModelsForecastPacfDetectOrders,tempModelsForecastPacfDetectOrders)}
}


ChosenModelsMape = NULL
for (i in 1:length(ChosenModelsForecast)){
  ChosenModelsMape[i]=accuracy(ChosenModelsForecast[[i]])[5]
}
write.csv(ChosenModelsMape,'G:/Porashuna/MS Thesis/ChosenModelsMape.csv')


defaultArimaModels=list()
for(j in 1:ncol(allSeries)){                                                
  tempModel <- auto.arima(allSeries[,j])
  defaultArimaModels[[j]]=tempModel}

defaultArimaModelsForecast = lapply(defaultArimaModels,ChosenModelsForecastFuntion)

defaultArimaModelsMape=NULL
for (i in 1:length(defaultArimaModelsForecast)){
  defaultArimaModelsMape[i]=accuracy(defaultArimaModelsForecast[[i]])[5]
}
  

forecastData=data.frame(Forecast=as.numeric(),Lo80=as.numeric(),Hi80=as.numeric(),Lo95=as.numeric(),
                        Hi95=as.numeric())
for (i in 1:length(allSeriesNames)){
  tempDataFrame=as.data.frame(ChosenModelsForecast[[i]])
  forecastData[i,]=tempDataFrame
}
forecastData[,'SeriesName']=allSeriesNames
#write.csv(forecastData,'G:/Porashuna/MS Thesis/Forecast Data.csv')

defaultArimaModelsOrder=NULL
for(i in 1:length(allSeriesNames)){
  tempModelOrder=arimaorder(defaultArimaModels[[i]])
  defaultArimaModelsOrder[i]=paste(as.character(tempModelOrder)[1],as.character(tempModelOrder)[2],as.character(tempModelOrder)[3],sep=',')
}


ChosenModelsOrder=NULL
for(i in 12:length(allSeriesNames)){
  tempModelOrder=ChosenModels[[i]]
  if(class(tempModelOrder)=="ets"){
    ChosenModelsOrder=append(ChosenModelsOrder,tempModelOrder$method)
  }
  else{
  tempModelOrder=arimaorder(tempModelOrder)
  ChosenModelsOrder=append(ChosenModelsOrder,(paste('ARIMA','(',as.character(tempModelOrder)[1],',',as.character(tempModelOrder)[2],',',as.character(tempModelOrder)[3],')')))}
}
#write.csv(ChosenModelsOrder,'G:/Porashuna/MS Thesis/ChosenModelsOrder.csv')



  