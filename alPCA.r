install.packages("nnfor")
install.packages("tsintermittent")

library(tsintermittent); library(nnfor)
library(Metrics); library(forecast); library(tseries);library(dplyr)
library(tidyverse); library(prophet)
library(forecTheta);library(tsfgrnn) 
library(yager);library(univariateML);
library(thief)

## Data: The average monthly data of CO2 mole fraction (co2-mm-mlo), which represents CO2 emissions,
# have been taken as an illustrative example. The series is obtained from daily averages and covers 
# from 01/01/2000 to 31/12/2019 (20 years and 240 observations). The data are from the U.S. 
#  Government's Earth System Research Laboratory, Global Monitoring Division, and are available 
# from the Trends in Atmospheric Carbon Dioxide website \url{https://datahub.io/core/co2-ppm#data}
x <-ts(c(377.04,375.52,380.69,379.21,377.12,378.38,379.44,376.07,377.81,380.13,376.79,378.04,
         377.29,378.17,378.92,378.14,382.36,377.56,377.15,375.64,379.56,380.49,374.41,373.46,
         373.28,378.14,373.93,376.37,371.84,376.39,374.60,372.07,386.03,379.18,374.71,370.40,
         367.73,369.39,370.86,371.96,368.79,369.70,372.56,391.80,373.76,376.22,368.45,366.91,
         366.21,362.38,365.18,367.40,367.93,369.41,371.67,372.61,373.95,374.46,375.49,367.54,
         368.85,371.64,369.79,372.24,378.29,373.13,369.34,363.30,365.49,360.64,368.58,369.97,
         371.85,363.79,361.74,362.77,370.70,370.94,365.05,363.13,364.46,376.14,367.91,364.12,
         369.22,366.28,364.47,365.04,364.03,361.05,369.17,373.33,368.98,362.31,375.79,364.99,
         367.59,365.94,368.99,377.53,375.39,366.16,361.43,362.31,368.38,355.80,357.82,359.44,
         362.91,363.81,361.53,362.08,363.65,366.78,366.29,369.21,363.68,371.12,374.31,365.75,
         361.91,367.21,367.00,367.72,366.40,362.70,359.46,358.72,357.61,362.07,368.30,369.43,
         364.64,372.31,370.78,377.23,374.36,367.08,368.31,366.20,365.76,374.41,360.65,361.80,
         359.16,370.93,372.75,366.94,360.86,372.01,367.12,365.39,369.99,367.26,364.88,365.76,
         365.68,369.29,366.17,365.95,368.49,367.96,367.21,369.93,362.48,367.20,370.66,369.98,
         366.37,379.25,377.22,375.00,369.47,366.44,362.60,378.62,374.53,375.82,376.83,374.03,
         366.10,373.90,373.61,374.91,377.30,374.89,372.10,368.31,365.65,367.99,378.41,382.06,
         374.89,379.22,375.20,375.87,384.68,376.52,370.63,380.63,379.66,374.57,380.81,383.32,
         382.35,378.56,374.25,376.10,374.73,372.95,374.84,383.67,381.89,377.06,373.13,373.89,
         373.78,374.40,379.05,382.79,374.48,374.39,376.83,376.50,386.90,382.62,386.48,388.31,
         384.71,376.27,383.93,389.19,395.03,391.24,382.23,382.43,379.89,382.62,386.48,388.31),
       frequency = 12, start = c(2000, 1))
datos=list(x=x)

## Functions
## Functions used in alPCA
# calculo1(): Sets each configuration of a model in period T1, and generates the predictions for period T2 for each configuration.
# calculo2():Calculates the loss functions in T1 and T2
# calculo3():Constructs the error matrix
# calculo4.m():Calculates the CPs and Manhattan distance to the method associated with the smallest error
# calculo5():Calculates the weights associated to each combination according to the set coordinates
# calculo6():Obtain the predictions for the period T1+T2


calculo1<-function(j){
  # We define the general parameters
  serie<-ensayo[j]
  n_p<-12  
  n<-length(datos[[serie]]) 
  estac<-12 #monthly seasonality
  comienzo <- start(datos[[serie]])
  
  # We select the period T1(x) T2(x)
  x  <- ts(datos[[serie]][1:(n-(2*n_p))],start = comienzo, frequency = estac) #T1
  if(end(x)[2]==estac) 
  {comienzo2<-c(end(x)[1]+1,1)
  }else 
  {comienzo2<-c(end(x)[1],end(x)[2]+1)}
  x1 <- ts(datos[[serie]][(n-(2*n_p)+1):(n-n_p)], start = comienzo2,frequency = estac) #T2
  
  set.seed(846) 
  metodos<-list()
  #Naïve
  metodos[[1]]<-naive(x,h=n_p)
  #Naïve with seasonality
  metodos[[2]]<-snaive(x,h=n_p)
  #Exponential smoothing
  metodos[[3]]<-ets(x,model = "ZZZ", opt.crit="lik", damped = NULL)
  metodos[[3]]<-predict(metodos[[3]],n_p)
  metodos[[4]]<-ets(x,model = "ZZZ", opt.crit="mse", damped = NULL)
  metodos[[4]]<-predict(metodos[[4]],n_p)
  metodos[[5]]<-ets(x,model = "ZZZ", opt.crit="amse", damped = NULL)
  metodos[[5]]<-predict(metodos[[5]],n_p)
  metodos[[6]]<-ets(x,model = "ZZZ", opt.crit="sigma", damped = NULL)
  metodos[[6]]<-predict(metodos[[6]],n_p)
  metodos[[7]]<-ets(x,model = "ZZZ", opt.crit="mae", damped = NULL)
  metodos[[7]]<-predict(metodos[[7]],n_p)
  metodos[[8]]<-ets(x,model = "ZZZ", opt.crit="lik", damped = TRUE)
  metodos[[8]]<-predict(metodos[[8]],n_p)
  metodos[[9]]<-ets(x,model = "ZZZ", opt.crit="mse", damped = TRUE)
  metodos[[9]]<-predict(metodos[[9]],n_p)
  metodos[[10]]<-ets(x,model = "ZZZ", opt.crit="amse", damped = TRUE)
  metodos[[10]]<-predict(metodos[[10]],n_p)
  metodos[[11]]<-ets(x,model = "ZZZ", opt.crit="sigma", damped = TRUE)
  metodos[[11]]<-predict(metodos[[11]],n_p)
  metodos[[12]]<-ets(x,model = "ZZZ", opt.crit="mae", damped = TRUE)
  metodos[[12]]<-predict(metodos[[12]],n_p)
  #ARIMA
  a<-NULL
  metodos[[13]]<-auto.arima(x, ic="aic",approximation=T, trace=FALSE, allowdrift=F)
  a<-predict(metodos[[13]],n_p)
  metodos[[13]]$mean<-a$pred
  metodos[[14]]<-auto.arima(x, ic="bic",approximation=T, trace=FALSE, allowdrift=F)
  a<-predict(metodos[[14]],n_p)
  metodos[[14]]$mean<-a$pred
  #Theta
  metodos[[15]]<-otm.arxiv(x, h=n_p, s= NULL, g= "sAPE")
  metodos[[16]]<-dotm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="Nelder-Mead") #Dynamic Optimised Theta Model
  metodos[[17]]<-dotm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="L-BFGS-B") #Dynamic Optimised Theta Model
  metodos[[18]]<-dotm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="SANN") #Dynamic Optimised Theta Model
  metodos[[19]]<-dstm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="Nelder-Mead") #Dynamic Standard Theta Model
  metodos[[20]]<-dstm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="L-BFGS-B") #Dynamic Standard Theta Model
  metodos[[21]]<-dstm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="SANN") #Dynamic Standard Theta Model
  metodos[[22]]<-otm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="Nelder-Mead")  #Optimised Theta Model
  metodos[[23]]<-otm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="L-BFGS-B")  #Optimised Theta Model
  metodos[[24]]<-otm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="SANN")  #Optimised Theta Model
  metodos[[25]]<-stm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="Nelder-Mead")  #Standard Theta Model (Fiorucci)
  metodos[[26]]<-stm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="L-BFGS-B")  #Standard Theta Model (Fiorucci)
  metodos[[27]]<-stm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="SANN")  #Standard Theta Model (Fiorucci)
  metodos[[28]]<-stheta(x, h=n_p, s= NULL) #Standard Theta Model (Assimakopoulos)
  #STL (sólo series estacionales)
  metodos[[29]]<-stlf(x,h=n_p, method=c("ets"),robust = TRUE,allow.multiplicative.trend = TRUE)
  metodos[[30]]<-stlf(x,h=n_p, method=c("arima"),robust = TRUE,allow.multiplicative.trend = TRUE)
  metodos[[31]]<-stlf(x,h=n_p, forecastfunction=thetaf,robust = TRUE,allow.multiplicative.trend = TRUE)
  #Croston
  metodos[[32]]<-crost(x,h=n_p,type = "croston",cost = "mar",init.opt=TRUE)
  metodos[[32]]$fitted<-ts(metodos[[32]]$frc.in, start =comienzo,frequency = estac) 
  metodos[[32]]$mean<-  ts(metodos[[32]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[33]]<-crost(x,h=n_p,type = "croston",cost = "msr",init.opt=TRUE)
  metodos[[33]]$fitted<-ts(metodos[[33]]$frc.in, start =comienzo,frequency = estac) 
  metodos[[33]]$mean<-  ts(metodos[[33]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[34]]<-crost(x,h=n_p,type = "croston",cost = "mae",init.opt=TRUE)
  metodos[[34]]$fitted<-ts(metodos[[34]]$frc.in, start =comienzo,frequency = estac)
  metodos[[34]]$mean<-  ts(metodos[[34]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[35]]<-crost(x,h=n_p,type = "croston",cost = "mse",init.opt=TRUE)
  metodos[[35]]$fitted<-ts(metodos[[35]]$frc.in, start =comienzo,frequency = estac)
  metodos[[35]]$mean<-  ts(metodos[[35]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[36]]<-crost(x,h=n_p,type = "sba",cost = "mar",init.opt=TRUE)
  metodos[[36]]$fitted<-ts(metodos[[36]]$frc.in, start =comienzo,frequency = estac)
  metodos[[36]]$mean<-  ts(metodos[[36]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[37]]<-crost(x,h=n_p,type = "sba",cost = "msr",init.opt=TRUE)
  metodos[[37]]$fitted<-ts(metodos[[37]]$frc.in, start =comienzo,frequency = estac)
  metodos[[37]]$mean<-  ts(metodos[[37]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[38]]<-crost(x,h=n_p,type = "sba",cost = "mae",init.opt=TRUE)
  metodos[[38]]$fitted<-ts(metodos[[38]]$frc.in, start =comienzo,frequency = estac)
  metodos[[38]]$mean<-  ts(metodos[[38]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[39]]<-crost(x,h=n_p,type = "sba",cost = "mse",init.opt=TRUE)
  metodos[[39]]$fitted<-ts(metodos[[39]]$frc.in, start =comienzo,frequency = estac)
  metodos[[39]]$mean<-  ts(metodos[[39]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[40]]<-crost(x,h=n_p,type = "sbj",cost = "mar",init.opt=TRUE)
  metodos[[40]]$fitted<-ts(metodos[[40]]$frc.in, start =comienzo,frequency = estac)
  metodos[[40]]$mean<-  ts(metodos[[40]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[41]]<-crost(x,h=n_p,type = "sbj",cost = "msr",init.opt=TRUE)
  metodos[[41]]$fitted<-ts(metodos[[41]]$frc.in, start =comienzo,frequency = estac)
  metodos[[41]]$mean<-  ts(metodos[[41]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[42]]<-crost(x,h=n_p,type = "sbj",cost = "mae",init.opt=TRUE)
  metodos[[42]]$fitted<-ts(metodos[[42]]$frc.in, start =comienzo,frequency = estac)
  metodos[[42]]$mean<-  ts(metodos[[42]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[43]]<-crost(x,h=n_p,type = "sbj",cost = "mse",init.opt=TRUE)
  metodos[[43]]$fitted<-ts(metodos[[43]]$frc.in, start =comienzo,frequency = estac)
  metodos[[43]]$mean<-  ts(metodos[[43]]$frc.out, start =comienzo2,frequency = estac) 
  #Prophet
  set.seed(50)
  df<-data.frame(ds = seq(as.Date(paste(start(x)[1],"-",start(x)[2],"-","01",sep="")),
                          as.Date(paste(end(x)[1],"-",end(x)[2],"-","01",sep="")),by="m"), y = x)
  metodos[[44]]<-prophet(df,yearly.seasonality=TRUE,
                         weekly.seasonality = FALSE,
                         daily.seasonality = FALSE,seasonality.mode = "additive")
  a<-make_future_dataframe(metodos[[44]],periods = 18,freq = "month")
  b<-predict(metodos[[44]],a)
  metodos[[44]]$fitted <- ts(b$yhat[1:length(x)],start = comienzo,frequency = estac)
  metodos[[44]]$mean <- ts(b$yhat[(length(x)+1):(length(x)+n_p)],start = comienzo2,frequency = estac)
  metodos[[45]]<-prophet(df,yearly.seasonality=TRUE,
                         weekly.seasonality = FALSE,
                         daily.seasonality = FALSE,seasonality.mode = "multiplicative")
  a<-make_future_dataframe(metodos[[45]],periods = 18,freq = "month")
  b<-predict(metodos[[45]],a)
  metodos[[45]]$fitted <- ts(b$yhat[1:length(x)],start = comienzo,frequency = estac)
  metodos[[45]]$mean <- ts(b$yhat[(length(x)+1):(length(x)+n_p)],start = comienzo2,frequency = estac)
  set.seed(50)
  metodos[[46]]<-nnetar(x,p=estac,scale.inputs=T)
  a<-forecast(metodos[[46]],h=n_p)
  metodos[[46]]$mean<-a$mean  
  #TBATS model (Exponential smoothing state space model with Box-Cox transformation, ARMA errors, Trend and Seasonal components)
  set.seed(50)
  metodos[[47]]<-tbats(x, biasadj = TRUE)
  metodos[[47]]$fitted<-metodos[[47]]$fitted.values
  a<-forecast(metodos[[47]],h=n_p)
  metodos[[47]]$mean<-a$mean
  #GRNN: general regression neural network 
  set.seed(50)
  metodos[[48]]<-grnn_forecasting(x, h = n_p, transform = "additive")
  a<-rolling_origin(metodos[[48]],h = length(x)-estac-1)
  metodos[[48]]$fitted<-a$predictions[,1]
  metodos[[48]]$mean<-metodos[[48]]$prediction
  metodos[[49]]<-grnn_forecasting(x, h = n_p, transform = "multiplicative")
  a<-rolling_origin(metodos[[49]],h = length(x)-estac-1)
  metodos[[49]]$fitted<-a$predictions[,1]
  metodos[[49]]$mean<-metodos[[49]]$prediction
  #MLP: multilayer perceptron
  set.seed(50)
  metodos[[50]]<-mlp.thief(x, h=12)
  metodos[[51]]<-mlp(x, m=12, comb = c("median"),hd = c(10,5))
  a<-forecast(metodos[[51]],h=12)
  metodos[[51]]$mean<-a$mean
  metodos[[52]]<-mlp(x, m=12, comb = c("mean"),hd = c(10,5))
  a<-forecast(metodos[[52]],h=12)
  metodos[[52]]$mean<-a$mean
  gc()
  return(metodos)
}

calculo2<- function(j) {
  # We define the general parameters
  serie<-ensayo[j]
  n_p<-12 
  n<-length(datos[[serie]]) #nº de datos
  estac<-12 
  comienzo <- start(datos[[serie]])
  
  # We select the period T1(x) T2(x)
  x  <- ts(datos[[serie]][1:(n-(2*n_p))],start = comienzo, frequency = estac) #T1
  if(end(x)[2]==estac) 
  {comienzo2<-c(end(x)[1]+1,1)
  }else 
  {comienzo2<-c(end(x)[1],end(x)[2]+1)}
  x1 <- ts(datos[[serie]][(n-(2*n_p)+1):(n-n_p)], start = comienzo2,frequency = estac) #T2
  
  #sMAPE in T1: In some cases it is straightforward, but in others we lose initial values.
  metodos.T1[[j]][[1]]$sMapeT1 <- 100*Metrics::smape(x[2:length(x)],metodos.T1[[j]][[1]]$fitted[2:length(x)])
  metodos.T1[[j]][[2]]$sMapeT1 <- 100*Metrics::smape(x[(estac+1):length(x)],metodos.T1[[j]][[2]]$fitted[(estac+1):length(x)])
  for (i in 1:24) {
    metodos.T1[[j]][[2+i]]$sMapeT1 <- 100*Metrics::smape(x,metodos.T1[[j]][[2+i]]$fitted)
  }
  for (i in 1:17) {
    metodos.T1[[j]][[26+i]]$sMapeT1 <- 100*Metrics::smape(x[2:length(x)],metodos.T1[[j]][[26+i]]$fitted[2:length(x)])
  }
  for (i in 1:2) {
    metodos.T1[[j]][[43+i]]$sMapeT1 <- 100*Metrics::smape(x,metodos.T1[[j]][[43+i]]$fitted)
  }
  metodos.T1[[j]][[46]]$sMapeT1 <- 100*Metrics::smape(x[(estac+1):length(x)],metodos.T1[[j]][[46]]$fitted[(estac+1):length(x)])
  metodos.T1[[j]][[47]]$sMapeT1 <- 100*Metrics::smape(x,metodos.T1[[j]][[47]]$fitted)
  for (i in 1:2) {
    metodos.T1[[j]][[47+i]]$sMapeT1 <- 100*Metrics::smape(x[1:(length(x)-estac-1)],metodos.T1[[j]][[47+i]]$fitted)
  }
  metodos.T1[[j]][[50]]$sMapeT1 <- 100*Metrics::smape(x[(estac+3):length(x)],
                                                      metodos.T1[[j]][[50]]$fitted[(estac+3):length(x)])
  metodos.T1[[j]][[51]]$sMapeT1 <- 100*Metrics::smape(x,metodos.T1[[j]][[51]]$fitted)
  metodos.T1[[j]][[52]]$sMapeT1 <- 100*Metrics::smape(x,metodos.T1[[j]][[52]]$fitted)
  #MASE in T1: In some cases it is straightforward, but in others we lose initial values.
  metodos.T1[[j]][[1]]$maseT1 <- Metrics::mase(x[2:length(x)],metodos.T1[[j]][[1]]$fitted[2:length(x)],step_size=12)
  metodos.T1[[j]][[2]]$maseT1 <- Metrics::mase(x[(estac+1):length(x)],metodos.T1[[j]][[2]]$fitted[(estac+1):length(x)],step_size=12)
  for (i in 1:24) {
    metodos.T1[[j]][[2+i]]$maseT1 <- Metrics::mase(x,metodos.T1[[j]][[2+i]]$fitted,step_size=12)
  }
  for (i in 1:17) {
    metodos.T1[[j]][[26+i]]$maseT1 <- Metrics::mase(x[2:length(x)],metodos.T1[[j]][[26+i]]$fitted[2:length(x)],step_size=12)
  }
  for (i in 1:2) {
    metodos.T1[[j]][[43+i]]$maseT1 <- Metrics::mase(x,metodos.T1[[j]][[43+i]]$fitted,step_size=12)
  }
  metodos.T1[[j]][[46]]$maseT1 <- Metrics::mase(x[(estac+1):length(x)],metodos.T1[[j]][[46]]$fitted[(estac+1):length(x)],step_size=12)
  metodos.T1[[j]][[47]]$maseT1 <- Metrics::mase(x,metodos.T1[[j]][[47]]$fitted,step_size=12)
  for (i in 1:2) {
    metodos.T1[[j]][[47+i]]$maseT1 <- Metrics::mase(x[1:(length(x)-estac-1)],metodos.T1[[j]][[47+i]]$fitted,step_size=12)
  }
  metodos.T1[[j]][[50]]$maseT1 <- Metrics::mase(x[(estac+3):length(x)],metodos.T1[[j]][[50]]$fitted[(estac+3):length(x)])
  metodos.T1[[j]][[51]]$maseT1 <- Metrics::mase(x, metodos.T1[[j]][[51]]$fitted)
  metodos.T1[[j]][[52]]$maseT1 <- Metrics::mase(x, metodos.T1[[j]][[52]]$fitted)
  #sMAPE in T2
  for (i in 1:length(metodos.T1[[1]])) {
    metodos.T1[[j]][[i]]$sMapeT2 <- 100*Metrics::smape(x1,metodos.T1[[j]][[i]]$mean)
  }
  #RMSE in T2
  for (i in 1:length(metodos.T1[[1]])) {
    metodos.T1[[j]][[i]]$rmseT2 <- Metrics::rmse(x1,metodos.T1[[j]][[i]]$mean)
  }
  #mase en T2
  for (i in 1:length(metodos.T1[[1]])) {
    metodos.T1[[j]][[i]]$maseT2 <- Metrics::mase(x1,metodos.T1[[j]][[i]]$mean,step_size=1)
  }
  #OWA in T2
  for (i in 1:length(metodos.T1[[1]])) {
    metodos.T1[[j]][[i]]$owaT2 <- ((metodos.T1[[j]][[i]]$maseT2 / metodos.T1[[j]][[2]]$maseT2)+(metodos.T1[[j]][[i]]$sMapeT2 / metodos.T1[[j]][[2]]$sMapeT2))/2
  }
  return(metodos.T1[[j]])
}  

calculo3<-function(j){
  erroresT2[[j]]<-matrix(NA,nrow = length(metodos.T1[[1]]),ncol = 13)
  colnames(erroresT2[[j]])<-c("id","smapeT1","maseT1","rmse","smape",
                              "mase","owa","CP1","CP2", "dist", "w1","w2","w3")
  for (i in 1:length(metodos.T1[[1]])) {
    erroresT2[[j]][i,1]<-i
    erroresT2[[j]][i,2]<-metodos.T1[[j]][[i]]$sMapeT1
    erroresT2[[j]][i,3]<-metodos.T1[[j]][[i]]$maseT1
    erroresT2[[j]][i,4]<-metodos.T1[[j]][[i]]$rmseT2
    erroresT2[[j]][i,5]<-metodos.T1[[j]][[i]]$sMapeT2
    erroresT2[[j]][i,6]<-metodos.T1[[j]][[i]]$maseT2
    erroresT2[[j]][i,7]<-metodos.T1[[j]][[i]]$owaT2
  }
  dat1<-as.data.frame(round(erroresT2[[j]],3))
  dat1<-dat1 %>% distinct(smapeT1, .keep_all = TRUE)
  dat1<-dat1 %>% distinct(rmse, .keep_all = TRUE)
  erroresT2[[j]]<-as.matrix(dat1)
  return(erroresT2[[j]])
}

calculo4<-function(j){
  pca[[j]] <- prcomp(erroresT2[[j]][,2:7], scale = TRUE)
  if(sum(pca[[j]][["rotation"]][,1])<0) {erroresT2[[j]][,8]<- round(-(pca[[j]][["x"]][,1]),3)} else {erroresT2[[j]][,8]<- round(pca[[j]][["x"]][,1],3)}
  erroresT2[[j]][,9]<- round(pca[[j]][["x"]][,2],3)
  for (i in 1:dim(erroresT2[[j]])[1]) {
    erroresT2[[j]][i,10]<- round(abs(erroresT2[[j]][i,8] - min(erroresT2[[j]][,8])),3)
  }
  resultado[[1]][[j]]<-erroresT2[[j]]
  resultado[[2]][[j]]<-pca[[j]]
  return(resultado)
}

calculo4.m<-function(j){
  pca[[j]] <- prcomp(erroresT2[[j]][,2:7], scale = TRUE)
  a<-summary(pca[[j]])
  erroresT2[[j]][,8:9]<-round(pca[[j]][["x"]][,1:2],3)
  if(a[["importance"]][2,1]>0.80) {
    if(sum(pca[[j]][["rotation"]][3:6,1])<0) {for (i in 1:nrow(erroresT2[[j]])) {erroresT2[[j]][i,10]<- round(abs(erroresT2[[j]][i,8] - max(erroresT2[[j]][,8])),3)} } 
    else { for (i in 1:nrow(erroresT2[[j]])) {erroresT2[[j]][i,10]<- round(abs(erroresT2[[j]][i,8] - min(erroresT2[[j]][,8])),3)} }
  } else {
    if(sum(pca[[j]][["rotation"]][3:6,1])>0 & sum(pca[[j]][["rotation"]][1:2,2])>0) {
      m1 <- min(erroresT2[[j]][,8])
      m2 <- min(erroresT2[[j]][,9])} else if(sum(pca[[j]][["rotation"]][3:6,1])>0 & sum(pca[[j]][["rotation"]][1:2,2])<0) {
        m1 <- min(erroresT2[[j]][,8])
        m2 <- max(erroresT2[[j]][,9])} else if(sum(pca[[j]][["rotation"]][3:6,1])<0 & sum(pca[[j]][["rotation"]][1:2,2])>0) {
          m1 <- max(erroresT2[[j]][,8])
          m2 <- min(erroresT2[[j]][,9])} else {
            m1 <- max(erroresT2[[j]][,8])
            m2 <- max(erroresT2[[j]][,9])}
    #Manhattan distance
    for (i in 1:nrow(erroresT2[[j]])) {
      erroresT2[[j]][i,10]<- round(abs(erroresT2[[j]][i,8] - m1)+abs(erroresT2[[j]][i,9] - m2),3)
    }
  }
  resultado[[1]][[j]]<-erroresT2[[j]]
  resultado[[2]][[j]]<-pca[[j]]
  #If they exist, distances 0 must be replaced by distances 0.001
  for (i in 1:nrow(resultado[[1]][[j]])) {
    resultado[[1]][[j]][i,10]<-ifelse(resultado[[1]][[j]][i,10]==0.000, 0.001,resultado[[1]][[j]][i,10])
  }
  comparacion_bic <- BIC(
    mlexp(resultado[[1]][[j]][,10]),
    mlinvgamma(resultado[[1]][[j]][,10]),
    mlgamma(resultado[[1]][[j]][,10])
  )
  resultado[[3]][[j]]<- comparacion_bic %>% rownames_to_column(var = "distribucion") %>% arrange(BIC)
  return(resultado)
}

calculo5<-function(j,i){
  if(dis.ajust[j] ==1) {obj<-mlexp(resultado[[1]][[j]][,10])} else
    if(dis.ajust[j]==2) {obj<-mlgamma(resultado[[1]][[j]][,10])} else {obj<-mlinvgamma(resultado[[1]][[j]][,10])}
  dat<-subset(resultado[[1]][[j]],resultado[[1]][[j]][,10]<=round(qml(cota[i], obj = obj),4) )
  if(nrow(dat)<1) {dat<-subset(resultado[[1]][[j]], resultado[[1]][[j]][,1]==1)}
  a<-NULL
  b<-NULL
  for (m in 1:nrow(dat)) {
    a[m]<-round(1/dat[m,2],4) 
    b[m]<-round(1/dat[m,5],4) 
  }
  for (m in 1:nrow(dat)) {
    dat[m,11]<- round(abs(dat[m,8])/sum(abs(dat[,8])),4)  #w1
    dat[m,12]<- round(a[m] / sum(a),4)  #w2
    dat[m,13]<- round(b[m] / sum(b),4)  #w3
  }
  resultado[[4]][[j]][[i]]<-dat
  return(resultado)
}

calculo6<-function(j){
  # We define the general parameters
  serie<-ensayo[j]
  n_p<-12 
  n<-length(datos[[serie]]) 
  estac<-12 
  comienzo <- start(datos[[serie]])
  # Select the period T1(x) + T2(x)
  x  <- ts(datos[[serie]][1:(n-n_p)],start = comienzo, frequency = estac) #T1
  if(end(x)[2]==estac) 
  {comienzo2<-c(end(x)[1]+1,1)
  }else 
  {comienzo2<-c(end(x)[1],end(x)[2]+1)}
  x2 <- ts(datos[[serie]][(n-n_p+1):n],start = comienzo2,frequency = estac) #T2
  set.seed(846) 
  metodos<-list()
  #Naïve
  metodos[[1]]<-naive(x,h=n_p)
  #Naïve with seasonality
  metodos[[2]]<-snaive(x,h=n_p)
  #Exponential smoothing
  metodos[[3]]<-ets(x,model = "ZZZ", opt.crit="lik", damped = NULL)
  metodos[[3]]<-predict(metodos[[3]],n_p)
  metodos[[4]]<-ets(x,model = "ZZZ", opt.crit="mse", damped = NULL)
  metodos[[4]]<-predict(metodos[[4]],n_p)
  metodos[[5]]<-ets(x,model = "ZZZ", opt.crit="amse", damped = NULL)
  metodos[[5]]<-predict(metodos[[5]],n_p)
  metodos[[6]]<-ets(x,model = "ZZZ", opt.crit="sigma", damped = NULL)
  metodos[[6]]<-predict(metodos[[6]],n_p)
  metodos[[7]]<-ets(x,model = "ZZZ", opt.crit="mae", damped = NULL)
  metodos[[7]]<-predict(metodos[[7]],n_p)
  metodos[[8]]<-ets(x,model = "ZZZ", opt.crit="lik", damped = TRUE)
  metodos[[8]]<-predict(metodos[[8]],n_p)
  metodos[[9]]<-ets(x,model = "ZZZ", opt.crit="mse", damped = TRUE)
  metodos[[9]]<-predict(metodos[[9]],n_p)
  metodos[[10]]<-ets(x,model = "ZZZ", opt.crit="amse", damped = TRUE)
  metodos[[10]]<-predict(metodos[[10]],n_p)
  metodos[[11]]<-ets(x,model = "ZZZ", opt.crit="sigma", damped = TRUE)
  metodos[[11]]<-predict(metodos[[11]],n_p)
  metodos[[12]]<-ets(x,model = "ZZZ", opt.crit="mae", damped = TRUE)
  metodos[[12]]<-predict(metodos[[12]],n_p)
  #ARIMA
  a<-NULL
  metodos[[13]]<-auto.arima(x, ic="aic",approximation=T, trace=FALSE, allowdrift=F)
  a<-predict(metodos[[13]],n_p)
  metodos[[13]]$mean<-a$pred
  metodos[[14]]<-auto.arima(x, ic="bic",approximation=T, trace=FALSE, allowdrift=F)
  a<-predict(metodos[[14]],n_p)
  metodos[[14]]$mean<-a$pred
  #Theta
  metodos[[15]]<-otm.arxiv(x, h=n_p, s= NULL, g= "sAPE")
  metodos[[16]]<-dotm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="Nelder-Mead") #Dynamic Optimised Theta Model
  metodos[[17]]<-dotm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="L-BFGS-B") #Dynamic Optimised Theta Model
  metodos[[18]]<-dotm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="SANN") #Dynamic Optimised Theta Model
  metodos[[19]]<-dstm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="Nelder-Mead") #Dynamic Standard Theta Model
  metodos[[20]]<-dstm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="L-BFGS-B") #Dynamic Standard Theta Model
  metodos[[21]]<-dstm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="SANN") #Dynamic Standard Theta Model
  metodos[[22]]<-otm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="Nelder-Mead")  #Optimised Theta Model
  metodos[[23]]<-otm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="L-BFGS-B")  #Optimised Theta Model
  metodos[[24]]<-otm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="SANN")  #Optimised Theta Model
  metodos[[25]]<-stm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="Nelder-Mead")  #Standard Theta Model (Fiorucci)
  metodos[[26]]<-stm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="L-BFGS-B")  #Standard Theta Model (Fiorucci)
  metodos[[27]]<-stm(x, h=n_p, s= NULL, estimation=TRUE, opt.method="SANN")  #Standard Theta Model (Fiorucci)
  metodos[[28]]<-stheta(x, h=n_p, s= NULL) #Standard Theta Model (Assimakopoulos)
  #STL (seasonal series only)
  metodos[[29]]<-stlf(x,h=n_p, method=c("ets"),robust = TRUE,allow.multiplicative.trend = TRUE)
  metodos[[30]]<-stlf(x,h=n_p, method=c("arima"),robust = TRUE,allow.multiplicative.trend = TRUE)
  metodos[[31]]<-stlf(x,h=n_p, forecastfunction=thetaf,robust = TRUE,allow.multiplicative.trend = TRUE)
  #Croston
  metodos[[32]]<-crost(x,h=n_p,type = "croston",cost = "mar",init.opt=TRUE)
  metodos[[32]]$fitted<-ts(metodos[[32]]$frc.in, start =comienzo,frequency = estac) 
  metodos[[32]]$mean<-  ts(metodos[[32]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[33]]<-crost(x,h=n_p,type = "croston",cost = "msr",init.opt=TRUE)
  metodos[[33]]$fitted<-ts(metodos[[33]]$frc.in, start =comienzo,frequency = estac) 
  metodos[[33]]$mean<-  ts(metodos[[33]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[34]]<-crost(x,h=n_p,type = "croston",cost = "mae",init.opt=TRUE)
  metodos[[34]]$fitted<-ts(metodos[[34]]$frc.in, start =comienzo,frequency = estac)
  metodos[[34]]$mean<-  ts(metodos[[34]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[35]]<-crost(x,h=n_p,type = "croston",cost = "mse",init.opt=TRUE)
  metodos[[35]]$fitted<-ts(metodos[[35]]$frc.in, start =comienzo,frequency = estac)
  metodos[[35]]$mean<-  ts(metodos[[35]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[36]]<-crost(x,h=n_p,type = "sba",cost = "mar",init.opt=TRUE)
  metodos[[36]]$fitted<-ts(metodos[[36]]$frc.in, start =comienzo,frequency = estac)
  metodos[[36]]$mean<-  ts(metodos[[36]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[37]]<-crost(x,h=n_p,type = "sba",cost = "msr",init.opt=TRUE)
  metodos[[37]]$fitted<-ts(metodos[[37]]$frc.in, start =comienzo,frequency = estac)
  metodos[[37]]$mean<-  ts(metodos[[37]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[38]]<-crost(x,h=n_p,type = "sba",cost = "mae",init.opt=TRUE)
  metodos[[38]]$fitted<-ts(metodos[[38]]$frc.in, start =comienzo,frequency = estac)
  metodos[[38]]$mean<-  ts(metodos[[38]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[39]]<-crost(x,h=n_p,type = "sba",cost = "mse",init.opt=TRUE)
  metodos[[39]]$fitted<-ts(metodos[[39]]$frc.in, start =comienzo,frequency = estac)
  metodos[[39]]$mean<-  ts(metodos[[39]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[40]]<-crost(x,h=n_p,type = "sbj",cost = "mar",init.opt=TRUE)
  metodos[[40]]$fitted<-ts(metodos[[40]]$frc.in, start =comienzo,frequency = estac)
  metodos[[40]]$mean<-  ts(metodos[[40]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[41]]<-crost(x,h=n_p,type = "sbj",cost = "msr",init.opt=TRUE)
  metodos[[41]]$fitted<-ts(metodos[[41]]$frc.in, start =comienzo,frequency = estac)
  metodos[[41]]$mean<-  ts(metodos[[41]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[42]]<-crost(x,h=n_p,type = "sbj",cost = "mae",init.opt=TRUE)
  metodos[[42]]$fitted<-ts(metodos[[42]]$frc.in, start =comienzo,frequency = estac)
  metodos[[42]]$mean<-  ts(metodos[[42]]$frc.out, start =comienzo2,frequency = estac) 
  metodos[[43]]<-crost(x,h=n_p,type = "sbj",cost = "mse",init.opt=TRUE)
  metodos[[43]]$fitted<-ts(metodos[[43]]$frc.in, start =comienzo,frequency = estac)
  metodos[[43]]$mean<-  ts(metodos[[43]]$frc.out, start =comienzo2,frequency = estac) 
  #Prophet
  set.seed(50)
  df<-data.frame(ds = seq(as.Date(paste(start(x)[1],"-",start(x)[2],"-","01",sep="")),
                          as.Date(paste(end(x)[1],"-",end(x)[2],"-","01",sep="")),by="m"), y = x)
  metodos[[44]]<-prophet(df,yearly.seasonality=TRUE,
                         weekly.seasonality = FALSE,
                         daily.seasonality = FALSE,seasonality.mode = "additive")
  a<-make_future_dataframe(metodos[[44]],periods = 18,freq = "month")
  b<-predict(metodos[[44]],a)
  metodos[[44]]$fitted <- ts(b$yhat[1:length(x)],start = comienzo,frequency = estac)
  metodos[[44]]$mean <- ts(b$yhat[(length(x)+1):(length(x)+n_p)],start = comienzo2,frequency = estac)
  metodos[[45]]<-prophet(df,yearly.seasonality=TRUE,
                         weekly.seasonality = FALSE,
                         daily.seasonality = FALSE,seasonality.mode = "multiplicative")
  a<-make_future_dataframe(metodos[[45]],periods = 18,freq = "month")
  b<-predict(metodos[[45]],a)
  metodos[[45]]$fitted <- ts(b$yhat[1:length(x)],start = comienzo,frequency = estac)
  metodos[[45]]$mean <- ts(b$yhat[(length(x)+1):(length(x)+n_p)],start = comienzo2,frequency = estac)
  #NNAR: neural network autoregression
  set.seed(50)
  metodos[[46]]<-nnetar(x,p=estac,scale.inputs=T)
  a<-forecast(metodos[[46]],h=n_p)
  metodos[[46]]$mean<-a$mean
  #TBATS model (Exponential smoothing state space model with Box-Cox transformation, ARMA errors, Trend and Seasonal components)
  set.seed(50)
  metodos[[47]]<-tbats(x, biasadj = TRUE)
  metodos[[47]]$fitted<-metodos[[47]]$fitted.values
  a<-forecast(metodos[[47]],h=n_p)
  metodos[[47]]$mean<-a$mean
  #GRNN: general regression neural network 
  set.seed(50)
  metodos[[48]]<-grnn_forecasting(x, h = n_p, transform = "additive")
  a<-rolling_origin(metodos[[48]],h = length(x)-estac-1)
  metodos[[48]]$fitted<-a$predictions[,1]
  metodos[[48]]$mean<-metodos[[48]]$prediction
  metodos[[49]]<-grnn_forecasting(x, h = n_p, transform = "multiplicative")
  a<-rolling_origin(metodos[[49]],h = length(x)-estac-1)
  metodos[[49]]$fitted<-a$predictions[,1]
  metodos[[49]]$mean<-metodos[[49]]$prediction
  #MLP: multilayer perceptron
  # En fitted del 14:n, las primeras n_p+2 son NA
  set.seed(50)
  metodos[[50]]<-mlp.thief(x, h=12)
  # En fitted del 5:n, las primeras 4 son vacías
  metodos[[51]]<-mlp(x, m=12, comb = c("median"),hd = c(10,5))
  a<-forecast(metodos[[51]],h=12)
  metodos[[51]]$mean<-a$mean
  metodos[[52]]<-mlp(x, m=12, comb = c("mean"),hd = c(10,5))
  a<-forecast(metodos[[52]],h=12)
  metodos[[52]]$mean<-a$mean
  
  return(metodos)
}


## The algorithm is designed to obtain the predictions of a monthly series
ensayo<-1:length(datos)

## Apply each method with all its configurations in period T1
metodos.T1<-list()
for (j in 1:length(ensayo)) {
  metodos.T1[[j]]<-calculo1(j)
}

## Obtains loss functions in T1 and T2
for (j in 1:length(ensayo)) {
  metodos.T1[[j]]<-calculo2(j) 
}
erroresT2<-list() 
for (j in 1:length(ensayo)) {
  erroresT2[[j]]<-calculo3(j)
}

## Calculates PCA and Manhattan distance to the method associated with the smallest error
resultado<-list()
pca<-resultado[[1]]<-resultado[[2]]<-resultado[[3]]<-list()
dis.ajust<-rep(NA,length(ensayo))
for (j in 1:length(ensayo)) {
  resultado<-calculo4.m(j) 
  if(str_detect(resultado[[3]][[j]][1,1], "mlexp")==TRUE) {dis.ajust[j]<-1} else 
    if(str_detect(resultado[[3]][[j]][1,1], "mlgamma")==TRUE) {dis.ajust[j]<-2} else {dis.ajust[j]<-3}
}

## alPCA establishes the percentiles of the adjusted distribution as cut-off points for method selection
cota<-c(0.05,0.50, 0.95) 
resultado[[4]]<-list()
for (j in 1:length(ensayo)) {  
  resultado[[4]][[j]]<-list()
  for (i in 1:length(cota)) {  
  resultado<-calculo5(j,i)
  }
}

## Obtain the predictions for the period T1+T2
metodos.T12<-list()
for (j in 1:length(ensayo)) {
  metodos.T12[[j]]<-calculo6(j)
}

##------------------------------------------------------------------------------------------
## Obtain the predictions for T3 according to the percentile chosen by the user
j=1;n_p=12;pred5=pred50=pred95=NULL;yu=1
b1<-rep(NA,n_p)
##for (j in 1:length(ensayo)) {
  for (i in 1:length(cota)) {
 
	  n_p<-12 #predicciones
	  dat<-resultado[[4]][[j]][[i]]
	  ##Matriz predicciones métodos para T2
	  pred_T2<-pred_T3<-matrix(NA,nrow=dim(dat)[1],ncol = n_p)
		  for (k in 1:dim(dat)[1]) {
			for (kk in 1:length(metodos.T1[[j]])) { 
			  #print(metodos.T12[[j]][[kk]]$mean)
			  if(dat[k,1]==kk){pred_T2[k,]<-metodos.T1[[j]][[kk]]$mean; pred_T3[k,]<-metodos.T12[[j]][[kk]]$mean}
			}
		  }
	  for (k in 1:n_p) {
		b1[k]<-pred_T3[,k]%*%dat[,11];if(i==1){pred5[k]=b1[k]}; if(i==2){pred50[k]=b1[k]}; if(i==3){pred95[k]=b1[k]}
	  }
  }
percentile.5=round(pred5, digits=2);percentile.50=round(pred50, digits=2);percentile.95=round(pred95, digits=2)
result_alPCA<-cbind.data.frame(percentile.5,percentile.50,percentile.95)

rm(dat);rm(erroresT2);rm(metodos.T1);rm(metodos.T12);rm(pca)
rm(pred_T2);rm(b1);
rm(cota);rm(dis.ajust);rm(ensayo);rm(i);rm(j);rm(k);rm(kk)
rm(n_p);rm(yu);rm(pred_T3);rm(resultado);rm(pred5);rm(pred50);rm(pred95)
rm(percentile.5);rm(percentile.50);rm(percentile.95)

## You can choose the percentile column and corresponding prediction
result_alPCA



