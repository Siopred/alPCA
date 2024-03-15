install.packages("nnfor")
install.packages("tsintermittent")

library(tsintermittent); library(nnfor)
library(Metrics); library(forecast); library(tseries);library(dplyr)
library(tidyverse); library(prophet)
library(forecTheta);library(tsfgrnn) 
library(yager);library(univariateML);
library(thief)

## Data upload
x <-readRDS("~/capsule/data/serieCO2_2000_2019.rds")
datos=list(x=x)

## Source file of functions
source("~/capsule/code/f_alPCA.monthly.r")

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
## Obtain the predictions for Obtains the predictions for T3 according to the percentile chosen by the user
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

#
##------------------------------------------------------------------------------------------

save(result_alPCA,file = "~/capsule/results/result_alPCA_monthly.Rdata")

