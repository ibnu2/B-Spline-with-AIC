#pengujian parameter
Uji_Parameter<-function(x,y,m,k){
  n<-length(y)
  knot<-c(k)
  knot<-as.matrix(knot)
  k1<-length(knot)
  orde = (m-1)
  bs1=bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  w = cbind(bs1)
  wtw = t(w) %*% w
  z = MPL(wtw)
  bet = (z %*% (t(w) %*% y))
  Beta = as.matrix(bet)
  cat("Nilai parameter adalah", "\n", "\n")
  print(Beta)
  S = (w %*% z %*% t(w))
  yhat = (w %*% bet)
  ybar <- mean(y)
  MSE = ((t(y-yhat)) %*% (y-yhat))/n
  MSE <- as.numeric(MSE)
  SSR <- sum((yhat-ybar)^2)
  MSR <- SSR/((orde+k1)-1)
  SSE <- sum((y-yhat)^2)
  MSER <- SSE/(n-(orde+k1))
  JKT = sum((y-ybar)^2)
  Fhit_1 = MSR/MSER
  cat("Nilai F hitung pengujian serentak adalah","\n", Fhit_1,"\n")
  cat("Kesimpulan hasil uji serentak","\n")
  cat("-----------------------------------","\n")
  cat("Analysis Of Variance (ANOVA)","\n")
  cat("===================================================","\n")
  cat("Sumber df   SS    MS      Fhitung","\n")
  cat("Regresi",((orde+k1)-1), "",SSR,"",MSR,"",Fhit_1,"\n")
  cat("Error", (n-(orde+k1)), "", SSE, "", MSER, "", "\n")
  cat("Total ", (n-1), "", JKT, "", "\n")
  print(MSE)
  cat("===================================================", "\n")
}


#uji parameter
Uji_Parameter(x,y,3, k=c(-0.1327273,-0.1185859,-0.1044444))
Uji_Parameter(x,y,4, k=c(-0.07616162,0.5319192,0.5460606,0.560202))
  

#pengujian parameter 2
Uji_tahap2<-function(x,y,m,k){
  knot<-c(k)
  knot<-as.matrix(knot)
  k<-length(knot)
  orde = (m-1)
  bs1=bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  w2 = cbind(bs1)
  w2 = as.matrix(w2)
  wtw2 = t(w2) %*% w2
  z2 = MPL(wtw2)
  bet2 = (z2 %*% (t(w2) %*% y))
  Beta2 = as.matrix(bet2)
  cat("Nilai parameter adalah", "\n", "\n")
  print(Beta2)
  
  # yhat = (w2 %*% bet2)
  # MSE = ((t(y-yhat)) %*% (y-yhat))/n
  # S = (w2 %*% z2 %*% t(w2))
  # ybar <- mean(y)
  # 
  # MSE <- as.numeric(MSE)
  # SSR <- sum((yhat-ybar)^2)
  # MSR <- SSR/((orde+k)-1)
  # SSE <- sum((y-yhat)^2)
  # MSER <- SSE/(n-(orde+k))
  # JKT = sum((y-ybar)^2)
  # Fhit_1 = MSR/MSER
  # cat("Nilai F hitung pengujian serentak adalah","\n", Fhit_1,"\n")
  # cat("Kesimpulan hasil uji serentak","\n")
  # cat("-----------------------------------","\n")
  # cat("Analysis Of Variance (ANOVA)","\n")
  # cat("===================================================","\n")
  # cat("Sumber df   SS    MS      Fhitung","\n")
  # cat("Regresi",((orde+k)-1), "",SSR,"",MSR,"",Fhit_1,"\n")
  # cat("Error", (n-(orde+k)), "", SSE, "", MSER, "", "\n")
  # cat("Total ", (n-1), "", JKT, "", "\n")
  # cat("===================================================", "\n")
}



#uji individu
uji_individu(x,y,2, k=c(1.252727))

#mencari nilai ftabel
var.test(x, y, conf.level = 0.95)
t.test(x, alternative = "greater" , conf.level = 0.95)
alpha = 0.05
qt(c(alpha/2, 1-(alpha/2)), df = 71)
qt(0.05, df = 71, lower.tail = FALSE)

#uji parameter tahap 2
Uji_tahap2(x,y,2,k=c(1.252727))

#predik
predik<-function(x,y,m,k){
  knot<-c(k)
  knot<-as.matrix(knot)
  k1<-length(knot)
  bs1=bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  w = cbind(bs1)
  w = as.matrix(w)
  wtw = t(w) %*% w
  z = MPL(wtw)
  bet = (z %*% (t(w) %*% y))
  Beta = as.matrix(bet)
  yhat = (w %*% Beta)
  
  predik = cbind(y,yhat)
  predik_w = as.matrix(predik)
  cat("Nilai prediksi adalah \n \n")
  print(predik_w)
  
}

predik(x,y,3, k=c(-0.1327273,-0.1185859,-0.1044444))
predik(x,y,4, k=c(-0.07616162,0.5319192,0.5460606,0.560202))

install.packages("MLmetrics")
library(MLmetrics)
#mape
yhat = predik(x,y,4, k=c(-0.07616162,0.5319192,0.5460606,0.560202))[,2]
# MAPE = mean(abs((y-yhat)/y))*100
# MAPE
MAPE(yhat, y)
plot(x,type = "o", col = "#1ECBE1")
lines(yhat, type = "l" ,col = "#E11ECB")
lines(yhat1, col = "#66710F")
legend(1, 1.15, legend=c("Data Aktual", "Data Prediksi AIC", "Data Prediksi GCV"), cex=1,
       fill = c("#1ECBE1","#E11ECB","#66710F")
)
# plot(yhat,y)
# plot(sapply(y, function(x) min(150, x)))

#mape
#mape
yhat1 = predik(x,y,3, k=c(-0.1327273,-0.1185859,-0.1044444))[,2]
# MAPE = mean(abs((y-yhat1)/y))*100
# MAPE
MAPE(yhat1, y)
plot(x,type = "l", col = "red")
lines(yhat1, col = "#CBE11E")
legend(2, 1, legend=c("Data Aktual", "Data Prediksi"), cex=0.5,
       fill = c("red","purple")
)
