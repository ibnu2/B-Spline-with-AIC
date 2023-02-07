#pengujian parameter
Uji_Parameter<-function(x,y,m,k){
  knot<-c(k)
  knot<-as.matrix(knot)
  k<-length(knot)
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
  MSR <- SSR/((orde+k)-1)
  SSE <- sum((y-yhat)^2)
  MSER <- SSE/(n-(orde+k))
  JKT = sum((y-ybar)^2)
  Fhit_1 = MSR/MSER
  cat("Nilai F hitung pengujian serentak adalah","\n", Fhit_1,"\n")
  cat("Kesimpulan hasil uji serentak","\n")
  cat("-----------------------------------","\n")
  cat("Analysis Of Variance (ANOVA)","\n")
  cat("===================================================","\n")
  cat("Sumber df   SS    MS      Fhitung","\n")
  cat("Regresi",((orde+k)-1), "",SSR,"",MSR,"",Fhit_1,"\n")
  cat("Error", (n-(orde+k)), "", SSE, "", MSER, "", "\n")
  cat("Total ", (n-1), "", JKT, "", "\n")
  cat("===================================================", "\n")
  
  #pengujian individu
  cat("Uji Individu \n")
  Sbett<-matrix(0,nrow = m+k+1, ncol = 1)
  Sbett<-sqrt(diag(MSE*(t(w)%*%w)))
  Sbett<-as.matrix(Sbett)
  thitung<-bet/Sbett 
  thitung
  # SEbeta<-matrix(0, nrow = (orde+k), ncol = (orde+k))
  # Ai<- solve(wtw)
  # SEbeta <- sqrt(diag((MSE*Ai)))
  # SEbeta <- as.matrix(SEbeta)
  # Thit<- Beta/SEbeta
  # print(Thit)
}


#uji model paramter
Uji_Parameter(x,y,2,k=c(1.252727))
Uji_Parameter(x,y,4,k=c(1.016768, 1.046263, 1.075758, 1.105253))  
  

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
  k<-length(knot)
  bs1=bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  w = cbind(bs1[,2:3])
  w = as.matrix(w)
  wtw = t(w) %*% w
  z = MPL(wtw)
  bet = (z %*% (t(w) %*% y))
  Beta = as.matrix(bet)
  yhat = (w %*% bet)
  
  predik = cbind(y,yhat)
  predik_w = as.matrix(predik)
  cat("Nilai prediksi adalah \n \n")
  print(predik_w)
  
}


#uji terakhir
predic <- function(x,y,m,k)
{
  knot<-c(k)
  knot<-as.matrix(knot)
  k<-length(knot)
  n<-length(x)
  bs1=bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  bs1n = bs1[,2:3]
  wn = cbind(bs1n)
  wn = as.matrix(wn)
  
  N <- matrix(0, ncol = n, nrow = n)
  diag(N) = (1/x)
  
  A <- t(wn) %*% N %*% wn
  A = as.matrix(A)
  B <- t(wn) %*% N %*% y
  Beta <- A %*% B
  Beta_w <- as.matrix(Beta)
  yhat_w <- wn %*% Beta_w
  
  pred = cbind(y,yhat_w)
  pred = as.matrix(pred)
  cat("Nilai prediksi adalah \n" )
  print(pred)
}

predic(x,y,2,k=c(1.252727))

#mape
yhat = predik(x,y,2,k = c(1.252727))[,2]
MAPE = mean(abs(y-yhat)/y)*100
MAPE
plot(y,type = "l", col = "red")
lines(yhat, col = "blue")
# plot(yhat,y)
# plot(sapply(y, function(x) min(150, x)))

#mape
#mape
yhat1 = predik(x,y,4,k = c(1.016768, 1.046263, 1.075758, 1.105253))[,2]
MAPE = mean(abs(y-yhat)/y)*100
MAPE
plot(y,type = "l", col = "red")
lines(yhat, col = "blue")
