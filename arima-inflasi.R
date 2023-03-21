library(forecast)
library(tseries)

x<-c(0.55,-0.05,0.15,0.10,0.08,0.46,0.56,-0.26,-0.11,0.13,0.46,0.57,0.42,-0.08,0.26,0.46,0.42,0.25,0.05,0.07,-0.07,0.18,0.31,0.46,0.27,0.40,0.07,-0.24,0.22,0.08,-0.08,-0.04,0.03,0.08,0.13,0.48,0.54,0.14,0.08,0.01,0.07,0.05,0.11,0.05,-0.17,0.24,0.45,0.71,0.59,0.05,0.77,1.14,0.75,0.52,0.47,-0.12,1.05,0.11,0.32,0.65)
y<-c(-0.05,0.15,0.10,0.08,0.46,0.56,-0.26,-0.11,0.13,0.46,0.57,0.42,-0.08,0.26,0.46,0.42,0.25,0.05,0.07,-0.07,0.18,0.31,0.46,0.27,0.40,0.07,-0.24,0.22,0.08,-0.08,-0.04,0.03,0.08,0.13,0.48,0.54,0.14,0.08,0.01,0.07,0.05,0.11,0.05,-0.17,0.24,0.45,0.71,0.59,0.05,0.77,1.14,0.75,0.52,0.47,-0.12,1.05,0.11,0.32,0.65,0.17)

inflasi = cbind(x,y)
inflasits = ts(inflasi)
adf.test(x)
plot(inflasits)

Arima1 = auto.arima(inflasits[,1])
Arima1
autoplot(Arima1)
autoplot(forecast(Arima1))

autoplot(fArima1)
a=forecast(Arima1)
plot(x,type = "o", col = "#1ECBE1")
lines(a)



acf(x)
pacf(x)
r<-diff(y)
summary(r.arma <- arma(r, order =c(0,2)))
r<-diff(x)
summary(r.arma <- arma(r, order =c(0,2)))
fit1=Arima(inflasi[1], order=c(0,1,2))
Arima2 = auto.arima(inflasi[,2])
Arima2

arima.predict(inflasi[1])
arima.forecase(inflasi[1])
autoplot(forecast(arimamodel1))
checkresiduals(arimamodel1)
autoplot(forecast(arimamodel1))
lines(arimamodel1, type = "l" ,col = "green")
accuracy(Arima1)
accuracy(Arima2)
a=forecast(arimamodel1)
autoplot(a[1])
arimamodel1=Arima(inflasi[,1], order=c(0,1,2))
arimamodel2=arima(inflasi[,2], order=c(0,1,2))
accuracy(arimamodel1)
accuracy(arimamodel2)
plot(arimamodel1)

#uji terakhir
predic <- function(x,y,m,k)
{
  knot<-c(k)
  knot<-as.matrix(knot)
  k<-length(knot)
  n<-length(x)
  bs1=bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  bs1n = bs1[,2:2]
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

predic(x,y,4, k=c(-0.07616162,0.5319192,0.5460606,0.560202))
