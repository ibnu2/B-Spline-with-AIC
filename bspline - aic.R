library(splines)
MPL<-function(x,eps=1e-009)
{
  x<-as.matrix(x)
  xsvd<-svd(x)
  diago<-xsvd$d[xsvd$d>eps]
  if(length(diago)==1)
  {
    xplus<-as.matrix(xsvd$v[,1]) %*% t(as.matrix(xsvd$u[,1])/diago)
  }
  else
  {
    xplus<-xsvd$v[,1:length(diago)] %*% diag(1/diago) %*% t(xsvd$u[,1:length(diago)])
  }
  return(xplus)
}

# x<-c(1.25,1.01,0.56,0.21,1.08,2.51,1.31,0.67,1.15,0.62,0.07,-0.11,0.09,0.32,0.18,-0.34,0.27,0.18,0.32,0.77,0.8,-0.03,0.09,0.24,0.57,0.31,0.13,0.25,0.14,1.26,1.4,0.43,1.06,0.28,0.62,0.72,0.84,0.1,0.21,-0.28,0.13,0.26,0.9,0.63,0.19,0.04,0.33,0.48,0.25,0.1,0.36,0.11,0.05,0.75,0.76,0.42,0.19,0.38,0.2,0.66,0.96,0.93,0.79,-0.3,-0.29,0.84,2.58,0.87,-0.24,0.61,0.2,0.17)
# y<-c(1.01,0.56,0.21,1.08,2.51,1.31,0.67,1.15,0.62,0.07,-0.11,0.09,0.32,0.18,-0.34,0.27,0.18,0.32,0.77,0.8,-0.03,0.09,0.24,0.57,0.31,0.13,0.25,0.14,1.26,1.4,0.43,1.06,0.28,0.62,0.72,0.84,0.1,0.21,-0.28,0.13,0.26,0.9,0.63,0.19,0.04,0.33,0.48,0.25,0.1,0.36,0.11,0.05,0.75,0.76,0.42,0.19,0.38,0.2,0.66,0.96,0.93,0.79,-0.3,-0.29,0.84,2.58,0.87,-0.24,0.61,0.2,0.17,1.05)
# x=x*100
# y=y*100

x<-c(0.55,-0.05,0.15,0.10,0.08,0.46,0.56,-0.26,-0.11,0.13,0.46,0.57,0.42,-0.08,0.26,0.46,0.42,0.25,0.05,0.07,-0.07,0.18,0.31,0.46,0.27,0.40,0.07,-0.24,0.22,0.08,-0.08,-0.04,0.03,0.08,0.13,0.48,0.54,0.14,0.08,0.01,0.07,0.05,0.11,0.05,-0.17,0.24,0.45,0.71,0.59,0.05,0.77,1.14,0.75,0.52,0.47,-0.12,1.05,0.11,0.32,0.65)
y<-c(-0.05,0.15,0.10,0.08,0.46,0.56,-0.26,-0.11,0.13,0.46,0.57,0.42,-0.08,0.26,0.46,0.42,0.25,0.05,0.07,-0.07,0.18,0.31,0.46,0.27,0.40,0.07,-0.24,0.22,0.08,-0.08,-0.04,0.03,0.08,0.13,0.48,0.54,0.14,0.08,0.01,0.07,0.05,0.11,0.05,-0.17,0.24,0.45,0.71,0.59,0.05,0.77,1.14,0.75,0.52,0.47,-0.12,1.05,0.11,0.32,0.65,0.17)
plot(x, type = "l")

# ynaive = x

gcv11<-function(x,y,a,b)
{
  n=length(y)
  print(n)
  gcv=10^10
  aic=10^10
  t1=seq(min(x),max(x),length.out=100)
  p=rep(0,(n-2))
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1]
  }
  comb1=combn(p,1,FUN = NULL)
  c1=t(comb1)
  for (m in 2:4)
  {
    k1=c1[,1]
    v1=length(k1)
    for (j1 in 1:v1) {
      bs1=bs(x, df=NULL, knot=k1[j1], degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      w = cbind(bs1)
      wtw = t(w) %*% w
      deter=det(wtw)
      z = MPL(wtw)
      beta= z %*% (t(w) %*% y)
      S= w %*% z %*% t(w)
      mu= w %*% beta
      MSE = t(y-mu) %*% (y-mu)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
        I[j,j]=1
      k = sum(diag(I-S))
      AIc = (n+((n*log(2*pi))+n*(log(MSE))+(2*(1+m))))
      GCV = MSE / (k/n)^2
      if(gcv>=GCV)
      {
        gcv=GCV
        knot1=k1[j1]
        det=deter
        orde1=m
      }
      if(aic>=AIc)
      {
        aic=AIc
        knot1=k1[j1]
        det=deter
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "GCV = ", GCV, "AIC = ", AIc, "determinan = ", deter,"\n")
    
    gcv=10^10
    aic=10^10
  }
}

aic1<-function(x,y,a,b)
{
  n=length(y)
  aic=10^10
  t1=seq(min(x),max(x),length.out=100)
  p=rep(0,(n-2))
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1]
  }
  comb1=combn(p,1,FUN = NULL)
  c1=t(comb1)
  for (m in 2:4)
  {
    k1=c1[,1]
    v1=length(k1)
    for (j1 in 1:v1) {
      bs1=bs(x, df=NULL, knot=k1[j1], degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      w = cbind(bs1)
      wtw = t(w) %*% w
      deter=det(wtw)
      z = MPL(wtw)
      beta= z %*% (t(w) %*% y)
      S= w %*% z %*% t(w)
      mu= w %*% beta
      MSE = t(y-mu) %*% (y-mu)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
        I[j,j]=1
      k = sum(diag(I-S))
      AIc = (n+((n*log(2*pi))+n*(log(MSE))+(2*(1+m))))
      if(gcv>=GCV)
      {
        gcv=GCV
        knot1=k1[j1]
        knot2=k2[j1]
        orde1=m
      }
      if(aic>=AIc)
      {
        aic=AIc
        knot1=k1[j1]
        det=deter
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "AIC = ", AIc, "determinan = ", deter,"\n")
    
    aic=10^10
  }
}

gcv22<-function(x,y,a,b)
{
  n=length(y)
  gcv=10^10
  aic=10^10
  t1=seq(min(x),max(x),length.out=100)
  p=rep(0,(n-2))
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1]
  }
  comb1=combn(p,2,FUN = NULL)
  c1=t(comb1)
  for (m in 2:4)
  {
    k1=c1[,1]
    k2=c1[,2]
    v1=length(k1)
    for (j1 in 1:v1) {
      bs1=bs(x, df=NULL, knot=c(k1[j1],k2[j1]), degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      w = cbind(bs1)
      wtw = t(w) %*% w
      deter=det(wtw)
      z = MPL(wtw)
      beta= z %*% (t(w) %*% y)
      S= w %*% z %*% t(w)
      mu= w %*% beta
      MSE = t(y-mu) %*% (y-mu)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
        I[j,j]=1
      k = sum(diag(I-S))
      AIc = (n+((n*log(2*pi))+n*(log(MSE))+(2*(2+m))))
      GCV = MSE / (k/n)^2
      if(gcv>=GCV)
      {
        gcv=GCV
        knot1=k1[j1]
        knot2=k2[j1]
        orde1=m
      }
      if(aic>=AIc)
      {
        aic=AIc
        knot1=k1[j1]
        knot2=k2[j1]
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "knot 2 = ", knot2, "GCV = ", GCV, "AIC = ", AIc, "determinan = ", deter,"\n")
    
    gcv=10^10
    aic=10^10
  }
}

gcv33<-function(x,y,a,b)
{
  n=length(y)
  gcv=10^10
  aic=10^10
  t1=seq(min(x),max(x),length.out=100)
  p=rep(0,(n-2))
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1]
  }
  comb1=combn(p,3,FUN = NULL)
  c1=t(comb1)
  for (m in 2:4)
  {
    k1=c1[,1]
    k2=c1[,2]
    k3=c1[,3]
    v1=length(k1)
    for (j1 in 1:v1) {
      bs1=bs(x, df=NULL, knot=c(k1[j1],k2[j1],k3[j1]), degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      w = cbind(bs1)
      wtw = t(w) %*% w
      deter=det(wtw)
      z = MPL(wtw)
      beta= z %*% (t(w) %*% y)
      S= w %*% z %*% t(w)
      mu= w %*% beta
      MSE = t(y-mu) %*% (y-mu)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
        I[j,j]=1
      k = sum(diag(I-S))
      AIc = (n+((n*log(2*pi))+n*(log(MSE))+(2*(3+m))))
      GCV = MSE / (k/n)^2
      if(gcv>=GCV)
      {
        gcv=GCV
        knot1=k1[j1]
        knot2=k2[j1]
        knot3=k3[j1]
        orde1=m
      }
      if(aic>=AIc)
      {
        aic=AIc
        knot1=k1[j1]
        knot2=k2[j1]
        knot3=k3[j1]
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "knot 2 = ", knot2, "knot 3 = ", knot3, "AIC = ", AIc, "GCV = ", GCV, "determinan = ", deter,"\n")
    
    gcv=10^10
    aic=10^10
  }
}

gcv44<-function(x,y,a,b)
{
  n=length(y)
  gcv=10^10
  aic = 10^10
  t1=seq(min(x),max(x),length.out=100)
  p=rep(0,(n-2))
  for(z in 1:(n-2))
  {
    p[z]=t1[z+1]
  }
  comb1=combn(p,4,FUN = NULL)
  c1=t(comb1)
  for (m in 2:4)
  {
    k1=c1[,1]
    k2=c1[,2]
    k3=c1[,3]
    k4=c1[,4]
    v1=length(k1)
    for (j1 in 1:v1) {
      bs1=bs(x, df=NULL, knot=c(k1[j1],k2[j1],k3[j1],k4[j1]), degre=m-1, intercept=TRUE, Boundary.knots=range(x))
      w = cbind(bs1)
      wtw = t(w) %*% w
      deter=det(wtw)
      z = MPL(wtw)
      beta= z %*% (t(w) %*% y)
      S= w %*% z %*% t(w)
      mu= w %*% beta
      MSE = t(y-mu) %*% (y-mu)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
        I[j,j]=1
      k = sum(diag(I-S))
      AIc = (n+((n*log(2*pi))+n*(log(MSE))+(2*(4+m))))
      GCV = MSE / (k/n)^2
      if(gcv>=GCV)
      {
        gcv=GCV
        knot1=k1[j1]
        knot2=k2[j1]
        knot3=k3[j1]
        knot4=k4[j1]
        orde1=m
      }
      if(aic>=AIc)
      {
        aic=AIc
        knot1=k1[j1]
        knot2=k2[j1]
        knot3=k3[j1]
        knot4=k4[j1]
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "knot 2 = ", knot2, "knot 3 = ", knot3, "knot 4 = ", knot4, "AIC = ", AIc, "GCV = ", GCV, "determinan = ", deter,"\n")
    
    gcv=10^10
    aic=10^10
  }
}

#estimasi parameter
bspline<-function(x,y,m,k){
  n<-length(y)
  print(n)
  knot<-c(k)
  knot<-as.matrix(knot)
  k1<-length(knot)
  bs1 = bs(x,df=NULL, knots=k, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
  w = cbind(bs1)
  wtw = t(w) %*% w
  z = MPL(wtw)
  bet = (z %*% (t(w) %*% y))
  cat("Nilai parameter adalah", "\n", bet, "\n")
  S = (w %*% z %*% t(w))
  yhat = (w %*% bet)
  MSE = ((t(y-yhat)) %*% (y-yhat))/n
  I <- matrix(0, ncol=n, nrow = n)
  for (j in 1:n) 
    I[j,j] = 1
  l = sum(diag(I-S))
  GCV=(MSE/(l/n)^2)
  AIc = (n+((n*log(2*pi))+n*(log(MSE))+(2*(3+m))))
  cat("NIlai GCV adalah ","\n",GCV,"\n")
  cat("NIlai AIC adalah ","\n",AIc,"\n")
}
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
n

#mencari knot dari orde 2 sampai 4 dan 1knot sampai 4knot
gcv11(x,y,a,b)
gcv22(x,y,a,b)
gcv33(x,y,a,b)
gcv44(x,y,a,b)

#menentukan parameter dari titik knot optimal
bspline(x,y,3,k=c(-0.1327273,-0.1185859,-0.1044444))
bspline(x,y,4,k=c(-0.07616162,0.5319192,0.5460606,0.560202))

#uji parameter
Uji_Parameter(x,y,3, k=c(-0.1327273,-0.1185859,-0.1044444))
Uji_Parameter(x,y,4, k=c(-0.07616162,0.5319192,0.5460606,0.560202))

#f table
qf(p = 0.05, df1=1, df2=70, lower.tail = FALSE)

#uji individu
SEbeta <- matrix(0, nrow = (orde+k), ncol = (orde+k))
Ai <- solve(wtw)
SEbeta <- sqrt(diag(MSE * Ai))
SEbeta <- as.matrix(SEbeta)
Thit <- bet/SEbeta
Thit
wtw

#tambah kondisi buat di orde jika m = 2 maka knot mulai dari 1 jika m = 3 maka mulai dari 2 jika m = 4 maka mulai dari 3
