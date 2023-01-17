library(splines)
gcv11<-function(x,y,a,b)
{
  n=length(y)
  gcv=10^10
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
      mu= w %*% beta
      MSE = t(y-mu) %*% (y-mu)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
      I[j,j]=1
      AIc = (n+((n*log(2*pi))+n*(log(MSE/n))+(2*(1+m))))
      GCV = MSE / (k/n)^2
      if(gcv>=AIc)
      {
        gcv=AIc
        knot1=k1[j1]
        det=deter
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "AIC = ", AIc, "GCV = ", GCV, "determinan = ", deter,"\n")
    
    gcv=10^10
  }
}

gcv22<-function(x,y,a,b)
{
  n=length(y)
  gcv=10^10
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
      mu= w %*% beta
      MSE = t(y-mu) %*% (y-mu)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
        I[j,j]=1
      AIc = (n+((n*log(2*pi))+n*(log(MSE/n))+(2*(2+m))))
      GCV = MSE / (k/n)^2
      if(gcv>=AIc)
      {
        gcv=AIc
        knot1=k1[j1]
        knot2=k2[j1]
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "knot 2 = ", knot2, "AIC = ", AIc, "GCV = ", GCV, "determinan = ", deter,"\n")
    
    gcv=10^10
  }
}

gcv33<-function(x,y,a,b)
{
  n=length(y)
  gcv=10^10
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
      mu= w %*% beta
      MSE = t(y-mu) %*% (y-mu)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
        I[j,j]=1
      AIc = (n+((n*log(2*pi))+n*(log(MSE/n))+(2*(3+m))))
      GCV = MSE / (k/n)^2
      if(gcv>=AIc)
      {
        gcv=AIc
        knot1=k1[j1]
        knot2=k2[j1]
        knot3=k3[j1]
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "knot 2 = ", knot2, "knot 3 = ", knot3, "AIC = ", AIc, "GCV = ", GCV, "determinan = ", deter,"\n")
    
    gcv=10^10
  }
}

gcv44<-function(x,y,a,b)
{
  n=length(y)
  gcv=10^10
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
      mu= w %*% beta
      MSE = t(y-mu) %*% (y-mu)/n
      I <- matrix(0, ncol = n, nrow = n)
      for (j in 1:n) 
        I[j,j]=1
      AIc = (n+((n*log(2*pi))+n*(log(MSE/n))+(2*(4+m))))
      GCV = MSE / (k/n)^2
      if(gcv>=AIc)
      {
        gcv=AIc
        knot1=k1[j1]
        knot2=k2[j1]
        knot3=k3[j1]
        knot4=k4[j1]
        orde1=m
      }
    }
    cat("orde = ", orde1, "knot 1 = ", knot1, "knot 2 = ", knot2, "knot 3 = ", knot3, "knot 4 = ", knot4, "AIC = ", AIc, "GCV = ", GCV, "determinan = ", deter,"\n")
    
    gcv=10^10
  }
}

bspline<-function(x,y,m,k=c(...)){
  n<-length(y)
  knot<-c(k)
  knot<-as.matrix(knot)
  k<-length(knot)
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
  k = sum(diag(I-S))
  GCV=(MSE/(k/n)^2)
  AIc = (n+((n*log(2*pi))+n*(log(MSE/n))+(2*(4+m))))
  cat("NIlai GCV adalah ","\n",GCV,"\n")
  cat("NIlai AIC adalah ","\n",AIc,"\n")
}

#pengujian parameter
Uji_Parameter<-function(x,y,m,k){
knot<-c(k1)
knot<-as.matrix(knot)
k<-length(knot)
m=2
orde<-(m-1)
bs1=bs(x,df=NULL, knots=k1, degree=m-1, intercept=TRUE, Boundary.knots=range(x))
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
}

#pengujian individu
SEbeta<-matrix(0, nrow = (orde+k), ncol = (orde+k))
Ai<- solve(wtw)
SEbeta <- sqrt(diag((MSE*Ai)))
SEbeta <- as.matrix(SEbeta)
Thit<- Beta/SEbeta


#mencari knot dari orde 2 sampai 4 dan 1knot sampai 4knot
gcv11(x,y,a,b)
gcv22(x,y,a,b)
gcv33(x,y,a,b)
gcv44(x,y,a,b)

#menentukan parameter dari titik knot optimal
bspline(x,y,2,k=c(1.252727))
bspline(x,y,4,k=c(1.016768, 1.046263, 1.075758, 1.105253))

#uji model paramter
Uji_Parameter(x,y,2,k=c(1.252727))
Uji_Parameter(x,y,4,k=c(1.016768, 1.046263, 1.075758, 1.105253))

#tambah kondisi buat di orde jika m = 2 maka knot mulai dari 1 jika m = 3 maka mulai dari 2 jika m = 4 maka mulai dari 3
