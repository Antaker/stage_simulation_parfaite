normalfromcauchy <- function(){ #changer le return pour obtenir une réalisation ou le temps mis à l'obtenir#
  C=0
  x=0
  n=0
  while(C!=1){
    x = rcauchy(1)
    C = rbinom(1,1,((sqrt(exp(1))/2)*(1+x^2)*exp(-(1/2)*x^2)))
    n=n+1
  }
  return(x)
}

s=0
for (i in 1:10^6) {   #changer le return pour obtenir le nombre moyen de lancers nécessaires
  s=s+normalfromcauchy()
}
s = s/10^6

t = matrix(0,1,10^5)
for (i in 1:10^5) {   #changer le return pour obtenir des réalisations
  t[i]=normalfromcauchy()

}
hist(t)
hist(sort(t))


AR_basic <- function(n,p,a){
  C=0
  t=0
  k = 0
  while (C != 1){
    
    t = rbinom(1,n,p)
    C = (t>=a)*1
    k=k+1
  }
  
  return(k)
  
}

s1=0
for (i in 1:10^6) {   #changer le return pour obtenir le nombre moyen de lancers nécessaires
  s1=s1+ AR_basic(10,0.1,5)
}
s1= s1/10^6


AR_chernoffs <- function(n,p,a){
  C=0
  gam = (a/n)*(1-p)/(p*(1-(a/n)))
  t=0
  k = 0
  while (C != 1){
    
    t = rbinom(n,1,a/n)
    x = sum(t)
    C = rbinom(1,1,(x>=a)*1*gam**(a-x))
    k=k+1
  }
  
  return(k)
  
}

s2=0
for (i in 1:10^5) {   #changer le return pour obtenir le nombre moyen de lancers nécessaires
  s2=s2 + AR_chernoffs(10,0.1,5)
}
s2 = s2/10^5