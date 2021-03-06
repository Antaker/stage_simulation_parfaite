normalfromcauchy <- function(){ #changer le return pour obtenir une r�alisation ou le temps mis � l'obtenir#
  
  ##On initialise les valeurs ##
  C=0
  x=0
  n=0
  
  ##Tant que la bernouilli n'est pas 1##
  while(C!=1){
    
    ##On g�n�re une Cauchy standard##
    x = rcauchy(1)
    
    ##On genere une bernouilli du param�tre indiqu�##
    C = rbinom(1,1,((sqrt(exp(1))/2)*(1+x^2)*exp(-(1/2)*x^2)))
    
    ##mise � jour du nombre de lancers##
    n=n+1
  }
  return(x)
}

##On calcule le nombre moyen de lancers n�cessaires pour atteindre la sortie de l'algorithme pr�c�dent##

##On initialise le nombre moyen##
s=0

##On lance 10^6 fois l'algorithme qui retourne le nombre de lancers effectu�##
for (i in 1:10^6) {
  s=s+normalfromcauchy()
}

##On moyenne sur les 10^6 lancers
s = s/10^6



##On cr�� un histogramme illustrant 10^5 lancers de l'algorithme rendant une r�alisation de la normale standard##

##On initialise le tableau qui contiendra les 10^5 lancers
t = matrix(0,1,10^5)

##On remplit le tableau en lan�ant 10^5 fois l'algorithme##
for (i in 1:10^5) {   #changer le return pour obtenir des r�alisations
  t[i]=normalfromcauchy()

}

##On produit l'histogramme du tableau t et on rajoute la courbe de la densit� de la loi normale standard##
hist(t,prob=T)

curve(dnorm(x, 0, 1), col="red",xlim=c(-4,4),add=T)



AR_basic <- function(n,p,a){##fonction permettant de g�n�rer une r�alisation de la loi binomiale de param�tres n et p, conditionnellement � �tre sup�rieure o� �gale � a en utilisant la m�thode de rejet basique##
  
  ##On initialise le crit�re d'arr�t##
  C=0
  
  ##On initialise la Binomiale de param�tres n et p##
  t = 0
  
  ##On initialise le nombre de lancers##
  k = 0
  
  ##Tant que le crit�re d'arr�t n'est pas 1##
  while (C != 1){
    
    ##On g�n�re une r�alisation de la Binomiale de param�tres n et p##
    t = rbinom(1,n,p)
    
    ##On met � jour le crit�re d'arr�t selon que la binomiale soit sup�rieure � a ou non##
    C = (t>=a)*1
    
    ##On met � jour le nombre de lancers##
    k=k+1
  }
  
  ##On rend soit le nombre de lancers, soit la r�alisation de la loi demand�e (ie, k ou t))##
  return(k)
  
}


##On calcule ici le nombre moyen de lancers n�cessaires pour sortir de l'algorithme pr�c�dent
s1=0
for (i in 1:10^6) {   #changer le return pour obtenir le nombre moyen de lancers n�cessaires
  s1=s1+ AR_basic(10,0.1,5)
}
s1= s1/10^6


AR_chernoffs <- function(n,p,a){ ##fonction permettant de g�n�rer une r�alisation de la loi binomiale de param�tres n et p, conditionnellement � �tre sup�rieure o� �gale � a en utilisant les in�galit�s de Chernoff##
  
  
  ##On initialise la Bernouilli##
  C=0
  
  ##On cr�� la valeur gamma pr�sent�e pr�c�demment##
  gam = (a/n)*(1-p)/(p*(1-(a/n)))
  
  ##On initialise la Binomiale ##
  t = 0
  
  ##On initialise le nombre de lancers##
  k = 0
  
  ##Tant que la Bernouilli n'est pas 1##
  while (C != 1){
    
    ##On g�n�re n Bernouillis de param�tres a/n## 
    t = rbinom(n,1,a/n)
    
    ##On effectue la somme pour obtenir la binomiale consid�r�e##
    x = sum(t)
    
    ##On g�n�re une r�alisation de la Bernouilli du param�tre indiqu�
    C = rbinom(1,1,(x>=a)*1*gam**(a-x))
    
    ##On met � jour le nombre de lancers effectu�s##
    k=k+1
  }
  
  ##On rend soit le nombre de lancers, soit la r�alisation de la loi demand�e (ie, k ou x))##
  return(k)
  
}


##On calcule ici le nombre moyen de lancers n�cessaires pour sortir de l'algorithme pr�c�dent

s2=0
for (i in 1:10^6) {   #changer le return pour obtenir le nombre moyen de lancers n�cessaires
  s2=s2 + AR_chernoffs(10,0.1,5)
}
s2 = s2/10^6


AR_boule_unite <- function(n){
  
  u = runif(n,-1,1)
  s = u**2
  s = sum(s)
  i=1
  while(s>1){
    u = runif(n,-1,1)
    s=u**2
    s = sum(s)
    i=i+1
  }
  return(u)
  
  
}


s3=0
for (i in 1:10^6) {   #changer le return pour obtenir le nombre moyen de lancers n�cessaires
  s3=s3+ AR_boule_unite(5)
}
s3= s3/10^6

t=c()
for(i in 1:100){
  
  t = c(t, (gamma(1+i/2)*2**i)/(pi**(i/2)))
  
}


u=c()
for(i in 1:10^5){
  
  u = rbind(u,AR_boule_unite(2))
  
}
u = u[,2]/u[,1]
u=as.matrix(u)
hist(u,breaks = 15)
curve((gamma(1+x/2)*2**x)/(pi**(x/2)),0,12, xlab = NULL)






