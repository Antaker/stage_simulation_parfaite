X = matrix(sample(c(-1,1),100,replace = TRUE),100,100)

IsingGibbsUpdate <- function(X,u,v,beta,mu){       ## X = graphe (rectangulaire ou carré), u = lancé uniforme sur [0,1], v = noeud de X ##
  
  ##On entoure la matrice de valeurs aberrantes pour la manipuler plus facilement##
  ##On change les indices pour s'y conformer##
  
  i = v[1] + 1
  j = v[2] + 1
  
  X = rbind(rep(0,dim(X)[2]),X)
  X = rbind(X,rep(0,dim(X)[2]))
  X = cbind(rep(0,dim(X)[1]),X)  
  X = cbind(X,rep(0,dim(X)[1]))
  
  if(u< exp(mu + beta*((X[i+1,j]==1)*1 + (X[i-1,j]==1)*1 +(X[i,j+1]==1)*1 + (X[i,j-1]==1)*1))/((exp(mu + beta*((X[i+1,j]==1)*1 + (X[i-1,j]==1)*1 +(X[i,j+1]==1)*1 + (X[i,j-1]==1)*1))) + exp(-mu + beta*((X[i+1,j]==-1)*1 + (X[i-1,j]==-1)*1 +(X[i,j+1]==-1)*1 + (X[i,j-1]==-1)*1)))){
    X[i,j] =  1
  }
  else{
    X[i,j] = -1
  }
  
  X=X[,-1]
  X=X[,-dim(X)[2]]
  X=X[-1,]
  X=X[-dim(X)[1],]
  
  X
} 


IsingGibbssteps <- function(X,beta,mu,n){ ##effectue n pas de la chaîne de Markov pour le sampler de Gibbs, pour un graphe rectangulaire X##  

  for (i in 1:n){
    k = sample(1:dim(X)[1],1)
    l = sample(1:dim(X)[2],1)
    u = runif(1)
    X=IsingGibbsUpdate(X,u,c(k,l),beta,mu)
  }

  X

}

IsingGibbsconditionnalstop <- function(X,beta,mu,epsilon){ ##effectue un nombre aléatoire de pas de la chaîne de Markov pour le sampler de Gibbs, pour un graphe rectangulaire X##  
  Y=c(mean(X==1))
  bol = TRUE
  i=1
  while (bol){
    k = sample(1:dim(X)[1],1)
    l = sample(1:dim(X)[2],1)
    u = runif(1)
    X=IsingGibbsUpdate(X,u,c(k,l),beta,mu)
    Y = c(Y,mean(X==1))
    ##if ((i%%((dim(X)[1]**2)/epsilon)==0)){image(X)}
    image(X)
    if ((i>=((log(dim(X)[1]**2)*dim(X)[1]**2)/epsilon)) && ((mean(Y[i+1-floor((dim(X)[1]**(1/2))/epsilon):i])-Y[i])<epsilon)){
      break()
    }
    i=i+1
  }
  
  return(list(X,i))
  
}

MetropolisHastingsUpdate <- function(X,u,v,beta){       ## X = graphe (rectangulaire ou carré), u = lancé uniforme sur [0,1], v = noeud de X ##
  
  ##On entoure la matrice de valeurs aberrantes pour la manipuler plus facilement##
  ##On change les indices pour s'y conformer##
  
  i = v[1] + 1
  j = v[2] + 1
  
  c=sample(c(-1,1),1,prob=c(1/2,1/2))
  
  X = rbind(rep(0,dim(X)[2]),X)
  X = rbind(X,rep(0,dim(X)[2]))
  X = cbind(rep(0,dim(X)[1]),X)  
  X = cbind(X,rep(0,dim(X)[1]))
  
  if(u< exp(beta*((X[i+1,j]==c)*1 + (X[i-1,j]==c)*1 +(X[i,j+1]==c)*1 + (X[i,j-1]==c)*1))/exp(beta*((X[i+1,j]==X[i,j])*1 + (X[i-1,j]==X[i,j])*1 +(X[i,j+1]==X[i,j])*1 + (X[i,j-1]==X[i,j])*1))){
    X[i,j] =  c
  }
  
  X=X[,-1]
  X=X[,-dim(X)[2]]
  X=X[-1,]
  X=X[-dim(X)[1],]
  
  X
} 

MetropolisHastingssteps <- function(X,beta,n){ ##effectue n pas de la chaîne de Markov pour le sampler de Gibbs, pour un graphe rectangulaire X##  
  
  for (i in 1:n){
    k = sample(1:dim(X)[1],1)
    l = sample(1:dim(X)[2],1)
    u = runif(1)
    X=MetropolisHastingsUpdate(X,u,c(k,l),beta)
  }
  
  X
  
}

MetropolisHastingsconditionnalstop <- function(X,beta,epsilon){ ##effectue un nombre aléatoire de pas de la chaîne de Markov pour le sampler de Gibbs, pour un graphe rectangulaire X##  
  Y=c(mean(X==1))
  bol = TRUE
  i=1
  while (bol){
    k = sample(1:dim(X)[1],1)
    l = sample(1:dim(X)[2],1)
    u = runif(1)
    X=MetropolisHastingsUpdate(X,u,c(k,l),beta)
    Y = c(Y,mean(X==1))
    if ((i%%((dim(X)[1]**2)/epsilon)==0)){image(X)}
    if ((i>=((log(dim(X)[1]**2)*dim(X)[1]**2)/epsilon)) && ((mean(Y[i+1-floor((dim(X)[1]**(1/2))/epsilon):i])-Y[i])<epsilon)){
      break()
    }
    i=i+1
  }
  
  return(list(X,i))
  
}
