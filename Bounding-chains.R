
Bounding_chain_voisins_1 <- function(Y,k,l){ ## Entrée : bounding chain, coordonnées d'un noeud du graphe. Sortie : nombre de voisins de label 1 du noeud choisi##
  
  
  ##Comme on ne vérifie que les 1, on a besoin seulement de la deuxième partie de la bounding chain, de par sa construction##
  X = Y
  
  ##On entoure la matrice de valeurs aberrant##
  X = rbind(rep(-99,dim(X)[2]),X)
  X = rbind(X,rep(-99,dim(X)[2]))
  X = cbind(rep(-99,dim(X)[1]),X)  
  X = cbind(X,rep(-99,dim(X)[1]))
  
  ##On change les coordonnées du noeud d'entrée en conséquence##
  k = k+1
  l = l+1
  
  ##On compte le nombre de 1 autour du noeud considéré##
  
  N1 = 1*(X[k-1,l]==1) + 1*(X[k,l+1]==1) + 1*(X[k+1,l]==1) + 1*(X[k,l-1]==1)
  
  return(N1)
  
}

Bounding_chain_voisins_01 <- function(Y,k,l){ ## Entrée : bounding chain, coordonnées d'un noeud du graphe. Sortie : nombre de voisins de label 1 du noeud choisi##
  
  X = Y

  
  ##On entoure les matrices de valeurs aberrantes##
  X = rbind(rep(-99,dim(X)[1]),X)
  X = rbind(X,rep(-99,dim(X)[2]))
  X = cbind(rep(-99,dim(X)[1]),X)  
  X = cbind(X,rep(-99,dim(X)[1]))
  
  ##On change les coordonnées du noeud d'entrée en conséquence##
  k = k+1
  l = l+1
  
  ##On compte le nombre de 2 autour du noeud considéré##
  
  N01 = 1*(X[k-1,l]==2) + 1*(X[k,l+1]==2) + 1*(X[k+1,l]==2) + 1*(X[k,l-1]==2)
  
  return(N01)
  
  
}



Bounding_chain_update_HCGM <- function(Y,u,lambda){##mise à jour d'un noeud pour la bounding chain du modèle hard core gas##
  
  for(i in 1:length(u)){
    
    ##On tire un noeud uniformément dans le graphe##
    k = sample(1:dim(Y)[1],1)
    l = sample(1:dim(Y)[2],1)
    
    ##On compte le nombre de voisins de label {1} et de label {0,1} = 2##
    
    N1 = Bounding_chain_voisins_1(Y,k,l)
    
    N01 = Bounding_chain_voisins_01(Y,k,l)
    
    ##mise à jour selon le pseudo code fourni par Huber##
    
    if((u[i] > lambda/(1+lambda)) | (N1>0)){
      
      Y[k,l] = 0
      
    }else if(N01 == 0){
      
      Y[k,l] = 1
      
    }else{
      
      Y[k,l] = 2
      
    }
  }
  
  return(Y)
  
}


Bounding_chain_cftp_HCGM <- function(t,lambda,Nrows,Ncols){##Effectue la méthode de bounding chains pour l'obtention d'un échantillon de la loi stationnaire pour le modèle hardcore gas de paramètre lambda et une mise à jour de type Gibbs pour un graphe de taille Nrows*Ncols##
  
  ##tirage des aléatoires##
  u = runif(t)
  
  ##état initial de la bounding chain##
  Y = matrix(2,Nrows,Ncols) 
  
  ##on met à jour la bounding chain selon les choix aléatoires u##
  for(i in 1:t){
    Y = Bounding_chain_update_HCGM(Y,u[i],lambda)
  }
  
  ##si la bounding chain a plus d'un état sur au moins un des noeuds##
  
  if(sum(Y==2)!=0){
    
    
    ##après récursions, X sera la bounding chain où tout noeud n'aura qu'un élément##
    X = Bounding_chain_cftp_HCGM(2*t,lambda,Nrows,Ncols)
    
    ##on représente la bounding chain par une matrice remplie de 0,1 et 2, où 2 représente le label {0,1}$##
    Y = X
    
    for (i in 1:t){
      
      Y = Bounding_chain_update_HCGM(Y,u[i],lambda)
      
    }
  
  }
  return(Y)
  
}



Bounding_chain_update_shift <- function(Y,u,lambda,pshift){##mise à jour d'un noeud pour la bounding chain du modèle hard core gas shift##
  
  
  for(i in 1:length(u)){
    
    ##On génère une réalisation de la Bernouilli de paramètre pshift##
    S = rbinom(1,1,pshift)
    
    ##On tire un noeud uniformément dans le graphe##
    k = sample(1:dim(Y)[1],1)
    l = sample(1:dim(Y)[2],1)
    
    ##On compte le nombre de voisins (du noeud de coordonnées k,l) de label {1} et de label {0,1}##
    
    N1 = Bounding_chain_voisins_1(Y,k,l)
    
    N01 = Bounding_chain_voisins_01(Y,k,l)
    
    ##On vérifie le tableau des cas donné dans le chapitre d'Huber associé##
    
    if(u[i]>lambda/(1 + lambda)){
      
      Y[k,l] = 0
      
    }else if((N1 == 0) && (N01 == 0) && (u[i] < lambda/(1+lambda)) ){
      
      Y[k,l] = 1
      
    }else if((N1 == 0) && (N01 == 1) && (u[i]<lambda/(1+lambda)) && (S == 0)){
      
      Y[k,l] = 2

    }else if ((N1 == 1) && (S==0) && u[i]<lambda/(1+lambda)){
      
      Y[k,l] = 0
      
    }else if(N1>=2 && u[i]<lambda/(1+lambda)){
      
      Y[k,l] = 0
      
    }else if(N1 == 0 && N01 == 1 && u[i]<lambda/(1+lambda) && S==1){
      
      k = k+1
      l = l+1
      
      X = Y
      
      X = rbind(rep(-99,dim(X)[1]),X)
      X = rbind(X,rep(-99,dim(X)[2]))
      X = cbind(rep(-99,dim(X)[1]),X)  
      X = cbind(X,rep(-99,dim(X)[1]))
      
      w = c(k-1,l)*(X[k-1,l] == 2) + c(k,l+1)*(X[k,l+1] == 2) + c(k+1,l)*(X[k+1,l] == 2) + c(k,l-1)*(X[k,l-1] == 2)
      
      w = w - c(1,1)
      k=k-1
      l=l-1
      Y[w] = 0
      Y[k,l] = 1
      
    }else if((N1==1) && (N01 == 0) && (u[i]<lambda/(1+lambda)) && (S == 1) ){
      
      k = k+1
      l = l+1
      
      X = Y
      
      X = rbind(rep(-99,dim(X)[1]),X)
      X = rbind(X,rep(-99,dim(X)[2]))
      X = cbind(rep(-99,dim(X)[1]),X)  
      X = cbind(X,rep(-99,dim(X)[1]))
      
      w = c(k-1,l)*(X[k-1,l] == 1) + c(k,l+1)*(X[k,l+1] == 1) + c(k+1,l)*(X[k+1,l] == 1) + c(k,l-1)*(X[k,l-1] == 1)
      
      w = w - c(1,1)
      k=k-1
      l=l-1
      Y[w] = 0
      Y[k,l] = 1
      
      
    }else if((N1==1) && (N01 >= 1) && (u[i]<lambda/(1+lambda)) && (S == 1)){
      
      k = k+1
      l = l+1
      
      X = Y
      
      X = rbind(rep(-99,dim(X)[1]),X)
      X = rbind(X,rep(-99,dim(X)[2]))
      X = cbind(rep(-99,dim(X)[1]),X)  
      X = cbind(X,rep(-99,dim(X)[1]))
      
      w = c(k-1,l)*(X[k-1,l] == 1) + c(k,l+1)*(X[k,l+1] == 1) + c(k+1,l)*(X[k+1,l] == 1) + c(k,l-1)*(X[k,l-1] == 1)
      
      w = w - c(1,1)
      k=k-1
      l=l-1
      Y[w] = 2
      Y[k,l] = 2
      
    }
    
    
  }
  
  return(Y)
  
}

Bounding_chain_cftp_shift <- function(t,lambda,pshift,Nrows,Ncols){##Effectue la méthode de bounding chains pour l'obtention d'un échantillon de la loi stationnaire pour le modèle hardcore gas shift de paramètre lambda et une mise à jour de type Gibbs pour un graphe de taille Nrows*Ncols##
  
  ##tirage des aléatoires##
  u = runif(t)
  
  ##état initial de la bounding chain##
  Y = matrix(2,Nrows,Ncols) 
  
  ##on met à jour la bounding chain selon les choix aléatoires u##
  for(i in 1:t){
    Y = Bounding_chain_update_shift(Y,u[i],lambda,pshift)
  }
  
  ##si la bounding chain a plus d'un état sur au moins un des noeuds##
  
  if(sum(Y==2)!=0){
    
    
    ##après récursions, X sera la bounding chain où tout noeud n'aura qu'un élément##
    X = Bounding_chain_cftp_shift(2*t,lambda,pshift,Nrows,Ncols)
    
    ##on représente la bounding chain par une matrice remplie de 0,1 et 2, où 2 représente le label {0,1}$##
    Y = X
    
    for (i in 1:t){
      
      Y = Bounding_chain_update_shift(Y,u[i],lambda,pshift)
      
    }
    
  }
  return(Y)
  
}



Bounding_chain_cftp_HCGM(1,0.5,5,5)
Bounding_chain_cftp_shift(1,0.5,0.5,5,5)


