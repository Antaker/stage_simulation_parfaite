
Bounding_chain_voisins_1 <- function(Y,k,l){ ## Entrée : bounding chain, coordonnées d'un noeud du graphe. Sortie : nombre de voisins de label 1 du noeud choisi##
  
  X = Y
  
  ##On entoure la matrice de valeurs aberrantes##
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
  
  a = 1
  ##tirage des aléatoires##
  u = runif(t)
  
  ##état initial de la bounding chain##
  Y = matrix(2,Nrows,Ncols) 
  
  ##on met à jour la bounding chain selon les choix aléatoires u##

  Y = Bounding_chain_update_HCGM(Y,u,lambda)

  
  ##si la bounding chain a plus d'un état sur au moins un des noeuds##
  
  if(sum(Y==2)!=0){
    
    
    ##après récursions, X sera la bounding chain où tout noeud n'aura qu'un élément##
    X = Bounding_chain_cftp_HCGM(2*t,lambda,Nrows,Ncols)
    a = 2*X[[2]]
    
    ##on représente la bounding chain par une matrice remplie de 0,1 et 2, où 2 représente le label {0,1}$##
    Y = X[[1]]
    
    Y = Bounding_chain_update_HCGM(Y,u,lambda)
  
  }
    return(list(Y,a))
  
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
    
    if(u[i]>=lambda/(1 + lambda)){
      
      Y[k,l] = 0
      
    }else if((N1 == 0) && (N01 == 0)){
      
      Y[k,l] = 1
      
    }else if((N1 == 0) && (N01 == 1) && (S == 0)){
      
      Y[k,l] = 2

    }else if ((N1 == 1) && (S==0)){
      
      Y[k,l] = 0
      
    }else if(N1>=2){
      
      Y[k,l] = 0
      
    }else if((N1 == 0) && (N01 == 1) && (S==1)){
      
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
      Y[w[1],w[2]] = 0
      Y[k,l] = 1
      
    }else if((N1==1) && (N01 == 0) && (S == 1)){
      
      k = k+1
      l = l+1
      
      X = Y
      
      X = rbind(rep(-99,dim(X)[1]),X)
      X = rbind(X,rep(-99,dim(X)[2]))
      X = cbind(rep(-99,dim(X)[1]),X)  
      X = cbind(X,rep(-99,dim(X)[1]))
      
      w = c(k-1,l)*(X[k-1,l] == 1) + c(k,l+1)*(X[k,l+1] == 1) + c(k+1,l)*(X[k+1,l] == 1) + c(k,l-1)*(X[k,l-1] == 1)
      
      w = w - c(1,1)
      k = k-1
      l = l-1
      Y[w[1],w[2]] = 0
      Y[k,l] = 1
      
      
    }else if((N1==1) && (N01 >= 1) && (S == 1)){
      
      k = k+1
      l = l+1
      
      X = Y
      
      X = rbind(rep(-99,dim(X)[1]),X)
      X = rbind(X,rep(-99,dim(X)[2]))
      X = cbind(rep(-99,dim(X)[1]),X)  
      X = cbind(X,rep(-99,dim(X)[1]))
      
      w = c(k-1,l)*(X[k-1,l] == 1) + c(k,l+1)*(X[k,l+1] == 1) + c(k+1,l)*(X[k+1,l] == 1) + c(k,l-1)*(X[k,l-1] == 1)
      
      w = w - c(1,1)
      k =k-1
      l =l-1
      Y[w[1],w[2]] = 2
      Y[k,l] = 2
      
    }
    
  }
  
  return(Y)
  
}

Bounding_chain_cftp_shift <- function(t,lambda,pshift,Nrows,Ncols){##Effectue la méthode de bounding chains pour l'obtention d'un échantillon de la loi stationnaire pour le modèle hardcore gas shift de paramètre lambda et une mise à jour de type Gibbs pour un graphe de taille Nrows*Ncols##
  a = 1
  ##tirage des aléatoires##
  u = runif(t)
  
  ##état initial de la bounding chain##
  Y = matrix(2,Nrows,Ncols) 
  
  ##on met à jour la bounding chain selon les choix aléatoires u##
  
  Y = Bounding_chain_update_shift(Y,u,lambda,pshift)
  
  
  ##si la bounding chain a plus d'un état sur au moins un des noeuds##
  
  if(sum(Y==2)!=0){
    
    
    ##après récursions, X sera la bounding chain où tout noeud n'aura qu'un élément##
    X = Bounding_chain_cftp_shift(2*t,lambda,pshift,Nrows,Ncols)
    a = 2*X[[2]]
    ##on représente la bounding chain par une matrice remplie de 0,1 et 2, où 2 représente le label {0,1}$##
    Y = X[[1]]
    
    Y = Bounding_chain_update_shift(Y,u,lambda,pshift)
      
  }
  return(list(Y,a))
  
}

##nombre d'uniformes utilisées = 2*a - 1#
##complexité = nombre de calculs fait pour chaque uniforme roulée##

t = c()
for (i in 1:10^4){t = c(t,Bounding_chain_cftp_shift(1,0.5,0.5,10,10)[[2]])
print(i)}
mean(t)
Bounding_chain_cftp_HCGM(1,0.5,10,10)[[2]]
X = Bounding_chain_cftp_shift(1,0.5,0.5,10,10)[[1]]
X
image(X)


Bounding_chain_voisins_1_bipartite <- function(Y,A,k){ ## Entrée : bounding chain, coordonnées d'un noeud du graphe. Sortie : nombre de voisins de label 1 du noeud choisi##
  
    ##On récupère des informations sur l'état actuel## 
    X = Y
    V1 = X[[1]]
    V2 = X[[2]]
    
    ##On initialise le compteur d'états 1##
    N1 = 0
    
    ##numvoisin contiendra les numéros des noeuds voisins de k##
    numvoisin = (A[k,] == 1)
    a = 1:dim(A)[2]
    numvoisin = a[numvoisin]
    
    ##On compte le nombre de 1 autour du noeud considéré selon qu'il soit dans une partie du graphe ou l'autre##
    
    if(k <= dim(V1)[2]){
      
      for (i in numvoisin){
        
        N1 = N1 + 1*(V2[i - length(V1)]==1)
        
      }
    }else{
      
      for (i in numvoisin){
        
        N1 = N1 + 1*(V1[i]==1)
        
      }
    }
    
    return(N1)
  
}



Bounding_chain_voisins_01_bipartite <- function(Y,A,k){ ## Entrée : bounding chain, coordonnées d'un noeud du graphe. Sortie : nombre de voisins de label 1 du noeud choisi##
  
    ##On récupère des informations sur l'état actuel##  
    X = Y
    V1 = X[[1]]
    V2 = X[[2]]
    
    ##On initialise le compteur d'états inconnus##
    N01 = 0
    
    ##numvoisin contiendra les numéros des noeuds voisins de k##
    numvoisin = (A[k,] == 1)
    a = 1:dim(A)[2]
    numvoisin = a[numvoisin]
    
    ##On compte le nombre de 1 autour du noeud considéré selon qu'il soit dans une partie du graphe ou l'autre##
    
    if(k <= dim(V1)[2]){
      
      for (i in numvoisin){
        
        N01 = N01 + 1*(V2[i - length(V1)]==2)
        
      }
    }else{
      
      for (i in numvoisin){
        
        N01 = N01 + 1*(V1[i]==2)
        
      }
    }
    
    return(N01)
    
}
  









Bounding_chain_update_HCGM_bipartite <- function(Y,A,u,lambda){
  
  for(i in 1:length(u)){
    
    ##On tire un noeud uniformément dans le graphe##
    ij = sample(1:dim(A)[1],1,replace = TRUE)
    
    ##On compte le nombre de voisins de label {1} et de label {0,1} = 2##
    
    N1 = Bounding_chain_voisins_1_bipartite(Y,A,ij)
    
    N01 = Bounding_chain_voisins_01_bipartite(Y,A,ij)
    
    ##mise à jour selon le pseudo code fourni par Huber##
    
    if((u[i] > lambda/(1+lambda)) | (N1>0)){
      
      if(ij <= length(Y[[1]])){Y[[1]][ij] = 0}else{Y[[2]][ij - length(Y[[1]])] = 0}
      
    }else if(N01 == 0){
      
      if(ij <= length(Y[[1]])){Y[[1]][ij] = 1}else{Y[[2]][ij - length(Y[[1]])] = 1}
      
    }else{
      
      if(ij <= length(Y[[1]])){Y[[1]][ij] = 2}else{Y[[2]][ij - length(Y[[1]])] = 2}
      
    }
  }
  
  return(Y)
  
}



Bounding_chain_HCGM_Bipartite <- function(t,A,lV1,lV2,lambda){ ##t = nombre de pas effectués, A = matrice d'adjacence du graphe biparti, V1 = ensemble des noeuds d'un côté du graphe biparti, V2 = ensemble des noeuds de l'autre côté du graphe biparti
  
  #compteur d'uniformes utilisées##
  a = 1
  
  ##Initialisation du graphe##
  V1 = matrix(2,1,lV1)
  V2 = matrix(2,1,lV2)
  
  ##tirage des aléatoires##
  u = runif(t)
  
  ##on représente le graphe par une liste de 2 vecteurs##
  Y = list(V1,V2)
  
  ##On met à jour le graphe selon les choix aléatoires##
  Y = Bounding_chain_update_HCGM_bipartite(Y,A,u,lambda)
  
  
  if (sum(Y[[1]]==2)>=1 | sum(Y[[2]]==2)>=1){ ##s'il existe encore au moins un noeud de label inconnu##
    
    ##la récursion s'opère##
    X = Bounding_chain_HCGM_Bipartite(2*t,A,lV1,lV2,lambda)
    
    ##On met à jour le compteur##
    a = 2*X[[3]]
    
    ##On met à jour la chaîne résultant de la récursion selon les choix aléatoires effectués##
    Y = Bounding_chain_update_HCGM_bipartite(X[1:2],A,u,lambda)
  }
  
  ##On rend le résultat##
  return(c(Y,a))
  
}


A1bcomptage2 = c()
for (i in 1:10^4){A1bcomptage2 = c(A1bcomptage2,Bounding_chain_HCGM_Bipartite(t,A1,lV1,lV2,2)[[3]])
print(i)}
mean(A1bcomptage2)
table(A1bcomptage2)
hist(A1bcomptage2)



Bounding_chain_update_ising <- function(Y,u,mu,beta){
  
  for(i in 1:length(u)){
    
    ##On tire un noeud uniformément dans le graphe##
    k = sample(1:dim(Y)[1],1)
    l = sample(1:dim(Y)[2],1)
    
    ##On compte le nombre de voisins de label {1} et de label {-1,1} = 2##
    
    N1 = Bounding_chain_voisins_1(Y,k,l)
    
    N01 = Bounding_chain_voisins_01(Y,k,l)
    
    
    ##On prépare l'entourage de la matrice par des valeurs aberrantes##
    k = k+1
    l = l+1
    
    X = Y
    
    X = rbind(rep(-99,dim(X)[1]),X)
    X = rbind(X,rep(-99,dim(X)[2]))
    X = cbind(rep(-99,dim(X)[1]),X)  
    X = cbind(X,rep(-99,dim(X)[1]))
    
    
    ##On calcule le degré du noeud choisi en comptant les valeurs aberrantes autour de celui-ci##
    w = 1*(X[k-1,l] == -99) + 1*(X[k,l+1] == -99) + 1*(X[k+1,l] == -99) + 1*(X[k,l-1] == -99)
    
    
    ##On met à jour selon le degré du noeud, le nombre de valeurs aberrantes et le nombre d'états inconnus##
    if(N01 == 0){
      
      if(u[i] < exp(mu + beta*N1)/(exp(mu + beta*N1) + exp(beta*(4 - N1 - w) - mu))){
        
        Y[k-1,l-1] = 1
        
      }else{
        
        Y[k-1,l-1] = -1
        
      }
      
      
    }else if(w == 0){
      
      if(u[i] < exp(beta*N1 + mu)/(exp(beta*N1 + mu) + exp(beta*(4-N1) - mu))){
        
        Y[k-1,l-1] = 1
        
      }else if(u[i] > exp(beta*(N1 + N01) + mu)/(exp(beta*(N1+N01) + mu) + exp(beta*(4-N1-N01) - mu))){
        
        Y[k-1,l-1] = -1
        
      }else{
        
        Y[k-1,l-1] = 2
        
      }
      
      
    }else if (w==1){
      
      if(u[i] < exp(beta*N1 + mu)/(exp(beta*N1 + mu) + exp(beta*(3-N1) - mu))){
        
        Y[k-1,l-1] = 1
        
      }else if(u[i] > exp(beta*(N1 + N01) + mu)/(exp(beta*(N1+N01) + mu) + exp(beta*(3-N1-N01) - mu))){
        
        Y[k-1,l-1] = -1
        
      }else{
        
        Y[k-1,l-1] = 2
        
      }
      
    }else if (w==2){
      
      if(u[i] < exp(beta*N1 + mu)/(exp(beta*N1 + mu) + exp(beta*(2-N1) - mu))){
        
        Y[k-1,l-1] = 1
        
      }else if(u[i] > exp(beta*(N1 + N01) + mu)/(exp(beta*(N1+N01) + mu) + exp(beta*(2-N1-N01) - mu))){
        
        Y[k-1,l-1] = -1
        
      }else{
        
        Y[k-1,l-1] = 2
        
      }
      
    }
  }
    
  return(Y)
}









Bounding_chain_cftp_ising <- function(t,mu,beta,Nrows,Ncols){##Effectue la méthode de bounding chains pour l'obtention d'un échantillon de la loi stationnaire pour le modèle d'Ising de paramètre mu et beta et une mise à jour de type Gibbs pour un graphe rectangulaire de taille Nrows*Ncols##
  
  a = 1
  
  ##tirage des aléatoires##
  u = runif(t)
  
  ##état initial de la bounding chain##
  Y = matrix(2,Nrows,Ncols) 
  
  ##on met à jour la bounding chain selon les choix aléatoires u##
  
  Y = Bounding_chain_update_ising(Y,u,mu,beta)
  
  
  ##si la bounding chain a plus d'un état sur au moins un des noeuds##
  
  if(sum(Y==2)!=0){
    print(Y)
    image(Y)
    ##après récursions, X sera la bounding chain où tout noeud n'aura qu'un élément##
    X = Bounding_chain_cftp_ising(2*t,mu,beta,Nrows,Ncols)
    a = 2*X[[2]]
    ##on représente la bounding chain par une matrice remplie de -1,1 et 2, où 2 représente le label {-1,1}##
    Y = X[[1]]
    
    Y = Bounding_chain_update_ising(Y,u,mu,beta)
    
  }
  return(list(Y,a))
  
}



