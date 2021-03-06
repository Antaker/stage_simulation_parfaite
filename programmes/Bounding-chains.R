
Bounding_chain_voisins_1 <- function(Y,k,l){ ## Entr�e : bounding chain, coordonn�es d'un noeud du graphe. Sortie : nombre de voisins de label 1 du noeud choisi##
  
  X = Y
  
  ##On entoure la matrice de valeurs aberrantes##
  X = rbind(rep(-99,dim(X)[2]),X)
  X = rbind(X,rep(-99,dim(X)[2]))
  X = cbind(rep(-99,dim(X)[1]),X)  
  X = cbind(X,rep(-99,dim(X)[1]))
  
  ##On change les coordonn�es du noeud d'entr�e en cons�quence##
  k = k+1
  l = l+1
  
  ##On compte le nombre de 1 autour du noeud consid�r�##
  
  N1 = 1*(X[k-1,l]==1) + 1*(X[k,l+1]==1) + 1*(X[k+1,l]==1) + 1*(X[k,l-1]==1)
  
  return(N1)
  
}

Bounding_chain_voisins_01 <- function(Y,k,l){ ## Entr�e : bounding chain, coordonn�es d'un noeud du graphe. Sortie : nombre de voisins de label 1 du noeud choisi##
  
  
  X = Y

  ##On entoure les matrices de valeurs aberrantes##
  X = rbind(rep(-99,dim(X)[1]),X)
  X = rbind(X,rep(-99,dim(X)[2]))
  X = cbind(rep(-99,dim(X)[1]),X)  
  X = cbind(X,rep(-99,dim(X)[1]))
  
  ##On change les coordonn�es du noeud d'entr�e en cons�quence##
  k = k+1
  l = l+1
  
  ##On compte le nombre de 2 autour du noeud consid�r�##
  
  N01 = 1*(X[k-1,l]==2) + 1*(X[k,l+1]==2) + 1*(X[k+1,l]==2) + 1*(X[k,l-1]==2)
  
  return(N01)
  
  
}



Bounding_chain_update_HCGM <- function(Y,u,lambda){##mise � jour d'un noeud pour la bounding chain du mod�le hard core gas##
  
  for(i in 1:length(u)){
    
    ##On tire un noeud uniform�ment dans le graphe##
    k = sample(1:dim(Y)[1],1)
    l = sample(1:dim(Y)[2],1)
    
    ##On compte le nombre de voisins de label {1} et de label {0,1} = 2##
    
    N1 = Bounding_chain_voisins_1(Y,k,l)
    
    N01 = Bounding_chain_voisins_01(Y,k,l)
    
    ##mise � jour selon le pseudo code fourni par Huber##
    
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


Bounding_chain_cftp_HCGM <- function(t,lambda,Nrows,Ncols){##Effectue la m�thode de bounding chains pour l'obtention d'un �chantillon de la loi stationnaire pour le mod�le hardcore gas de param�tre lambda et une mise � jour de type Gibbs pour un graphe de taille Nrows*Ncols##
  
  a = 1
  ##tirage des al�atoires##
  u = runif(t)
  
  ##�tat initial de la bounding chain##
  Y = matrix(2,Nrows,Ncols) 
  
  ##on met � jour la bounding chain selon les choix al�atoires u##

  Y = Bounding_chain_update_HCGM(Y,u,lambda)

  
  ##si la bounding chain a plus d'un �tat sur au moins un des noeuds##
  
  if(sum(Y==2)!=0){
    
    
    ##apr�s r�cursions, X sera la bounding chain o� tout noeud n'aura qu'un �l�ment##
    X = Bounding_chain_cftp_HCGM(2*t,lambda,Nrows,Ncols)
    a = 2*X[[2]]
    
    ##on repr�sente la bounding chain par une matrice remplie de 0,1 et 2, o� 2 repr�sente le label {0,1}$##
    Y = X[[1]]
    
    Y = Bounding_chain_update_HCGM(Y,u,lambda)
  
  }
    return(list(Y,a))
  
}



Bounding_chain_update_shift <- function(Y,u,lambda,pshift){##mise � jour d'un noeud pour la bounding chain du mod�le hard core gas shift##
  
  
  for(i in 1:length(u)){
    
    ##On g�n�re une r�alisation de la Bernouilli de param�tre pshift##
    S = rbinom(1,1,pshift)
    
    ##On tire un noeud uniform�ment dans le graphe##
    k = sample(1:dim(Y)[1],1)
    l = sample(1:dim(Y)[2],1)
    
    ##On compte le nombre de voisins (du noeud de coordonn�es k,l) de label {1} et de label {0,1}##
    
    N1 = Bounding_chain_voisins_1(Y,k,l)
    
    N01 = Bounding_chain_voisins_01(Y,k,l)
    
    ##On v�rifie le tableau des cas donn� dans le chapitre d'Huber associ�##
    
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

Bounding_chain_cftp_shift <- function(t,lambda,pshift,Nrows,Ncols){##Effectue la m�thode de bounding chains pour l'obtention d'un �chantillon de la loi stationnaire pour le mod�le hardcore gas shift de param�tre lambda et une mise � jour de type Gibbs pour un graphe de taille Nrows*Ncols##
  a = 1
  ##tirage des al�atoires##
  u = runif(t)
  
  ##�tat initial de la bounding chain##
  Y = matrix(2,Nrows,Ncols) 
  
  ##on met � jour la bounding chain selon les choix al�atoires u##
  
  Y = Bounding_chain_update_shift(Y,u,lambda,pshift)
  
  
  ##si la bounding chain a plus d'un �tat sur au moins un des noeuds##
  
  if(sum(Y==2)!=0){
    
    
    ##apr�s r�cursions, X sera la bounding chain o� tout noeud n'aura qu'un �l�ment##
    X = Bounding_chain_cftp_shift(2*t,lambda,pshift,Nrows,Ncols)
    a = 2*X[[2]]
    ##on repr�sente la bounding chain par une matrice remplie de 0,1 et 2, o� 2 repr�sente le label {0,1}$##
    Y = X[[1]]
    
    Y = Bounding_chain_update_shift(Y,u,lambda,pshift)
      
  }
  return(list(Y,a))
  
}

##nombre d'uniformes utilis�es = 2*a - 1#
##complexit� = nombre de calculs fait pour chaque uniforme roul�e##

t = c()
for (i in 1:10^4){t = c(t,Bounding_chain_cftp_shift(1,0.5,0.5,10,10)[[2]])
print(i)}
mean(t)
Bounding_chain_cftp_HCGM(1,0.5,10,10)[[2]]
X = Bounding_chain_cftp_shift(1,0.5,0.5,10,10)[[1]]
X
image(X)


Bounding_chain_voisins_1_bipartite <- function(Y,A,k){ ## Entr�e : bounding chain, coordonn�es d'un noeud du graphe. Sortie : nombre de voisins de label 1 du noeud choisi##
  
    ##On r�cup�re des informations sur l'�tat actuel## 
    X = Y
    V1 = X[[1]]
    V2 = X[[2]]
    
    ##On initialise le compteur d'�tats 1##
    N1 = 0
    
    ##numvoisin contiendra les num�ros des noeuds voisins de k##
    numvoisin = (A[k,] == 1)
    a = 1:dim(A)[2]
    numvoisin = a[numvoisin]
    
    ##On compte le nombre de 1 autour du noeud consid�r� selon qu'il soit dans une partie du graphe ou l'autre##
    
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



Bounding_chain_voisins_01_bipartite <- function(Y,A,k){ ## Entr�e : bounding chain, coordonn�es d'un noeud du graphe. Sortie : nombre de voisins de label 1 du noeud choisi##
  
    ##On r�cup�re des informations sur l'�tat actuel##  
    X = Y
    V1 = X[[1]]
    V2 = X[[2]]
    
    ##On initialise le compteur d'�tats inconnus##
    N01 = 0
    
    ##numvoisin contiendra les num�ros des noeuds voisins de k##
    numvoisin = (A[k,] == 1)
    a = 1:dim(A)[2]
    numvoisin = a[numvoisin]
    
    ##On compte le nombre de 1 autour du noeud consid�r� selon qu'il soit dans une partie du graphe ou l'autre##
    
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
    
    ##On tire un noeud uniform�ment dans le graphe##
    ij = sample(1:dim(A)[1],1,replace = TRUE)
    
    ##On compte le nombre de voisins de label {1} et de label {0,1} = 2##
    
    N1 = Bounding_chain_voisins_1_bipartite(Y,A,ij)
    
    N01 = Bounding_chain_voisins_01_bipartite(Y,A,ij)
    
    ##mise � jour selon le pseudo code fourni par Huber##
    
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



Bounding_chain_HCGM_Bipartite <- function(t,A,lV1,lV2,lambda){ ##t = nombre de pas effectu�s, A = matrice d'adjacence du graphe biparti, V1 = ensemble des noeuds d'un c�t� du graphe biparti, V2 = ensemble des noeuds de l'autre c�t� du graphe biparti
  
  #compteur d'uniformes utilis�es##
  a = 1
  
  ##Initialisation du graphe##
  V1 = matrix(2,1,lV1)
  V2 = matrix(2,1,lV2)
  
  ##tirage des al�atoires##
  u = runif(t)
  
  ##on repr�sente le graphe par une liste de 2 vecteurs##
  Y = list(V1,V2)
  
  ##On met � jour le graphe selon les choix al�atoires##
  Y = Bounding_chain_update_HCGM_bipartite(Y,A,u,lambda)
  
  
  if (sum(Y[[1]]==2)>=1 | sum(Y[[2]]==2)>=1){ ##s'il existe encore au moins un noeud de label inconnu##
    
    ##la r�cursion s'op�re##
    X = Bounding_chain_HCGM_Bipartite(2*t,A,lV1,lV2,lambda)
    
    ##On met � jour le compteur##
    a = 2*X[[3]]
    
    ##On met � jour la cha�ne r�sultant de la r�cursion selon les choix al�atoires effectu�s##
    Y = Bounding_chain_update_HCGM_bipartite(X[1:2],A,u,lambda)
  }
  
  ##On rend le r�sultat##
  return(c(Y,a))
  
}

##permet l'obtention des histogrammes obtenus pour la comparaison de mod�les##
A1bcomptage2 = c()
for (i in 1:10^4){A1bcomptage2 = c(A1bcomptage2,Bounding_chain_HCGM_Bipartite(t,A1,lV1,lV2,2)[[3]])
print(i)}
mean(A1bcomptage2)
table(A1bcomptage2)
hist(A1bcomptage2)



Bounding_chain_update_ising <- function(Y,u,mu,beta){
  
  for(i in 1:length(u)){
    
    ##On tire un noeud uniform�ment dans le graphe##
    k = sample(1:dim(Y)[1],1)
    l = sample(1:dim(Y)[2],1)
    
    ##On compte le nombre de voisins de label {1} et de label {-1,1} = 2##
    
    N1 = Bounding_chain_voisins_1(Y,k,l)
    
    N01 = Bounding_chain_voisins_01(Y,k,l)
    
    
    ##On pr�pare l'entourage de la matrice par des valeurs aberrantes##
    k = k+1
    l = l+1
    
    X = Y
    
    X = rbind(rep(-99,dim(X)[1]),X)
    X = rbind(X,rep(-99,dim(X)[2]))
    X = cbind(rep(-99,dim(X)[1]),X)  
    X = cbind(X,rep(-99,dim(X)[1]))
    
    
    ##On calcule le degr� du noeud choisi en comptant les valeurs aberrantes autour de celui-ci##
    w = 1*(X[k-1,l] == -99) + 1*(X[k,l+1] == -99) + 1*(X[k+1,l] == -99) + 1*(X[k,l-1] == -99)
    
    
    ##On met � jour selon le degr� du noeud, le nombre de valeurs aberrantes et le nombre d'�tats inconnus##
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

Bounding_chain_cftp_ising <- function(t,mu,beta,Nrows,Ncols){##Effectue la m�thode de bounding chains pour l'obtention d'un �chantillon de la loi stationnaire pour le mod�le d'Ising de param�tre mu et beta et une mise � jour de type Gibbs pour un graphe rectangulaire de taille Nrows*Ncols##
  
  a = 1
  
  ##tirage des al�atoires##
  u = runif(t)
  
  ##�tat initial de la bounding chain##
  Y = matrix(2,Nrows,Ncols) 
  
  ##on met � jour la bounding chain selon les choix al�atoires u##
  
  Y = Bounding_chain_update_ising(Y,u,mu,beta)
  
  
  ##si la bounding chain a plus d'un �tat sur au moins un des noeuds##
  
  if(sum(Y==2)!=0){
    ##apr�s r�cursions, X sera la bounding chain o� tout noeud n'aura qu'un �l�ment##
    X = Bounding_chain_cftp_ising(2*t,mu,beta,Nrows,Ncols)
    a = 2*X[[2]]
    ##on repr�sente la bounding chain par une matrice remplie de -1,1 et 2, o� 2 repr�sente le label {-1,1}##
    Y = X[[1]]
    
    Y = Bounding_chain_update_ising(Y,u,mu,beta)
    
  }
  return(list(Y,a))
  
}

##permet l'obtention des histogrammes obtenus pour la comparaison de m�thodes##
Isingbcomptage881 = NULL
for(i in 1:10^2){
  
  X = *Bounding_chain_cftp_ising(1,0,0.5,8,8)[[1]]  
  Isingbcomptage881 = c(Isingbcomptage881,X)
  print(i)
}

table(Isingbcomptage881)
hist(Isingbcomptage881)



