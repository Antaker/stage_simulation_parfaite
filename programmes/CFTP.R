simple_update <- function(x,u,t){##fonction d'update présentée comme exemple par Huber##
  
  for (i in 1:t){
      
    if (u[i]<=0.5 && x>0){
      x = x-1
    }
    else if(u[i]>0.5 && x<2){
      x = x+1
    }
    
  }
  
  return(x)

}

Coupling_from_the_past_counting <- function(i=2){##effectue le CFTP pour la fonction d'update donnée plus haut, par la méthode donnée par Huber##
  
  
  ##Cet algorithme rend le nombre d'appels à la fonction d'update nécessaires à l'obtention d'une réalisation de CFTP##
  
  u = runif(2)
  
  if (u[1]<0.5 && u[2]<0.5){
    
    return(i)
    
  }
  else{
    
    i = Coupling_from_the_past_counting(i+2)
    return(i)
  }
  
}


Coupling_from_the_past_Huber <- function(){  ##effectue le CFTP pour la fonction d'update donnée plus haut, par la méthode donnée par Huber##
  
  ##Cet algorithme rend une réalisation de CFTP##
  
  u = runif(2)
  
  if (u[1]<0.5 && u[2]<0.5){
    
    return(0)
    
  }
  else{
    
    x = Coupling_from_the_past_Huber()
    return(simple_update(x,u,2))
  }
  
}

##On calcule le nombre moyen de choix aléatoires utilisés (ou on rend un tableau de sorties de l'algorithme)
t=NULL
for (i in 1:(10**6)){
  
  t=c(t, Coupling_from_the_past_counting())
  
}

hist(t,freq = F)


Coupling_from_the_past_Bacc <- function(){  ##effectue le CFTP pour la fonction d'update donnée plus haut, par la méthode donnée dans Baccelli et Bremaud##
  
  ##changer le return pour obtenir une réalisation ou le nombre de pas utilisé (ie, a ou dim(u)[2])##
  
  
  u = NULL  ##on initalise les choix aléatoires##
  
  
  while (TRUE) {
    
    u = cbind(t(t(runif(3))),u) ##on ajoute de nouveaux choix alétoires##
    
    ## on réinitialise les états de départ ##
    a=0
    b=1
    c=2
    
    
    ##on met à jour les états de départ selon les choix aléatoires dans u et la fonction d'update##
    for (i in 1:dim(u)[2]){
     
      a = simple_update(a,u[a+1,i],1)
      b = simple_update(b,u[b+1,i],1)
      c = simple_update(c,u[c+1,i],1)
      
    }
   
    ##si les trajectoires se sont rejointes en un même point, on rend ce point (ou le nombre de choix aléatoires pour l'obtenir##
    if(a==b && b==c){
      return(dim(u)[2])
    }
    
  }
  
}


##On calcule le nombre moyen de choix aléatoires utilisés (ou on rend un tableau de sorties de l'algorithme)
t=NULL
for (i in 1:10^4){
  
  t=c(t, Coupling_from_the_past_Bacc())
  
}
hist(t,prob = T)




IsingGibbsUpdate <- function(X,u,v,beta){       ## X = graphe (rectangulaire ou carré), u = lancé uniforme sur [0,1], v = noeud de X ##
  
  ##On entoure la matrice de valeurs aberrantes pour la manipuler plus facilement##
  ##On change les indices pour s'y conformer##
  
  i = v[1] + 1
  j = v[2] + 1
  
  X = rbind(rep(0,dim(X)[2]),X)
  X = rbind(X,rep(0,dim(X)[2]))
  X = cbind(rep(0,dim(X)[1]),X)  
  X = cbind(X,rep(0,dim(X)[1]))
  
  ##On applique la fonction de mise à jour selon le sampler de Gibbs, donnée par Huber##
  if(u< exp(beta*((X[i+1,j]==1)*1 + (X[i-1,j]==1)*1 +(X[i,j+1]==1)*1 + (X[i,j-1]==1)*1))/((exp(beta*((X[i+1,j]==1)*1 + (X[i-1,j]==1)*1 +(X[i,j+1]==1)*1 + (X[i,j-1]==1)*1))) + exp(beta*((X[i+1,j]==-1)*1 + (X[i-1,j]==-1)*1 +(X[i,j+1]==-1)*1 + (X[i,j-1]==-1)*1)))){
    X[i,j] =  1
  }
  else{
    X[i,j] = -1
  }
  
  ##On retire les valeurs aberrantes puis on rend la matrice##
  X=X[,-1]
  X=X[,-dim(X)[2]]
  X=X[-1,]
  X=X[-dim(X)[1],]
  
  return(X)
} 



Monotonic_Ising_Gibbs <- function(t,a,b,bet){ ##t = nombre de pas effectués, a = dim(X)[1], b = dim(X)[2]
  
  
  ##On tire t choix aléatoires##
  u = runif(t)
  
  ##On tire t noeuds du graphe rectangulaire##
  i = sample(1:a,t,replace = TRUE)
  j = sample(1:b,t,replace = TRUE)
  
  ##On initialise les états extremaux##
  xmin = matrix(-1,a,b)
  xmax = list(matrix(1,a,b),t)
  
  ##On fait évoluer les états extrémaux selon les mêmes choix aléatoires##
  for (k in 1:t){
    xmax[[1]] = IsingGibbsUpdate(xmax[[1]],u[k],c(i[k],j[k]),bet)
    xmin = IsingGibbsUpdate(xmin,u[k],c(i[k],j[k]),bet)
  }
  ##si les fins de trajectoires résultant des mises à jour des états extrémaux sont différentes, on fait une récursion##
  if (sum(xmax[[1]] != xmin)>=1){     
    xmax = Monotonic_Ising_Gibbs(2*t,a,b,bet)
    
    ##après récursion, xmax est un état tel que les états extrémaux se sont rejoints, on le met alors à jour selon les choix aléatoires effectués pré-récursion##
    for (m in 1:t){
    
      xmax[[1]] = IsingGibbsUpdate(xmax[[1]],u[m],c(i[m],j[m]),bet)  
      
    }
    
  }  
  
  return(xmax)
  
}


Isingcomptage8815 = NULL
for(i in 1:25){
  
  X = Monotonic_Ising_Gibbs(1,8,8,1.5)[[2]]  
  Isingcomptage8815 = c(Isingcomptage8815,X)
  print(i)
}

table(Isingcomptage8815)
hist(Isingcomptage8815)

voisins_1 <- function(Y,A,k){ ## Entrée : état du modèle HCGM,matrice d'adjacence,numéro du noeud du graphe. Sortie : nombre de voisins de label 1 du noeud choisi##
  
  X = Y
  V1 = X[[1]]
  V2 = X[[2]]
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

HCGM_Gibbs_update <- function(X,A,lambda,u){ ##met à jour un état X de la chaine pour le modèle HCGM biparti de matrice d'adjacence A
  
  for(i in 1: length(u)){
    ##tirage d'un noeud aléatoire uniformément##
    ij = sample(1:dim(A)[1],1,replace = TRUE)
    
    ##On compte les voisins de label 1##
    N1 = voisins_1(X,A,ij)
    
    ##On effectue la mise à jour pour le modèle HCGM donnée par Huber##
    if((u[i] < lambda/(lambda +1)) && (N1==0)){
      if(ij <= length(X[[1]])){X[[1]][ij] = 1}else{X[[2]][ij - length(X[[1]])] = 1}
    }else{if(ij <= length(X[[1]])){X[[1]][ij] = 0}else{X[[2]][ij - length(X[[1]])] = 0}}
  }
  return(X)
}

Monotonic_HCGM_Bipartite_Gibbs <- function(t,A,lV1,lV2,lambda){ ##t = nombre de pas effectués, A = matrice d'adjacence du graphe biparti, V1 = ensemble des noeuds d'un côté du graphe biparti, V2 = ensemble des noeuds de l'autre côté du graphe biparti
  
  ##On effectue les choix aléatoires##
  u = runif(t)
  
  ##On initialise les états extremaux##
  V1 = matrix(1,1,lV1)
  V2 = matrix(0,1,lV2)
  V3 = matrix(0,1,lV1)
  V4 = matrix(1,1,lV2)

  xmin = list(V1,V2)
  xmax = list(V3,V4)
  
  ##On fait évoluer les trajectoires selon les choix aléatoires##
    xmin = HCGM_Gibbs_update(xmin,A,lambda,u)
    xmax = HCGM_Gibbs_update(xmax,A,lambda,u)
    
  ##si les fins de trajectoires résultant des mises à jour des états extrémaux sont différentes, on fait une récursion##
  if ((sum(xmin[[1]]!=xmax[[1]])>=1 | sum(xmin[[2]]!=xmax[[2]])>=1) && t<8192){  
    
    xmax = Monotonic_HCGM_Bipartite_Gibbs(2*t,A,lV1,lV2,lambda)
      
    ##après récursion, on met à jour la sortie selon les choix aléatoires effectués pré-récursion##
    xmax[1:2] = HCGM_Gibbs_update(xmax,A,lambda,u)  
      
  
  }
  return(c(xmax,t,(sum(xmin[[1]]!=xmax[[1]])>=1 | sum(xmin[[2]]!=xmax[[2]])>=1)))
  
}

##On compte le nombre de choix aléatoires nécessaires à la sortie d'un algorithme##

A1comptage2 = NULL
for(i in 1:10^3){
  
  X = Monotonic_HCGM_Bipartite_Gibbs(t,A1,lV1,lV2,2)  
  if(X[[4]] == TRUE){A1comptage2 = c(A1comptage2,8192)}
  else{A1comptage2= c(A1comptage2,X[[3]]) }
  print(i)
}

table(A1comptage2)
hist(A1comptage2)
