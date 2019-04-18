IsingGibbsUpdate <- function(X,u,v,beta){       ## X = graphe (rectangulaire ou carré), u = lancé uniforme sur [0,1], v = noeud de X ##
  
      i = v[1]
      j = v[2]
      
      if (v[1] == 1){
        
        
        if (v[2] == 1){
          if (u<exp(beta*((X[i,j+1]==1) + (X[i+1,j]==1)))/(exp(beta*((X[i,j+1]==1) + (X[i+1,j]==1))) + exp(beta*((X[i,j+1]==-1) + (X[i+1,j]==-1))))){
            X[i,j] = 1
          }
          else{
            X[i,j] = -1
          }
         }
        
        
        else if (v[2] == dim(X)[2]){
          if (u<exp(beta*((X[i,j-1]==1) + (X[i+1,j]==1)))/(exp(beta*((X[i,j-1]==1) + (X[i+1,j]==1))) + exp(beta*((X[i,j-1]==-1) + (X[i+1,j]==-1))))){
            X[i,j] = 1
          }
          else{
            X[i,j] = -1
          }
        }
        
        
        else{
          if (u<exp(beta*((X[i,j+1]==1) + (X[i+1,j]==1) + (X[i,j-1]==1)))/(exp(beta*((X[i,j+1]==1) + (X[i+1,j]==1) + (X[i,j-1]==1))) + exp(beta*((X[i,j+1]==-1) + (X[i+1,j]==-1) + (X[i,j-1]==-1))))){
            X[i,j] = 1
          }
          else{
            X[i,j] = -1
          }
        }
        
        
        
        if (v[1] == dim(X)[1]){
          
          
          if(v[2] == 1){
            if (u<exp(beta*((X[i,j+1]==1) + (X[i-1,j]==1)))/(exp(beta*((X[i,j+1]==1) + (X[i-1,j]==1))) + exp(beta*((X[i,j+1]==-1) + (X[i-1,j]==-1))))){
              X[i,j] = 1
            }
            else{
              X[i,j] = -1
            }
          }
          
          
          else if (v[2] == dim(X)[2]){
            if (u<exp(beta*((X[i,j-1]==1) + (X[i-1,j]==1)))/(exp(beta*((X[i,j-1]==1) + (X[i-1,j]==1))) + exp(beta*((X[i,j-1]==-1) + (X[i-1,j]==-1))))){
              X[i,j] = 1
            }
            else{
              X[i,j] = -1
            }
          }
          
          
          
          else{
            if (u<exp(beta*((X[i,j+1]==1) + (X[i-1,j]==1) + (X[i,j-1]==1)))/(exp(beta*((X[i,j+1]==1) + (X[i-1,j]==1) + (X[i,j-1]==1))) + exp(beta*((X[i,j+1]==-1) + (X[i-1,j]==-1) + (X[i,j-1]==-1))))){
              X[i,j] = 1
            }
            else{
              X[i,j] = -1
            }
          }
      
      }else{
        
        
        if (v(2)==1){
          if (u<exp(beta*((X[i,j+1]==1) + (X[i-1,j]==1) + (X[i+1,j]==1)))/(exp(beta*((X[i,j+1]==1) + (X[i-1,j]==1) + (X[i+1,j]==1))) + exp(beta*((X[i,j+1]==-1) + (X[i-1,j]==-1) + (X[i+1,j]==-1))))){
            X[i,j] = 1
          }
          else{
            X[i,j] = -1
          }
        }
        
        
        else if (v(2)==dim(X)[2]){
          if (u<exp(beta*((X[i,j-1]==1) + (X[i-1,j]==1) + (X[i+1,j]==1)))/(exp(beta*((X[i,j-1]==1) + (X[i-1,j]==1) + (X[i+1,j]==1))) + exp(beta*((X[i,j-1]==-1) + (X[i-1,j]==-1) + (X[i+1,j]==-1))))){
            X[i,j] = 1
          }
          else{
            X[i,j] = -1
          }
        }
        
        else{
          if (u<exp(beta*((X[i,j-1]==1) + (X[i-1,j]==1) + (X[i+1,j]==1) + (X[i+1,j+1]==1)))/(exp(beta*((X[i,j-1]==1) + (X[i-1,j]==1) + (X[i+1,j]==1) + (X[i+1,j+1]==1))) + exp(beta*((X[i,j-1]==-1) + (X[i-1,j]==-1) + (X[i+1,j]==-1) + (X[i+1,j+1]==-1))))){
            X[i,j] = 1
          }
          else{
            X[i,j] = -1
          }
        }
          
          
      }
        
        
        
    }
  
      return(X)
}   
