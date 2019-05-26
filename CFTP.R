simple_update <- function(x,u,t){
  
  for (i in 1:t){
      
    if (u[i]<=0.5 && x>0){
      x = x-1
    }
    else if(u[i]>0.5 && x<2){
      x=x+1
    }
    
  }
  
  return(x)

}

Coupling_from_the_past <- function(i=0){
  
  u = runif(2)
  
  if (u[1]<0.5 && u[2]<0.5){
    
    return(i+2)
    
  }
  else{
    
    i=i+2
    x = Coupling_from_the_past(i)
    return(simple_update(x,u,2))
    
  }
  
}

t=NULL
for (i in 1:10^5){
  
  t=c(t, Coupling_from_the_past())
  
}



