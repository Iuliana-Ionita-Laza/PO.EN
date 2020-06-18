
prob.fun<-function(a){
  if(is.infinite(exp(a))){
    if(a>0){
      return(1)
    }else{
      return(0)
    }
  }else{
    return(exp(a)/(1+exp(a)))
  }
}
