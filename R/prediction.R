#' PO-EN predicting function
#'
#' A prediction function using the linear predictor of PO-EN fitting results.
#'@param X Input design matrix. Should not include the intercept vector.
#'@param beta A coefficients vector from the PO-EN fitting function.
#'@export
#'@examples
#'PO.EN.cv<-cv.PO.EN(x_train,y_train,input.pi=seq(0.01,0.4,length.out=10))
#'PO.EN.beta<-PO.EN(x_train,y_train,lambda=PO.EN.cv$lambda.min,
#'            true.prob=PO.EN.cv$pi,beta_start=rep(0,ncol(x_train)+1))
#'predictions<-PO.EN.predict(x_test,PO.EN.beta)
#'roc(y_test~predictions)
PO.EN.predict<-function(X,beta){
  X<-as.matrix(X)
  X<-cbind(1,X)
  mu<-X%*%beta
  Y.predict=exp(mu)/(1+exp(mu))
  Y.predict[which(is.nan(Y.predict))]=1

  return(Y.predict)
}
