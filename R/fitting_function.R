#' A robust presence-only model with Elastic Net penalty
#'
#' Fit a logistic regression with presence-only response via penalized maximum likelihood.
#' The regularization path is computed for the elastic-net penalty at a pair values of lambda and the prevalence parameter.
#'@usage PO.EN(x,y,o.iter=5, i.iter=5, lambda=.01,alpha=.5,
#'true.prob=0.5,beta_start,epsilon=1e-4, gram.input=F,XtX.input=0,
#'ytx.input=0,XtX_reduce.input)
#'@param x Input design matrix. Should not include the intercept vector.
#'@param y Response variable. Should be a binary vector, such that \eqn{0} represents background observations and \eqn{1} represents presence observations.
#'@param o.iter Number of outer loop iteration.
#'@param i.iter Number of inner loop iteration.
#'@param lambda A user supplied Elastic Net penalty parameter.
#'@param alpha The elastic net mixing parameter, where \eqn{0\le}\code{alpha}\eqn{\le 1}.
#'@param true.prob The  prevalence parameter, should be provided by users. Can be tuned in the cross-validation function.
#'@param epsilon The threshold for stopping the coordinate descent algorithm.
#'@param gram.input The function allows users to feed the gram matrix for fasting computation. The default setting is False, and the function
#'compute the gram matrix for computation.
#'@return \tabular{ll}{
#'\tab \cr
#'\code{beta} \tab The fitting vector of the coefficients, the intercept included.\cr
#'}
#'@details The function fits a presence-only model with an elastic net penalty.
#'@examples
#'data(example.data) # example datasets, including training dataset and testing dataset
#'train_data<-example.data$train.data
#'y_train=train_data$response;x_train=train_data[,-1]  # response and design matrix of training data
#'test_data<-example.data$test.data
#'y_test=test_data$response;x_test=test_data[,-1]  # response and design matrix of testing data
#'PO.EN.cv<-cv.PO.EN(x_train,y_train,input.pi=seq(0.01,0.4,length.out=10))
#'PO.EN.beta<-PO.EN(x_train,y_train,lambda=PO.EN.cv$lambda.min,
#'            true.prob=PO.EN.cv$pi,beta_start=rep(0,ncol(x_train)+1))
#'predictions<-PO.EN.predict(x_test,PO.EN.beta)
#'roc(y_test~predictions)
PO.EN<-function(x,y,o.iter=5, i.iter=5, lambda=.01,alpha=.5,
                               true.prob=0.5,beta_start,epsilon=1e-4, gram.input=F,XtX.input=0,ytx.input=0,XtX_reduce.input){
  x<-as.matrix(cbind(1,x))
  p<-dim(x)[2]
  N=dim(x)[1]
  n.l=sum(y==1)
  n.u=sum(y!=1)
  beta<-beta_start
  beta00<-beta

  if(gram.input==T){
    ytx=ytx.input
    XtX_reduce=XtX_reduce.input
    XtX=XtX.input
    xtx=diag(XtX)
  }else{
    ytx=t(y)%*%x
    XtX=t(x)%*%x
    xtx=diag(XtX)
    XtX_reduce<-c()
    for(o in 1:ncol(XtX)){
      XtX_reduce<-cbind(XtX_reduce,XtX[o,-o])
    }
  }

  ############ generate working response
  XTbeta<-x%*%as.matrix(beta)
  prob<-apply(as.matrix(XTbeta),1,prob.fun)
  working_response=XTbeta+4*(y-prob)


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
  expectation_func<-function(y,prob){
    if(y==1){
      return(1)
    }else{
      return(prob)
    }
  }
  expectation_func<-Vectorize(expectation_func,vectorize.args = c('y','prob'))

  y_working<-expectation_func(y,prob)
  adjust.prob.func<-function(a,correction_input){
    if(is.infinite(exp(a))){
      if(a>0){
        return(1)
      }else{
        return(0)
      }
    }else{
      return(exp(a+correction_input )/(1+exp(a+correction_input)))
    }
  }
  adjust.prob.func<-Vectorize(adjust.prob.func,vectorize.args = c('a'))
  for (s in 1:o.iter){
    XTbeta<-x%*%beta
    beta1=beta
    prob<-apply(as.matrix(XTbeta),1,prob.fun)
    y_working<-expectation_func(y,prob)

      correc<- log((n.l+true.prob*n.u)/(true.prob*n.u))
      prob.corr<-adjust.prob.func(a=XTbeta,correction_input = correc)
      working_response=XTbeta+4*(y_working-prob.corr)


    ytx=t(working_response)%*%x
    a<-sum_compute_single_rcpp(working_response,beta,ytx,XtX_reduce,lambda,alpha,xtx,i.iter)
    beta<-a$beta
    if((sum(abs(beta-c(beta1)))<epsilon)){
      beta.f<-beta*(1+lambda*(1-alpha)/(N))
      return(beta*(1+lambda*(1-alpha)/(N)))
    }

  }
  beta.f<-beta*(1+lambda*(1-alpha)/(N))
  return(beta*(1+lambda*(1-alpha)/(N)))

}
PO.EN<-Vectorize(PO.EN,vectorize.args = 'lambda')
