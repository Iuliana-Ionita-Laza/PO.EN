library(PO.EN)
library(PRROC)
library(mvtnorm)
library(data.table)
library(glmnet)
library(pROC)

rm(list=ls()); 
f.fun<-function(classier,true.y,a=1){
  retrieved <- sum(classier)
  
  precision <- sum(classier & true.y) / retrieved
  
  precision <- ifelse(is.nan(precision),0,precision)
  
  recall <- sum(classier & true.y) / sum(true.y)
  
  Fmeasure <- (1+a^2) * precision * recall / (a^2*precision + recall)
  if(sum(precision + recall)==0){
    Fmeasure=0
  }
  
  stand.pre<-sum(true.y)/length(true.y)
  stand.recall<-1
  F.stand<-2*stand.pre/(1+stand.pre)
  if(Fmeasure==F.stand){
    Fmeasure=0
  }
  return(Fmeasure)
}
aupr.fun<-function(pre,true){
  s0<-pre[true==1]
  s1<-pre[true==0]
  return(pr.curve(s0,s1))
}

   
    tissue='mutasea_E123';  # Name of training dataset, other datasets are 'mutasea_E086' ,'mutasea_E118','suresea_E118','suresea_E123'
    simu=5;     #How many replicates
    true.pre=seq(0.01,0.45,length.out = 10);# The candidate for pi
    cv.fold=10       # number of folds of c.v.
    type.measure='auc' # Choice of measurement for c.v. other option is 'deviance'
    imba.ratio<-5      # Imbalance ratio for training dataset
    posi.length=10;   # The number of positive responses included
    train.threshold<-1e-5 # Threshold for positive response based on p-values
    seed=1
    alpha=0.5 # Mixing parameter
    set.seed(seed)
  
    
  label.data<-fread(paste0('https://github.com/Iuliana-Ionita-Laza/PO.EN/raw/master/training_data/labeled_data_',tissue,'.csv.gz'))
  #######create imba.dataset
  Response<-as.numeric(label.data[,4]<=train.threshold)
  posi.index<-which(Response==1)
  nega.index<-which(Response==0)
  
  

  
  
  result<-data.frame(PO.EN=rep(NA,simu),EN=rep(NA,simu), FUNLDA=rep(NA,simu),Sig=rep(NA,simu),eQTL=rep(NA,simu),GWAS=rep(NA,simu),HGMD=rep(NA,simu))
  
  
  for(i in 1:simu){
    # Randomly select training and testing datasets from the dataset
    imba.index<-c(sample(posi.index,posi.length),sample(nega.index,posi.length*imba.ratio)) 
    train.data<-label.data[imba.index,]            
    y.train<-as.numeric(train.data[,4]<=train.threshold)  
    x<-train.data[,(ncol(label.data)-918):ncol(label.data)]
    m=floor(nrow(x)/2) 
    total.number<-length(y.train)
    
    
    fun.lda<-train.data$`fun-lda`     
    sig_scores<-train.data$func_sig
    eQTL<-train.data$eQTL
    GWAS<-train.data$GWAS
    HGMD<-train.data$HGMD
    
    ####### Formulate training dataset and testing dataset
    pos.num<-floor(posi.length*1/2)
    train.index<-c(sample(which(y.train==1 ) ,pos.num),sample(which(y.train==0 ) ,m-pos.num))
    print(length(train.index))
    test.index<-c(1:length(y.train)[1])[-train.index]
    
    Y.training<-y.train[train.index];
    Y.test<-y.train[test.index]
    
   
    X.training<-(x[train.index,]);X.test<-(x[test.index,])
    
    X.training<-as.matrix(X.training)
    ### Elastic Net 
    EN.cv<-cv.glmnet(X.training,Y.training,alpha=alpha,standardize = F,intercept = T,family='binomial',type.measure ='deviance')
    
    EN.fit<-glmnet(X.training,Y.training,alpha=alpha,standardize = F,intercept = T,family='binomial',lambda=EN.cv$lambda.min)
    EN.beta<-as.matrix(c(EN.fit$a0,as.numeric(EN.fit$beta)))
    EN.pre<-plogis(as.matrix(cbind(1,X.test))%*%EN.beta)
    ### PO-EN model
      
      P.O.lambda<-cv.PO.EN(X.training,Y.training,o.iter=10,i.iter=10,alpha=0.5,
                           epsilon=1e-2,nfolds=10,type.measure=type.measure,input.pi = true.pre)
      
      print(paste0('lambda.min: ',P.O.lambda$lambda.min,'; lambad.1se: ',P.O.lambda$lambda.1se,'; pi: ',P.O.lambda$pi))
    
      presence_beta<-PO.EN(X.training,Y.training,o.iter=10, i.iter=20, lambda=P.O.lambda$lambda.min,alpha=alpha,
                           true.prob = P.O.lambda$pi,beta_start=rep(0,ncol(X.training)+1))
      presence_pre<-plogis(as.matrix(cbind(1,X.test))%*%as.matrix(presence_beta))
      
      ### Functional scores of other models
      funlda.pre<-fun.lda[test.index]
      sig.pre<-sig_scores[test.index]
      eQTL.pre<-eQTL[test.index]
      GWAS.pre<-GWAS[test.index]
      HGMD.pre<-HGMD[test.index]
      ###Compute results
      PO.c<-apply(presence_pre,1,function(x){ifelse(x>0.5,1,0)})
      EN.c<-apply(EN.pre,1,function(x){ifelse(x>0.5,1,0)})
      funlda.c<-apply(as.matrix(funlda.pre),1,function(x){ifelse(x>0.5,1,0)})
      eQTL.c<-as.numeric(eQTL.pre>0.5)
      GWAS.c<-as.numeric(GWAS.pre>0.5)
      HGMD.c<-as.numeric(HGMD.pre>0.5)
      
      fb<-f.fun(PO.c,Y.test)
      fd<-f.fun(EN.c,Y.test)
      ff<-f.fun(funlda.c,Y.test)
      fe<-f.fun(eQTL.c,Y.test)
      fg<-f.fun(GWAS.c,Y.test)
      fh<-f.fun(HGMD.c,Y.test)
      
      result.f<-c(fb,fd,ff,fe,fg,fh)
      
      result[i,1]=as.numeric(roc(Y.test~as.numeric( presence_pre) )$auc)
      result[i,2]=as.numeric(roc(Y.test~as.numeric(EN.pre))$auc)
      result[i,3]=as.numeric(roc(Y.test~as.numeric(funlda.pre))$auc)
      ####deep sea
      result[i,4]=as.numeric(roc(Y.test~as.numeric(sig.pre))$auc)
      result[i,5]=as.numeric(roc(Y.test~as.numeric(eQTL.pre))$auc)
      result[i,6]=as.numeric(roc(Y.test~as.numeric(GWAS.pre))$auc)
      result[i,7]=as.numeric(roc(Y.test~as.numeric(HGMD.pre))$auc)
      
      po.aupr<-aupr.fun(presence_pre,Y.test)
      en.aupr<-aupr.fun(EN.pre,Y.test)
      fun.aupr<-aupr.fun(funlda.pre,Y.test)
     
      ####deepsea
      sig.aupr<-aupr.fun(sig.pre,Y.test)
      eQTL.aupr<-aupr.fun(eQTL.pre,Y.test)
      GWAS.aupr<-aupr.fun(GWAS.pre,Y.test)
      HGMD.aupr<-aupr.fun(HGMD.pre,Y.test)
      
      result.aupr<-c(po.aupr$auc.davis.goadrich,en.aupr$auc.davis.goadrich,fun.aupr$auc.davis.goadrich,
                     sig.aupr$auc.davis.goadrich,eQTL.aupr$auc.davis.goadrich,GWAS.aupr$auc.davis.goadrich,HGMD.aupr$auc.davis.goadrich)
      
      
      
    
      write.table(result[i,],file=paste0('result_table_auroc_',tissue,'_' , as.character(posi.length),'_',type.measure, '_imba_size.txt'),
                  row.names=F,append=T,col.names=F)
      
      write.table(t(result.aupr),file=paste0('result_table_','aupr','_',tissue,'_' , as.character(posi.length),'_',type.measure, '_imba_size.txt'),
                  row.names=F,append=T,col.names=F)
      
      write.table(t(result.f),file=paste0('result_table_','f','_',tissue,'_' , as.character(posi.length),'_',type.measure, '_imba_size.txt'),
                  row.names=F,append=T,col.names=F)
    
    
    
    
  }
  
  

