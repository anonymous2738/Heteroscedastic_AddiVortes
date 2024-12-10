if (!requireNamespace("invgamma", quietly = TRUE)) {
  install.packages("invgamma")
}

if (!requireNamespace("FNN", quietly = TRUE)) {
  install.packages("FNN")
}
if (!requireNamespace("energy", quietly = TRUE)) {
  install.packages("energy")
}

library('invgamma')
library('FNN')
library('energy')
library('scoringRules')

AddiVortes_Algorithm<-function(y,x,m = 200, m_var = 40 ,max_iter = 1200,burn_in= 200,nu = 6,q =0.85,k = 3 ,sd = 0.8 ,Omega = 3,lambda_rate = 25,YTest,XTest,IntialSigma = "Linear",thinning=1){

  #Scaling x and y
  yScaled=(y-(max(y)+min(y))/2)/(max(y)-min(y))
  xScaled=x;
  for (i in 1:length(x[1,])){
    xScaled[,i]=(x[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }
  
  
  for (i in 1:length(XTest[1,])){
    XTest[,i]=(XTest[,i]-(max(x[,i])+min(x[,i]))/2)/(max(x[,i])-min(x[,i]));
  }
  
  #Initialize:
  #Prediction Set (A list of vectors with the output values for each tessellation), 
  #Dimension set (A list of vectors with the covariates included in the tessellaions);
  #and Tessellation Set (A list of matrices that give the coordinates of the centers in the tessellations)
  
  Pred<-rep(list(matrix(mean(yScaled)/m)),m)
  Dim=vector(length = m)
  Tess=vector(length = m)
  for (i in 1:m){
    Dim[i]<-list(sample(1:length(x[1,]), 1))
    Tess[i]<-(list(matrix(rnorm(1,0,sd))))
  }
  
  #Prepare some variables used in the backfitting algorithm
  SumOfAllTess=rep(mean(yScaled),length(yScaled))
  SigmaSquaredMu=(0.5/(k*sqrt(m)))^2
  LastTessPred=matrix
  
  #Matrices that will hold the samples from the poseterior distribution for the training samples and test samples.
  posterior_samples<-floor((max_iter - burn_in) / thinning)
  PredictionMatrix<-array(dim=c(length(y),posterior_samples))
  TestMatrix<-array(dim=c(length(YTest),posterior_samples))
  
  #finding lambda
  if (IntialSigma=="Naive"){ # Usually used if p is greater then n. Uses Standard deviation of y to predict sigma.
    SigmaSquaredHat=var(yScaled)
  }
  else{  # Default method using residual standard deviation from a least-squared linear regression of y, to predict sigma.
    MultiLinear<-lm(yScaled ~ xScaled)
    SigmaSquaredHat=sum(MultiLinear$residuals^2)/(length(yScaled)-length(xScaled[1,])-1)
  }
  
  #Find lambda
  lambda=1;
  lambda <- optim(par = 1,
                  fitting_function,
                  method = "Brent",
                  lower = 0.001,
                  upper = 100,
                  q=q , nu=nu, sigmaSquared_hat=SigmaSquaredHat)$par

  lambda<-lambda^(1/m_var)
  nu<-2/(1-(1-2/nu)^(1/m_var))
  
  SigmaSquared_inital=SigmaSquaredCalculation_initial(yScaled,nu,lambda,SumOfAllTess)
  #SigmaSquared_inital=var(yScaled)
    
  ## variance intialistation ##
  
  Pred_var<-rep(list(matrix((SigmaSquared_inital)^(1/m_var))),m_var)
  Dim_var=vector(length = m)
  Tess_var=vector(length = m)
  for (i in 1:m){
    Dim_var[i]<-list(sample(1:length(x[1,]), 1))
    Tess_var[i]<-(list(matrix(rnorm(1,0,sd))))
  }
  
  Numerator_e<-vector(length = length(x[,1]))
  sigmaSquared<-vector(length = length(xScaled[1,]))
  
  #Prepare some variables used in the backfitting algorithm
  ProdOfAllTess=rep(SigmaSquared_inital,length(yScaled))
  LastTessPred_var=matrix
  Sigma_squared_stored<-matrix(ncol=max_iter-burn_in,nrow= length(yScaled))
  Sigma_squared_test<-array(dim=c(length(YTest),posterior_samples))
  
  for (i in 1:max_iter){
    
    #Sample Sigma squared using all tessellations to predict the outcome variables
    Numerator_e<-(yScaled-SumOfAllTess)^2
    sigmaSquared_model=SigmaSquaredCalculation(xScaled, yScaled,Tess_var,Dim_var,Pred_var, nu,lambda,ProdOfAllTess,Numerator_e,Omega,lambda_rate,m_var,sd)
    sigmaSquared<-sigmaSquared_model[[1]]
    Tess_var<-sigmaSquared_model[[2]]
    Dim_var<-sigmaSquared_model[[3]]
    Pred_var<-sigmaSquared_model[[4]]
    ProdOfAllTess<-sigmaSquared_model[[1]]
    ##sigmaSquared=rep(1/((max(y)-min(y))^2),length(yScaled))
    ##print('start')
    ##print(mean(sigmaSquared*((max(y)-min(y))^2)/1))
    if (i>burn_in ){ 
      Sigma_squared_stored[,i-burn_in]<-sigmaSquared*((max(y)-min(y))^2)
      Sigma_squared_test[,i-burn_in]<-TestPrediction_var(XTest,m_var,Tess_var,Dim_var,Pred_var)*((max(y)-min(y))^2)
    }
    for (j in 1:m){
      NewTessOutput<-NewTess(xScaled,j,Tess,Dim,sd) #Propose new Tessellation 
      TessStar<-NewTessOutput[[1]]  
      DimStar<-NewTessOutput[[2]]
      Modification<-NewTessOutput[[3]]
      
      ResidualsOutput<-CalculateResiduals(yScaled,xScaled,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred,sigmaSquared) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
      R_ijOld<-ResidualsOutput[[1]]   #Old and New refer to the original and proposed tessellations
      SigmaSquared_ijOld<-ResidualsOutput[[2]]
      R_ijNew<-ResidualsOutput[[3]]
      SigmaSquared_ijNew<-ResidualsOutput[[4]]
      SumOfAllTess<-ResidualsOutput[[5]] #Keeps track of the prediction for all tessellations to help sample sigma squared.
      IndexesStar<-ResidualsOutput[[6]] #Gives the row of each observation for the cell it falls in for the proposed tessellation.
      Indexes<-ResidualsOutput[[7]]  #Gives the row of each observation for the cell it falls in for the original tessellation.
      if (!any(SigmaSquared_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.
        
        LOGAcceptenceProb=AlphaCalculation(xScaled,TessStar,DimStar,j,R_ijOld,SigmaSquared_ijOld,R_ijNew,SigmaSquared_ijNew,Modification,SigmaSquaredMu,Omega,lambda_rate) #Gives the log of the acceptence probability.

        if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation. 
          Tess=TessStar
          Dim=DimStar
          Pred[[j]]=NewPredSet(j,TessStar,R_ijNew,SigmaSquaredMu,SigmaSquared_ijNew)
          LastTessPred=Pred[[j]][IndexesStar]
        }
        else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
          Pred[[j]]=NewPredSet(j,Tess,R_ijOld,SigmaSquaredMu,SigmaSquared_ijOld);
          LastTessPred=Pred[[j]][Indexes];
        }
      }
      else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
        Pred[[j]]=NewPredSet(j,Tess,R_ijOld,SigmaSquaredMu,SigmaSquared_ijOld);
        LastTessPred=Pred[[j]][Indexes];
      }
      if (j==m){ #If j equals m then adds the last tessellation output values to give a prediction.
        SumOfAllTess=SumOfAllTess+LastTessPred;
      }
    }
    
    if (i %% 100 == 0){
      cat(sprintf("Iteration %d out of %d", i, max_iter), "\n")
    }
    
    if (i>burn_in & (i-burn_in) %% thinning ==0 ){ #vectors that hold the predictions for each iteration after burn in.
      PredictionMatrix[,(i-burn_in)/thinning]=SumOfAllTess;
      TestMatrix[,(i-burn_in)/thinning]=TestPrediction(XTest,m,Tess,Dim,Pred);
    }
  }
  
  #finding the mean of the predition over the iterations and then unscaling the predictions.
  mean_yhat=(rowSums(PredictionMatrix)/(posterior_samples))*(max(y)-min(y))+((max(y)+min(y))/2)
  mean_yhat_Test=(rowSums(TestMatrix)/(posterior_samples))*(max(y)-min(y))+((max(y)+min(y))/2)
  
  Sigma_squared_sorted<-Sigma_squared_test
  #print(TestMatrix)
  LowerConfidenceTRAINValue<-vector(length=length(mean_yhat_Test))
  UpperConfidenceTRAINValue<-vector(length=length(mean_yhat_Test))

  for (i in 1:length(mean_yhat_Test)){
    Sigma_squared_sorted[i,]<-sort(sqrt(Sigma_squared_test[i,]))
    
    if ((((max_iter-burn_in+1)*0.05))== round((max_iter-burn_in+1)*0.05)){
      LowerConfidenceTRAINValue[i]<-(Sigma_squared_sorted[i,(max_iter-burn_in+1)*0.05])
      UpperConfidenceTRAINValue[i]<-(Sigma_squared_sorted[i,(max_iter-burn_in+1)*0.95])
    }
    else{
      LowerConfidenceTRAINValue[i]<-((Sigma_squared_sorted[i,trunc((max_iter-burn_in+1)*0.05)]+Sigma_squared_sorted[i,trunc(((max_iter-burn_in+1)*0.05)+1)])/2)
      UpperConfidenceTRAINValue[i]<-((Sigma_squared_sorted[i,trunc((max_iter-burn_in+1)*0.95)]+Sigma_squared_sorted[i,trunc(((max_iter-burn_in+1)*0.95)+1)])/2)
    }
  }
  
  LowerConfidenceTESTValue<-vector(length=length(mean_yhat_Test))
  UpperConfidenceTESTValue<-vector(length=length(mean_yhat_Test))
  
  for (i in 1:length(mean_yhat_Test)){
    TestMatrix[i,]<-sort(TestMatrix[i,])
  
    if ((((max_iter-burn_in+1)*0.05))== round((max_iter-burn_in+1)*0.05)){
      LowerConfidenceTESTValue[i]<-(TestMatrix[i,(max_iter-burn_in+1)*0.05])*(max(y)-min(y))+((max(y)+min(y))/2)
      UpperConfidenceTESTValue[i]<-(TestMatrix[i,(max_iter-burn_in+1)*0.95])*(max(y)-min(y))+((max(y)+min(y))/2)
    }
    else{
      LowerConfidenceTESTValue[i]<-((TestMatrix[i,trunc((max_iter-burn_in+1)*0.05)]+TestMatrix[i,trunc(((max_iter-burn_in+1)*0.05)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
      UpperConfidenceTESTValue[i]<-((TestMatrix[i,trunc((max_iter-burn_in+1)*0.95)]+TestMatrix[i,trunc(((max_iter-burn_in+1)*0.95)+1)])/2)*(max(y)-min(y))+((max(y)+min(y))/2)
    }
  }
  
  ##dims
  nd = nrow(t(TestMatrix))
  np = ncol(t(TestMatrix))

  ##draw from predictive
  pdraw = t(TestMatrix*(max(y)-min(y))+((max(y)+min(y))/2)) + t(sqrt(Sigma_squared_test)) * matrix(rnorm(nd*np),nrow=nd)

  ## predictive qqplot
  #yquantile = qsamp(f(X[TestSet,]),pdraw)
  yquantile = qsamp(YTest,pdraw)
  unifdr = runif(1000)
  plot(qqplot(yquantile,unifdr,plot.it = FALSE)$x,qqplot(yquantile,unifdr,plot.it = FALSE)$y,pch=20,cex=1.2)
  #points(qqplot(yquantile,unifdr,plot.it = FALSE)$x,qqplot(yquantile,unifdr,plot.it = FALSE)$y,pch=20,cex=1.2)
  abline(0,1,col='red',lwd=2)
  qvec = qsamp(YTest,pdraw)
  test_result<-edist(matrix(c(qvec,runif(1000)),ncol=1),c(np,1000))
  
 
  #combined_samples <- as.matrix(rbind(as.matrix(yquantile), as.matrix(unifdr)))
  #test_result <- eqdist.e(combined_samples, sizes = c(length(YTest), length(unifdr)), distance = FALSE)
  print(test_result)
  #print(pdraw)
  pdraw<-t(pdraw)
  # Calculate CRPS for each test observation 
  crps_values <- sapply(1:nrow(pdraw), function(i) { 
    crps_sample(y = YTest[i], dat = pdraw[i, ]) 
    }) 
  # Calculate the mean CRPS for the test set 
  mean_crps <- mean(crps_values) 
  # Print result 
  cat("Mean CRPS for test set:", mean_crps, "\n")
  
  return( #Returns the RMSE value for the test samples.
    ##test_result)
    list(
      In_sample_RMSE = sqrt(mean((y-mean_yhat)^2)),
      Out_of_sample_RMSE = sqrt(mean((YTest-mean_yhat_Test)^2)) ,
      Sigma_squared_stored = rowMeans(Sigma_squared_test),
      Sigma_squared_stored_col = colMeans(Sigma_squared_stored),
      LowerConfidenceTRAINValue_stored = LowerConfidenceTRAINValue,
      UpperConfidenceTRAINValue_stored = UpperConfidenceTRAINValue,
      mean_yhat_Test_stored = mean_yhat_Test,
      LowerConfidenceTESTValue_stored= LowerConfidenceTESTValue,
      UpperConfidenceTESTValue_stored =UpperConfidenceTESTValue,
      e_stat = test_result,
      pdraws= pdraw
    )
  )
}

qinsamp = function(y,ysamp) { ###get quantile of y in sample ysamp
  n=length(ysamp)
  return(which.min(abs(y-sort(ysamp)))/n)
}

##helper function
qsamp = function(y,yd) { ###get quantile of yi in ith column of yd
  nd=nrow(yd)
  n=ncol(yd)
  qvec=rep(0,n)
  for(i in 1:n) {
    qvec[i]=qinsamp(y[i],yd[,i])
  }
  return(qvec)
}

calculate_percentile <- function(y_true, y_pred_samples) {
  # Calculate the proportion of predictions less than or equal to the true y
  mean(y_pred_samples <= y_true)
}

NewTess<-function(x,j,Tess,Dim,sd){ #Propose a new tessellation
  
  p=runif(1,0,1) #Randomly sample p to decide the proposed modification to tessellation.
  
  DimStar=Dim # Let proposed dimension matrix equal original dimension matrix.
  TessStar=Tess #Similar for the tessellation matrix.
  
  if (p<0.2 & length(Dim[[j]])!=length(x[1,]) | length(Dim[[j]])==1 & p<0.4){ #Add a dimension if p is less then 0.2 or if p is less then 0.4 when there is only one dimension in the Tessellation due to adjustments (Supplementary Material).
    NumberOfCovariates=1:length(x[1,]) #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
    NumberOfCovariates=NumberOfCovariates[-Dim[[j]]] #Remove all values that for the covariates that are already in the tessellation.
    DimStar[[j]]<-c(Dim[[j]],sample(NumberOfCovariates,1)) # Uniformly sample a new covariate and add it to the dimension matrix.
    TessStar[[j]]=cbind(Tess[[j]],rnorm(length(Tess[[j]][,1]),0,sd)) # Sample new coordinates from Normal distribution for the new dimension and add it to the Tessellation matrix.
    Modification="AD"}
  else if (p<0.4){ #Remove a dimension if p is less then 0.4.
    RemovedDim=sample(1:length(Dim[[j]]),1) #Uniformly sample the dimension to be removed.
    DimStar[[j]]=DimStar[[j]][-RemovedDim] #Remove the dimension from the dimesion Matrix.
    TessStar[[j]]=matrix(TessStar[[j]][,-RemovedDim],ncol=length(DimStar[[j]])) #Remove the coordinates in the Tessellation matrix corresponding to the dimension removed.
    Modification="RD"}
  else if (p<0.6 || p<0.8 & length(Tess[[j]][,1])==1){ #Add a centre if p is less then 0.6 or if p is less then 0.4 when there is only one center in the Tessellation due to adjustments (Supplementary Material).
    TessStar[[j]]=rbind(Tess[[j]],rnorm(length(Dim[[j]]),0,sd)) #Add a new row of coordinates, sampled from a normal distribution, to the Tessellation matrix to add a center.
    Modification="AC"}
  else if (p<0.8){ #Add a centre if p is less then 0.8. 
    CenterRemoved=sample(1:length(TessStar[[j]][,1]),1) #Sample a row.
    TessStar[[j]]=matrix(TessStar[[j]][-CenterRemoved,],ncol=length(Dim[[j]])) #Remove row sampled.
    Modification="RC"}
  else if (p<0.9 || length(Dim[[j]])==length(x[1,])){ #Change a center if p is less then 0.9 or if the all the covariates are in the tessellation.
    TessStar[[j]][sample(1:length(TessStar[[j]][,1]),1),]=rnorm(length(Dim[[j]]),0,sd) # Sample a row in the tessellaion matrix and change the coordinates of the centre by sampling from a normal distribution.
    Modification="Change"}
  else{ #Swop a dimension.
    NumberOfCovariates=1:length(x[1,])  #Let NumberOfCovariates be a vector from 1 to the number of covariates considered.
    NumberOfCovariates=NumberOfCovariates[-Dim[[j]]]  #Remove all values that for the covariates that are already in the tessellation.
    DimToChange=sample(1:length(Dim[[j]]),1) #Uniformly sample a dimension to change.
    DimStar[[j]][DimToChange]=sample(NumberOfCovariates,1) #Replace the Dimension to a new uniforly sampled covariate that is not already in the tessellaion.
    TessStar[[j]][,DimToChange]=rnorm(length(Tess[[j]][,1]),0,sd) #Add new normally sampled coordinates new dimension added.
    Modification="Swop"}
  
  TessStar[[j]]<-matrix(TessStar[[j]],ncol=length(DimStar[[j]])) #Ensure the the Tessellation matrix is a "matrix" type.
  
  return(list(TessStar,DimStar,Modification)) #Return new proposed tessellation.
}

fitting_function<- function(lambda,q,nu,sigmaSquared_hat){ #function that calculates the squared difference between sigma squared hat and the inverse gamma function
  return((sigmaSquared_hat- qinvgamma(q, shape=nu/2, rate=nu*lambda/2))^2)
}

Indexes<-function(x,Tess,Dim){ #Gives the row (the center) of the tessellation that each obseravtion falls within.
  if (length(Tess[,1])==1){ #only 1 centre
    CellsForGivenTess=rep(1,length(x[,1]))
  }
  else{ #multiple
    CellsForGivenTess=knnx.index(Tess,matrix(x[,Dim],ncol = length(Dim)),1)
  }
  return(CellsForGivenTess)
}

AlphaCalculation<-function(x,Tess,Dim,j,R_ijOld, SigmaSquared_old, R_ijNew,SigmaSquared_new,Modification,SigmaSquaredMu,Omega,lambda_rate){ #Calculates the acceptence rate of the proposed tessellation.
  
  d=length(Dim[[j]]);
  NumCovariates=length(x[1,]);
  cStar=length(Tess[[j]][,1]);
  
  #The Log Likelihood Ratio in the acceptence ratio
  LOGlikelihoodRatio=0.5*(log(prod(1+SigmaSquaredMu*SigmaSquared_old))-log(prod(1+SigmaSquaredMu*SigmaSquared_new)))+((SigmaSquaredMu/2)*(-sum((R_ijOld^2)/(1+SigmaSquaredMu*SigmaSquared_old))+sum((R_ijNew^2)/(1+SigmaSquaredMu*SigmaSquared_new))))

    #Calculating the acceptence probablity for "AD"=Adding a dimension, "RD"=Removing a dimension, "AC"=Adding a center, "RC"=Removing a center, "Change"=Changing the coordinates of a center and Swopping a dimension.
  if (Modification == "AD"){ 
    TessStructure=(dbinom(d-1,NumCovariates-1,Omega/NumCovariates))/(dbinom(d-2,NumCovariates-1,Omega/NumCovariates)*(NumCovariates-d+1))
    TransitionRatio=(NumCovariates-d+1)/d;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (length(Dim[[j]])==1){
      AcceptenceProb=AcceptenceProb+log(1/2)
    }
    else if (length(Dim[[j]])==NumCovariates-1){
      AcceptenceProb=AcceptenceProb+log(2)}
  }
  else if (Modification == "RD"){
    TessStructure=(dbinom(d-1,NumCovariates,Omega/NumCovariates)*(NumCovariates-d))/(dbinom(d,NumCovariates,Omega/NumCovariates))
    TransitionRatio=(d+1)/(NumCovariates-d)
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (length(Dim[[j]])==NumCovariates){
      AcceptenceProb=AcceptenceProb+log(1/2)
    }
    else if (length(Dim[[j]])==2){
      AcceptenceProb=AcceptenceProb+log(2);
    }
  }
  else if(Modification == "AC"){
    TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar-2,lambda_rate)
    TransitionRatio=1/cStar;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (cStar==1){
      AcceptenceProb=AcceptenceProb+log(1/2);
    }
  }
  else if (Modification == "RC"){
    TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar,lambda_rate);
    TransitionRatio=cStar+1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (cStar==2){
      AcceptenceProb=AcceptenceProb+log(2);
    }
  }
  else if (Modification == "Change"){
    TessStructure=1;
    TransitionRatio=1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
  }
  else {
    TessStructure=1;
    TransitionRatio=1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
  }
  return(AcceptenceProb)
}

CalculateResiduals<-function(y,x,j,SumOfAllTess,Tess,Dim,Pred,TessStar,DimStar,LastTessPred, sigma_squared_vector){ #A function that calculates the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation.
  if (j==1){
    indexes=Indexes(x,Tess[[j]],Dim[[j]]);
    CurrentTessPred<-Pred[[j]][indexes]
    SumOfAllTess=SumOfAllTess-CurrentTessPred}
  else{
    indexes=Indexes(x,Tess[[j]],Dim[[j]]);
    CurrentTessPred<-Pred[[j]][indexes]
    SumOfAllTess=SumOfAllTess+LastTessPred-CurrentTessPred;
  }
  
  IndexesStar=Indexes(x,TessStar[[j]],DimStar[[j]]);
  R_j<-y-SumOfAllTess
  
  #Initializing Sizes
  
  R_ijOld=rep(0,length(Pred[[j]]))
  Sigma_squared_ijOld=rep(0,length(Pred[[j]]))
  
  for (i in 1:length(Pred[[j]])){
    R_ijOld[i]<-sum(R_j[indexes==i]/sigma_squared_vector[indexes==i])
    Sigma_squared_ijOld[i]<-sum(1/sigma_squared_vector[indexes==i])
  }
  
  R_ijNew=rep(0,length(TessStar[[j]][,1]))
  Sigma_squared_ijNew=rep(0,length(TessStar[[j]][,1]))
  
  for (i in 1:length(TessStar[[j]][,1])){
    R_ijNew[i]<-sum(R_j[IndexesStar==i]/sigma_squared_vector[IndexesStar==i])
    Sigma_squared_ijNew[i]<-sum(1/sigma_squared_vector[IndexesStar==i])
  }
  
  return(list(R_ijOld,Sigma_squared_ijOld,R_ijNew,Sigma_squared_ijNew,SumOfAllTess,IndexesStar,indexes))
}

NewPredSet<-function(j,Tess,R_ijNew,sigmaSquaredMu,SigmaSquared){ #Sampling the new output values for the new tessellation.
  PredSet=rep(0,length(Tess[[j]][,1]))
  for (i in 1:length(Tess[[j]][,1])){
    PredSet[i]=rnorm(1,sigmaSquaredMu*R_ijNew[i]/(sigmaSquaredMu*SigmaSquared[i]+1),((sigmaSquaredMu)/(1+sigmaSquaredMu*SigmaSquared[i]))^0.5);
  }
  return(PredSet)
}

### Variance Functions ###

SigmaSquaredCalculation<-function(xScaled, yScaled,Tess_var,Dim_var,Pred_var, nu,lambda,ProdOfAllTess,Numerator_e,Omega,lambda_rate,m_var,sd){ #Sample sigma squared from inverse gamma distribution
  
  for (j in 1:m_var){
    NewTessOutput<-NewTess(xScaled,j,Tess_var,Dim_var,sd) #Propose new Tessellation 
    TessStar_var<-NewTessOutput[[1]]  
    DimStar_var<-NewTessOutput[[2]]
    Modification_var<-NewTessOutput[[3]]
    
    ResidualsOutput<-Calculate_variance_Residuals(yScaled,xScaled,j,ProdOfAllTess,Tess_var,Dim_var,Pred_var,TessStar_var,DimStar_var,LastTessPred_var,Numerator_e) #Calculate the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation and the number of observations in each cell.
    e_squared_ijOld<-ResidualsOutput[[1]]
    Resid_by_s_ijOld<-ResidualsOutput[[2]]
    n_ijOld<-ResidualsOutput[[3]]
    e_squared_ijNew<-ResidualsOutput[[4]]
    Resid_by_s_ijNew<-ResidualsOutput[[5]]
    n_ijNew<-ResidualsOutput[[6]]
    ProdOfAllTess<-ResidualsOutput[[7]]
    IndexesStar<-ResidualsOutput[[8]]
    indexes<-ResidualsOutput[[9]]

    if (!any(e_squared_ijNew==0)){ #automatically rejects proposed tessellation if there exists a cell with no observations in.

      LOGAcceptenceProb=AlphaCalculation_Variance(xScaled,TessStar_var,DimStar_var,j,e_squared_ijOld,n_ijOld,e_squared_ijNew,n_ijNew,Modification_var,nu,lambda,Omega,lambda_rate) #Gives the log of the acceptence probability.
      
      if (log(runif(n=1, min=0, max=1))<LOGAcceptenceProb){ #Accepts the proposed tessellation is accepted then calculates the new output values for the new tessellation. 
        Tess_var=TessStar_var
        Dim_var=DimStar_var
        Pred_var[[j]]=New_Variance_set(j,Tess_var,e_squared_ijNew,nu,lambda,n_ijNew)
        LastTessPred_var=Pred_var[[j]][IndexesStar]
        #if(Modification_var=='AC'){
          #print('accept')
          #print(Modification_var)
        #}
      }
      else { #Rejects the proposed tesellation then calculates new output values for the original tessellation.
        Pred_var[[j]]=New_Variance_set(j,Tess_var,e_squared_ijOld,nu,lambda,n_ijOld)
        LastTessPred_var=Pred_var[[j]][indexes];
        #print('reject1')
      }
    }
    else{ #Rejects the proposed tesellation then calculates new output values for the original tessellation.
      Pred_var[[j]]=New_Variance_set(j,Tess_var,e_squared_ijOld,nu,lambda,n_ijOld)
      LastTessPred_var=Pred_var[[j]][indexes];
      #print('reject2')
    }
    if (j==m_var){ #If j equals m then adds the last tessellation output values to give a prediction.
      ProdOfAllTess=ProdOfAllTess*LastTessPred_var;
    }
  }
  #print(Pred_var)
  return(list(ProdOfAllTess,Tess_var,Dim_var,Pred_var))
}

SigmaSquaredCalculation_initial<-function(yScaled,nu,lambda,SumOfAllTess){ #Sample sigma squared from inverse gamma distribution
  
  n=length(yScaled)
  SigmaSquared<-rinvgamma(1,shape=(nu+n)/2,rate=(nu*lambda+sum((yScaled-SumOfAllTess)^2))/2)
  
  return(SigmaSquared)
}

Calculate_variance_Residuals<-function(y,x,j,ProdOfAllTess,Tess_var,Dim_var,Pred_var,TessStar_var,DimStar_var,LastTessPred_var, Numerator_e){ #A function that calculates the n-vector of partial residuals derived from a fitting process that excludes the jth tessellation.
  if (j==1){
    indexes=Indexes(x,Tess_var[[j]],Dim_var[[j]]);
    CurrentTessPred<-Pred_var[[j]][indexes]
    ProdOfAllTess=ProdOfAllTess/CurrentTessPred}
  else{
    indexes=Indexes(x,Tess_var[[j]],Dim_var[[j]]);
    CurrentTessPred<-Pred_var[[j]][indexes]
    ProdOfAllTess=ProdOfAllTess*LastTessPred_var/CurrentTessPred;
  }

  IndexesStar=Indexes(x,TessStar_var[[j]],DimStar_var[[j]]);

  #Initializing Sizes
  e_squared_ijOld=rep(0,length(Pred_var[[j]]))
  Resid_by_s_ijOld=rep(0,length(Pred_var[[j]]))
  n_ijOld=rep(0,length(Pred_var[[j]]))
  
  for (i in 1:length(Pred_var[[j]])){
    e_squared_ijOld[i]<-sum(Numerator_e[indexes==i]/ProdOfAllTess[indexes==i])
    Resid_by_s_ijOld[i]<-sum(Numerator_e[indexes==i]/(ProdOfAllTess[indexes==i])^2)
    n_ijOld[i]<-sum(indexes==i)
  }
  
  e_squared_ijNew=rep(0,length(TessStar_var[[j]][,1]))
  Resid_by_s_ijNew=rep(0,length(TessStar_var[[j]][,1]))
  n_ijNew=rep(0,length(TessStar_var[[j]][,1]))
  
  for (i in 1:length(TessStar_var[[j]][,1])){
    e_squared_ijNew[i]<-sum(Numerator_e[IndexesStar==i]/ProdOfAllTess[IndexesStar==i])
    Resid_by_s_ijNew[i]<-sum(Numerator_e[IndexesStar==i]/(ProdOfAllTess[IndexesStar==i])^2)
    n_ijNew[i]<-sum(IndexesStar==i)

    
  }

  return(list(e_squared_ijOld,Resid_by_s_ijOld,n_ijOld,e_squared_ijNew,Resid_by_s_ijNew,n_ijNew,ProdOfAllTess,IndexesStar,indexes))
}

AlphaCalculation_Variance<-function(x,Tess_var,Dim_var,j,e_squared_ijOld,n_ijOld,e_squared_ijNew,n_ijNew,Modification,nu,lambda,Omega,lambda_rate){ #Calculates the acceptence rate of the proposed tessellation.
  
  n<-length(x[,1])
  d=length(Dim_var[[j]]);
  NumCovariates=length(x[1,]);
  cStar=length(Tess_var[[j]][,1]);
  
  #The Log Likelihood Ratio in the acceptence ratio
  LOGlikelihoodRatio=sum(lgamma((nu+n_ijNew)/2))-sum(lgamma((nu+n_ijOld)/2))-sum(((nu+n_ijNew)/2)*log((nu*lambda+e_squared_ijNew)/2))+sum(((nu+n_ijOld)/2)*log((nu*lambda+e_squared_ijOld)/2))

  ##print('Log likihood')
  ##print(Modification)
  ##print(sum(lgamma((nu+n_ijNew)/2))-sum(lgamma((nu+n_ijOld)/2)))
  ##print(-sum(((nu+n_ijNew)/2)*log(nu*lambda+e_squared_ijNew))+sum(((nu+n_ijOld)/2)*log(nu*lambda+e_squared_ijOld)))
  ##print((nu/2)*log(nu*lambda/2)-log(lgamma(nu/2)))
  
  #Calculating the acceptence probablity for "AD"=Adding a dimension, "RD"=Removing a dimension, "AC"=Adding a center, "RC"=Removing a center, "Change"=Changing the coordinates of a center and Swopping a dimension.
  if (Modification == "AD"){ 
    TessStructure=(dbinom(d-1,NumCovariates-1,Omega/NumCovariates))/(dbinom(d-2,NumCovariates-1,Omega/NumCovariates)*(NumCovariates-d+1))
    TransitionRatio=(NumCovariates-d+1)/d;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (length(Dim_var[[j]])==1){
      AcceptenceProb=AcceptenceProb+log(1/2)
    }
    else if (length(Dim_var[[j]])==NumCovariates-1){
      AcceptenceProb=AcceptenceProb+log(2)}
  }
  else if (Modification == "RD"){
    TessStructure=(dbinom(d-1,NumCovariates,Omega/NumCovariates)*(NumCovariates-d))/(dbinom(d,NumCovariates,Omega/NumCovariates))
    TransitionRatio=(d+1)/(NumCovariates-d)
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)
    
    #Adjustments.
    if (length(Dim_var[[j]])==NumCovariates){
      AcceptenceProb=AcceptenceProb+log(1/2)
    }
    else if (length(Dim_var[[j]])==2){
      AcceptenceProb=AcceptenceProb+log(2);
    }
  }
  else if(Modification == "AC"){
    TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar-2,lambda_rate)
    TransitionRatio=1/cStar;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)+(nu/2)*log(nu*lambda/2)-lgamma(nu/2)

    #Adjustments.
    if (cStar==1){
      AcceptenceProb=AcceptenceProb+log(1/2);
    }
  }
  else if (Modification == "RC"){
    TessStructure=dpois(cStar-1,lambda_rate)/dpois(cStar,lambda_rate);
    TransitionRatio=cStar+1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio)-(nu/2)*log(nu*lambda/2)+lgamma(nu/2)
    
    #Adjustments.
    if (cStar==2){
      AcceptenceProb=AcceptenceProb+log(2);
    }
  }
  else if (Modification == "Change"){
    TessStructure=1;
    TransitionRatio=1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
  }
  else {
    TessStructure=1;
    TransitionRatio=1;
    AcceptenceProb=LOGlikelihoodRatio+log(TessStructure)+log(TransitionRatio);
  }
  return(AcceptenceProb)
}

New_Variance_set<-function(j,Tess_var,Resid_by_s,nu,lambda,n){ #Sampling the new output values for the new tessellation.
  Pred_var_Set=rep(0,length(Tess_var[[j]][,1]))
  for (i in 1:length(Tess_var[[j]][,1])){
    Pred_var_Set[i]<-rinvgamma(1,shape=(nu+n[i])/2,rate=(nu*lambda+(Resid_by_s[i]))/2)
  }

  return(Pred_var_Set)
  
}

TestPrediction<-function(x,m,Tess,Dim,Pred){ #A function that derives a prediction for the output variable give an observation, for one iteration in AddiVortes.
  Prediction=rep(0,length(x[,1]));
  for (j in 1:m){
    NewTessIndexes=Indexes(x,Tess[[j]],Dim[[j]]);
    Prediction=Prediction+Pred[[j]][NewTessIndexes]
  }
  return(Prediction)
}

TestPrediction_var<-function(x,m,Tess_var,Dim_var,Pred_var){ #A function that derives a prediction for the output variable give an observation, for one iteration in AddiVortes.
  Prediction=rep(1,length(x[,1]));
  for (j in 1:m){
    NewTessIndexes=Indexes(x,Tess_var[[j]],Dim_var[[j]]);
    Prediction=Prediction*Pred_var[[j]][NewTessIndexes]
  }

  return(Prediction)
}
