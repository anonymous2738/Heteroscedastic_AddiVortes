par(mfrow=c(1,2))
par(mgp = c(2.5, 1, 0))

set.seed(6)
n=length(Y)
TrainSet=sort(sample.int(n,3*n/5));
TestSet=1:n;
TestSet=TestSet[! TestSet %in% TrainSet];

Original_AddiVortes_non<-Homo_AddiVortes_Algorithm(Y[TrainSet],as.matrix(X[TrainSet,]),m=200 ,max_iter = 1200,burn_in = 300,6,0.85,0.5,0.8,5,25,Y[TestSet],as.matrix(X[TestSet,]))

Original_AddiVortes<-AddiVortes_Algorithm(Y[TrainSet],as.matrix(X[TrainSet,]),m=200,m_var = 40 ,max_iter = 1200,burn_in = 600,3,0.75,1,0.8,3,25,Y[TestSet],as.matrix(X[TestSet,]))

plotCI(sqrt(Original_AddiVortes$Sigma_squared_stored),sqrt(Original_AddiVortes$Sigma_squared_stored), Original_AddiVortes$UpperConfidenceTRAINValue -sqrt(Original_AddiVortes$Sigma_squared_stored),sqrt(Original_AddiVortes$Sigma_squared_stored) - Original_AddiVortes$LowerConfidenceTRAINValue,sfrac = 0, scol = 'grey',ylab = expression(hat(s) * (bold(x)) ~ "Posterior Intervals"), xlab = expression(hat(s)* (bold(x))), cex.lab = 1,pch = NA)
abline(mean(sqrt(Original_AddiVortes_non$sigma_squared_test)),0,lwd=3)
abline(quantile(sqrt(Original_AddiVortes_non$sigma_squared_test),0.95),0, lty =2,lwd=3)
abline(quantile(sqrt(Original_AddiVortes_non$sigma_squared_test),0.05),0,lty = 2,lwd=3)
