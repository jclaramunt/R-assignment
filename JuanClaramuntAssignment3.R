#####################################
#Part 1                             #
#####################################



#Define the data
GroupExp<-rnorm(50,160,15)#group which receive the treatment
GroupCnt<-rnorm(50,150,15)#control group

#Run a t test using the in-biult t.test() function
ttest<-t.test(GroupExp,GroupCnt)
ttpvalue<-ttest$p.value#store the pvalue
ttpvalue#print the p value
#The p-value is lower than 0.05 so it seems that there is a significant difference between the means of the two groups
#which in fact should be as the mean of the control group is 150 and the mean of the treatment group is 160
#This means that there is a probability of less than 5% that both means are equal.  

#compute the theoretical power of the previous test
ttpower<-power.t.test(n=50,delta=10,sd=15,sig.level = 0.05)



#Approximate the power by simulation. 

#Make a function which computes the power of a t test function in which the two samples are normally distributed
TtestfunctionPower<- function(n,n2,mean,mean2,sd,sd2,nrep,alpha){#the inputs are the number of subjects in each sample, their means and standard deviations, the number of times the procedure is repeated and alpha level
  tstat<-matrix(rep(0,nrep),nrep,1)#generate the matrix to store the t-statistic values
  tpval<-matrix(rep(0,nrep),nrep,1)#generate the matrix to store the p-values values
    for (i in 1:nrep) {
    DataExp<-rnorm(n,mean,sd)#generate the treatment data
    DataCont<-rnorm(n2,mean2,sd2)#generate the control data
    tstat[i]<-t.test(DataExp,DataCont)$statistic#obtain the t-statistic
    tpval[i]<-t.test(DataExp,DataCont)$p.value#obtain the p-value
    }
  #power is prob of rejecting H0 given H1 is true
  Power<-sum(tpval<alpha)/length(tpval)#compute the empirical power
  #plot a histogram with the distribution of the p-values and the cut-off value
  hist(tpval,breaks=100,main='Distribution of the p-values',xlab = "p-values")
  abline(v=alpha,col="blue")
  theoretical.tpower<-power.t.test(n=n,delta=(mean-mean2),sd=sd,sig.level = alpha)$power #compute the theoretical power
  Error<-theoretical.tpower-Power #compute the error between the empirical and the theoretical power
  output<-list(tStatistic=tstat,tPvalues=tpval,Power=Power,Theoretical.Power=theoretical.tpower,Error=Error) #list with the outputs
} 
#Test the results
resultttest<-TtestfunctionPower(50,50,160,150,15,15,1000,0.05)
resultttest$Power#ask for the empirical power
resultttest$Error#ask for the difference between the empirical and the theoretical power
#the greater the amount of repetitions the smaller the error.
#With 1000 simulations the error is in general lower than 0.01

#####################################
#Part 2                             #
#####################################


#vectors with the values of the example to test the functions
Data1<-matrix(c(2,5,5,6,6,7,8,9),nrow=8,ncol=1)
Data2<-matrix(c(1,1,2,3,3,4,5,7,7,8),nrow=10,ncol=1)


#We use bootstrap as a default method as it considers less assumptions than permutation and then it is able to deal with more inputs.  Then we think Bootstrap is the best choice for this exercise.
#However, as there are cases in which Permutation can be better we also allow the user to use the more convenient
MyTtestresamp<- function(Sample1,Sample2,resampling=c("Bootstrap","Permutation"),nres){#the imputs are both samples, the desired method and the number of resamples
Samp1<- na.omit(Sample1) #listwise deletion of missing inputs
Samp2<- na.omit(Sample2) #listwise deletion of missing inputs
SampleA<-as.matrix(Samp1) #transform the data into matrix form
SampleB<-as.matrix(Samp2) #transform the data into matrix form 
Sampletotal<-c(SampleA,SampleB) #joint the data
sizeA<-dim(SampleA)[1]#compute the number of subjects in Sample1 after deletion of missing values
sizeB<-dim(SampleB)[1]#compute the number of subjects in Sample2 after deletion of missing values
Samptval<-t.test(SampleA, SampleB)$statistic#compute the t-statistic of the original data
tstat<-matrix(rep(0,nres),nrow=nres,ncol=1)#generate the matrix to store the t-statistic values

if (resampling=="Bootstrap"){#To introduce the resampling method desired by the user, in this case botstrap

  for (i in 1:nres) {#for loop to do bootstrapping
    Dataresamp<-sample(Sampletotal,replace=TRUE)#sample with replacement
    DataA<-Dataresamp[1:sizeA]#assign the first resampled elements to the first group
    DataB<-Dataresamp[(sizeA+1):(sizeA+sizeB)]#assign the rest to the second group
    tstat[i]<-t.test(DataA, DataB)$statistic#compute the t-statistic
  }
}else if (resampling=="Permutation") { #To use permutation as resampling method
  
  for (i in 1:nres) {#for loop to do Permutation
    Dataresamp<-sample(Sampletotal,replace=FALSE)#sample without replacement
    DataA<-Dataresamp[1:sizeA]#assign the first resampled elements to the first group
    DataB<-Dataresamp[(sizeA+1):(sizeA+sizeB)]#assign the rest to the second group
    tstat[i]<-t.test(DataA, DataB)$statistic#compute the t-statistic
  }
}


tval.larger.samptval<- sum(abs(tstat)>abs(Samptval))#compute the number of resamples whose t-statistic is larger in absolute value than the original t-statistic
pvalresamp<-tval.larger.samptval/nres#compute the p value as the proportion of larger t-values over the total number of t values
#plot a histogram with the Distribution of the t-statistic
hist(tstat,main='Distribution of the t-statistic along the samples',xlab = "t-statistic")
abline(v=Samptval,col="blue")#add a vertical lines with the original t-value 
CI<-quantile(tstat,c(0.025,0.975))#compute the CI
abline(v=CI,col="red")#add vertical lines with the CI

output<-list(sample.tvalue=Samptval,tvalues=tstat,pvalue=pvalresamp, CI=CI)  #list with the outputs
}
#test the function
result1<-MyTtestresamp(Data1,Data2,resampling = "Bootstrap",nres=1000)
# We can compare it with permutation by uncomment and run the following line
#result1.2<-MyTtestresamp(Data1,Data2,resampling = "Permutation",B=1000)
#We obtain a p-value greater than 0.05, so we consider that the mean of both samples are equal.



MyTtestresamp2<- function(Sample1,Sample2,resampling=c("Bootstrap","Permutation"),nres){#the imputs are both samples, the desired method and the number of resamples
  Samp1<- na.omit(Sample1) #listwise deletion of missing inputs
  Samp2<- na.omit(Sample2) #listwise deletion of missing inputs
  SampleA<-as.matrix(Samp1) #transform the data into matrix form
  SampleB<-as.matrix(Samp2) #transform the data into matrix form
  Sampletotal<-c(SampleA,SampleB) #combine the data
  sizeA<-dim(SampleA)[1]#compute the number of subjects in Sample1 after deletion of missing values
  sizeB<-dim(SampleB)[1]#compute the number of subjects in Sample2 after deletion of missing values
  sizeT<-sizeA+sizeB#compute the total size
  Samptval<-t.test(SampleA, SampleB)$statistic#compute the t-statistic of the original data
  tstat<-matrix(rep(0,nres),nrow=nres,ncol=1)#generate the matrix to store the t-statistic values
  
  
  if (resampling=="Bootstrap"){ #To introduce the resampling method desired by the user, in this case botstrap 
  for (i in 1:nres) {#for loop to do bootstrapping
    Dataresamp<-matrix(rep(0,sizeT),nrow=sizeT,ncol=1)#generate the matrix to store the samples
    for (j in 1:sizeT) {#for to create each sample
      samp1<-runif(n=1,min=0,max=sizeT)#obtain a random number from 0 to 16
      samp2<-ceiling(samp1)#round it such that we obtain integer numbers between 1 and 16
      Dataresamp[j]=Sampletotal[samp2]#store the number of the total sample corresponding to the samp2-th position
    }
    DataA<-Dataresamp[1:sizeA]#assign the first resampled elements to the first group
    DataB<-Dataresamp[(sizeA+1):sizeT]#assign the rest to the second group
    tstat[i]<-t.test(DataA, DataB)$statistic#compute the t-statistic
  }
  }else if (resampling=="Permutation") { #To use permutation as resampling method 
    for (i in 1:nres) {#for loop to do Permutation
      Dataresamp<-matrix(rep(0,sizeT),nrow=sizeT,ncol=1)#generate the matrix to store the samples
      
      cont<-0 #variable to count the numbers in the elements for each sample
      while (cont<sizeT) {#while loop to sample until the sample is complete
        samp1<-runif(n=1,min=0,max=sizeT)#obtain a random number from 0 to 16
        samp2<-ceiling(samp1)#round it such that we obtain integer numbers between 1 and 16
        if (!(samp2 %in% Dataresamp)){#only add a number if it is not in the sample to avoid sample with replacement
        Dataresamp[(cont+1)]=Sampletotal[samp2]#store the number of the total sample corresponding to the samp2-th position
        cont<-cont+1#add one to the variable which counts elements as we have one more element in the sample
        }
      }
      
      
      DataA<-Dataresamp[1:sizeA]#assign the first resampled elements to the first group
      DataB<-Dataresamp[(sizeA+1):sizeT]#assign the rest to the second group
      tstat[i]<-t.test(DataA, DataB)$statistic#compute the t-statistic
    }
    
  }
  
  tval.larger.samptval<- sum(abs(tstat)>abs(Samptval))#compute the number of resamples whose t-statistic is larger in absolute value than the original t-statistic
  pvalresamp<-tval.larger.samptval/nres#compute the p value as the proportion of larger t-values over the total number of t values
  #plot a histogram with the Distribution of the t-statistic
  hist(tstat,main='Distribution of the t-statistic along the samples',xlab = "t-statistic")
  abline(v=Samptval,col="blue")#add a vertical lines with the original t-value 
  CI<-quantile(tstat,c(0.025,0.975))#compute the CI
  abline(v=CI,col="red")#add vertical lines with the CI
  output<-list(sample.tvalue=Samptval,tvalue=tstat,pvalue=pvalresamp, CI=CI)  #list with the outputs
}

#test the function
Result2<-MyTtestresamp2(Data1,Data2,"Bootstrap",1000)
# We can compare it with permutation by uncomment and run the following line
#Result2.1<-MyTtestresamp2(Data1,Data2,"Permutation",1000)
#Again, we obtain a p-value greater than 0.05, so we consider that the mean of both samples are equal.

#If we want to compare the results we set a number and run both functions:
set.seed(1000)
Comp1<-MyTtestresamp(Data1,Data2,"Bootstrap",1000)
set.seed(1000)
Comp2<-MyTtestresamp2(Data1,Data2,"Bootstrap",1000)
#We observe that both functions lead to the same results


#We use bootstrap again because the datasets are small, and consequently Jackknife, which is the other alternative, is not so stable
#Furthermore, we use jacknife after bootstrap to obtain estimates of bias and SE
MyTtestresamp3<- function(Sample1,Sample2,nres){ #inputs are two independent samples and the number of resamples desired
  Samp1<- na.omit(Sample1) #listwise deletion of missing inputs
  Samp2<- na.omit(Sample2) #listwise deletion of missing inputs
  SampleA<-as.matrix(Samp1) #transform the data into matrix form
  SampleB<-as.matrix(Samp2) #transform the data into matrix form 
  sizeA<-dim(SampleA)[1] #compute the number of elements in Sample1 after deletion of missing values
  sizeB<-dim(SampleB)[1] #compute the number of elements in Sample2 after deletion of missing values
  sizeT<-sizeA+sizeB #compute the total number of elements in both samples
  Smeandiff<-mean(SampleA)-mean(SampleB)#compute the mean difference in the original samples
  meandiffresamp<-matrix(rep(0,nres),nrow=nres,ncol=1) #generate the matrix to store the  mean difference values obtained by resampling
  indices<-matrix(0, nrow = nres, ncol = sizeT) #generate a matrix to store the sampled indices (to do jacknife after bootstrap) 
  for (i in 1:nres) {#for loop to resample
    indicesdataA<-sample(1:sizeA,replace=TRUE)#sample with replacement the indices for sampleA
    DataAestim<-SampleA[indicesdataA]#Obtain the data of sampleA for the sampled indices
    indicesdataB<-sample((sizeA+1):sizeT,replace=TRUE)#sample with replacement the indices for sampleB
    DataBestim<-SampleB[(indicesdataB-sizeA)]#Obtain the data of sampleB for the sampled indices
    meandiffresamp[i]<-mean(DataAestim)-mean(DataBestim)#compute the mean difference
    indices[i, ] <- c(indicesdataA,indicesdataB)#store the indices
  }
  sd.jack <- numeric(sizeT) #generate vector to store jackknife standard error
  for (i in 1:sizeT) { #for loop to run jackknife
    keep <- (1:nres)[apply(indices, MARGIN = 1,
                             FUN = function(k) {!any(k == i)})] #keep the indices different from i
    sd.jack[i] <- sd(meandiffresamp[keep]) #compute the standard deviation of the mean differences of the samples which not contain the i-th index
  }
  #plot a histogram of the distribution of the mean differences obtained by resampling
  hist(meandiffresamp,main='Distribution of the estimated mean difference along the samples',xlab = "Mean difference")
  abline(v=Smeandiff,col="blue")#add the mean difference obtained from the original data
  CI<-quantile(meandiffresamp,c(0.025,0.975))#compute the confidence intervals 
  abline(v=CI,col="red")#add the CI to the histogram
  mmdiffs<-mean(meandiffresamp)#compute the mean of the parameter estimated by resampling 
  Bias=mean(mmdiffs)-Smeandiff#compute the bias
  sd=sqrt((1/(nres-1))*sum((mmdiffs-(mmdiffs/nres))^2))#compute the standard error
  sdjack<-(sqrt((sizeT-1)*mean((sd.jack-mean(sd.jack))^2)))#compute the standard error using jacknife after bootstrap obtaining the SE of the standard error and bias
  output<-list(mean.differences.resamp=meandiffresamp, CI.estimation=CI, Bias=Bias, standar.error=sd, standard.error.jacknife=sd.jack) #list with the outputs  
}
#Test the function
Result3<-MyTtestresamp3(Data1,Data2,1000)
#The resampled mean differences are centered around the mean difference of the original samples








