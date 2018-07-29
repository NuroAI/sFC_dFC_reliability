
# sliding window metrices icc


library(R.matlab)
library(lme4)
library(stringr)

icc=function(data2){
  
  icc=array(0, dim=c(12720,1))
  
  for (i in 1:12720){
    
    data3=as.vector(data2[,,i])
    G=data.frame(subject=gl(820,4,820*4), time=gl(4,1,820), score=data3)
    
    # a random-intercept model (equivalent to the one-way ANOVA)
    a=VarCorr(lmer(score ~ 1 + (1|subject), data=G))   # extract std dev of random effects
    sd1=as.numeric(attributes(a$subject)$stddev)
    sd2=attributes(a)$sc
    icc[i]=sd1^2/(sd1^2+sd2^2)
    
    print(i)
  }
  
  m = matrix(0, 160, 160)
  m[lower.tri(m)] = 1
  m1=as.vector(m)
  ind=which(m1>0)
  m1[ind]=icc
  m1=matrix(m1,160,160)
  m1=m1+t(m1)
  return(m1)
}




icc=array(0, dim=c(12720,1))

for (i in 1:12720){
  
  data3=as.vector(data2[,,i])
  G=data.frame(subject=gl(820,4,820*4), time=gl(4,1,820), score=data3)
  
  # a random-intercept model (equivalent to the one-way ANOVA)
  a=VarCorr(lmer(score ~ 1 + (1|subject), data=G))   # extract std dev of random effects
  sd1=as.numeric(attributes(a$subject)$stddev)
  sd2=attributes(a)$sc
  icc[i]=sd1^2/(sd1^2+sd2^2)
  
  print(i)
}

m = matrix(0, 160, 160)
m[lower.tri(m)] = 1
m1=as.vector(m)
ind=which(m1>0)
m1[ind]=icc
m1=matrix(m1,160,160)
m1=m1+t(m1)





# ICA time series test
for (j in 2:7){
  
  sw=str_pad(j*10, 3, pad = "0")
  
  data1 = readMat(paste('E:/imp_Chao/projects/07_TRT/dynamic_820_ICAd100_slidingwindow_bandpass_highpass_win',sw,'_leonardi2015.mat',sep=''))
  #d1 = data1$sw.d1       
  #d2 = data1$sw.d2
  #d3 = data1$sw.d3
  d4 = data1$sw.d4
  
  #icc1=icc(d1)
  #icc2=icc(d2)
  #icc3=icc(d3)
  icc4=icc(d4)
  
  writeMat(paste('E:/imp_Chao/projects/07_TRT/icc_dynamic_820_ICAd100_slidingwindow_bandpass_highpass_win',sw,'_leonardi2015.mat',sep=''),icc4=icc4)
  print(j)
}











for (j in 1:20){
  
  sw=str_pad(j*10, 3, pad = "0")
  
  data1 = readMat(paste('E:/imp_Chao/projects/07_TRT/dynamic_820_slidingwindow_bandpass_highpass_win',sw,'_leonardi2015.mat',sep=''))
  d1 = data1$sw.d1       # 2*20*12720 instead of 4*820*12720
  d2 = data1$sw.d2
  d3 = data1$sw.d3
  #d4 = data1$sw.d4
  
  icc1=icc(d1)
  icc2=icc(d2)
  icc3=icc(d3)
  #icc4=icc(d4)
  
  writeMat(paste('E:/imp_Chao/projects/07_TRT/dynamic_820_slidingwindow_bandpass_highpass_win',sw,'_leonardi2015.mat',sep=''),icc1=icc1,icc2=icc2,icc3=icc3)
  print(j)
}

rm(list=ls()) 


for (j in 1:20){
  
  sw=str_pad(j*10, 3, pad = "0")
  
  data1 = readMat(paste('E:/imp_Chao/projects/DCC/dynamic_820_slidingwindow_Nobandpass_Nohighpass_win',sw,'_shift025_std.mat',sep=''))
  #d1 = data1$sw.d1       # 2*20*12720 instead of 4*820*12720
  d2 = data1$sw.d2
  #d3 = data1$sw.d3
  #d4 = data1$sw.d4
  
  #icc1=icc(d1)
  icc2=icc(d2)
  #icc3=icc(d3)
  #icc4=icc(d4)
  
  writeMat(paste('E:/imp_Chao/projects/DCC/icc_dynamic_820_slidingwindow_Nobandpass_Nohighpass_win',sw,'_shift025_std.mat',sep=''),icc2=icc2)
  print(j)
}



data1 = readMat('E:/imp_Chao/projects/DCC/static_fc_bandpass_DOS160_HCP20_run1run2.mat')
d1 = data1$fc 
icc1=icc(d1)
writeMat('E:/imp_Chao/projects/DCC/icc_static_fc_bandpass_DOS160_HCP20_run1run2.mat',icc1=icc1)






data1 = readMat(paste('E:/imp_Chao/projects/DCC/std_','tc.mat',sep=''))
#d1 = data1$sw.d1       # 2*20*12720 instead of 4*820*12720
d2 = data1$std.tc
icc2=icc(d2)
writeMat(paste('E:/imp_Chao/projects/DCC/icc_std','_tc.mat',sep=''),icc2=icc2)






