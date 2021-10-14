#############################################
##### POPULATION SCENARIO TO TEST ###########
#############################################

if (scenario==1){
#Scenario 1 

covid_exp1<-function(t){
  result=rep(0,length(t))
  result[t<=0.35]<-10
  result[t>0.35 & t<0.4]<-9.999987e+14*exp(-92.1034*t[t>0.35 & t<0.4]) #50 times smaller
  result[t>=0.4]<-0.1
  return(result)
}


covid_exp2<-function(t){
  t <- t
  result=rep(0,length(t))
  result[t<=0.35]<-10
  result[t>0.35 & t<0.4]<-9.999987e+14*exp(-92.1034*t[t>0.35 & t<0.4]) #50 times smaller
  result[t>=0.4]<-0.1
  return(0.5*result^1.5)
}

} else if (scenario==2){
#Scenario 2 

covid_exp1<-function(t){
  result=rep(0,length(t))
  result[t<=0.35]<-10
  result[t>0.35 & t<0.4]<-9.999987e+14*exp(-92.1034*t[t>0.35 & t<0.4]) #50 times smaller
  result[t>=0.4]<-0.1
  return(result)
}

covid_exp2<-function(t){
  t <- t
  result=rep(0.5,length(t))
  return(result)
}


} else if (scenario==3){
#Scenario 3


covid_exp1<-function(t){
  result=rep(0,length(t))
  result[t<=0.3]<-5
  result[t>0.3 & t<0.4]<- 625001.9*exp(-39.12023*t[t>0.3 & t<0.4]) #50 times smaller
  result[t>=0.4]<-0.1
  return(result)
}
covid_exp2<-function(t){
  t <- t +0.1
  result=rep(0,length(t))
  result[t<=0.3]<-5
  result[t>0.3 & t<0.4]<-625001.9*exp(-39.12023*t[t>0.3 & t<0.4]) #50 times smaller
  result[t>=0.4]<-0.1
  return(result)
}


} else if (scenario==4){
#Scenario 4 


covid_exp1<-function(t){
  result=rep(0,length(t))
  result[t<=0.3]<-5
  result[t>0.3 & t<0.4]<- 625001.9*exp(-39.12023*t[t>0.3 & t<0.4]) #50 times smaller
  result[t>=0.4]<-0.1
  return(result)
}
covid_exp2<-function(t){
  t <- t +0.3
  result=rep(0,length(t))
  result[t<=0.3]<-5
  result[t>0.3 & t<0.4]<-625001.9*exp(-39.12023*t[t>0.3 & t<0.4]) #50 times smaller
  result[t>=0.4]<-0.1
  return(result)
}


} else if (scenario==5){
#Scenario 5 

covid_exp1<-function(t){
  result=rep(0,length(t))
  result[t<0.3]<- exp(13.04008*t[t<0.3]-4.60517) #50 times smaller
  result[t>=0.3]<-0.5
  return(result)
}
covid_exp2<-function(t){
  result=rep(0,length(t))
  result[t<0.3]<-4*exp(-12.2962*t[t<0.3]) #50 times smaller
  result[t>=0.3]<-0.1
  return(result)
}

} else if (scenario==6){
  
  covid_exp1<-function(t){
    result=rep(0,length(t))
    result[t<=0.35]<-10
    result[t>0.35 & t<0.4]<-9.999987e+14*exp(-92.1034*t[t>0.35 & t<0.4]) #50 times smaller
    result[t>=0.4]<-0.1
    return(result)
  }
  
  covid_exp2<-function(t){
    result=rep(0,length(t))
    result[t<=0.35]<-10
    result[t>0.35 & t<0.4]<-9.999987e+14*exp(-92.1034*t[t>0.35 & t<0.4]) #50 times smaller
    result[t>=0.4]<-0.1
    return(result)
  }
  
}
