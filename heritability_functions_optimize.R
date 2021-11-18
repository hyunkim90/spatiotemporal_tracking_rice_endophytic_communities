##heritbility functions

#library(coxme)
library(HLMdiag)
library(lme4)

## ---- H2lmekin (multimatrix)

## H2lmekin2=function(model, g="id"){
##     varID=model$vcoef[names(model$vcoef)%in%g]
##     varexp=model$vcoef
##     varexp=varexp[names(varexp)!=g]
##     H2=sum(unlist(varID))/sum(unlist(varID), unlist(varexp),model$sigma^2)
##     return(H2)
## }

##
## H2lmekin=function(model, g="id"){
##     varID=model$vcoef[g]
##     varexp=model$vcoef
##     varexp=varexp[names(varexp)!=g]
##     H2=unlist(varID)/sum(unlist(varID), unlist(varexp),model$sigma^2)
##     return(H2)
## }


## ---- end-of-H2lmekin


## ---- varid

var_ID=function(m){
    require(HLMdiag)
    vc=varcomp.mer(m)
    H2=vc[2]/sum(vc)
    return(as.numeric(H2))
}
## ---- end-of-varid

## ---- H2

## H2=function(Y, K, id){
##     data=data.frame(Y=Y, id=id)
##     data$id=as.factor(data$id)
##     model=lmekin(Y~ (1|id) ,data, varlist=K, method="REML")
##     return(H2lmekin(model, g="id"))
## }

## ---- end-of-H2

H2nokin=function(Y,id,cov=NULL){
    if(is.null(cov)){
        data=data.frame(Y=Y, id=id)
        m=lmer(Y~ (1|id) ,data)}else{
            data=data.frame(Y=Y, id=id, cov=cov)
            m=lmer(Y~  cov + (1|id) ,data)}
    print(summary(m))
    return(var_ID(m))
    }

H2nokin.multi=function(Y,id1,id2,cov=NULL,rep){
    if(is.null(cov)){
        data=data.frame(Y=Y, id1=id1, id2=id2)
        m=lmer(Y~ id1+id2 + (1|rep) ,data)}else{
            data=data.frame(Y=Y,id1=id1, id2=id2, cov=cov)
            m=lmer(Y~ id1 + id2 + cov + (1|rep) ,data)}
    print(summary(m))
       return(var_ID(m))}


H2bootNokin=function(Y, id, cov=NULL,rep, nsim=100, ncpus=10){
    if(is.null(cov)){
        data=data.frame(Y=Y, id=id)
        m=lmer(Y~ id + (1|rep) ,data)
        b=bootMer(m, nsim=nsim, ncpus=ncpus, FUN=var_ID, parallel="multicore")}else{
            data=data.frame(Y=Y, id=id, cov=cov)
            m=lmer(Y~ id + cov + (1|rep) ,data)
            b=bootMer(m, nsim=nsim, ncpus=ncpus, FUN=var_ID, parallel="multicore")}
    return(b$t)
}


H2bootNokin.ori=function(Y, id, cov=NULL,rep, nsim=100, ncpus=10){
    if(is.null(cov)){
        data=data.frame(Y=Y, id=id)
        m=lmer(Y~ (1|id) ,data)
        b=bootMer(m, nsim=nsim, ncpus=ncpus, FUN=var_ID, parallel="multicore")}else{
            data=data.frame(Y=Y, id=id, cov=cov)
            m=lmer(Y~ cov + (1|id) ,data)
            b=bootMer(m, nsim=nsim, ncpus=ncpus, FUN=var_ID, parallel="multicore")}
    return(b$t)
}
