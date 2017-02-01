getStanPValue = function(obj){
pnorm(-abs(getStanMean(obj)/getStanSD(obj)))*2
}
fit.stan = sampling(stanDso, data = lData, pars = c('beta', 'iSize', 'betaSigma'),
iter = 5000, chains = 4)
m = as.matrix(fit.stan)
print(fit.stan)
getStanSummary = function(obj){
return(apply(extract(obj)$beta, 2, sd))
}
getStanSummary(fit.stan)
getStanSD = function(obj){
return(apply(extract(obj)$beta, 2, sd))
}
getStanMean = function(obj){
return(apply(extract(obj)$beta, 2, mean))
}
getStanPValue = function(obj){
pnorm(-abs(getStanMean(obj)/getStanSD(obj)))*2
}


### prepare data for stan
lStanData = list(Ntotal=nrow(dfData), Ncols=2, modmatrix= model.matrix(resp ~ pred, data=dfData), y=dfData$resp,
                 betaShape=lData$betaParam['shape'], betaRate=lData$betaParam['rate'])

library(rstan)

stanDso.old = rstan::stan_model(file='nb_glm.stan')
fit.stan = sampling(stanDso.old, data = lStanData, pars = c('beta', 'iSize', 'betaSigma'),
                    iter = 5000, chains = 4)

print(fit.stan)

