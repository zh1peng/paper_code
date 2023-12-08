library(rjags)
library(HDInterval)
library(MCMCvis)

#################################
# functions
#################################
cbind.fill <- function(...){ # can use list instead
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

# get all mode for posterior distribution
get_all_mode <- function(zm_post){
  mt=as.matrix(zm_post)
  all_mode=apply(mt, 2, get_mode)
  return(all_mode)
}

get_mode <- function(paramSampleVec){
  mcmcDensity = density(paramSampleVec)
  mo = mcmcDensity$x[which.max(mcmcDensity$y)]
  return(mo)
}

#######################################
# load observations: effect size 
#                    simulated variance
#######################################
setwd('F:\\Google Drive\\post-doc\\Bayesian_Project\\new_model\\revision\\gamma_wider')
all_site_es=read.csv('all_site_es.csv')
all_site_sd=read.csv('all_site_sd.csv')


#######################################
# MCMC model
#######################################
# MCMC:  Gelman-Rubin statistic
# the value of 1 indicate the chains are fully converged
# gelman.diag
# 
# ESS as how stable and accuarte the chain is
# ESS were set to above 10000 so that the sampler is stable

regions_list=list(L_ct=colnames(all_site_es)[grep('^L_.*\\_thickavg$',colnames(all_site_es))],
                 R_ct=colnames(all_site_es)[grep('^R_.*\\_thickavg$',colnames(all_site_es))])

for (regions_i in names(regions_list)){
# loop through all regions
col2test=regions_list[[regions_i]]
# an alternative way is to save them as list
all_mean=data.frame()
all_mode=data.frame()
all_rhat=data.frame()
all_ess=data.frame()
all_lower=data.frame()
all_upper=data.frame()
all_post=list()
n=0
for (col in col2test){
  n=n+1
  print(sprintf('%d--%s',n, col))
# regional data
site_col_es=all_site_es[,col]
site_col_sd=all_site_sd[,col]
site_w=site_col_sd
site_obv=site_col_es


# MCMC model
# site_mean-->�i
# mu-->M
# itau-->??
# site_w-->??i
# common_sigma-->??
# thes-->??i*??

model.sat.text<-"
  model {
for(i in 1:N) {
site_mean[i] ~ dnorm(mu,itau)
thes[i] <- 1/pow(site_w[i]*common_sigma,2)
site_obv[i] ~ dnorm(site_mean[i],thes[i])
}

mu ~ dnorm(0,10)
itau   ~ dgamma(1.01005, 0.01005) # mode 1.0, sd 100
tau <- pow(1/itau,1/2)
icommon_sigma ~ dgamma(1.01005, 0.01005)  # mode 1.0, sd 100
common_sigma <- pow(1/icommon_sigma,1/2)
}
"
model.sat.spec<-textConnection(model.sat.text)


zm1 <- jags.model(model.sat.spec,
                  data=list('site_w'=site_w,
                            'site_obv'=site_obv,
                            'N'=length(site_obv)
                  ),
                  n.adapt = 1000,
                  n.chains = 4)

zm_post=coda.samples(zm1,
                     c('mu','tau','site_mean','common_sigma'),
                     n.iter=100000)

# save all posteriors for sharing and plotting
all_post[[col]]=as.matrix(zm_post)

# use MCMCsummary to obtain some information
results=MCMCsummary(zm_post)
col_mean=subset(results,select='mean')
colnames(col_mean)=col
col_rhat=subset(results,select='Rhat')
colnames(col_rhat)=col
col_ess=subset(results,select='n.eff')
colnames(col_ess)=col

# get mean hdi get site hdi
hdi_results=hdi(zm_post, method = "HDI")
col_lower=data.frame(hdi_results['lower',])
colnames(col_lower)=col
col_upper=data.frame(hdi_results['upper',])
colnames(col_upper)=col

# get mode 
col_mode=data.frame(get_all_mode(zm_post))
colnames(col_mode)=col

all_mean=cbind.fill(all_mean,col_mean)
all_mode=cbind.fill(all_mode,col_mode)
all_rhat=cbind.fill(all_rhat,col_rhat)
all_ess=cbind.fill(all_ess,col_ess)
all_lower=cbind.fill(all_lower,col_lower)
all_upper=cbind.fill(all_upper,col_upper)
}
saveRDS(all_post, file = sprintf("Posterior_%s.rds",regions_i))
rm(all_post)
save(all_mean,all_mode,all_rhat,all_ess,all_lower,all_upper, file = sprintf("Posterior_extraction_%s.RData",regions_i))
rm(all_mean,all_mode,all_rhat,all_ess,all_lower,all_upper)
gc()
}








