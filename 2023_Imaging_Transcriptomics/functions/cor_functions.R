

  


cor2p <- function(r,n=34){
  #n is 34 here but this is how to find pairwise n
  #n <- t(!is.na(x)) %*% (!is.na(y)) # same as count.pairwise(x,y) from psych package
  # a matrix version of cor.test https://stackoverflow.com/questions/13112238/a-matrix-version-of-cor-test
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  return(p)
}

  
corr_brain_gene <- function(gene_data,
                            brain_data,
                            method=c('pearson',
                                     'spearman'),
                            r2z=T){
  #make sure the regions are in a same order
  method=match.arg(method)
  if (!identical(rownames(brain_data),rownames(gene_data))){
    stop('Regions in brain_data and gene_data are not matched')
  }
  geneList=cor(gene_data,brain_data,method = method)
  
  if (r2z){
    geneList=apply(geneList,2,function(r){0.5*(log((1+r)/(1-r)))})
  }
  attr(geneList,'is_fisherz')=r2z
  attr(geneList,'n.region')=dim(gene_data)[1]
  return(geneList) # use the term from clusterprofiler
}
  
  
aggregate_geneSet <- function(geneList,# geneList and correlation
                              geneSet, # one geneSet of interest
                              method = c('mean', 
                                         'median', 
                                         'meanabs', 
                                         'meansqr',
                                         'maxmean', 
                                         'sig_n',
                                         'sign_test',
                                         'rank_sum',
                                         'ks_orig',
                                         'ks_weighted',
                                         'ks_sum',
                                         #'sig_n_chisqr', NA issue
                                         'locfdr'),
                              FUN=NULL) {
  if (!is.matrix(geneList)){
    stop('geneList should be a matrix of m*n, m is the genes, n is the number of (null) brains')
  }
  
  if (missing(method) &
      missing(FUN)) {
    stop('Either method or FUN is required')
  }
  if (missing(method) &
      !is.function(FUN)) {
    stop("Argument FUN is not a function")
  }
  
  if (method == 'mean') {
    FUN = function(genelist, geneSet) {
      geneSet = intersect(geneSet, names(genelist))
      hits = names(genelist) %in% geneSet # this is calculated for each function because the hit is not always the same for es
      return(mean(genelist[hits]))
    }
  } else if (method == 'median') {
    FUN = function(genelist, geneSet) {
      geneSet = intersect(geneSet, names(genelist))
      hits = names(genelist) %in% geneSet ## logical
      return(median(genelist[hits]))
    }
  } else if (method == 'meanabs') {
    FUN = function(genelist, geneSet) {
      geneSet = intersect(geneSet, names(genelist))
      hits = names(genelist) %in% geneSet 
      return(mean(abs(genelist[hits])))
    }
  } else if (method == 'meansqr') {
    FUN = function(genelist, geneSet) {
      geneSet = intersect(geneSet, names(genelist))
      hits = names(genelist) %in% geneSet
      return(mean((genelist[hits]) ^ 2))
    }
  } else if (method == 'maxmean') {
    FUN = function(genelist, geneSet) {
      geneSet = intersect(geneSet, names(genelist))
      hits = names(genelist) %in% geneSet 
      X = genelist[hits]
      # Mean of positive numbers
      pos.mean = mean(X[X > 0])
      # Mean of negative numbers
      neg.mean = mean(X[X < 0])
      
      if (is.na(pos.mean)){pos.mean=0} # when all value are neg
      if (is.na(neg.mean)){neg.mean=0} # when all value are pos
      
      maxmean = ifelse(pos.mean > abs(neg.mean), pos.mean, neg.mean)
      return(maxmean)
    }
  } else if (method == 'ks_orig') {
    FUN = function(genelist, geneSet) {
      #function taken from clusterprofiler
      genelist = sort(genelist, decreasing = TRUE)
      geneSet = intersect(geneSet, names(genelist))
      N = length(genelist)
      Nh = length(geneSet)
      Phit = Pmiss = numeric(N)
      hits = names(genelist) %in% geneSet
      Phit[hits] = abs(genelist[hits]) ^ 0 #raw rank
      NR = sum(Phit)
      Phit = cumsum(Phit / NR)
      Pmiss[!hits] = 1 / (N - Nh)
      Pmiss = cumsum(Pmiss)
      runningES = Phit - Pmiss
      max.ES = max(runningES)
      min.ES = min(runningES)
      if (abs(max.ES) > abs(min.ES)) {
        ES = max.ES
      } else {
        ES = min.ES
      }
      return(ES)
      }
  } else if (method == 'ks_weighted') {
    FUN = function(genelist, geneSet) {
      #function taken from clusterprofiler
      genelist = sort(genelist, decreasing = TRUE)
      geneSet = intersect(geneSet, names(genelist))
      N = length(genelist)
      Nh = length(geneSet)
      Phit = Pmiss = numeric(N)
      hits = names(genelist) %in% geneSet
      Phit[hits] = abs(genelist[hits]) ^ 1 #weighted rank
      NR = sum(Phit)
      Phit = cumsum(Phit / NR)
      Pmiss[!hits] = 1 / (N - Nh)
      Pmiss = cumsum(Pmiss)
      runningES = Phit - Pmiss
      max.ES = max(runningES)
      min.ES = min(runningES)
      if (abs(max.ES) > abs(min.ES)) {
        ES = max.ES
      } else {
        ES = min.ES
      }
      return(ES)
    }
  } else if (method == 'ks_sum_pos_neg') {
    FUN = function(genelist, geneSet) {
      #function taken from clusterprofiler
      genelist = sort(genelist, decreasing = TRUE)
      geneSet = intersect(geneSet, names(genelist))
      N = length(genelist)
      Nh = length(geneSet)
      Phit = Pmiss = numeric(N)
      hits = names(genelist) %in% geneSet
      Phit[hits] = abs(genelist[hits]) ^ 1 #weighted rank
      NR = sum(Phit)
      Phit = cumsum(Phit / NR)
      Pmiss[!hits] = 1 / (N - Nh)
      Pmiss = cumsum(Pmiss)
      runningES = Phit - Pmiss
      max.ES = max(runningES)
      min.ES = min(runningES)
      # if (abs(max.ES) > abs(min.ES)) {
      #   ES = max.ES
      # } else {
      #   ES = min.ES
      # }
      ES=max.ES+min.ES
      return(ES)
    }
  }else if (method == 'ks_pos') {
    FUN = function(genelist, geneSet) {
      #function taken from clusterprofiler
      genelist = sort(genelist, decreasing = TRUE)
      geneSet = intersect(geneSet, names(genelist))
      N = length(genelist)
      Nh = length(geneSet)
      Phit = Pmiss = numeric(N)
      hits = names(genelist) %in% geneSet
      Phit[hits] = abs(genelist[hits]) ^ 1 #weighted rank
      NR = sum(Phit)
      Phit = cumsum(Phit / NR)
      Pmiss[!hits] = 1 / (N - Nh)
      Pmiss = cumsum(Pmiss)
      runningES = Phit - Pmiss
      max.ES = max(runningES)
      min.ES = min(runningES)
      # if (abs(max.ES) > abs(min.ES)) {
      #   ES = max.ES
      # } else {
      #   ES = min.ES
      # }
      ES=max.ES
      return(ES)
    }
  }else if (method == 'ks_neg') {
    FUN = function(genelist, geneSet) {
      #function taken from clusterprofiler
      genelist = sort(genelist, decreasing = TRUE)
      geneSet = intersect(geneSet, names(genelist))
      N = length(genelist)
      Nh = length(geneSet)
      Phit = Pmiss = numeric(N)
      hits = names(genelist) %in% geneSet
      Phit[hits] = abs(genelist[hits]) ^ 1 #weighted rank
      NR = sum(Phit)
      Phit = cumsum(Phit / NR)
      Pmiss[!hits] = 1 / (N - Nh)
      Pmiss = cumsum(Pmiss)
      runningES = Phit - Pmiss
      max.ES = max(runningES)
      min.ES = min(runningES)
      # if (abs(max.ES) > abs(min.ES)) {
      #   ES = max.ES
      # } else {
      #   ES = min.ES
      # }
      ES=min.ES
      return(ES)
    }
  }else if (method=='sign_test'){
    FUN=function(genelist,geneSet){
      geneSet = intersect(geneSet, names(genelist))
      hits = names(genelist) %in% geneSet 
      X = genelist[hits]
      n.pos=sum(X>0)
      n.neg=sum(X<0)
      n.smaller=ifelse(n.pos<n.neg,n.pos,n.neg)
      n.total=n.pos+n.neg
      z=ifelse(n.total<=25,n.smaller,((n.smaller+0.5)-(n.total/2))/(sqrt(n.total)/2)) 
      return(z)
    }
  }else if (method=='locfdr'){
    FUN=function(genelist, geneSet){
    geneSet = intersect(geneSet, names(genelist))
    hits = names(genelist) %in% geneSet ## logical
    df=data.frame(vals=as.numeric(genelist),hits=as.numeric(hits))
    
    fit.pos=glm(hits~vals,family = binomial(link = "logit"),data=df[df$vals>=0,])
    fit.neg=glm(hits~vals,family = binomial(link = "logit"),data=df[df$vals<0,])
    S.pos=coef(summary(fit.pos))['vals',"z value"]
    S.neg=coef(summary(fit.neg))['vals',"z value"]
    S=ifelse(abs(S.pos)>abs(S.neg),S.pos,S.neg)
    return(S)
        }
  }else if (method=='sig_n'){
    FUN=function(genelist,geneSet,
                 n=attr(geneList,'n.region'),
                 z2r=attr(geneList,'is_fisherz')){
      geneSet = intersect(geneSet, names(genelist))
      hits = names(genelist) %in% geneSet 
      r = genelist[hits]
      if(z2r){r=tanh(r)}
      t = (r*sqrt(n-2))/sqrt(1-r^2)
      p = 2*(1 - pt(abs(t),(n-2)))
      adj_p=p.adjust(p,'fdr')
      sig_n=sum(adj_p<0.05)
      return(sig_n)
    }
  } else if (method=='rank_sum'){
    FUN=function(genelist,geneSet){
      abslist.sorted = sort(abs(genelist), decreasing = F) # sort by abs value
      ranklist=c(1:length(genelist)) # create rank
      names(ranklist)=names(abslist.sorted) # set the gene names
      
      geneSet = intersect(geneSet, names(genelist))
      hits = names(genelist) %in% geneSet
      pos = genelist>0
      neg = genelist<0
      pos.name=names(genelist)[hits&pos]
      neg.name=names(genelist)[hits&neg]
      
      pos.ranksum=sum(ranklist[pos.name])
      neg.ranksum=sum(ranklist[neg.name])
      W=ifelse(pos.ranksum<neg.ranksum,pos.ranksum,neg.ranksum)
      return(W)
    }
  }
    
  # }else if (method=='sig_n_chisqr'){
  #   FUN=function(genelist,geneSet,
  #                n=attr(geneList,'n.region'),
  #                z2r=attr(geneList,'is_fisherz')){
  #     geneSet = intersect(geneSet, names(genelist))
  #     hits = names(genelist) %in% geneSet 
  #     r = genelist
  #     if(z2r){r=tanh(r)}
  #     t = (r*sqrt(n-2))/sqrt(1-r^2)
  #     p = 2*(1 - pt(abs(t),(n-2)))
  #     
  #     hits.n=sum(hits)
  #     hits.sig.n=sum(p.adjust(p[hits],'fdr')<0.05)
  #     
  #     all.n=length(genelist)
  #     all.sig.n=sum(p.adjust(p,'fdr')<0.05)
  #     
  #     mx=matrix(c(hits.sig.n,hits.n-hits.sig.n,all.sig.n, all.n-all.sig.n),2,2)
  #     chisqr=as.numeric(chisq.test(mx)['statistic']) # it tests if the freq observed in geneSet is independent of that observed in bg-genes
  #     return(chisqr)
  #   }
  
  
  if (!is.function(FUN)){stop('wrong method name')}
  gs.score=sapply(c(1:dim(geneList)[2]),
                  function(x){FUN(genelist=geneList[,x], # each col represents a null model # this will give a named num
                                  geneSet=geneSet)})
 
  return(gs.score)
}


aggregate_geneSetList <- function(geneSetList, # mutiple geneSets of interest after filtering
                                  geneList,...){
  allgs.scores=lapply(geneSetList, function(gs)
    aggregate_geneSet(geneList=geneList,
               geneSet=gs,...))
  return(allgs.scores)
}


caculate_pvals <- function(statList.true, statList.null) {
  if (!identical(names(statList.true),names(statList.null))){
    stop('statList.true and statList.null are not matched')
  }
  pvals=purrr::map2(statList.true, statList.null, function(true_stat, null_stat) {
    ifelse(true_stat >= 0,
           (sum(null_stat >= true_stat) + 1) / (sum(null_stat >= 0) + 1),
           (sum(null_stat <= true_stat) + 1) / (sum(null_stat < 0) + 1))
  })
  return(pvals)
}




