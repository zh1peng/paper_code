
library(pls)
corr_brain_gene <- function(gene_data, brain_data, method = c('pearson', 'spearman', 'pls1','loo','pls1w'), r2z = TRUE) {
  # Function to calculate the correlation between gene data and brain data
  
  # Arguments:
  #   gene_data: Data frame or matrix of gene expression data
  #   brain_data: Data frame or matrix of brain data
  #   method: Correlation method to use, can be 'pearson', 'spearman', 'pls1','loo'
  #           'pearson': Pearson correlation
  #           'spearman': Spearman correlation
  #            'pls1': PLS1 regression
  #            'loo': leave one region out approach (requested by reviewer)  
  #                   The correlation is calculated with each region left out
  #                   and then averaged across iterations.
  #   r2z: Logical indicating whether to convert correlation values to Fisher's z scores

  # Returns:
  #   A named matrix of correlations of coefficient between gene data and brain data
  #   This works for true model (when brain_data contains one columne) and null models (when brain_data contains multiple columns)
  
  method <- match.arg(method)
  
  # Make sure the regions are in the same order
  if (!identical(rownames(brain_data), rownames(gene_data))) {
    stop('Regions in brain_data and gene_data are not matched')
  }
  
  if (method %in% c('pearson', 'spearman')) {
    # Calculate correlation using Pearson or Spearman method
    geneList <- cor(gene_data, brain_data, method = method)
    
    # Set up for the p-value-based approach
    if (r2z) {
      geneList <- apply(geneList, 2, function(r) 0.5 * (log((1 + r) / (1 - r))))
    }
    
    # Set attributes
    attr(geneList, 'is_fisherz') <- r2z
    attr(geneList, 'n.region') <- dim(gene_data)[1]
    attr(geneList, 'cor_type') <- method
  }  else if (method == 'pls1') {
    # Calculate correlation using PLS1 regression
    # The coefficents are extracted for the following analysis (aggregate_geneSet)
    X <- apply(gene_data, 2, scale)
    y <- apply(brain_data, 2, scale)
      # do pls1 for each column
      pls1 <- apply(y,2,function(y_col) pls::plsr(y_col ~ X, 1))
      # extract the coefficients for the first component
      coef_list <- lapply(pls1,function(pls_i) drop(pls_i$coefficients))
      geneList=do.call(cbind,coef_list)
    # Set attributes
    attr(geneList, 'is_fisherz') <- NULL
    attr(geneList, 'n.region') <- dim(gene_data)[1]
    attr(geneList, 'cor_type') <- method
  } else if (method == 'pls1w') {
    # Calculate correlation using PLS1 regression
    # The coefficents are extracted for the following analysis (aggregate_geneSet)
    X <- apply(gene_data, 2, scale)
    y <- apply(brain_data, 2, scale)
      # do pls1 for each column
      pls1 <- apply(y,2,function(y_col) plsr(y_col ~ X, 1))
      # extract the coefficients for the first component
      weights_list <- lapply(pls1,function(pls_i) drop(loading.weights(pls_i)))
      geneList=do.call(cbind,weights_list)
    # Set attributes
    attr(geneList, 'is_fisherz') <- NULL
    attr(geneList, 'n.region') <- dim(gene_data)[1]
    attr(geneList, 'cor_type') <- method
  } else if (method == 'loo'){
    geneList <- loo_pearson(gene_data, brain_data)
    
    # Set attributes
    attr(geneList, 'is_fisherz') <- TRUE
    attr(geneList, 'n.region') <- dim(gene_data)[1]
    attr(geneList, 'cor_type') <- method

  }
  
  return(geneList)
}


loo_pearson <- function(gene_data, brain_data) {
  # Function to calculate the correlation between gene data and brain data
  # using leave one region out approach
  # the correlation is calculated with each region left out
  # and then averaged across iterations
  geneList_all <- list()
  
  for (region_i in 1:nrow(brain_data)) {
    # leave one row (region) out and do the correlation
    brain_data_loo <- brain_data[-region_i, ]
    gene_data_loo <- gene_data[-region_i, ]
  
    geneList_loo <- cor(gene_data_loo, brain_data_loo, method = 'pearson')
    geneList_loo <- apply(geneList_loo, 2, function(r) 0.5 * (log((1 + r) / (1 - r))))
    geneList_all[[region_i]] <- geneList_loo
  }
  
  geneList_avg <- average_matrices(geneList_all)
  return(geneList_avg)
}

average_matrices <- function(matrix_list) {
  # Calculate the sum of matrices using Reduce and +
  sum_matrix <- Reduce(`+`, matrix_list)
  
  # Divide the sum_matrix by the number of matrices
  avg_matrix <- sum_matrix / length(matrix_list)
  
  return(avg_matrix)
}
  
aggregate_geneSet <- function(geneList,# named correlation/coefficient(pls1) matrix
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

  # Function to aggregate geneList based on geneSet (one geneSet at a time)

  # Arguments:
  #   geneList: named correlation/coefficient(pls1) matrix; each col represents the true model or null models 
  #             geneList should be a matrix of genes*models (true model or null models)
  #   geneSet: A geneSet of interest
  #   method:  Method to use for aggregation
  #   FUN:     Function to use for aggregation, if specified, method will be ignored


  if (!is.matrix(geneList)){
    stop('geneList should be a matrix of genes*models (true model or null models)')
  }
  
  method = match.arg(method)

  if (missing(method) &
      missing(FUN)) {
    stop('Either method or FUN is required')
  }


  if (missing(method) &
      !is.function(FUN)) {
    stop("Argument FUN is not a function")
  }

  if (!missing(method) &
      !missing(FUN)) {
    stop("Both method and FUN are specified")
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

  if (!is.function(FUN)){stop('wrong method name')}
  gs.score=sapply(c(1:dim(geneList)[2]),
                  function(x){FUN(genelist=geneList[,x], # each col represents a null model # this will give a named num
                                  geneSet=geneSet)})
 
  return(gs.score)
}


aggregate_geneSetList <- function(geneSetList, # mutiple geneSets of interest after filtering
                                  geneList, # can be either true geneList or null geneList
                                  ...){
  allgs.scores=lapply(geneSetList, function(gs)
    aggregate_geneSet(geneList=geneList,
               geneSet=gs,...))
  return(allgs.scores)
}


aggregate_geneSetList_with_constrain <- function(geneSetList, # mutiple geneSets of interest after filtering
                                                  sampled_geneSetList,
                                                  geneList.true, # should be true geneList only!
                                                  ...){


  # if geneList.true has more than 2 columns, it is a null geneList
if (!is.matrix(geneList.true)){
  stop('geneList.true should be a m x 1 vector/matrix. Pls include drop=F when subsetting')
}else{
    if (dim(geneList.true)[2]>1){
    stop('geneList.true should be a m x 1 vector/matrix.')}
}

# make sure the order in geneSetList and sampled_geneSetList are the same
if (!identical(names(geneSetList),names(sampled_geneSetList))){
  stop('geneSetList and sampled_geneSetList are not matched')
}

# lapply does not work well here
allgs.scores=mapply(function(gs, sampled_gs){
  geneList.null=swap_geneList(geneList.true=geneList.true,
                                   orig_gs=gs,
                                   sampled_gs = sampled_gs)
  gs.score=aggregate_geneSet(geneList=geneList.null,
                        geneSet=gs,
                        ...)
 return(gs.score)},
 geneSetList,sampled_geneSetList,SIMPLIFY = FALSE)

  return(allgs.scores)
}





caculate_pvals <- function(statList.true, statList.null,method=c('standard','split_pos_neg')) {
  
  method=match.arg(method)
  if (!identical(names(statList.true),names(statList.null))){
    stop('statList.true and statList.null are not matched')
  }

  if (method=='standard'){
  pvals=purrr::map2(statList.true, statList.null, function(true_stat, null_stat) {
    sum(abs(null_stat) >= abs(true_stat)) / (length(null_stat) + 1)
  })} else if (method=='split_pos_neg'){
  pvals=purrr::map2(statList.true, statList.null, function(true_stat, null_stat) {
    ifelse(true_stat >= 0,
           (sum(null_stat >= true_stat) + 1) / (sum(null_stat >= 0) + 1),
           (sum(null_stat <= true_stat) + 1) / (sum(null_stat < 0) + 1))
  })  
  }
  return(pvals)
}




cor2p <- function(r,n){
  t <- (r*sqrt(n-2))/sqrt(1-r^2)
  p <- 2*(1 - pt(abs(t),(n-2)))
  return(p)
}

