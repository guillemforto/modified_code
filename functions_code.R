
#' wrapper
#' The main wrapper function
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param strataname  Name of the variable to use for stratification
#' @param varname  Name(s) of the variable(s) of interest
#' @param gn  Population size
#' @param gnh  Population size in each stratum
#' @param method  Sampling design: si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#' @param esttype  Type of estimator
#' 
#' @return Returns two dataframes: the first one gives the robust estimator and acts as a summary of the winsorisation. The second one details the weight changes of every unit.
#' @export
#' 

"wrapper" <- function(data,
                      varname = NULL,
                      strataname = NULL,
                      gn,
                      gnh = NULL,
                      est_type = c("BHR", "standard", "DT"),
                      method = c("si", "poisson", "rejective"),
                      pii = NULL,
                      di = NULL,
                      id = NULL) {
  # Conditions
  if (missing(pii) & missing(di)) {
    stop("the vector of probabilities is missing\n")
  } else if (missing(pii) & !missing(di)) {
    pii <- 1 / di
  } else if (!missing(pii) & !missing(di)) {
    warning("Warning: di is redundant. Only pii is being used\n")
  }
  if (missing(id)) {
    warning("Warning: no column is specified as an identifier, one is added automatically (id)\n")
    id <- "id"
    identifier <- c(1:nrow(data))
  } else {
    if (!(id %in% colnames(data))) {
      stop("the specified identifier is not a column name\n")
    }
    identifier <- data[[id]]
  }
  
  # initialisation data frames
  if (!missing(strataname) & !missing(gnh)) {
    df <- data.frame(stratum=numeric(),
                     est_type=character(),
                     var=character(),
                     RHT=numeric(),
                     tuning_const=numeric(),
                     HT=numeric(),
                     rel_diff=numeric(),
                     modif_weights=integer())
  } else {
    df <- data.frame(est_type=character(),
                     var=character(),
                     RHT=numeric(),
                     tuning_const=numeric(),
                     HT=numeric(),
                     rel_diff=numeric(),
                     modif_weights=integer())
  }
  
  df2 <- data.frame(matrix(0, ncol=0, nrow=nrow(data)))
  df2[,id] <- identifier
  df2$init_weight <- 1/pii
  if (!missing(strataname) & !missing(gnh)) {
    df2$stratum <- data[,strataname]
  }
  
  for (var in varname) {
    if (!missing(strataname) & !missing(gnh)) {
      # conditional bias
      bi <- strata_HTcondbiasest(data, strataname, var, gnh, method, pii, remerge=FALSE)
      # robust estimator
      RHT <- strata_robustest(data, strataname, var, gnh, method, pii)[[1]]
    } else {
      # conditional bias
      bi <- HTcondbiasest(data, var, gn, method, pii, id="none", remerge=FALSE)
      # robust estimator
      RHT <- robustest(data, var, gn, method, pii)[[1]]
    }
    
    # filling df2
    df2[,var] <- data[,var]
    df2[,paste("condbias", var, sep="_")] <- bi$condbias
    
    for (t in est_type) {
      # tuning constant
      if (t == "BHR") {
        tun_const <- tuningconst(bi$condbias)
      } else if (t == "standard") {
        tun_const <- determinconstws(pii, data[,var], bi$condbias)
      } else if (t == "DT") {
        tun_const <- determinconstwDT(pii, data[,var], bi$condbias)
      }
      # weights
      if (!missing(strataname) & !missing(gnh)) {
        new_weights <- strata_robustweights.r(data, strataname, var, gnh, method, pii, typewin=t, remerge=F)[,]
      } else {
        new_weights <- robustweights.r(data, var, gn, method, ech$piks, typewin=t, remerge=F)[,]
      }
      # HT estimator
      if (!missing(strataname) & !missing(gnh)) {
        HT <- c()
        for (i in unique(data[,strataname])) {
          HT[i] <- crossprod(data[data[,strataname]==i, var], 1/pii[data[,strataname]==i])
        }
      } else {
        HT <- crossprod(data[,var], 1/pii)
      }
      
      # relative difference + modif_weights + nb_modif_weights
      rel_diff <- round((RHT - HT) / HT * 100, 2)
      modif_weights <- diag(outer(1/pii, new_weights, Vectorize(all.equal)))
      modif_weights[modif_weights != TRUE] <- FALSE
      modif_weights <- !as.logical(modif_weights)
      if (!missing(strataname) & !missing(gnh)) {
        nb_modif_weights <- as.vector(table(modif_weights, data[,strataname])["TRUE",])
      } else {
        nb_modif_weights <- sum(modif_weights)
      }
      
      # filling df
      if (!missing(strataname) & !missing(gnh)) {
        for (i in unique(data[,strataname])) {
          df <- rbind(df, list(i, t, var, RHT[i], tun_const, HT[i], rel_diff[i], nb_modif_weights[i]), stringsAsFactors=FALSE)
        }
      } else {
        df <- rbind(df, list(t, var, RHT, tun_const, HT, rel_diff, nb_modif_weights), stringsAsFactors=FALSE)
      }
      
      # filling df2
      df2[,paste("new_weights", var, t, sep="_")] <- new_weights
      df2[,paste("modifed", var, t, sep="_")] <- modif_weights
    }
  }
  if (!missing(strataname) & !missing(gnh)) {
    colnames(df) <- c("stratum", "est_type", "var", "RHT", "tuning_const", "HT", "rel_diff", "nb_modif_weights")
  } else {
    colnames(df) <- c("est_type", "var", "RHT", "tuning_const", "HT", "rel_diff", "nb_modif_weights")
  }
  
  return(list(df, df2))
}


#' HTcondbiasest
#' The function to estimate the conditional biases for the HT estimator of a total, for the non-stratiﬁed sampling designs. For one or several given variables of interest as input, this function returns a vector with the estimates of the conditional bias of the HT estimator where each row corresponds to a sample unit.
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param varname  Name(s) of the variable(s) of interest
#' @param gn  Population size
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#' @param di  Inverse of the first order inclusion probabilities
#' @param id  Primary key to keep as an identifier of every unit
#' @param remerge  True/False to remerge the resulting conditional biases with the original data set
#'
#' @return Returns a data set with the conditional bias of each sampled unit
#' @export
#'

"HTcondbiasest" <- function (data,
                             varname = NULL,
                             gn,
                             method = c("si", "poisson", "rejective"),
                             pii = NULL,
                             di = NULL,
                             id = NULL,
                             remerge = T) {
  if (nrow(data) <= 1) {
    stop("your sample must contain at least 2 units\n")
  }
  if (sum(is.na(data[,varname])) != 0) {
    stop("at least one missing value in the variable(s) of interest\n")
  }
  if (any(data[,varname] < 0)) {
    warning("Warning: your variable(s) of interest contain negative values\n")
  }
  if (missing(method) | length(method) > 1) {
    warning("Warning: the method is not specified or multiple methods are specified.\nBy default, the method is 'si'\n")
    method = "si"
  }
  if (!(method %in% c("si", "poisson", "rejective"))) {
    stop("the name of the method is wrong\n")
  }
  if (missing(pii) & missing(di)) {
    stop("the vector of probabilities is missing\n")
  } else if (missing(pii) & !missing(di)) {
    pii <- 1 / di
  } else if (!missing(pii) & !missing(di)) {
    warning("Warning: di is redundant. Only pii is being used\n")
  }
  if (missing(gn)) {
    stop("the population size is missing\n")
  }
  if (method == "si" & !all.equal(gn, sum(1/pii))) {
    warning("Warning: the sum of the inclusion probabilities is not equal to N\n")
    print("N:, ", gn)
    print("sum(1/pii):", sum(1/pii))
  }
  if (remerge==F) {
    if (missing(id)) {
      warning("Warning: no column is specified as an identifier, one is added automatically (id)\n")
      id <- "id"
      identifier <- c(1:nrow(data))
    } else {
      if (id != "none") {
        if (!(id %in% colnames(data))) {
          stop("the specified identifier is not a column name\n")
        }
        identifier <- data[[id]]
      }
    }
  }
  
  data = data.frame(data)
  pn = nrow(data)
  index = 1:nrow(data)
  m = match(varname, colnames(data))
  if (any(is.na(m))) {
    stop("the name of the variable is wrong\n")
  }
  data2 = cbind.data.frame(data[, m], index)
  colnames(data2) = c(varname, "index")
  if (method == "si") {
    if (length(m)==1){
      bc = (pn/(pn-1))*(gn/pn-1)*(data[,m]-mean(data[,m]))
    } else {
      bc = (pn/(pn-1))*(gn/pn-1)*(data[,m]-matrix(data=colMeans(data[,m]),nrow=pn,ncol=length(m), byrow =T))
    }
  }
  if (method == "poisson") {
    bc = (1/pii-1) * data[,m]
  }
  if (method == "rejective") {
    gd = sum(1-pii)
    gb = t(1/pii-1) %*% as.matrix(data[,m]) / gd
    print(paste0("D: ", round(gd, 3)))
    print(paste0("N/D: ", round(gn/gd, 3)))
    warning("Warning: please make sure that D is large enough and N/D is bounded\n")
    bc = (1/pii - 1) * (data[,m] - t(t(gb) %*% c(pii)))
  }
  if (remerge == T) {
    result = cbind.data.frame(data, bc)
    colnames(result) = c(colnames(data), paste0("condbias", colnames(bc)))
  } else {
    result = data.frame(bc)
    if (length(m) == 1) {
      colnames(result) = c("condbias")
    } else {
      colnames(result) = c(paste0("condbias", colnames(bc)))
    }
    if (id != "none") {
      result <- cbind.data.frame(identifier, result)
      colnames(result)[1] <- id
    }
  }
  result
}




#' strata_HTcondbiasest
#' The function to estimate the conditional biases for the HT estimator of a total, for the stratified sampling designs. For one or several given variables of interest as input, this function returns a vector with the estimates of the conditional bias of the HT estimator where each row corresponds to a sample unit.
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param strataname  Name of the variable to use for stratification
#' @param varname  Name(s) of the variable(s) of interest
#' @param gnh  Population size in each stratum
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#' @param di  Inverse of the first order inclusion probabilities
#' @param id  Primary key to keep as an identifier of every unit
#' @param remerge  True/False to remerge the resulting conditional biases with the original data set
#'
#' @return Returns a dataframe with the conditional bias of each sampled unit
#' @export


"strata_HTcondbiasest" <- function(data,
                                   strataname = NULL,
                                   varname = NULL,
                                   gnh,
                                   method = c("si","poisson","rejective"),
                                   pii = NULL,
                                   di = NULL,
                                   id = NULL,
                                   remerge = T) {
  if (missing(pii) & missing(di)) {
    stop("the vector of probabilities is missing\n")
  } else if (missing(pii) & !missing(di)) {
    pii <- 1 / di
  } else if (!missing(pii) & !missing(di)) {
    warning("Warning: di is redundant. Only pii is being used\n")
  }
  if (missing(gnh)) {
    stop("the population size vector is missing\n")
  }
  if (missing(strataname) | is.null(strataname)) {
    stop("no variable name to use for stratification has been specified\n")
  }
  if (any(table(data[,strataname]) <= 1)) {
    stop("each stratum must contain at least 2 units\n")
  }
  if (remerge==F) {
    if (missing(id)) {
      warning("Warning: no column is specified as an identifier, one is added automatically (id)\n")
      id <- "id"
      identifier <- c(1:nrow(data))
    } else {
      if (!(id %in% colnames(data))) {
        stop("the specified identifier is not a column name\n")
      }
      identifier <- data[[id]]
    }
  }

  data = data.frame(data)
  index = 1:nrow(data)
  m = match(varname, colnames(data))
  if (any(is.na(m)))
    stop("the name of the variable is wrong\n")
  ms = match(strataname, colnames(data))
  if (any(is.na(ms)))
    stop("the name of the strata is wrong\n")
  x1 = data.frame(unique(data[,ms]))
  bc = matrix(0, nrow=nrow(data), ncol=length(m))
  cgn = 0

  for (i in 1:nrow(x1)) {
    print(paste("Stratum ", i, ":", sep=""))
    datastr = as.data.frame(data[(data[,ms]==i),m])
    colnames(datastr) = colnames(data)[m]
    nh = nrow(as.data.frame(datastr))
    piisrt = pii[(data[,ms]==i)]
    resint = HTcondbiasest(data=datastr, varname=varname, gn=gnh[i], method=method, pii=piisrt, id="none", remerge=F)
    bc[(cgn+1):(cgn+nh),] = as.matrix(resint)
    # bc[(cgn+1):(cgn+nh)]=nh[i]/(nh[i]-1)*(gnh[i]/nh-1)*(y-mean(y))
    cgn = cgn+nh
  }
  if (remerge==T) {
    result = cbind.data.frame(data, bc)
    colnames(result) = c(colnames(data), colnames(resint))
  } else {
    result = data.frame(bc)
    if (length(m) == 1) {
      colnames(result) = c("condbias")
    } else {
      colnames(result) = c(paste0("condbias", varname))
    }
    result <- cbind.data.frame(identifier, result)
    colnames(result)[1] <- id
  }
  print("Done!")
  result
}



#' robustest
#'
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param varname  Name(s) of the variable(s) of interest
#' @param gn  Population size
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#'
#' @return Computes the robust estimator of Beaumont et al.(2013) using the conditional bias and the minmax criterion to compute the tuning constant
#' @export
#'

"robustest" <- function (data, 
                         varname = NULL, 
                         gn,
                         method = c("si","poisson","rejective"), 
                         pii) {
  if (method == "si" & !all.equal(gn, sum(1/pii))) {
    warning("Warning: the sum of the inclusion probabilities is not equal to N\n")
    print("N:, ", gn)
    print("sum(1/pii):", sum(1/pii))
  }
  if (missing(method) | length(method) > 1) {
    warning("Warning: the method is not specified or multiple methods are specified.\nBy default, the method is 'si'\n")
    method = "si"
  }
  if (!(method %in% c("si", "poisson", "rejective")))
    stop("the name of the method is wrong\n")
  if (method %in% c("poisson", "rejective") & missing(pii))
    stop("the vector of probabilities is missing\n")
  if  (missing(gn))
    stop("the population size is missing\n")
  data = data.frame(data)
  pn = nrow(data)
  index = 1:nrow(data)
  m = match(varname, colnames(data))
  if (any(is.na(m)))
    stop("the name of the variable is wrong\n")
  data2 = cbind.data.frame(data[, m], index)
  colnames(data2) = c(varname, "index")
  if (any(is.na(data[,m])))
    stop("Missing values for some y-variables\n")
  if (length(m) == 1) {
    htestim = crossprod(data[,m],1/pii)
    htbc = HTcondbiasest(data=data, varname = varname, gn=gn, method=method, pii=pii, id="none")[,c(seq(1+length(data),length(data)+length(varname)))]
    result = htestim - (min(htbc)+max(htbc))/2
  } else {
    htestim = (1/pii) %*% as.matrix(data[,m])
    htbc = HTcondbiasest(data=data, varname=varname, gn=gn, method=method, pii=pii, id="none")[,c(seq(1+length(data),length(data)+length(varname)))]
    result = htestim-(apply(X=htbc,MARGIN = 2,min )+apply(X =htbc,MARGIN = 2,max ))/2
  }
  colnames(result) = c(paste0("RHT_", colnames(htbc)))
  result
}




#' strata_robustest
#'
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param strataname  Name of the variable to use for stratification
#' @param varname  Name(s) of the variable(s) of interest
#' @param gnh  Population size in each stratum
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#'
#' @return Computes the robust estimator of Beaumont et al.(2013) by stratum using the conditional bias and the minmax criterion to compute the tuning constant
#' @export

"strata_robustest" <- function(data, 
                               strataname = NULL, 
                               varname = NULL, 
                               gnh,
                               method = c("si","poisson","rejective"), 
                               pii) {
  if (missing(gnh)) {
    stop("the population size vector is missing\n")
  }
  if (missing(strataname) | is.null(strataname)) {
    stop("no variable name to use for stratification has been specified\n")
  }
  data = data.frame(data)
  index = 1:nrow(data)
  m = match(varname, colnames(data))
  if (any(is.na(m))) {
    stop("the name of the variable is wrong\n")
  }
  ms = match(strataname, colnames(data))
  if (any(is.na(ms))) {
    stop("the name of the strata is wrong\n")
  }
  x1 = data.frame(unique(data[,ms]))
  rht = matrix(0, nrow=nrow(x1), ncol=length(m))

  for (i in 1:nrow(x1)) {
    print(paste("Stratum ", i, ":", sep=""))
    datastr = as.data.frame(data[(data[,ms]==i),m])
    colnames(datastr) = colnames(data)[m]
    nh = nrow(as.data.frame(datastr))
    piisrt = pii[(data[,ms]==i)]
    rht[i,] = robustest(data=datastr, varname=varname, gn=gnh[i], method=method, pii=piisrt)
  }
  result = cbind.data.frame(rht)
  colnames(result) = c(paste0("RHT_", colnames(datastr)))
  rownames(result) = c(paste("Stratum", 1:nrow(x1)))
  print("Done!")
  result
}




#' tuningconst
#'
#' @param bi  Conditional biases of the variable of interest
#' @importFrom stats optimize
#' @return Computes the tuning constant associated to the Beaumont et al (2013) estimator
#' @export
#'

tuningconst = function(bi) {
  tuningconstBHR = function(c, bi) {
    res = sapply(c, hub.psi, x=bi)
    return(abs(colSums(res-bi) + 0.5*(min(bi)+max(bi))))
  }
  resultat <- optimize(tuningconstBHR, interval=c(min(abs(bi)), max(abs(bi))), bi=bi, maximum=FALSE, tol=.Machine$double.eps^0.9)
  return(resultat$minimum)
}


#' robustweights
#'
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param varname  Name(s) of the variable(s) of interest
#' @param gn  Population size
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#' @param typewin  Winsorized estimator : Beaumont et al., Standard or Dalen-Tambay
#' @param maxit  Maximum number of iterations for the research of the minimum
#' @param remerge True/False to remerge the conditional bias with the original data set
#'
#' @return Computes the robust weights associated to the winsorized estimator
#' @export

"robustweights.r"<- function (data, 
                              varname = NULL,
                              gn,
                              method = "si",
                              pii,
                              typewin = "BHR",
                              maxit = 10000,
                              remerge = T) {
  if (missing(method) | length(method) > 1) {
      warning("Warning: the method is not specified or multiple methods are specified.\nBy default, the method is 'si'\n")
      method = "si"
  }
  if (!(method %in% c("si", "poisson", "rejective"))) {
    stop("the name of the method is wrong\n")
  }
  if (missing(typewin)) {
    warning("Warning: the type of winsorization is not specified; by default, the method is BHR\n")
    typewin = "BHR"
  }
  if (!(typewin %in% c("BHR", "standard", "DT"))) {
    stop("the name of the type of winsorization  is wrong\n")
  }
  if (missing(pii)) {
    stop("the vector of probabilities is missing\n")
  }
  data = data.frame(data)
  pn = nrow(data)
  index = 1:nrow(data)
  m = match(varname, colnames(data))
  if (any(is.na(m))) {
    stop("the name of the variable is wrong\n")
  }
  data2 = cbind.data.frame(data[, m], index)
  colnames(data2) = c(varname, "index")
  bc = HTcondbiasest(data, varname, gn, method, pii, id="none", remerge=F)
  if (length(m) != ncol(bc)) {
    stop("the number of conditional bias is different from the number of variable of study\n")
  }
  if (typewin == "BHR") {
    if (length(m) == 1){
      tc = tuningconst(bc[,1])
      ditilde = (1/pii) - (bc[,1]-hub.psi(bc[,1],b=tc))/data[,m]
      ditilde[is.nan(ditilde)] = 1/pii[is.nan(ditilde)]
    } else {
      tc = apply(bc, MARGIN = 2, tuningconst)
      ditilde = (1/pii) - (bc-mapply(x=bc,b=tc,hub.psi))/data[,m]
      ##Some weight might be Nan, since some data are equals to 0
      ##Since they don't contribute to the total, an arbitrary can be add, by default the orginal weight
      for (j in 1:length(bc)) {
        ditilde[is.nan(ditilde[,j]),j] = 1/pii[is.nan(ditilde[,j])]
      }
    }
  }
  if (typewin == "standard") {
    if (length(m)==1){
      ctws = determinconstws(pii=pii, x=data[,m], bi=bc, maxit)
      ditilde = (1/pii) * apply(cbind(data[,m],ctws*pii),MARGIN=1,min)/data[,m]
      ditilde[is.nan(ditilde)] = 1/pii[is.nan(ditilde)]
    } else {
      tc = mapply(determinconstws, x=data[,m], bi=bc, MoreArgs = list(pii=pii, maxit=maxit))
      df = data[,m]
      df[data[,m] - t(c(tc)%*%t(c(pii)))>0] = t(c(tc)%*%t(c(pii)))[data[,m]-t(c(tc)%*%t(c(pii)))>0]
      ditilde = (1/pii)*df/data[,m]
      ##Some weight might be Nan, since some data are equals to 0
      ##Since they don't contribute to the total, an arbitrary can be add, by default the orginal weight
      for (j in 1:length(bc)) {
        ditilde[is.nan(ditilde[,j]),j] = 1/pii[is.nan(ditilde[,j])]
      }
    }
  }
  if (typewin == "DT") {
    if (length(m)==1){
      copt = determinconstwDT(pii=pii, x=data[,m], bi=bc, maxit)
      ditilde = 1+(1/pii-1) * apply(cbind(data[,m],copt*pii),MARGIN=1,min)/data[,m]
      ditilde[is.nan(ditilde)] = 1/pii[is.nan(ditilde)]
    } else {
      tc = mapply(determinconstwDT, x=data[,m], bi=bc, MoreArgs=list(pii=pii, maxit=maxit))
      df = data[,m]
      df[data[,m] - t(c(tc)%*%t(c(pii)))>0] = t(c(tc)%*%t(c(pii)))[data[,m]-t(c(tc)%*%t(c(pii)))>0]
      ditilde = 1+(1/pii-1)*df/data[,m]
      ##Some weight might be Nan, since some data are equals to 0
      ##Since they don't contribute to the total, an arbitrary can be add, by default the orginal weight
      for (j in 1:length(bc)) {
        ditilde[is.nan(ditilde[,j]),j] = 1/pii[is.nan(ditilde[,j])]
      }
    }
  }
  if (remerge) {
    result = cbind.data.frame(data,ditilde)
    colnames(result)=c(colnames(data),paste0("weights",typewin,colnames(ditilde)))
  } else {
    result =data.frame(ditilde)
    if (length(m)==1) {
      colnames(result) = c(paste0("weights",typewin))
    } else {
      colnames(result) = c(paste0("weights",typewin,colnames(ditilde)))
    }
  }
  result
}




#' strata_robustweights
#'
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param strataname  Name of the variable to use for stratification
#' @param varname  Name(s) of the variable(s) of interest
#' @param gnh  Population size in each stratum
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#' @param typewin  Winsorized estimator : Beaumont et al., Standard or Dalen-Tambay
#' @param maxit  Maximum number of iterations for the research of the minimum
#' @param remerge  True/False to remerge the conditional bias with the original data set
#'
#'
#' @return Computes the robust weights associated to the winsorized estimator
#' @export

"strata_robustweights.r" <- function (data, 
                                      strataname = NULL,
                                      varname = NULL, 
                                      gnh, 
                                      method = c("si", "poisson", "rejective"), 
                                      pii,
                                      typewin = "BHR",
                                      maxit = 10000,
                                      remerge = T) {
  if (missing(gnh)) {
    stop("the population size vector is missing\n")
  }
  if (missing(strataname) | is.null(strataname)) {
    stop("no variable name to use for stratification has been specified\n")
  }
  data = data.frame(data)
  index = 1:nrow(data)
  m = match(varname, colnames(data))
  if (any(is.na(m))) {
    stop("the name of the variable is wrong\n")
  }
  ms = match(strataname, colnames(data))
  if (any(is.na(ms))) {
    stop("the name of the strata is wrong\n")
  }
  if (missing(typewin)) {
    warning("Warning: the type of winsorization is not specified; by default, the method is BHR\n")
    typewin = "BHR"
  }
  if (!(typewin %in% c("BHR", "standard", "DT"))) {
    stop("the name of the type of winsorization  is wrong\n")
  }
  x1 = data.frame(unique(data[,ms]))
  matw = matrix(0,nrow=nrow(data),ncol=length(m))
  cgn = 0
  for (i in 1:nrow(x1)) {
    datastr = as.data.frame(data[(data[,ms]==i),m])
    colnames(datastr) = colnames(data)[m]
    nh = nrow(as.data.frame(datastr))
    piisrt = pii[(data[,ms]==i)]
    resint = robustweights.r(data=datastr, varname=varname, gn=gnh[i], method=method , pii=piisrt, typewin=typewin, maxit=maxit, remerge=F)
    matw[(cgn+1):(cgn+nh),] = as.matrix(resint)
    # bc[(cgn+1):(cgn+nh)]=nh[i]/(nh[i]-1)*(gnh[i]/nh-1)*(y-mean(y))
    cgn=cgn+nh
  }
  if (remerge) {
    result = cbind.data.frame(data,matw)
    colnames(result) = c(colnames(data), colnames(resint))
    result
  } else {
    result = cbind.data.frame(matw)
    colnames(result) = c(colnames(resint))
    result
  }
}




#' determinconstws
#'
#' @param pii  First order inclusion probabilities
#' @param x  Variable(s) of interest
#' @param bi  Conditional bias
#' @param maxit  Maximum number of iterations for the research of the minimum
#' @importFrom stats uniroot
#' @return Computes the winsorisation constant associated to the standard winsorized estimator
#' @export
#'

determinconstws = function(pii, x, bi, maxit=10000) {
  if (max(bi) < -min(bi)) {
    stop("the condition for unicity is not satisfied\n")
  }
  di <- 1 / pii
  functiontws = function(a, di, x, bi) {
    testpos = function(x, tconst=1.345){
      res = rep(0, length(x))
      res[(x-tconst) > 0] = (x-tconst)[(x-tconst) > 0]
      return(res)
    }
    rest = -colSums(sapply(a, testpos, x=as.matrix(di*x)))
    copt = rest + 0.5*(min(bi)+max(bi))
    return(copt)
  }
  resultat <- uniroot(functiontws, c(0,max(di*x)), check.conv=FALSE, tol=.Machine$double.eps^10, maxiter=maxit, trace=0, di=di, x=x, bi=bi)
  return(resultat$root)
}




#' determinconstwDT
#'
#' @param pii  First order inclusion probabilities
#' @param x  Variable(s) of interest
#' @param bi  Conditional bias
#' @param maxit  Maximum number of iterations for the research of the minimum
#' @importFrom stats uniroot
#' @return Computes the winsorisation constant associated to the Dalen-Tambay winsorized estimator
#' @export
#'

determinconstwDT = function(pii, x, bi, maxit=10000) {
  if (max(bi) < -min(bi)) {
    stop("the condition for unicity is not satisfied\n")
  }
  di <- 1 / pii
  functiontDT = function(a, di, x, bi){
    testpos = function(x, tconst=1.345){
      res = rep(0, length(x))
      res[(x-tconst)>0] = (x-tconst)[(x-tconst)>0]
      return(res)
    }
    rest = -colSums(((di-1)/di) * sapply(a, testpos, x=as.matrix(di*x)))
    copt = rest + 0.5*(min(bi)+max(bi))
    return(copt)
  }
  resultat <- uniroot(functiontDT, c(0, max(di*x)), check.conv=FALSE, tol=.Machine$double.eps^10, maxiter=maxit, trace=0, di=di, x=x, bi=bi)
  return(resultat$root)
}



hub.psi <- function(x, b = 1.345) {
  psi <- ifelse(abs(x) <= b, x, sign(x) * b)
  return(psi = psi)
}



