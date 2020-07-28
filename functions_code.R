

#' HTcondbiasest
#' The function to estimate the conditional biases for the HT estimator of a total, for the non-stratiﬁed sampling designs. For one or several given variables of interest as input, this function returns a vector with the estimates of the conditional bias of the HT estimator where each row corresponds to a sample unit.
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param varname  Name(s) of the variable(s) of interest
#' @param gn  Population size
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#' @param di  Inverse of the first order inclusion probabilities
#' @param pkey  Primary key to keep as an identifier of every unit
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
                             pkey = NULL,
                             remerge = T) {
  if (nrow(data) <= 1) {
    stop("your sample must contain at least 2 units\n")
  }
  if (missing(method)) {
    warning("Warning: the method is not specified; by default, the method is si\n")
    method = "si"
  }
  if (!(method %in% c("si", "poisson", "rejective"))) {
    stop("the name of the method is wrong\n")
  }
  if (missing(pii) & missing(di)) {
    stop("the vector of probabilities is missing\n")
  } else if (missing(pii) & !missing(di)) {
    pii <- 1 / di
  }
  if (missing(gn)) {
    stop("the population size is missing\n")
  }
  if (method == "si" & !all.equal(gn, sum(1/pii))) {
    warning("Warning: the sum of the probabilities is not equal to N\n")
  }
  if (remerge==F) {
    if (missing(pkey)) {
      warning("Warning: no column is specified as a primary key, one is added automatically (id)\n")
      id <- c(1:nrow(data))
      pkey <- "id"
    } else {
      if (!(pkey %in% colnames(data))) {
        stop("the specified primary key is not a column name\n")
      }
      id <- data[[pkey]]
    }
  }

  data = data.frame(data)
  pn = nrow(data)
  index = 1:nrow(data)
  m = match(varname, colnames(data))
  if (any(is.na(m)))
    stop("the name of the variable is wrong\n")
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
    bc = (1/pii -1) * data[,m]
  }
  if (method == "rejective") {
    gd = sum(1-pii)
    gb = t(1/pii-1)%*%as.matrix(data[,m])/gd
    bc = (1/pii -1) * (data[,m]-t(t(gb)%*%c(pii)))
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
    result <- cbind.data.frame(id, result)
    colnames(result)[1] <- pkey
  }
  result
}




#' strata_HTcondbiasest
#' The function to estimate the conditional biases for the HT estimator of a total, for the stratified sampling designs. For one or several given variables of interest as input, this function returns a vector with the estimates of the conditional bias of the HT estimator where each row corresponds to a sample unit.
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param strataname  Name of the variable to use for stratification
#' @param varname  Name(s) of the variable(s) of interest
#' @param gnh  Lopulation size in each stratum
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#' @param di  Inverse of the first order inclusion probabilities
#' @param pkey  Primary key to keep as an identifier of every unit
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
                                   pkey = NULL,
                                   remerge = T) {
  if (missing(pii) & missing(di)) {
    stop("the vector of probabilities is missing\n")
  } else if (missing(pii) & !missing(di)) {
    pii <- 1 / di
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
    if (missing(pkey)) {
      warning("Warning: no column is specified as a primary key, one is added automatically (id)\n")
      id <- c(1:nrow(data))
      pkey <- "id"
    } else {
      if (!(pkey %in% colnames(data))) {
        stop("the specified primary key is not a column name\n")
      }
      id <- data[[pkey]]
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
    datastr = as.data.frame(data[(data[,ms]==i),m])
    colnames(datastr) = colnames(data)[m]
    datastr[pkey] <- c(1:nrow(datastr))
    nh = nrow(as.data.frame(datastr))
    piisrt = pii[(data[,ms]==i)]
    resint = HTcondbiasest(data=datastr, varname=varname, gn=gnh[i], method=method, pii=piisrt, pkey=pkey, remerge=F)
    resint[,pkey] <- NULL
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
    result <- cbind.data.frame(id, result)
    colnames(result)[1] <- pkey
  }
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
#' @return Computes the robust estimator of Beaumont and al.(2013) using the conditional bias and the minmax criterion to compute the tuning constant
#' @export
#'

"robustest" <- function (data, varname=NULL, gn,
                         method=c("si","poisson","rejective"), pii) {
  if (gn != sum(1/pii)) {
    warning("the sum of the inclusion probabilities is not equal to N\n")
  }
  if (missing(method)) {
    warning("the method is not specified; by default, the method is si\n")
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
    htbc = HTcondbiasest(data=data, varname = varname, gn=gn, method = method, pii=pii)[,c(seq(1+length(data),length(data)+length(varname)))]
    result = htestim-(min(htbc)+max(htbc))/2
  } else {
    htestim = (1/pii) %*% as.matrix(data[,m])
    htbc = HTcondbiasest(data=data, varname = varname, gn=gn, method = method, pii=pii)[,c(seq(1+length(data),length(data)+length(varname)))]
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
#' @return Computes the robust estimator of Beaumont and al.(2013) using the conditional bias and the minmax criterion to compute the tuning constant
#' @export

"strata_robustest" <- function(data, strataname=NULL, varname=NULL, gnh,
                               method=c("si","poisson","rejective"), pii) {
  if (missing(gnh))
    stop("the population size vector is missing\n")
  if (missing(strataname) | is.null(strataname))
    stop("no variable name to use for stratification has been specified\n")
  data = data.frame(data)
  index = 1:nrow(data)
  m = match(varname, colnames(data))
  if (any(is.na(m)))
    stop("the name of the variable is wrong\n")
  ms = match(strataname, colnames(data))
  if (any(is.na(ms)))
    stop("the name of the strata is wrong\n")
  x1 = data.frame(unique(data[,ms]))
  rht = matrix(0, nrow=nrow(x1), ncol=length(m))

  for (i in 1:nrow(x1)) {
    datastr = as.data.frame(data[(data[,ms]==i),m])
    colnames(datastr) = colnames(data)[m]
    nh = nrow(as.data.frame(datastr))
    piisrt = pii[(data[,ms]==i)]
    rht[i,] = robustest(data=datastr, varname=varname, gn=gnh[i], method=method, pii=piisrt)
  }
  result = cbind.data.frame(rht)
  colnames(result) = c(paste0("RHT_", colnames(datastr)))
  result
}




#' tuningconst
#'
#' @param bi  Conditional bias
#' @param tailleseq  Maximum number of iterations for the research of the minimum
#' @importFrom stats optimize
#' @return Computes the robust weight associated to the Beaumont et al (2013) estimator
#' @export
#'

tuningconst = function(bi, tailleseq=1000) {
  tuningconstBDR = function(c, bi, tailleseq=1000) {
    res = sapply(c, hub.psi, x=bi)
    return(abs(colSums(res-bi) + 0.5*(min(bi)+max(bi))))
  }
  resultat <- optimize(tuningconstBDR, interval=c(min(abs(bi)),max(abs(bi))), bi=bi, maximum = FALSE, tol = .Machine$double.eps^0.9)
  return(resultat$minimum)
}




#' weightswin
#'
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param varname  Name(s) of the variable(s) of interest
#' @param gn  Population size
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#' @param typewin  Winsorized estimator : Beaumont et al., Standard or Dalen-Tambay
#' @param tailleseq  Maximum number of iterations for the research of the minimum
#' @param remerge True/False to remerge the conditional bias with the original data set
#'
#' @return Computes the robust weight associated to the winsorized estimator
#' @export

"weightswin.r"<-
  function (data, varname = NULL, gn, method="si", pii, typewin="BHR",
            tailleseq=10000, remerge=T) {
    if (missing(method)) {
      warning("the method is not specified; by default, the method is si\n")
      method = "si"
    }
    if (!(method %in% c("si", "poisson", "rejective")))
      stop("the name of the method is wrong\n")
    if (missing(typewin)) {
      warning("the type of winsorization is not specified; by default, the method is BHR\n")
      typewin = "BHR"
    }
    if (!(typewin %in% c("BHR", "standard", "DT")))
      stop("the name of the type of winsorization  is wrong\n")
    if (missing(pii))
      stop("the vector of probabilities is missing\n")
    data = data.frame(data)
    pn = nrow(data)
    index = 1:nrow(data)
    m = match(varname, colnames(data))
    if (any(is.na(m)))
      stop("the name of the variable is wrong\n")
    data2 = cbind.data.frame(data[, m], index)
    colnames(data2) = c(varname, "index")
    bc = HTcondbiasest(data, varname, gn, method, pii,remerge=F)
    if(length(m) != ncol(bc))
      stop("the number of conditionnal bias is different from the number of variable of study\n")
    if (typewin == "BHR") {
      if (length(m) == 1){
        tc = tuningconst(bc[,1],tailleseq=tailleseq)
        ditilde = (1/pii)-(bc[,1]-hub.psi(bc[,1],b=tc))/data[,m]
        ditilde[is.nan(ditilde)] = 1/pii[is.nan(ditilde)]
      } else {
        tc = apply(bc, MARGIN = 2, tuningconst, tailleseq=tailleseq)
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
        ctws = determinconstws(di=1/pii,x=data[,m],bc,tailleseq)
        ditilde = (1/pii)*apply(cbind(data[,m],ctws*pii),MARGIN=1,min)/data[,m]
        ditilde[is.nan(ditilde)] = 1/pii[is.nan(ditilde)]
      } else {
        tc = mapply(determinconstws,x=data[,m],bc=bc, MoreArgs = list(di = 1/pii,tailleseq=tailleseq))
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
        copt = determinconstwDT(di=1/pii,x=data[,m],bc,tailleseq)
        ditilde = 1+(1/pii-1)*apply(cbind(data[,m],copt*pii),MARGIN=1,min)/data[,m]
        ditilde[is.nan(ditilde)] = 1/pii[is.nan(ditilde)]
      } else {
        tc = mapply(determinconstwDT,x=data[,m],bc=bc, MoreArgs = list(di = 1/pii,tailleseq=tailleseq))
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




#' strata_weightswin
#'
#' @param data  Original data set with the sample, the variable(s) of interest
#' @param strataname  Name of the variable to use for stratification
#' @param varname  Name(s) of the variable(s) of interest
#' @param gnh  Population size in each stratum
#' @param method  Sampling design : si for simple random sampling, poisson, rejective
#' @param pii  First order inclusion probabilities
#' @param typewin  Winsorized estimator : Beaumont et al., Standard or Dalen-Tambay
#' @param tailleseq  Maximum number of iterations for the research of the minimum
#' @param remerge  True/False to remerge the conditional bias with the original data set
#'
#'
#' @return  Returns the robust weights associated to the winsorized estimator
#' @export

"strata_weightswin.r" <-
  function (data, strataname = NULL,varname = NULL, gnh, method = c("si",
                                                                    "poisson", "rejective"), pii,typewin="BHR",tailleseq=10000,remerge=T)
  {
    if  (missing(gnh))
      stop("the population size vector is missing\n")
    if (missing(strataname) | is.null(strataname))
      stop("no variable name to use for stratification has been specified\n")
    data = data.frame(data)
    index = 1:nrow(data)
    m = match(varname, colnames(data))
    if (any(is.na(m)))
      stop("the name of the variable is wrong\n")
    ms = match(strataname, colnames(data))
    if (any(is.na(ms)))
      stop("the name of the strata is wrong\n")
    if (missing(typewin)) {
      warning("the type of winsorization is not specified; by default, the method is BHR\n")
      typewin = "BHR"
    }
    if (!(typewin %in% c("BHR", "standard", "DT")))
      stop("the name of the type of winsorization  is wrong\n")
    x1 = data.frame(unique(data[,ms]))
    matw=matrix(0,nrow=nrow(data),ncol=length(m))
    cgn=0
    for (i in 1:nrow(x1)) {
      datastr=as.data.frame(data[(data[,ms]==i),m])
      colnames(datastr)=colnames(data)[m]
      nh=nrow(as.data.frame(datastr))
      piisrt=pii[(data[,ms]==i)]
      resint=weightswin.r(data=datastr, varname =varname , gn=gnh[i], method =method , pii=piisrt,typewin=typewin,tailleseq=tailleseq,remerge = F)
      matw[(cgn+1):(cgn+nh),]=as.matrix(resint)

      # bc[(cgn+1):(cgn+nh)]=nh[i]/(nh[i]-1)*(gnh[i]/nh-1)*(y-mean(y))
      cgn=cgn+nh

    }
    if(remerge){
      result = cbind.data.frame(data,matw)
      colnames(result)=c(colnames(data),colnames(resint))
      result
    }else{
      result = cbind.data.frame(matw)
      colnames(result)=c(colnames(resint))
      result
    }
  }




#' determinconstws
#'
#' @param pii  First order inclusion probabilities
#' @param x  Variable(s) of interest
#' @param bi  Conditional bias
#' @param tailleseq  Maximum number of iterations for the research of the minimum
#' @importFrom stats uniroot
#' @return Computes the robust weights associated to the standard winsorized estimator
#' @export
#'

determinconstws = function(pii, x, bi, tailleseq) {
  if (max(bi) < -min(bi)) {
    stop("the condition for existence is not satisfied\n")
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
  resultat <- uniroot(functiontws, c(0,max(di*x)), check.conv=FALSE, tol=.Machine$double.eps^10, maxiter=tailleseq, trace=0, di=di, x=x, bi=bi)
  return(resultat$root)
}




#' determinconstwDT
#'
#' @param pii  First order inclusion probabilities
#' @param x  Variable(s) of interest
#' @param bi  Conditional bias
#' @param tailleseq  Maximum number of iterations for the research of the minimum
#' @importFrom stats uniroot
#' @return Computes the robust weight associated to the Dalen-Tambay winsorized estimator
#' @export
#'

determinconstwDT = function(pii, x, bi, tailleseq=1000) {
  if (max(bi) < -min(bi)) {
    stop("the condition for existence is not satisfied\n")
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
  resultat <- uniroot(functiontDT, c(0, max(di*x)), check.conv=FALSE, tol=.Machine$double.eps^10, maxiter=tailleseq, trace=0, di=di, x=x, bi=bi)
  return(resultat$root)
}




hub.psi <- function(x, b = 1.345) {
  psi <- ifelse(abs(x) <= b, x, sign(x) * b)
  der.psi <- ifelse(abs(x) <= b, 1, 0)
  return(psi = psi)
}
