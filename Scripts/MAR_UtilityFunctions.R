

estMARcomm_lv_overdisp <- function(Y, nr_replicates, X = NULL, n_d = 2, n_lv=3, nu=1,
                          family=NULL, phiStart=NULL, sigmaStart=NULL,
                          silent = TRUE, get_joint_precision=FALSE, 
                          a_i_Random=TRUE, c_i_Random=TRUE, # random effect on a_i and c_i?
                          a_i_All_Equal=FALSE, c_i_All_Equal=FALSE){ # a_i and/or c_i equal for all species
  
  # familyCode:
  # 1: Poisson
  # 2: Negative binomial quadratic
  #    Var = mu(1+mu/phi) = mu + mu^2/phi (glmmTMB:nbinom2)
  #    R: dnbinom(yobs, mu =  mu, size = phi)
  # 3: linear parameterization, Var = mu(1+phi) (glmmTMB:nbinom1)
  #    R: dnbinom(yobs, mu =  mu, size = mu/phi)
  if (is.null(family)) stop("'family' not specified")
  
  if (a_i_Random && a_i_All_Equal) stop("a_i : a_i_Random==TRUE and a_i_All_Equal==TRUE ???")
  if (c_i_Random && c_i_All_Equal) stop("a_i : c_i_Random==TRUE and c_i_All_Equal==TRUE ???")
  if (ncol(Y) != sum(nr_replicates)) stop("mismatch between number of columns of Y and sum of nr_replicates")
  
  familyCode <- NULL
  if (family=="Poisson") familyCode = 1
  if (family=="NB_quad") familyCode = 2
  if (family=="NB_lin") familyCode = 3
  if (is.null(familyCode)) stop("Unknown family")
  
  n_s <- nrow(Y)
  n_t <- length(nr_replicates)
  
  prms <- list()
  map <- list()
  
  #a_i_hat <- matrix(apply(log(obs_N+1)-log(nu), 1, function(x) ar.yw(x, order.max = 1, aic=FALSE)$ar), ncol=n_s, nrow=1)
  a_i_hat <- rep(0.5, n_s)
  c_i_hat <- as.numeric(rowMeans(log(Y+1), na.rm = TRUE)*(1-a_i_hat) -log(nu))
  
  # random effects 
  if (a_i_Random){
    prms$a_i_bar <- mean(a_i_hat)
    prms$U_a_i <- a_i_hat - mean(a_i_hat)
    prms$log_sigma_a <- log(0.5)
  }
  if (c_i_Random){
    prms$c_i_bar <- mean(c_i_hat)
    prms$U_c_i = c_i_hat - mean(c_i_hat)
    prms$log_sigma_c <- log(sd(c_i_hat))
  }
  
  # no random effects
  if (!a_i_Random){
    prms$a_i_bar <- a_i_hat
    prms$U_a_i <- rep(0, n_s)
    prms$log_sigma_a <- 1
    map$U_a_i <- factor(rep(NA, length(prms$U_a_i)))
    map$log_sigma_a <- factor(rep(NA, length(prms$log_sigma_a)))
    # one common parameter across species
    if (a_i_All_Equal){
      prms$a_i_bar <- mean(a_i_hat)
      map$a_i_bar <- factor(rep(1, length(prms$a_i_bar)))
    }
  }
  if (!c_i_Random){
    prms$c_i_bar <- c_i_hat
    prms$U_c_i = rep(0, ns)
    prms$log_sigma_c <- 1
    map$U_c_i <- factor(rep(NA, length(prms$U_c_i)))
    map$log_sigma_c <- factor(rep(NA, length(prms$log_sigma_c)))
    # one common parameter across species
    if (c_i_All_Equal){
      prms$c_i_bar <- mean(c_i_hat)
      map$c_i_bar <- factor(rep(1, length(prms$c_i_bar)))
    }
  }
  
  if (is.null(sigmaStart)){
    prms$log_sigma <- rep(0.3, n_s)
  }else{
    prms$log_sigma <- rep(sigmaStart, n_s)
  }
  
  year_index <- rep.int(1:length(nr_replicates), times = nr_replicates)
  Xstart <- t(aggregate(t(Y), list(year_index), mean)[,-1])
  Xstart[is.na(Xstart)] <- 0
  Xstart <- log(Xstart+1 - log(nu))
  prms$X <- Xstart
  
  # first fit a model with independent env effects among species
  prms$bU = rep(0, n_s*n_lv - sum(((1:n_lv)-1)))
  prms$U = matrix(0, nrow=n_lv, ncol=n_t-1)
  
  prms$log_phi <- ifelse(is.null(phiStart), 2, phiStart)
  
  data_list <- list(
    obs_N=Y,
    n_r=nr_replicates,
    nu=nu,
    familyCode = as.integer(familyCode),
    yearOK = as.integer(!apply(Y, 2, function(x) any(is.na(x))))
  )
  
  if (family=="Poisson"){
    map$log_phi <- factor(rep(NA, length(prms$log_phi)))
  }
  
  map$U <- factor(rep(NA, length(prms$U)))
  map$bU <- factor(rep(NA, length(prms$bU)))
  #map$log_sigma_c <- factor(rep(NA, length(prms$log_sigma_c)))
  #map$log_sigma <- factor(rep(NA, length(prms$log_sigma)))
  #map$log_sigma <- rep(as.factor(1, length(prms$log_sigma))) ### OBS OBS
  
  cat("Model data:", "\n")
  print(str(data_list))
  #obj <- TMB::MakeADFun(data_list, prms, DLL = "tmb_MAR1_comm_lv",
  #                      random= c("X", "U_a_i", "U_c_i", "U"), silent = FALSE)
  
  library(TMB)
  dll_loaded <- ("tmb_MAR1_comm_lv_replicates" %in% names(getLoadedDLLs()))
  if (!dll_loaded) dyn.load("tmb_MAR1_comm_lv_replicates")
  
  ################################################################################
  #                               STEP 1                                         #
  ################################################################################
  #prms$a_i_bar <- 0
  #prms$c_i_bar <- 0
  #prms$log_sigma_a <- 2
  #prms$log_sigma_c <- 2
  #prms$U_a_i <- rep(0, n_s)
  #prms$U_c_i <- rep(0, n_s)
  #prms$log_sigma <- rep(1, n_s)
  sta <- proc.time()
  cat("Step 1 start: ", "\t")
  obj <- TMB::MakeADFun(data_list, prms, DLL = "tmb_MAR1_comm_lv_replicates",
                        random= c("X", "U_a_i", "U_c_i"), silent = silent, map = map)
  
  opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr,
                control=list(iter.max=1000, eval.max=1500))
  sto <- proc.time()
  tim <- sto-sta
  cat("Time used: ", round(tim[3],2), "secs", "\t")
  if (opt$convergence==0){
    cat(" Model convergence", "\n")
  }else{
    stop("Model did not converge: stopping")
  }
  
  
  ################################################################################
  #                               STEP 2                                         #
  ################################################################################
  
  cat("Step 2 start: ", "\t")
  
  map_all_params <- function(parameters){
    fill_vals <- function(x,vals){rep(as.factor(vals), length(x))}
    map_use <- list()
    for (i in 1:length(parameters)) {
      map_use[[i]] <- fill_vals(parameters[[i]],NA)
      names(map_use)[[i]] <- names(parameters)[[i]]
    }
    map_use
  }
  
  #resids <- obj$env$report()$eps
  #fa <- factanal(resids, factors=3)
  #fa <- try(factanal(resids, factors = n_lv, 
  #                   scores = "regression"), silent = T)
  #U.hat <- fa$loadings
  #bu.hat <- fa$scores
  
  map_use <- map_all_params(prms)
  map_use$U <- NULL
  map_use$bU <- NULL
  map_use$log_sigma <- factor(rep(1, length(prms$log_sigma)))
  params_use <- obj$env$parList(opt$par)
  params_use$log_sigma = rep(log(0.05), length(prms$log_sigma))
  params_use$bU = rnorm(length(params_use$bU), 0,0.001)
  params_use$U = matrix(rnorm(prod(dim(params_use$U)), 0,1), nrow=nrow(params_use$U), ncol=ncol(params_use$U))
  # Fit the model
  sta2 <- proc.time()
  obj2 <- TMB::MakeADFun(data_list,params_use,random=c("U"),DLL="tmb_MAR1_comm_lv_replicates",map=map_use, silent = silent)
  #obj <- TMB::MakeADFun(data,params_use,DLL=DLL_use,map=map_use)
  #TMB::newtonOption(obj2,smartsearch=FALSE)
  opt2 <- nlminb(obj2$par,obj2$fn,obj2$gr, control=list(iter.max=1000, eval.max=1500))
  #rep2 <- TMB::sdreport(obj2)
  sto <- proc.time()
  tim <- sto-sta2
  cat("Time used: ", round(tim[3],2), "secs", "\t")
  if (opt2$convergence==0){
    cat(" Model convergence", "\n")
  }else{
    stop("Model did not converge: stopping")
  }
  
  
  ################################################################################
  #                               STEP 3                                         #
  ################################################################################
  
  cat("Step 3 start: ", "\t")
  map_use <- list()
  map_use$log_phi <- factor(rep(NA, length(prms$log_phi)))
  map_use$log_sigma <- factor(rep(1, length(prms$log_sigma)))
  #map_use$log_sigma <- factor(rep(NA, 1))
  params_use <- obj2$env$parList(opt2$par)
  #params_use$log_sigma = log(0.01)
  sta3 <- proc.time()
  obj3 <- TMB::MakeADFun(data_list,params_use,random=c("U", "X", "U_a_i", "U_c_i"),
                         DLL="tmb_MAR1_comm_lv_replicates", map=map_use, silent = silent)
  #obj <- TMB::MakeADFun(data,params_use,DLL=DLL_use,map=map_use)
  #TMB::newtonOption(obj3,smartsearch=FALSE)
  opt3 <- nlminb(obj3$par,obj3$fn,obj3$gr, control=list(iter.max=1000, eval.max=1500))
  #rep3 <- TMB::sdreport(obj3)
  sto <- proc.time()
  tim <- sto-sta3
  cat("Time used: ", round(tim[3],2), "secs", "\t")
  tmb_obj <- obj3
  tmb_opt <- opt3
  
  nlminb_loops <- 5
  newton_loops <- 2
  
  
  sta5 <- proc.time()
  
  if(tmb_opt$convergence!=0){
    cat("Model did not converge:", opt$message, "\n")
    cat("Loops start: ", "\t")
    #sta5 <- proc.time()
    if (nlminb_loops > 1) {
      if (!silent) cat("running extra nlminb loops\n")
      for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
        temp <- tmb_opt[c("iterations", "evaluations")]
        tmb_opt <- stats::nlminb(
          start = tmb_opt$par, objective = tmb_obj$fn, gradient = tmb_obj$gr)#,
        #control = .control, lower = lim$lower, upper = lim$upper)
        tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + temp[["iterations"]]
        tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + temp[["evaluations"]]
        if (tmb_opt$convergence==0) break
      }
    }
    
    
    if (newton_loops > 0) {
      if (!silent) cat("running newtonsteps\n")
      for (i in seq_len(newton_loops)) {
        g <- as.numeric(tmb_obj$gr(tmb_opt$par))
        h <- stats::optimHess(tmb_opt$par, fn = tmb_obj$fn, gr = tmb_obj$gr)
        tmb_opt$par <- tmb_opt$par - solve(h, g)
        tmb_opt$objective <- tmb_obj$fn(tmb_opt$par)
        if (tmb_opt$convergence==0) break
      }
    }
    
    if (nlminb_loops > 1) {
      if (!silent) cat("running extra nlminb loops\n")
      for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
        temp <- tmb_opt[c("iterations", "evaluations")]
        tmb_opt <- stats::nlminb(
          start = tmb_opt$par, objective = tmb_obj$fn, gradient = tmb_obj$gr)#,
        #control = .control, lower = lim$lower, upper = lim$upper)
        tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + temp[["iterations"]]
        tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + temp[["evaluations"]]
        if (tmb_opt$convergence==0) break
      }
    }
  }
  sto <- proc.time()
  tim <- sto-sta5
  cat("Time used: ", tim[3], "secs")
  
  opt3 <- tmb_opt
  obj3 <- tmb_obj 
  
  if (opt3$convergence==0){
    cat(" Model convergence", "\n")
  }else{
    cat("Model did not converge:", opt$message, "stopping", "\n")
    cat("Will still proceed to return parameter estimates: ", "\n")
  }
  
  ################################################################################
  #                       ADREPORT, RETURN RESULTS                               #
  ################################################################################
  
  
  
  out_structure <- structure(list(
    data       = data_list))
  
  sta4 <- proc.time()
  cat("ADREP start: ", "\t")
  rep <- TMB::sdreport(obj3, getJointPrecision = get_joint_precision)
  sto <- proc.time()
  tim <- sto-sta4
  cat("Time used: ", tim[3], "secs", "\n")
  
  
  rep.rep <- summary(rep, "report", p.value=TRUE)
  rep.ran <- summary(rep, "random", p.value=TRUE)
  rep.fix <- summary(rep, "fixed", p.value=TRUE)
  eps <- obj3$report()
  
  conv <- get_convergence_diagnostics(rep)
  
  out <- c(out_structure, list(
    model      = opt3,
    sd_report  = rep,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    pos_def_hessian = rep$pdHess))
  cat("\n")
  cat("SANITY CHECKS:", "\n")
  sanityCheck <- sanity(out, se_ratio = 10000, gradient_thresh = 0.001)
  out <- c(out, sanity=sanityCheck)
  cat("\n")
  
  TMBAIC=function(opt, p=2, n=Inf){
    k = length(opt[["par"]])
    if( all(c("par","objective") %in% names(opt)) ) negloglike = opt[["objective"]]
    if( all(c("par","value") %in% names(opt)) ) negloglike = opt[["value"]]
    Return = p*k + 2*negloglike + 2*k*(k+1)/(n-k-1)
    return( Return )
  }
  
  aic_corr_model <- TMBAIC(opt3)
  aic_ind_noise_model <- TMBAIC(opt)
  aicd <- data.frame(aic_ind_noise_model, aic_corr_model)
  
  tim <- sto-sta
  cat("Total time used: ", round(tim[3],2), "secs", "\n")
  cat("\n")
  print(aicd)
  return(list(data=data_list, out = out, n_lv = n_lv, aic=aicd, eps=eps,
              rep.rep=rep.rep, rep.ran=rep.ran, rep.fix=rep.fix, opt=opt3))
  
}



estMARcomm_lv <- function(Y, X = NULL, n_d = 2, n_lv=3, nu=1,
                          family=NULL, phiStart=NULL, sigmaStart=NULL,
                          silent = TRUE, get_joint_precision=FALSE, 
                          a_i_Random=TRUE, c_i_Random=TRUE, # random effect on a_i and c_i?
                          a_i_All_Equal=FALSE, c_i_All_Equal=FALSE){ # a_i and/or c_i equal for all species
  
  # familyCode:
  # 1: Poisson
  # 2: Negative binomial quadratic
  #    Var = mu(1+mu/phi) = mu + mu^2/phi (glmmTMB:nbinom2)
  #    R: dnbinom(yobs, mu =  mu, size = phi)
  # 3: linear parameterization, Var = mu(1+phi) (glmmTMB:nbinom1)
  #    R: dnbinom(yobs, mu =  mu, size = mu/phi)
  if (is.null(family)) stop("'family' not specified")
  familyCode <- NULL
  if (family=="Poisson") familyCode = 1
  if (family=="NB_quad") familyCode = 2
  if (family=="NB_lin") familyCode = 3
  if (is.null(familyCode)) stop("Unknown family")
  
  n_s <- nrow(Y)
  n_t <- ncol(Y)
  
  prms <- list()
  map <- list()
  
  a_i_hat <- matrix(apply(log(obs_N+1)-log(nu), 1, function(x) ar.yw(x, order.max = 1, aic=FALSE)$ar), ncol=n_s, nrow=1)
  c_i_hat <- as.numeric(rowMeans(log(obs_N+1))*(1-a_i_hat) -log(nu))
  
  prms$a_i_bar <- mean(a_i_hat)
  prms$c_i_bar <- mean(c_i_hat)
  prms$log_sigma_a <- log(sd(a_i_hat))
  prms$log_sigma_c <- log(sd(c_i_hat))
  if (is.null(sigmaStart)){
    prms$log_sigma <- rep(0.3, n_s)
  }else{
    prms$log_sigma <- rep(sigmaStart, n_s)
  }
  
  
  
  prms$U_a_i <- a_i_hat - mean(a_i_hat)
  prms$X = log(obs_N+1)-log(nu)
  # first fit a model with independent env effects among species
  prms$bU = rep(0, n_s*n_lv - sum(((1:n_lv)-1)))
  prms$U = matrix(0, nrow=n_lv, ncol=n_t-1)
  
  prms$U_c_i = c_i_hat - mean(c_i_hat)
  prms$log_phi <- ifelse(is.null(phiStart), 2, phiStart)
  
  data_list <- list(
    obs_N=obs_N,
    nu=nu,
    familyCode = as.integer(familyCode),
    yearOK = as.integer(!apply(obs_N, 2, function(x) any(is.na(x))))
  )
  
  if (family=="Poisson"){
    map$log_phi <- factor(rep(NA, length(prms$log_phi)))
  }
  
  map$U <- factor(rep(NA, length(prms$U)))
  map$bU <- factor(rep(NA, length(prms$bU)))
  #map$log_sigma_c <- factor(rep(NA, length(prms$log_sigma_c)))
  #map$log_sigma <- factor(rep(NA, length(prms$log_sigma)))
  #map$log_sigma <- rep(as.factor(1, length(prms$log_sigma))) ### OBS OBS

  #cat("Model data:", "\n")
  #print(str(data_list))
  #obj <- TMB::MakeADFun(data_list, prms, DLL = "tmb_MAR1_comm_lv",
  #                      random= c("X", "U_a_i", "U_c_i", "U"), silent = FALSE)
  
  library(TMB)
  dll_loaded <- ("tmb_MAR1_comm_lv" %in% names(getLoadedDLLs()))
  if (!dll_loaded) dyn.load("tmb_MAR1_comm_lv")
  
  ################################################################################
  #                               STEP 1                                         #
  ################################################################################
  prms$a_i_bar <- 0
  prms$c_i_bar <- 0
  prms$log_sigma_a <- 2
  prms$log_sigma_c <- 2
  prms$U_a_i <- rep(0, n_s)
  prms$U_c_i <- rep(0, n_s)
  prms$log_sigma <- rep(1, n_s)
  sta <- proc.time()
  cat("Step 1 start: ", "\t")
  obj <- TMB::MakeADFun(data_list, prms, DLL = "tmb_MAR1_comm_lv",
                        random= c("X", "U_a_i", "U_c_i"), silent = silent, map = map)
  
  opt <- nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr,
                control=list(iter.max=1000, eval.max=1500))
  sto <- proc.time()
  tim <- sto-sta
  cat("Time used: ", round(tim[3],2), "secs", "\t")
  if (opt$convergence==0){
    cat(" Model convergence", "\n")
  }else{
    stop("Model did not converge: stopping")
  }
  
  
  ################################################################################
  #                               STEP 2                                         #
  ################################################################################
  
  cat("Step 2 start: ", "\t")
  
  map_all_params <- function(parameters){
    fill_vals <- function(x,vals){rep(as.factor(vals), length(x))}
    map_use <- list()
    for (i in 1:length(parameters)) {
      map_use[[i]] <- fill_vals(parameters[[i]],NA)
      names(map_use)[[i]] <- names(parameters)[[i]]
    }
    map_use
  }
  
  #resids <- obj$env$report()$eps
  #fa <- factanal(resids, factors=3)
  #fa <- try(factanal(resids, factors = n_lv, 
  #                   scores = "regression"), silent = T)
  #U.hat <- fa$loadings
  #bu.hat <- fa$scores
  
  map_use <- map_all_params(prms)
  map_use$U <- NULL
  map_use$bU <- NULL
  map_use$log_sigma <- factor(rep(1, length(prms$log_sigma)))
  params_use <- obj$env$parList(opt$par)
  params_use$log_sigma = rep(log(0.05), length(prms$log_sigma))
  params_use$bU = rnorm(length(params_use$bU), 0,0.001)
  params_use$U = matrix(rnorm(prod(dim(params_use$U)), 0,1), nrow=nrow(params_use$U), ncol=ncol(params_use$U))
  # Fit the model
  sta2 <- proc.time()
  obj2 <- TMB::MakeADFun(data_list,params_use,random=c("U"),DLL="tmb_MAR1_comm_lv",map=map_use, silent = silent)
  #obj <- TMB::MakeADFun(data,params_use,DLL=DLL_use,map=map_use)
  #TMB::newtonOption(obj2,smartsearch=FALSE)
  opt2 <- nlminb(obj2$par,obj2$fn,obj2$gr, control=list(iter.max=1000, eval.max=1500))
  #rep2 <- TMB::sdreport(obj2)
  sto <- proc.time()
  tim <- sto-sta2
  cat("Time used: ", round(tim[3],2), "secs", "\t")
  if (opt2$convergence==0){
    cat(" Model convergence", "\n")
  }else{
    stop("Model did not converge: stopping")
  }
  
  
  ################################################################################
  #                               STEP 3                                         #
  ################################################################################
  
  cat("Step 3 start: ", "\t")
  map_use <- list()
  map_use$log_phi <- factor(rep(NA, length(prms$log_phi)))
  map_use$log_sigma <- factor(rep(1, length(prms$log_sigma)))
  #map_use$log_sigma <- factor(rep(NA, 1))
  params_use <- obj2$env$parList(opt2$par)
  #params_use$log_sigma = log(0.01)
  sta3 <- proc.time()
  obj3 <- TMB::MakeADFun(data_list,params_use,random=c("U", "X", "U_a_i", "U_c_i"),
                         DLL="tmb_MAR1_comm_lv", map=map_use, silent = silent)
  #obj <- TMB::MakeADFun(data,params_use,DLL=DLL_use,map=map_use)
  #TMB::newtonOption(obj3,smartsearch=FALSE)
  opt3 <- nlminb(obj3$par,obj3$fn,obj3$gr, control=list(iter.max=1000, eval.max=1500))
  #rep3 <- TMB::sdreport(obj3)
  sto <- proc.time()
  tim <- sto-sta3
  cat("Time used: ", round(tim[3],2), "secs", "\t")
  tmb_obj <- obj3
  tmb_opt <- opt3
  
  nlminb_loops <- 5
  newton_loops <- 2
  
  
  sta5 <- proc.time()
  
  if(tmb_opt$convergence!=0){
    cat("Model did not converge:", opt$message, "\n")
    cat("Loops start: ", "\t")
    #sta5 <- proc.time()
    if (nlminb_loops > 1) {
      if (!silent) cat("running extra nlminb loops\n")
      for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
        temp <- tmb_opt[c("iterations", "evaluations")]
        tmb_opt <- stats::nlminb(
          start = tmb_opt$par, objective = tmb_obj$fn, gradient = tmb_obj$gr)#,
        #control = .control, lower = lim$lower, upper = lim$upper)
        tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + temp[["iterations"]]
        tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + temp[["evaluations"]]
        if (tmb_opt$convergence==0) break
      }
    }
    
    
    if (newton_loops > 0) {
      if (!silent) cat("running newtonsteps\n")
      for (i in seq_len(newton_loops)) {
        g <- as.numeric(tmb_obj$gr(tmb_opt$par))
        h <- stats::optimHess(tmb_opt$par, fn = tmb_obj$fn, gr = tmb_obj$gr)
        tmb_opt$par <- tmb_opt$par - solve(h, g)
        tmb_opt$objective <- tmb_obj$fn(tmb_opt$par)
        if (tmb_opt$convergence==0) break
      }
    }
    
    if (nlminb_loops > 1) {
      if (!silent) cat("running extra nlminb loops\n")
      for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
        temp <- tmb_opt[c("iterations", "evaluations")]
        tmb_opt <- stats::nlminb(
          start = tmb_opt$par, objective = tmb_obj$fn, gradient = tmb_obj$gr)#,
        #control = .control, lower = lim$lower, upper = lim$upper)
        tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + temp[["iterations"]]
        tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + temp[["evaluations"]]
        if (tmb_opt$convergence==0) break
      }
    }
  }
  sto <- proc.time()
  tim <- sto-sta5
  cat("Time used: ", tim[3], "secs")

  opt3 <- tmb_opt
  obj3 <- tmb_obj 

  if (opt3$convergence==0){
    cat(" Model convergence", "\n")
  }else{
    cat("Model did not converge:", opt$message, "stopping", "\n")
    cat("Will still proceed to return parameter estimates: ", "\n")
  }

################################################################################
#                       ADREPORT, RETURN RESULTS                               #
################################################################################



  out_structure <- structure(list(
    data       = data_list))
  
  sta4 <- proc.time()
  cat("ADREP start: ", "\t")
  rep <- TMB::sdreport(obj3, getJointPrecision = get_joint_precision)
  sto <- proc.time()
  tim <- sto-sta4
  cat("Time used: ", tim[3], "secs", "\n")
  
  
  rep.rep <- summary(rep, "report", p.value=TRUE)
  rep.ran <- summary(rep, "random", p.value=TRUE)
  rep.fix <- summary(rep, "fixed", p.value=TRUE)
  eps <- obj3$report()
  
  conv <- get_convergence_diagnostics(rep)
  
  out <- c(out_structure, list(
    model      = opt3,
    sd_report  = rep,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    pos_def_hessian = rep$pdHess))
  cat("\n")
  cat("SANITY CHECKS:", "\n")
  sanityCheck <- sanity(out, se_ratio = 10000, gradient_thresh = 0.001)
  out <- c(out, sanity=sanityCheck)
  cat("\n")
  
  TMBAIC=function(opt, p=2, n=Inf){
    k = length(opt[["par"]])
    if( all(c("par","objective") %in% names(opt)) ) negloglike = opt[["objective"]]
    if( all(c("par","value") %in% names(opt)) ) negloglike = opt[["value"]]
    Return = p*k + 2*negloglike + 2*k*(k+1)/(n-k-1)
    return( Return )
  }
  
  aic_corr_model <- TMBAIC(opt3)
  aic_ind_noise_model <- TMBAIC(opt)
  aicd <- data.frame(aic_ind_noise_model, aic_corr_model)
  
  tim <- sto-sta
  cat("Total time used: ", round(tim[3],2), "secs", "\n")
  cat("\n")
  print(aicd)
return(list(data=data_list, out = out, n_lv = n_lv, aic=aicd, eps=eps,
            rep.rep=rep.rep, rep.ran=rep.ran, rep.fix=rep.fix, opt=opt3))

}



#' Sanity check of model
#'
#' @param fit Fitted model from [sdmTMB()]
#' @param se_ratio SE ratio to abs(parameter values) to issue warning
#' @param gradient_thresh Gradient threshold to issue warning
#'
#' @return An invisible named list of checks
#' }

sanity <- function(fit, se_ratio = 10, gradient_thresh = 0.001) {
  
  hessian_ok <- eigen_values_ok <- gradients_ok <- se_magnitude_ok <- FALSE
  nlminb_ok <- FALSE
  
  simplify_msg <- "Try simplifying the model or adding priors"
  
  if (identical(fit$model$convergence, 0L)) {
    msg <- "Non-linear minimizer suggests successful convergence"
    cli::cli_alert_success(msg)
    nlminb_ok <- TRUE
  } else {
    msg <- "Non-linear minimizer did not converge: do not trust this model!"
    cli::cli_alert_danger(msg)
    cli::cli_alert_info(simplify_msg)
    cat("\n")
  }
  
  if (isFALSE(fit$pos_def_hessian)) {
    msg <- "Non-positive-definite Hessian matrix: model may not have converged"
    cli::cli_alert_danger(msg)
    cli::cli_alert_info(simplify_msg)
    cat("\n")
  } else {
    msg <- "Hessian matrix is positive definite"
    cli::cli_alert_success(msg)
    hessian_ok <- TRUE
  }
  
  if (isTRUE(fit$bad_eig)) {
    msg <- "Extreme or very small eigen values detected: model may not have converged"
    cli::cli_alert_danger(msg)
    cli::cli_alert_info(simplify_msg)
    cat("\n")
  } else {
    msg <- "No extreme or very small eigen values detected"
    cli::cli_alert_success(msg)
    eigen_values_ok <- TRUE
  }
  
  g <- fit$gradients
  np <- names(fit$model$par)
  for (i in seq_along(g)) {
    if (g[i] > gradient_thresh) {
      cli::cli_alert_danger(c(
        "`", np[i],
        paste0("` gradient > ", gradient_thresh)
      ))
      msg <- "See `?run_extra_optimization()`"
      cli::cli_alert_info(msg)
      msg <- "Or refit with `control = control(newton_loops = 1)`"
      cli::cli_alert_info(msg)
      cat("\n")
    }
  }
  
  if (all(g <= gradient_thresh)) {
    msg <- "No gradients with respect to fixed effects are >= "
    cli::cli_alert_success(paste0(msg, gradient_thresh))
    gradients_ok <- TRUE
  }
  
  obj <- fit$tmb_obj
  random <- unique(names(obj$env$par[obj$env$random]))
  s <- summary(fit$sd_report)
  se <- s[,"Std. Error"]
  fixed_se <- !names(se) %in% random
  se <- se[fixed_se]
  np <- names(se)
  se_na_ok <- TRUE
  for (i in seq_along(se)) {
    if (is.na(se[i])) {
      cli::cli_alert_danger(c("`", np[i], paste0("` standard error is NA")))
      #par_message(np[i])
      #cli::cli_alert_info(simplify_msg)
      cat("\n")
      se_na_ok <- FALSE
    }
  }
  if (se_na_ok) {
    msg <- "No fixed-effect standard errors are NA"
    cli::cli_alert_success(msg)
  }
  
  est <- as.list(fit$sd_report, "Estimate")
  se <- as.list(fit$sd_report, "Std. Error")
  fixed <- !(names(est) %in% random)
  est <- est[fixed]
  se <- se[fixed]
  too_big <- function(est, se) {
    if (any(!is.na(se))) {
      ratio <- se[!is.na(se)] / abs(est[!is.na(se)])
      if (any(ratio > se_ratio)) return(TRUE)
    }
  }
  se_big <- mapply(too_big, est, se)
  for (i in seq_along(se_big)) {
    if (isTRUE(se_big[[i]])) {
      msg <- paste0(
        "` standard error may be large (> ",
        se_ratio,
        "x parameter estimate)"
      )
      cli::cli_alert_danger(c("`", names(se_big)[i], msg))
      #par_message(names(se_big)[i])
      #cli::cli_alert_info(simplify_msg)
      cat("\n")
    }
  }
  if (all(unlist(lapply(se_big, is.null)))) {
    msg <- "No fixed-effect standard errors look unreasonably large"
    #cli::cli_alert_success(msg)
    se_magnitude_ok <- TRUE
  }
  
  # b <- tidy(fit, conf.int = TRUE)
  # b2 <- tidy(fit, "ran_pars", conf.int = TRUE)
  # b <- rbind(b, b2)
  # s <- grep("sigma", b$term)
  # sigmas_ok <- TRUE
  # if (length(s)) {
  #   for (i in s) {
  #     if (b$estimate[i] < 1e-3) {
  #       msg <- "` is smaller than 0.001"
  #       cli::cli_alert_danger(c("`", b$term[i], msg))
  #       par_message(b$term[i])
  #       msg <- "Consider omitting this part of the model"
  #       cli::cli_alert_info(msg)
  #       cat("\n")
  #       sigmas_ok <- FALSE
  #     }
  #   }
  # }
  # if (sigmas_ok) {
  #   msg <- "No sigma parameters are < 0.001"
  #   cli::cli_alert_success(msg)
  # }
  
  # r1 <- diff(range(fit$data[[fit$mesh$xy_cols[1]]]))
  # r2 <- diff(range(fit$data[[fit$mesh$xy_cols[2]]]))
  # r <- max(r1, r2)
  # range_ok <- TRUE
  # if ("range" %in% b$term) {
  #   if (max(b$estimate[b$term == "range"]) > r) {
  #     msg <- "A `range` parameter looks fairly large (> greatest distance in data)"
  #     cli::cli_alert_danger(msg)
  #     cli::cli_alert_info(simplify_msg)
  #     cli::cli_alert_info("Also make sure your spatial coordinates are not too big or small (e.g., work in UTM km instead of UTM m)", wrap = TRUE)
  #     cat("\n")
  #     range_ok <- FALSE
  #   } else {
  #     nr <- length(grep("range", b$term))
  #     if (nr == 1L) msg <- "Range parameter doesn't look unreasonably large"
  #     if (nr > 1L) msg <- "Range parameters don't look unreasonably large"
  #     cli::cli_alert_success(msg)
  #   }
  # }
  
  ret <- list(
    hessian_ok=hessian_ok, eigen_values_ok=eigen_values_ok, nlminb_ok=nlminb_ok, #range_ok,
    gradients_ok=gradients_ok, se_magnitude_ok=se_magnitude_ok, se_na_ok=se_na_ok#, sigmas_ok
  )
  all_ok <- all(unlist(ret))
  ret <- c(ret, all_ok = all_ok)
  invisible(ret)
}


get_convergence_diagnostics <- function(sd_report) {
  final_grads <- sd_report$gradient.fixed
  bad_eig <- FALSE
  if (!is.null(sd_report$pdHess)) {
    if (!sd_report$pdHess) {
      warning("The model may not have converged: ",
              "non-positive-definite Hessian matrix.", call. = FALSE)
    } else {
      eigval <- try(1 / eigen(sd_report$cov.fixed)$values, silent = TRUE)
      if (is(eigval, "try-error") || (min(eigval) < .Machine$double.eps * 10)) {
        warning("The model may not have converged: ",
                "extreme or very small eigen values detected.", call. = FALSE)
        bad_eig <- TRUE
      }
      if (any(final_grads > 0.01))
        warning("The model may not have converged. ",
                "Maximum final gradient: ", max(final_grads), ".", call. = FALSE)
    }
  }
  pdHess <- isTRUE(sd_report$pdHess)
  invisible(list(final_grads=final_grads, bad_eig=bad_eig, pdHess=pdHess))
}


extractEstimates <- function(m,sp_names=NULL){
  if (is.null(sp_names)){
    sp_names <- rownames(m$data)
  }
  rho_x <- m$rep.rep[grep("rho_x", rownames(m$rep.rep)),1]
  rho_e <- m$rep.rep[grep("rho_e", rownames(m$rep.rep)),1]
  DC <- m$rep.rep[grep("DC", rownames(m$rep.rep)),1]
  A <- m$rep.rep[grep("a_i", rownames(m$rep.rep)),1]
  # change sd_corr_U : bad wording
  sigma2 <- (m$rep.rep[grep("sd_corr_U", rownames(m$rep.rep)),1]) + (m$rep.rep[grep("sigma", rownames(m$rep.rep)),1])^2
  rho_x_mat <- matrix(NA, nrow=length(A), ncol=length(A), dimnames = list(sp_names,sp_names))
  rho_e_mat <- DC_mat <- rho_x_mat
  ltri <- lower.tri(rho_x_mat)
  rho_e_mat[ltri] <- rho_e
  rho_x_mat[ltri] <- rho_x
  DC_mat[ltri] <- DC
  c_i <- m$rep.rep[grep("c_i", rownames(m$rep.rep))]
  AR1_pars <- matrix(c(c_i, A, sigma2), ncol=3, dimnames = list(sp_names, c("c_i", "A", "sigma2")))
  return(list(AR1_pars = AR1_pars, rho_x=rho_x, rho_e=rho_e, DC=DC,
              rho_x_mat=rho_x_mat, rho_e_mat=rho_e_mat, DC_mat=DC_mat))
}

