kappa_irr  <- function (ratings, exact = FALSE, detail = FALSE) 
{
  ratings <- as.matrix(na.omit(ratings))
  ns <- nrow(ratings)
  nr <- ncol(ratings)
  lev <- levels(as.factor(ratings))
  for (i in 1:ns) {
    frow <- factor(ratings[i, ], levels = lev)
    if (i == 1) 
      ttab <- as.numeric(table(frow))
    else ttab <- rbind(ttab, as.numeric(table(frow)))
  }
  ttab <- matrix(ttab, nrow = ns)
  agreeP <- sum((apply(ttab^2, 1, sum) - nr)/(nr * (nr - 1))/ns)
  if (!exact) {
    method <- "Fleiss' Kappa for m Raters"
    chanceP <- sum(apply(ttab, 2, sum)^2)/(ns * nr)^2
  }
  else {
    method <- "Fleiss' Kappa for m Raters (exact value)"
    for (i in 1:nr) {
      rcol <- factor(ratings[, i], levels = lev)
      if (i == 1) 
        rtab <- as.numeric(table(rcol))
      else rtab <- rbind(rtab, as.numeric(table(rcol)))
    }
    rtab <- rtab/ns
    chanceP <- sum(apply(ttab, 2, sum)^2)/(ns * nr)^2 - sum(apply(rtab, 
                                                                  2, var) * (nr - 1)/nr)/(nr - 1)
  }
  value <- (agreeP - chanceP)/(1 - chanceP)
  if (!exact) {
    pj <- apply(ttab, 2, sum)/(ns * nr)
    qj <- 1 - pj
    varkappa <- (2/(sum(pj * qj)^2 * (ns * nr * (nr - 1)))) * 
      (sum(pj * qj)^2 - sum(pj * qj * (qj - pj)))
    SEkappa <- sqrt(varkappa)
    u <- value/SEkappa
    p.value <- 2 * (1 - pnorm(abs(u)))
    if (detail) {
      pj <- apply(ttab, 2, sum)/(ns * nr)
      pjk <- (apply(ttab^2, 2, sum) - ns * nr * pj)/(ns * 
                                                       nr * (nr - 1) * pj)
      kappaK <- (pjk - pj)/(1 - pj)
      varkappaK <- 2/(ns * nr * (nr - 1))
      SEkappaK <- sqrt(varkappaK)
      uK <- kappaK/SEkappaK
      p.valueK <- 2 * (1 - pnorm(abs(uK)))
      tableK <- as.table(round(cbind(kappaK, uK, p.valueK,pj,pjk), 
                               digits = 3))
      rownames(tableK) <- lev
      colnames(tableK) <- c("Kappa", "z", "p.value","po","pc")
    }
  }
  if (!exact) {
    if (!detail) {
      rval <- list(method = method, subjects = ns, raters = nr, 
                   irr.name = "Kappa", value = value, po = agreeP, pc=chanceP)
    }
    else {
      rval <- list(method = method, subjects = ns, raters = nr, 
                   irr.name = "Kappa", value = value, detail = tableK, po = agreeP, pc=chanceP)
    }
    rval <- c(rval, stat.name = "z", statistic = u, p.value = p.value, po = agreeP, pc=chanceP)
  }
  else {
    rval <- list(method = method, subjects = ns, raters = nr, 
                 irr.name = "Kappa", value = value, po = agreeP, pc=chanceP)
  }
  class(rval) <- "irrlist"
  if(!exact)  return(list(P.o = agreeP, P.c = chanceP, Kappa = (agreeP-chanceP)/(1-chanceP)))
  if(exact)  return(list(P.o = agreeP, P.c = chanceP, Kappa = (agreeP-chanceP)/(1-chanceP)))
}
kappa_DescTools <- function (x, method = c("Fleiss", "Conger", "Light"), conf.level = NA) 
{
  if (is.matrix(x)) 
    x <- as.data.frame(x)
  x <- na.omit(x)
  ns <- nrow(x)
  nr <- ncol(x)
  lev <- levels(factor(unlist(x)))
  xx <- do.call(cbind, lapply(x, factor, levels = lev))
  ttab <- apply(Abind(lapply(as.data.frame(xx), function(z) Dummy(z, 
                                                                  method = "full", levels = seq_along(lev))), along = 3), 
                c(1, 2), sum)
  agreeP <- sum((rowSums(ttab^2) - nr)/(nr * (nr - 1))/ns)
  switch(match.arg(method, choices = c("Fleiss", "Conger", 
                                       "Light")), Fleiss = {
                                         chanceP <- sum(colSums(ttab)^2)/(ns * nr)^2
                                         value <- (agreeP - chanceP)/(1 - chanceP)
                                         pj <- colSums(ttab)/(ns * nr)
                                         qj <- 1 - pj
                                         varkappa <- (2/(sum(pj * qj)^2 * (ns * nr * (nr - 1)))) * 
                                           (sum(pj * qj)^2 - sum(pj * qj * (qj - pj)))
                                         SEkappa <- sqrt(varkappa)
                                         ci <- value + c(1, -1) * qnorm((1 - conf.level)/2) * 
                                           SEkappa
                                       }, Conger = {
                                         rtab <- apply(Abind(lapply(as.data.frame(t(xx)), function(z) Dummy(z, 
                                                                                                            method = "full", levels = seq_along(lev))), along = 3), 
                                                       c(1, 2), sum)
                                         rtab <- rtab/ns
                                         chanceP <- sum(colSums(ttab)^2)/(ns * nr)^2 - sum(apply(rtab, 
                                                                                                 2, var) * (nr - 1)/nr)/(nr - 1)
                                         value <- (agreeP - chanceP)/(1 - chanceP)
                                         ci <- c(NA, NA)
                                       }, Light = {
                                         m <- DescTools::PairApply(x, DescTools::CohenKappa, symmetric = TRUE)
                                         value <- mean(m[upper.tri(m)])
                                         levlen <- length(lev)
                                         for (nri in 1:(nr - 1)) for (nrj in (nri + 1):nr) {
                                           for (i in 1:levlen) for (j in 1:levlen) {
                                             if (i != j) {
                                               r1i <- sum(x[, nri] == lev[i])
                                               r2j <- sum(x[, nrj] == lev[j])
                                               if (!exists("dis")) dis <- r1i * r2j else dis <- c(dis, 
                                                                                                  r1i * r2j)
                                             }
                                           }
                                           if (!exists("disrater")) disrater <- sum(dis) else disrater <- c(disrater, 
                                                                                                            sum(dis))
                                           rm(dis)
                                         }
                                         B <- length(disrater) * prod(disrater)
                                         chanceP <- 1 - B/ns^(choose(nr, 2) * 2)
                                         varkappa <- chanceP/(ns * (1 - chanceP))
                                         SEkappa <- sqrt(varkappa)
                                         ci <- value + c(1, -1) * qnorm((1 - conf.level)/2) * 
                                           SEkappa
                                       })
  if (is.na(conf.level)) {
    res <- list(P.o = agreeP, P.c=chanceP,
                Kappa = value)
  }
  else {
    res <- list(P.o = agreeP, P.c=chanceP,
                Kappa = value, lwr.ci = ci[1], upr.ci = ci[2])
  }
  return(res)
} 
kappa_Conger_irr <- function(ratings){
  ratings <- as.matrix(na.omit(ratings))
  N <- nrow(ratings)
  M <- ncol(ratings)
  for (i in 1:(M - 1)) for (j in (i + 1):M) {
    if ((i == 1) & (j == (i + 1))){
      res = kappa2_HB(ratings[, c(i, j)], weight = "u")
      P.os = res$P.o
      P.cs = res$P.c
      Kappas <- res$Kappa
    }
    
    else {
      res = kappa2_HB(ratings[, c(i, j)], weight = "u")
      P.os = c(P.os,res$P.o)
      P.cs = c(P.cs,res$P.c)
      Kappas = c(Kappas,res$Kappa)
    }
  }
  Kappas = mean(Kappas)
  P.o = mean(P.os)
  P.c = mean(P.cs)
  return(list(P.o = P.o, P.c=P.c,Kappa =(P.o-P.c)/(1-P.c)))
} #Exact Kappa (Conger) from irr package
kappa_Fleiss_irr <- function(ratings){
  N = nrow(ratings)
  M = ncol(ratings)
  mysum=sum(ratings)
  if(mysum==0  || mysum==N*M){
    outp = list(P.o = 1, P.c=1,Kappa = 1)
  }
  else{
    n.mat = t(apply(ratings,1,function(x){n0=sum(x==0); n1=sum(x==1); c(n0,n1)}))
    p.mat = apply(n.mat,2,sum)/N/M
    P.mat = apply(n.mat,1,function(x){x[1]*(x[1]-1)+x[2]*(x[2]-1)})/M/(M-1)
    
    P.o = mean(P.mat)
    P.c = sum(p.mat^2)
    outp=list(P.o = P.o, P.c=P.c,Kappa = (P.o-P.c)/(1-P.c))
  }
  outp
} #Non-exact Kappa (Fleiss) from irr package
kappa_Light_Desc <- function(ratings) {kappa_DescTools(ratings,method="Light",conf.level = NA)}
kappa_HB_exact <- function(ratings){
  N = nrow(ratings)
  M = ncol(ratings)
  mysum=sum(ratings)
  if(mysum==0  || mysum==N*M){
    outp = list(P.o = 1, P.c=1,Kappa = 1)
  }
  else{
    n.mat = t(apply(ratings,1,function(x){n0=sum(x==0); n1=sum(x==1); c(n0,n1)}))
    m.mat = t(apply(ratings,2,function(x){n0=sum(x==0); n1=sum(x==1); c(n0,n1)}))
    
    
    m.mat = m.mat/N
    P.c <- sum(apply(n.mat, 2, sum)^2)/(N * M)^2 - sum(apply(m.mat, 
                                                             2, var) * (M - 1)/M)/(M - 1)
    
    p.mat = apply(n.mat,2,sum)/N/M
    P.mat = apply(n.mat,1,function(x){x[1]*(x[1]-1)+x[2]*(x[2]-1)})/M/(M-1)
    
    P.o = mean(P.mat)
    outp=list(P.o = P.o, P.c=P.c,Kappa = (P.o-P.c)/(1-P.c))
  }
  outp
} #faster than exact
kappa2_HB <- function(ratings, weight = c("unweighted", "equal", "squared"),sort.levels=FALSE){
  ratings <- as.matrix(na.omit(ratings))
  if (is.character(weight)) 
    weight = match.arg(weight)
  ns <- nrow(ratings)
  nr <- ncol(ratings)
  if (nr > 2) {
    stop("Number of raters exeeds 2. Try kappam.fleiss or kappam.light.")
  }
  r1 <- ratings[, 1]
  r2 <- ratings[, 2]
  if ((is.numeric(r1)) | (is.numeric(r2))) 
    sort.levels <- TRUE
  if (!is.factor(r1)) 
    r1 <- factor(r1)
  if (!is.factor(r2)) 
    r2 <- factor(r2)
  if (length(levels(r1)) >= length(levels(r2))) {
    lev <- c(levels(r1), levels(r2))
  }
  else {
    lev <- c(levels(r2), levels(r1))
  }
  if (sort.levels) 
    lev <- sort(lev)
  lev <- lev[!duplicated(lev)]
  r1 <- factor(ratings[, 1], levels = lev)
  r2 <- factor(ratings[, 2], levels = lev)
  ttab <- table(r1, r2)
  nc <- ncol(ttab)
  if (is.numeric(weight)) 
    w <- 1 - (weight - min(weight))/(max(weight) - min(weight))
  else if (weight == "equal") 
    w <- (nc - 1):0/(nc - 1)
  else if (weight == "squared") 
    w <- 1 - (0:(nc - 1))^2/(nc - 1)^2
  else w <- c(1, rep(0, nc - 1))
  wvec <- c(sort(w, decreasing = FALSE), w[2:length(w)])
  nw <- length(w)
  weighttab <- matrix(0, nrow = nw, ncol = nw)
  for (i in 1:nw) {
    weighttab[i, ] <- wvec[(nw - (i - 1)):(2 * nw - i)]
  }
  agreeP <- sum(ttab * weighttab)/ns
  tm1 <- apply(ttab, 1, sum)
  tm2 <- apply(ttab, 2, sum)
  eij <- outer(tm1, tm2)/ns
  chanceP <- sum(eij * weighttab)/ns
  return(list(P.o = agreeP, P.c=chanceP,Kappa = (agreeP-chanceP)/(1-chanceP)))
} # 2raters only
kappa_HB_old <- function(ratings){
  N = dim(ratings)[1]
  M = dim(ratings)[2]
  mysum=sum(ratings)
  if(mysum==0  || mysum==N*M){
    outp = list(P.o = 1, P.c=1,Kappa = 1)
  }
  else{
    n.mat = t(apply(ratings,1,function(x){n0=sum(x==0); n1=sum(x==1); c(n0,n1)}))
    m.mat = t(apply(ratings,2,function(x){n0=sum(x==0); n1=sum(x==1); c(n0,n1)}))
    
    N0 = sum(ratings==0)
    N1 = sum(ratings==1)
    
    N0 = apply(n.mat,2,sum)[1]
    N1 = apply(n.mat,2,sum)[2]
    
    num = choose(N0,2) +  choose(N1,2) - sum(apply(n.mat,1,function(x){choose(x[1],2)+choose(x[2],2)}))- 
      sum(apply(m.mat,1,function(x){choose(x[1],2)+choose(x[2],2)}))
    
    den = N*(N-1)*M*(M-1)/2
    P.c = num/den
    P.o = mean(apply(ratings,1,function(x){n=length(x); n0=sum(x==0);(choose(n0,2)+choose(n-n0,2))/choose(n,2)}))
    
    outp=list(P.o= P.o, P.c=P.c, Kappa = (P.o-P.c)/(1-P.c))
  }
  outp
}
kappa_HB_new <- function(ratings){
  N = dim(ratings)[1]
  M = dim(ratings)[2]
  mysum=sum(ratings)
  if(mysum==0  || mysum==N*M){
    outp = list(P.o = 1, P.c=1,Kappa = 1)
  }
  else{
    n.mat = t(apply(ratings,1,function(x){n0=sum(x==0); n1=sum(x==1); c(n0,n1)}))
    m.mat = t(apply(ratings,2,function(x){n0=sum(x==0); n1=sum(x==1); c(n0,n1)}))
    
    N0 = sum(ratings==0)
    N1 = sum(ratings==1)
    
    # N0 = apply(n.mat,2,sum)[1]
    # N1 = apply(n.mat,2,sum)[2]
    
    num = choose(N0,2) +  choose(N1,2) 
    
    den = choose(N*M,2)
    P.c = num/den
    P.o = mean(apply(ratings,1,function(x){n=length(x); n0=sum(x==0);(choose(n0,2)+choose(n-n0,2))/choose(n,2)}))
    
    outp=list(P.o= P.o, P.c=P.c, Kappa = (P.o-P.c)/(1-P.c))
  }
  outp
}
kappa_ordinal <- function(ratings){
  N = dim(ratings)[1]
  M = dim(ratings)[2]
  mysum=sum(ratings)
  if(mysum==0  || mysum==N*M){
    outp = list(P.o = 1, P.c=1,Kappa = 1)
  }
  else{
    outputnames <- c("esteta","estseeta","estsigma2u","estsesigma2u","estsigma2v","estsesigma2v","estkappam","estsekappam")
    dat = melt(as.matrix(ratings))
    colnames(dat)=c("case","rater","cat")
    m1 <- clmm(as.factor(cat) ~ -1 + (1|case) + (1|rater), link = "probit", threshold = "flexible", data=dat)
    
    ### save number of obs, raters, and items
    numobs <- m1$dims$n
    numitems <- as.numeric(m1$dims$nlev.re[1])
    numraters <- as.numeric(m1$dims$nlev.re[2])
    
    ### save estimate and standard error of the intercept
    esteta <- -summary(m1)$coef[1]
    estseeta <- summary(m1)$coef[2]
    
    ### Save variance
    estsigma2u <- VarCorr(m1)$case[1,1]
    estsigma2v <- VarCorr(m1)$rater[1,1]
    
    ### hessian gives the variance of the standard deviation -- use delta method to convert
    thehessian <- summary(m1)$Hessian
    estsesigma2u <- sqrt(4 *  estsigma2u * diag(solve(thehessian))[2])
    estsesigma2v <- sqrt(4 *  estsigma2v * diag(solve(thehessian))[3])
    
    ### kappa
    estrho <- estsigma2u/(estsigma2u + estsigma2v + 1)
    estetastar <- esteta/sqrt(estsigma2u + estsigma2v + 1)
    
    integrand=function(z){
      pnorm(z*sqrt(estrho)/sqrt(1-estrho)) * (1-pnorm(z*sqrt(estrho)/sqrt(1-estrho))) * dnorm(z)
    }
    result <- integrate(integrand,lower=-100,upper=100)
    estkappam <- 1-4*result$value
    
    integrand_po=function(z){
      pnorm((z*sqrt(estrho)+estetastar)/sqrt(1-estrho)) * (1-pnorm((z*sqrt(estrho)+estetastar)/sqrt(1-estrho))) * dnorm(z)
    }
    result <- integrate(integrand_po,lower=-100,upper=100)
    est_po <- 1-2*result$value
    est_pc <- pnorm(estetastar)
    
    outp=list(P.o= est_po, P.c=est_pc, Kappa = (est_po-est_pc)/(1-est_pc))
  }
  outp
}
get_ns <- function(df){
  sapply(df,function(x){length(unique(x))})
}
df_to_ratings <- function(df,inter = TRUE, subset = FALSE){  ### verify stacking is done on remaining columns also
  #change df (subject,rater,time,y) to long format for Kappa calculation
  #df contains subject, rater, time, response
  #inter=FALSE for intrarater
  n = sapply(df,function(x){length(unique(x))})
  if(inter & subset) ratings = df %>% 
      pivot_wider(names_from=rater,values_from=y) %>%
      dplyr::select((ncol(df)-1):(ncol(df)+n["rater"]-2))
  if(inter & !subset) ratings = df %>% 
      pivot_wider(names_from=rater,values_from=y)
  if(!inter & subset) ratings = df %>% 
      pivot_wider(names_from=time,values_from=y) %>%
      dplyr::select((ncol(df)-1):(ncol(df)+n["time"]-2))
  if (!inter & !subset) ratings = df %>% 
      pivot_wider(names_from=time,values_from=y)
  return(ratings)
} 
ratings_to_df <- function(ratings, inter = TRUE) {
  #ratings should contain columns subject, time, remaining K columns for raters/time{ 
  if(inter) df = melt(ratings,id=c("subject","time")) %>% 
      dplyr::rename(rater=variable,y=value)
  if(!inter) df = melt(ratings,id=c("subject","rater")) %>% 
      dplyr::rename(time=variable,y=value)
  return(df)
}
f.sim.ind = function(beta0=4, I=20, J=5, K=2, mean_sub=0, sigma_sub=2, mean_rater = 0,
                     sigma_rater = 2, mean_time=0, sigma_time=1){
  mu_subject = rnorm(I,mean_sub,sigma_sub) #subject effect on probit scale
  mu_rater = rnorm(J,mean_rater,sigma_rater) #rater effect on probit scale
  mu_time = rnorm(K,mean_time,sigma_time) #time effect on probit scale
  
  mu_overall = array(dim=c(I,J,K))
  prob = array(dim=c(I,J,K))
  y = array(dim=c(I,J,K))
  
  
  for(i in 1:I){
    for(j in 1:J){
      for(k in 1:K){
        mu_overall[i,j,k] = beta0 + mu_subject[i] + mu_rater[j] + mu_time[k]
        prob[i,j,k] = pnorm(mu_overall[i,j,k])
        y[i,j,k] = rbinom(1,1,prob[i,j,k])
      }
    }
  }
  df = melt(y)
  colnames(df) = c("subject","rater","time","y")
  return(list(df=df, prob=prob))
}
f.sim.partial = function(beta0=4, I=20, J=5, K=2, mean_sub=0, sigma_sub=2, mean_rater = 0,
                         sigma_rater = 2, rho_rater=0.3){
  Rho_rater <- diag(1-rho_rater,K,K)+matrix(rho_rater,K,K)
  Sigma_rater <- diag(sigma_rater,K) %*% Rho_rater %*% diag(sigma_rater,K)
  
  mu_subject = rnorm(I,mean_sub,sigma_sub) #subject effect on probit scale
  mu_rater = mvrnorm(J,rep(mean_rater,K),Sigma_rater)  #rater effect on probit scale
  
  
  mu_overall = array(dim=c(I,J,K))
  prob = array(dim=c(I,J,K))
  y = array(dim=c(I,J,K))
  
  
  for(i in 1:I){
    for(j in 1:J){
      for(k in 1:K){
        mu_overall[i,j,k] = beta0 + mu_subject[i] + mu_rater[j,k]
        prob[i,j,k] = pnorm(mu_overall[i,j,k])
        y[i,j,k] = rbinom(1,1,prob[i,j,k])
      }
    }
  }
  df = melt(y)
  colnames(df) = c("subject","rater","time","y")
  return(df)
}
f.sim.full = function(beta0=4, I=20, J=5, K=2, mean_sub=0, sigma_sub=2, mean_rater = 0,
                      sigma_rater = 2, rho_rater=0.5, mean_time=0, sigma_time = 2,
                      rho_time=0.9){
  
  Rho_rater <- diag(1-rho_rater,J,J)+matrix(rho_rater,J,J)
  Sigma_rater <- diag(sigma_rater,J) %*% Rho_rater %*% diag(sigma_rater,J)
  
  Rho_time <- diag(1-rho_time,K,K)+matrix(rho_time,K,K)
  Sigma_time <- diag(sigma_time,K) %*% Rho_time %*% diag(sigma_time,K)
  
  mu_subject = rnorm(I,mean_sub,sigma_sub) #subject effect on probit scale
  mu_rater = mvrnorm(I,rep(mean_rater,J),Sigma_rater)  #rater effect on probit scale
  
  mu_time = array(dim=c(I,J,K))
  for(i in 1:I){
    mu_time[i,,]=mvrnorm(J,rep(mean_time,K),Sigma_time)}
  
  mu_overall = array(dim=c(I,J,K))
  prob = array(dim=c(I,J,K))
  y = array(dim=c(I,J,K))
  
  
  for(i in 1:I){
    for(j in 1:J){
      for(k in 1:K){
        mu_overall[i,j,k] = beta0 + mu_subject[i] + mu_rater[j,k] + mu_time[i,j,k]
        prob[i,j,k] = pnorm(mu_overall[i,j,k])
        y[i,j,k] = rbinom(1,1,prob[i,j,k])
      }
    }
  }
  df = melt(y)
  colnames(df) = c("subject","rater","time","y")
  return(df)
}
f.sim.ind.gen = function(beta=c(rep(0.1,3),rep(-0.1,3)), X = matrix(rep(1,6),ncol=6), I=32, J=3, K=2, mean_sub=0, 
                         sigma_sub=2, mean_rater = 0,sigma_rater = 2, mean_time=0, sigma_time=1){
  mu_subject = rnorm(I,mean_sub,sigma_sub) #subject effect on probit scale
  mu_rater = rnorm(J,mean_rater,sigma_rater) #rater effect on probit scale
  mu_time = rnorm(K,mean_time,sigma_time) #time effect on probit scale
  N=nrow(X)
  mu_overall = rep(NA,N)
  prob = rep(NA,N)
  y = rep(NA,N)
  
  for(n in 1:N){
    myX = as.matrix(X[n,])
    mu_overall[n] = myX%*% beta + mu_subject[myX[,"subject"]] + mu_rater[myX[,"rater"]] + mu_time[myX[,"time"]] 
    prob[n]  =  pnorm(mu_overall[n])
    y[n] = rbinom(1,1,prob[n])
  } 
  df = cbind(y,X[,-1])
  return(list(df=df, prob=prob))
}
f.sim.partial.gen = function(beta=c(rep(0.1,3),rep(-0.1,3)), X = matrix(rep(1,6),ncol=6), I=32, J=3, K=2, mean_sub=0, 
                             sigma_sub=2, mean_rater = 0, sigma_rater = 2, rho_rater=0.3){
  Rho_rater <- diag(1-rho_rater,K,K)+matrix(rho_rater,K,K)
  Sigma_rater <- diag(sigma_rater,K) %*% Rho_rater %*% diag(sigma_rater,K)
  
  mu_subject = rnorm(I,mean_sub,sigma_sub) #subject effect on probit scale
  mu_rater = mvrnorm(J,rep(mean_rater,K),Sigma_rater)  #rater effect on probit scale
  
  N=nrow(X)
  mu_overall = rep(NA,N)
  prob = rep(NA,N)
  y = rep(NA,N)
  
  for(n in 1:N){
    myX = as.matrix(X[n,])
    mu_overall[n] = myX%*% beta + mu_subject[myX[,"subject"]] + mu_rater[myX[,"rater"],myX[,"time"]] 
    prob[n]  =  pnorm(mu_overall[n])
    y[n] = rbinom(1,1,prob[n])
  } 
  df = cbind(y,X[,-1])
  return(list(df=df, prob=prob))
}
f.sim.full.gen = function(beta=rep(1,6), X = matrix(rep(1,6),ncol=6), I=20, J=5, K=2, mean_sub=0, sigma_sub=2, mean_rater = 0,
                          sigma_rater = 2, rho_rater=0.5, mean_time=0, sigma_time = 2,
                          rho_time=0.9){
  
  Rho_rater <- diag(1-rho_rater,J,J)+matrix(rho_rater,J,J)
  Sigma_rater <- diag(sigma_rater,J) %*% Rho_rater %*% diag(sigma_rater,J)
  
  Rho_time <- diag(1-rho_time,K,K)+matrix(rho_time,K,K)
  Sigma_time <- diag(sigma_time,K) %*% Rho_time %*% diag(sigma_time,K)
  
  mu_subject = rnorm(I,mean_sub,sigma_sub) #subject effect on probit scale
  mu_rater = mvrnorm(I,rep(mean_rater,J),Sigma_rater)  #rater effect on probit scale
  
  mu_time = array(dim=c(I,J,K))
  for(i in 1:I){
    mu_time[i,,]=mvrnorm(J,rep(mean_time,K),Sigma_time)}
  
  N=nrow(X)
  mu_overall = rep(NA,N)
  prob = rep(NA,N)
  y = rep(NA,N)
  
  for(n in 1:N){
    myX = as.matrix(X[n,])
    mu_overall[n] = myX%*% beta + mu_subject[myX[,"subject"]] + mu_rater[myX[,"rater"],myX[,"time"]] +
      mu_time[myX[,"subject"],myX[,"rater"],myX[,"time"]]
    prob[n]  =  pnorm(mu_overall[n])
    y[n] = rbinom(1,1,prob[n])
  } 
  df = cbind(y,X[,-1])
  return(list(df=df, prob=prob))
}


kappa_sim <- function(df, inter = TRUE, N=1000){
  I=ns["subject"]
  J=ns["rater"]
  K=ns["time"]
  po.ind = rep(NA,N)
  pc.ind = rep(NA,N)
  if(inter){
    for(n in 1:N){
      i=sample(1:I,1)
      j=sample(1:J,2,replace = FALSE)
      k=sample(1:K,1)
      po.ind[n]=mean(df[df$subject==i & df$rater==j[1] & df$time==k,"y"]==
                       df[df$subject==i & df$rater==j[2] & df$time==k,"y"])
      
      i=sample(1:I,2)
      j=sample(1:J,2,replace = FALSE)
      k=sample(1:K,1)
      pc.ind[n]=mean(df[df$subject==i[1] & df$rater==j[1] & df$time==k,"y"]==
                       df[df$subject==i[2] & df$rater==j[2] & df$time==k,"y"])
    }
    
    if(!inter){
      for(n in 1:N){
        i=sample(1:I,1)
        k=sample(1:K,2,replace = FALSE)
        j=sample(1:J,1)
        po.ind[n]=mean(df[df$subject==i & df$rater==j & df$time==k[1],"y"]==
                         df[df$subject==i & df$rater==j & df$time==k[2],"y"])
        
        i=sample(1:I,2)
        k=sample(1:K,2,replace = FALSE)
        j=sample(1:K,1)
        pc.ind[n]=mean(df[df$subject==i[1] & df$rater==j & df$time==k[1],"y"]==
                         df[df$subject==i[2] & df$rater==j & df$time==k[2],"y"])
      }
      
      po = mean(po.ind)
      pc = mean(pc.ind)
      kappa = (po-pc)/(1-pc)
    }
    return(kappa)
  }
} #fixed effects par
model_BFN <- function(df,cov_T_str = "unstructured",
                      cov_R_str = "unstructured", addFixedEff = TRUE,
                      fixed_eff_intercept = TRUE, beta_a = 5, 
                      beta_b = 5, gamma_a = 3, gamma_b = 1.5, 
                      beta_mean=0, beta_sigma=1/0.3, betak_mean=0,
                      betak_sigma=1/0.3, rho_T_eta=1,rho_R_eta=1,
                      niters = 1000, nwarmup = 200, nchains = 1){
  load("compiled_full_new-ab.Rdata")
  
  I=length(unique(df$subject))
  J=length(unique(df$rater))
  K=length(unique(df$time))
  
  fixedeff = df%>%dplyr::select(-c(y,subject,rater,time))
  n_fixedeff = sapply(fixedeff,function(x){length(unique(x))})
  if(fixed_eff_intercept){
    X = cbind(rep(1,nrow(df)),df) %>%
      dplyr::rename("intercept"="rep(1, nrow(df))") %>%
      dplyr::select(-y)
  }
  if(!fixed_eff_intercept){
    X = df %>% dplyr::select(-y)
  }
  
  cov_T_str = ifelse(cov_T_str=="unstructured",0,ifelse(cov_T_str=="separate",1,2)) #0 for unstructured each term. 
  #1 for separate Corr matrices, 2 for common
  cov_R_str = ifelse(cov_R_str=="unstructured",0,ifelse(cov_R_str=="separate",1,2))
  
  
  N=dim(df)[1]
  p=ncol(X)
  nlevels = get_ns(as.data.frame(X)) %>% as.integer()
  beta_a = beta_a
  beta_b = beta_b
  gamma_a = gamma_a
  gamma_b = gamma_b
  beta_mean = rep(beta_mean,p)
  beta_sigma_vec = ifelse(p>1,list(diag(rep(beta_sigma,p))),beta_sigma)  #use this  when p=1
  beta_sigma_vec2 = matrix(unlist(beta_sigma_vec),nrow=p,ncol=p)
  if(p==1) {
    X=as.vector(X)
    beta_sigma = as.numeric(beta_sigma)
  }
  df_list = append(as.list(df), list(X=X,N=N,I=I,J=J,K=K,p=p,nlevels=nlevels,cov_T_str=cov_T_str,cov_R_str=cov_R_str,
                                     beta_a=beta_a,beta_b=beta_b,gamma_a=gamma_a,gamma_b=gamma_b,
                                     beta_mean=beta_mean,beta_sigma=beta_sigma_vec2,betak_mean=betak_mean,
                                     betak_sigma=betak_sigma,rho_T_eta = rho_T_eta,rho_R_eta=rho_R_eta))
  myinit = function(){list(beta_parm=rep(0,p),tau_S=0.5, tau_R=0.5, tau_T=0.5,rho_R=0.8,rho_T=rep(0.8,J))}
  nactual = niters-nwarmup
  
  m1 <- sampling(compiled,data=df_list,init=myinit,warmup = nwarmup,iter = niters,chains=nchains)
  return (m1)
}

measures_BFN <- function(df, m1, inter = TRUE, subset = FALSE, KappaMethod = "Exact",
                         MargCorr = TRUE, cov_T_str = "unstructured",
                         cov_R_str = "unstructured", Nsims=1000){
  extracted.full <- rstan::extract(m1)
  yhat.full =  extracted.full$y_hat  #(nwarmup-niters) x N
  mu_hat.full = extracted.full$mu_hat
  ns = get_ns(df)
  I = ns["subject"] %>%as.integer()
  J = ns["rater"] %>%as.integer()
  K = ns["time"] %>%as.integer()
  
  niter = m1@sim$iter
  nwarmup = m1@sim$warmup
  nactual = niter-nwarmup
  
  Kappa.est = rep(NA,nactual)
  if(cov_R_str == "unstructured") {
    rhoR.est = array(dim=c(nactual,J^2))
    CorrR.est = array(dim=c(nactual,J^2))}
  
  if(cov_R_str != "unstructured") rhoR.est = array(dim=nactual)
  if(cov_T_str == "unstructured") rhoT.est = array(dim=c(nactual,J*K*2))
  if(cov_T_str != "unstructured") rhoT.est = array(dim=nactual)
  
  
  rho_R_parm = ifelse(cov_R_str == "unstructured",
                      "Rho_R_unstr","Rho_R")
  rho_T_parm = ifelse(cov_T_str == "unstructured",
                      "Rho_T_unstr","Rho_T")
  for(i in 1:nactual){
    if(i%%10==0) print(paste0(round(i/nactual*100),"%"))
    if(!subset) posterior.data.full = data.frame(y_old = df$y,
                                                 y_new = yhat.full[i,],
                                                 mu_new = mu_hat.full[i,],
                                                 subject = df$subject,
                                                 rater = df$rater,
                                                 time = df$time,
                                                 foot = df$foot,
                                                 location = df$location)
    
    if(subset) posterior.data.full = data.frame(y_old = df$y,
                                                y_new = yhat.full[i,],
                                                mu_new = mu_hat.full[i,],
                                                subject = df$subject,
                                                rater = df$rater,
                                                time = df$time)
    
    ns = get_ns(df)
    
    df.inter.full = posterior.data.full %>% dplyr::select(-y_old,-mu_new) %>% 
      dplyr::rename(y=y_new) %>% df_to_ratings(subset=TRUE,inter=TRUE)
    df.intra.full = posterior.data.full %>% dplyr::select(-y_old,-mu_new) %>% 
      dplyr::rename(y=y_new) %>% df_to_ratings(subset=TRUE,inter=FALSE)
    
    if(inter){
      if(KappaMethod=="Fleiss") Kappa.est[i] = kappa_Fleiss_irr(df.inter.full)$Kappa
      if(KappaMethod=="Exact") Kappa.est[i] = kappa_Conger_irr(df.inter.full)$Kappa
      if(KappaMethod=="Light") Kappa.est[i] = kappa_Light_Desc(df.inter.full)$Kappa
      if(KappaMethod=="HB") Kappa.est[i] = kappa_HB_new(df.inter.full)$Kappa
      
    }
    if(!inter){
      if(KappaMethod=="Fleiss") Kappa.est[i] = kappa_Fleiss_irr(df.intra.full)$Kappa
      if(KappaMethod=="Exact") Kappa.est[i] = kappa_Conger_irr(df.intra.full)$Kappa
      if(KappaMethod=="Light") Kappa.est[i] = kappa_Light_Desc(df.intra.full)$Kappa
      if(KappaMethod=="HB") Kappa.est[i] = kappa_HB_new(df.intra.full)$Kappa
    }
  } #Kappa calc
  
  
  return(list(Kappa = mean(Kappa.est),
              rhoR = get_posterior_mean(m1,rho_R_parm),
              rhoT = get_posterior_mean(m1,rho_T_parm),
              CorrR = get_posterior_mean(m1,"Corr_R"),
              CorrT = get_posterior_mean(m1,"Corr_T")))
}



model_BIN <- function(df,addFixedEff = TRUE,
                      fixed_eff_intercept = TRUE, beta_a = 5, 
                      beta_b = 5, gamma_a = 3, gamma_b = 1.5, 
                      beta_mean=0, beta_sigma=1/0.3,
                      niters = 1000, nwarmup = 200, nchains = 1){
  load("compiled_ind_new-ab.Rdata")
  
  I=length(unique(df$subject))
  J=length(unique(df$rater))
  K=length(unique(df$time))
  
  # X = matrix(rep(1,nrow(df)))
  
  fixedeff = df%>%dplyr::select(-c(y,subject,rater,time))
  n_fixedeff = sapply(fixedeff,function(x){length(unique(x))})
  
  if(fixed_eff_intercept){
    X = cbind(rep(1,nrow(df)),df) %>%
      dplyr::rename("intercept"="rep(1, nrow(df))") %>%
      dplyr::select(-y)
  }
  if(!fixed_eff_intercept){
    X = df %>% dplyr::select(-y)
  }
  
  
  
  
  N=nrow(df)
  p=ncol(X)
  nlevels = get_ns(as.data.frame(X)) %>% as.integer()
  beta_a = beta_a
  beta_b = beta_b
  gamma_a = gamma_a
  gamma_b = gamma_b
  beta_mean_vec = rep(beta_mean,p)
  beta_sigma_vec = ifelse(p>1,list(diag(rep(beta_sigma,p))),beta_sigma)  #use this  when p=1
  beta_sigma_vec2 = matrix(unlist(beta_sigma_vec),nrow=p,ncol=p)
  if(p==1) {
    X=as.vector(X)
    beta_sigma = as.numeric(beta_sigma)
  }
  
  nactual = niters-nwarmup
  df_list = append(as.list(df), list(X=X,N=N,I=I,J=J,K=K,p=p,nlevels=nlevels,
                                     beta_a=beta_a,beta_b=beta_b,gamma_a=gamma_a,gamma_b=gamma_b,
                                     beta_mean=beta_mean_vec,beta_sigma=beta_sigma_vec2))
  # df_list = append(as.list(df), list(X=X,N=N,I=I,J=J,K=K,p=p,nlevels=nlevels,
  #                                    beta_a=beta_a,beta_b=beta_b,gamma_a=gamma_a,gamma_b=gamma_b,
  #                                    beta_mean=beta_mean, beta_mean_vec = beta_mean_vec,
  #                                    beta_sigma=beta_sigma, beta_sigma_vec2 = beta_sigma_vec2))
  myinit = function(){list(beta_parm=rep(0,p),tau_S=0.5, tau_R=0.5, tau_T=0.5)}
  
  m1 <- sampling(compiled,data=df_list,init=myinit,warmup = nwarmup,iter = niters,chains=nchains)
  return (m1)
}


measures_BIN <- function(df, m1, inter = TRUE, subset = FALSE, KappaMethod = "Exact",
                         MargCorr = TRUE){
  extracted.full <- rstan::extract(m1)
  yhat.full =  extracted.full$y_hat  #(nwarmup-niters) x N
  mu_hat.full = extracted.full$mu_hat
  ns = get_ns(df)
  I = ns["subject"] %>%as.integer()
  J = ns["rater"] %>%as.integer()
  K = ns["time"] %>%as.integer()
  
  niter = m1@sim$iter
  nwarmup = m1@sim$warmup
  nactual = niter-nwarmup
  
  Kappa.est = rep(NA,nactual)
  
  
  for(i in 1:nactual){
    if(i%%10==0) print(paste0(round(i/nactual*100),"%"))
    if(!subset) posterior.data.ind = data.frame(y_old = df$y,
                                                y_new = yhat.full[i,],
                                                mu_new = mu_hat.full[i,],
                                                subject = df$subject,
                                                rater = df$rater,
                                                time = df$time,
                                                foot = df$foot,
                                                location = df$location)
    
    if(subset) posterior.data.ind = data.frame(y_old = df$y,
                                               y_new = yhat.full[i,],
                                               mu_new = mu_hat.full[i,],
                                               subject = df$subject,
                                               rater = df$rater,
                                               time = df$time)
    
    
    df.inter.ind = posterior.data.ind %>% dplyr::select(-y_old,-mu_new) %>% 
      dplyr::rename(y=y_new) %>% df_to_ratings(subset=TRUE,inter=TRUE)
    df.intra.ind = posterior.data.ind %>% dplyr::select(-y_old,-mu_new) %>% 
      dplyr::rename(y=y_new) %>% df_to_ratings(subset=TRUE,inter=FALSE)
    
    if(inter){
      if(KappaMethod=="Fleiss") Kappa.est[i] = kappa_Fleiss_irr(df.inter.ind)$Kappa
      if(KappaMethod=="Exact") Kappa.est[i] = kappa_Conger_irr(df.inter.ind)$Kappa
      if(KappaMethod=="Light") Kappa.est[i] = kappa_Light_Desc(df.inter.ind)$Kappa
      if(KappaMethod=="HB") Kappa.est[i] = kappa_HB_new(df.inter.ind)$Kappa
      
    }
    if(!inter){
      if(KappaMethod=="Fleiss") Kappa.est[i] = kappa_Fleiss_irr(df.intra.ind)$Kappa
      if(KappaMethod=="Exact") Kappa.est[i] = kappa_Conger_irr(df.intra.ind)$Kappa
      if(KappaMethod=="Light") Kappa.est[i] = kappa_Light_Desc(df.intra.ind)$Kappa
      if(KappaMethod=="HB") Kappa.est[i] = kappa_HB_new(df.intra.ind)$Kappa
    }
  }
  
  return(list(Kappa = mean(Kappa.est),
              CorrR = get_posterior_mean(m1,"Corr_R"),
              CorrT = get_posterior_mean(m1,"Corr_T")))
}


model_BPN <- function(df,cov_T_str = "unstructured",addFixedEff = TRUE,
                      fixed_eff_intercept = TRUE, beta_a = 5, 
                      beta_b = 5, gamma_a = 3, gamma_b = 1.5, 
                      beta_mean=0, beta_sigma=1/0.3, betak_mean=0,
                      betak_sigma=1/0.3, rho_T_eta=1,rho_R_eta=1,
                      niters = 1000, nwarmup = 200, nchains = 1){
  
  load("compiled_partial_new-ab.Rdata")
  
  I=length(unique(df$subject))
  J=length(unique(df$rater))
  K=length(unique(df$time))
  
  fixedeff = df%>%dplyr::select(-c(y,subject,rater,time))
  n_fixedeff = sapply(fixedeff,function(x){length(unique(x))})
  
  if(fixed_eff_intercept){
    X = cbind(rep(1,nrow(df)),df) %>%
      dplyr::rename("intercept"="rep(1, nrow(df))") %>%
      dplyr::select(-y)
  }
  if(!fixed_eff_intercept){
    X = df %>% dplyr::select(-y)
  }
  
  cov_T_str = ifelse(cov_T_str=="unstructured",0,ifelse(cov_T_str=="separate",1,2)) #0 for unstructured each term. 
  
  
  N=nrow(df)
  p=ncol(X)
  nlevels = get_ns(as.data.frame(X)) %>% as.integer()
  beta_a = beta_a
  beta_b = beta_b
  gamma_a = gamma_a
  gamma_b = gamma_b
  beta_mean_vec = rep(beta_mean,p)
  beta_sigma_vec = ifelse(p>1,list(diag(rep(beta_sigma,p))),beta_sigma)  #use this  when p=1
  beta_sigma_vec2 = matrix(unlist(beta_sigma_vec),nrow=p,ncol=p)
  if(p==1) {
    X=as.vector(X)
    beta_sigma = as.numeric(beta_sigma)
  }
  Xk = as.vector(model.matrix(~time,as.data.frame(sapply(df, factor),stringsAsFactors=TRUE))[,2])
  nactual = niters-nwarmup
  df_list = append(as.list(df), list(X=X,N=N,I=I,J=J,K=K,p=p,nlevels=nlevels,cov_T_str=cov_T_str,
                                     beta_a=beta_a,beta_b=beta_b,gamma_a=gamma_a,gamma_b=gamma_b,
                                     beta_mean=beta_mean_vec,beta_sigma=beta_sigma_vec2,betak_mean=betak_mean,
                                     betak_sigma=betak_sigma,Xk=Xk,rho_T_eta = rho_T_eta))
  myinit = function(){list(beta_parm=rep(0,p),tau_S=0.5, tau_R=0.5, tau_T=0.5, rho_R=0.8,rho_T=rep(0.8,J))}
  
  m1 <- sampling(compiled,data=df_list,init=myinit,warmup = nwarmup,iter = niters,chains=nchains)
  return (m1)
}


measures_BPN <- function(df, m1, inter = TRUE, subset = FALSE, KappaMethod = "Exact",
                         MargCorr = TRUE){
  extracted.full <- rstan::extract(m1)
  yhat.full =  extracted.full$y_hat  #(nwarmup-niters) x N
  mu_hat.full = extracted.full$mu_hat
  ns = get_ns(df)
  I = ns["subject"] %>%as.integer()
  J = ns["rater"] %>%as.integer()
  K = ns["time"] %>%as.integer()
  
  niter = m1@sim$iter
  nwarmup = m1@sim$warmup
  nactual = niter-nwarmup
  
  Kappa.est = rep(NA,nactual)
  
  
  for(i in 1:nactual){
    if(i%%10==0) print(paste0(round(i/nactual*100),"%"))
    if(!subset) posterior.data.ind = data.frame(y_old = df$y,
                                                y_new = yhat.full[i,],
                                                mu_new = mu_hat.full[i,],
                                                subject = df$subject,
                                                rater = df$rater,
                                                time = df$time,
                                                foot = df$foot,
                                                location = df$location)
    
    if(subset) posterior.data.ind = data.frame(y_old = df$y,
                                               y_new = yhat.full[i,],
                                               mu_new = mu_hat.full[i,],
                                               subject = df$subject,
                                               rater = df$rater,
                                               time = df$time)
    
    
    df.inter.ind = posterior.data.ind %>% dplyr::select(-y_old,-mu_new) %>% 
      dplyr::rename(y=y_new) %>% df_to_ratings(subset=TRUE,inter=TRUE)
    df.intra.ind = posterior.data.ind %>% dplyr::select(-y_old,-mu_new) %>% 
      dplyr::rename(y=y_new) %>% df_to_ratings(subset=TRUE,inter=FALSE)
    
    if(inter){
      if(KappaMethod=="Fleiss") Kappa.est[i] = kappa_Fleiss_irr(df.inter.ind)$Kappa
      if(KappaMethod=="Exact") Kappa.est[i] = kappa_Conger_irr(df.inter.ind)$Kappa
      if(KappaMethod=="Light") Kappa.est[i] = kappa_Light_Desc(df.inter.ind)$Kappa
      if(KappaMethod=="HB") Kappa.est[i] = kappa_HB_new(df.inter.ind)$Kappa
      
    }
    if(!inter){
      if(KappaMethod=="Fleiss") Kappa.est[i] = kappa_Fleiss_irr(df.intra.ind)$Kappa
      if(KappaMethod=="Exact") Kappa.est[i] = kappa_Conger_irr(df.intra.ind)$Kappa
      if(KappaMethod=="Light") Kappa.est[i] = kappa_Light_Desc(df.intra.ind)$Kappa
      if(KappaMethod=="HB") Kappa.est[i] = kappa_HB_new(df.intra.ind)$Kappa
    }
  }
  
  return(list(Kappa = mean(Kappa.est),
              CorrT = get_posterior_mean(m1,"Corr_T")))
}

get_posteriors_general <- function(dataset = "drone", model = m1, modeltype = "IN"){
  mysum = summary(model)$summary
  myextracted = as.data.frame(extract(model))
  
  if(dataset == "drone"){
    betas_mysum = c("beta_parm[1]","beta_parm[2]","beta_parm[3]",
                    "beta_parm[4]","beta_parm[5]","beta_parm[6]")
    betas_myext = c("beta_parm.1","beta_parm.2","beta_parm.3",
                    "beta_parm.4","beta_parm.5","beta_parm.6")
  }
  
  
  if(dataset == "radio"){
    betas_mysum = c("beta_parm[1]","beta_parm[2]","beta_parm[3]",
                    "beta_parm[4]")
    betas_myext = c("beta_parm.1","beta_parm.2","beta_parm.3",
                    "beta_parm.4")
  }
  
  
  if(modeltype=="IN") {
    myparms = mysum[c(betas_mysum,"sigma_S", "sigma_R","sigma_T","Corr_R", "Corr_T" ),"mean"]
    mysamples = myextracted[,c(betas_myext,"sigma_S", "sigma_R","sigma_T","Corr_R", "Corr_T" )]
  }
  if(modeltype=="PN") {
    myparms = mysum[c(betas_mysum,"beta_k","rho_T[1]","rho_T[2]","rho_T[3]", "rho_T_common",
                      "sigma_S","sigma_T","Corr_T[1]","Corr_T[2]","Sigma_T[1]",
                      "Sigma_T[2]", "Corr_T[3]"),"mean"]
    mysamples = myextracted[,c(betas_myext,"beta_k","rho_T.1","rho_T.2","rho_T.3", "rho_T_common",
                               "sigma_S","sigma_T","Corr_T.1","Corr_T.2","Sigma_T.1",
                               "Sigma_T.2", "Corr_T.3")]
  }
  
  if(modeltype=="FN"){
    myparms = mysum[c(betas_mysum,"rho_R","rho_T[1]","rho_T[2]","rho_T[3]", "rho_T_common",
                      "sigma_S","sigma_R","sigma_T","Corr_R","Corr_T[1]",
                      "Corr_T[2]", "Corr_T[3]"),"mean"]
    mysamples = myextracted[,c(betas_myext,"rho_R","rho_T.1","rho_T.2","rho_T.3", "rho_T_common",
                               "sigma_S","sigma_R","sigma_T","Corr_R","Corr_T.1",
                               "Corr_T.2", "Corr_T.3")]
  }
  
  return(list(myparms = myparms, mysamples = mysamples))
}

get_Marg_Corr = function(model = "IN", sigma_S, sigma_R, sigma_T, rho_R, rho_T){
  if(model=="IN"){
    Corr_R = (sigma_S^2+sigma_T^2)/(sigma_S^2+sigma_R^2+sigma_T^2)
    Corr_T = (sigma_S^2+sigma_R^2)/(sigma_S^2+sigma_R^2+sigma_T^2)
  }
  if(model=="FN"){
    Corr_R = (sigma_S^2+rho_R*sigma_R^2)/(sigma_S^2+sigma_R^2+sigma_T^2)
    Corr_T = (sigma_S^2+rho_T*sigma_T^2)/(sigma_S^2+sigma_R^2+sigma_T^2)
  }
  if(model=="PN"){
    Corr_R = (sigma_S^2)/(sigma_S^2+sigma_T^2)
    Corr_T = (sigma_S^2+rho_T*sigma_T^2)/(sigma_S^2+sigma_T^2)
  }
}