#' @title Accommodate partial ordering in optimization
#' @description Accommodate partial ordering and the output will be used in the order constrained optimization program
#' @param order_list A list of strata order allowing for partial ordering grouped in a vector within a list
#' Example 1: \code{list(1,2,3)} means 3 stratum are in full ordering 1<2<3
#' Example 2: \code{list(1,c(2,3))} means 3 stratum are in partial ordering 1<2, 1<3 and there is no ordering between stratum 2 and 3
#' @return A matrix with columns stratum number and rows the stratum-wise comparison with 1 the larger one and -1 the smaller one.
#' @examples
#' library(constrselect)
#' partial_order(list(1,c(2,3),4,5))
#' @rdname partial_order
#' @export
partial_order <- function(order_list){
  ## partial ordering according to order list
  if(class(order_list)!="list") stop("The input or order_list must be a list")

  order_output <- NULL
  for (i in 1:(length(order_list)-1)){
    raw_output <- rep(0,length(unlist(order_list)))
    raw_output[order_list[[i]]]<--1
    raw_output[order_list[[i+1]]]<-1

    order_output <- rbind(order_output,raw_output)
  }
  order_matrix <- NULL
  for (i in 1:nrow(order_output)){
    if (sum(order_output[i,]==1)>1 | sum(order_output[i,]==-1)>1){
      raw_output <- rep(0,length(unlist(order_list)))

      id_neg <- which(order_output[i,]==-1)
      id_pos <- which(order_output[i,]==1)
      for (j in id_neg){
        raw_output[j] <- -1
        for (k in id_pos){
          raw_output[k] <- 1
          order_matrix <- rbind(order_matrix,raw_output)
          raw_output[k] <- 0
        }
        raw_output[j] <- 0

      }

    } else {
      order_matrix <- rbind(order_matrix,order_output[i,])
    }
  }
  rownames(order_matrix) <- NULL
  return(order_matrix)
}


#' @title Order constrained binary outcome estimation
#' @description For each stratum, given the order list, estimate the
#' order constrained binary outcome estimator using quadratic programming optimization
#' @param r A vector of responses for all strata
#' @param n A vector of total sample sizes for all strata
#' @param order_list A list of strata order allowing for partial ordering grouped in a vector within a list
#' @return A vector containing the solution of the quadratic programming problem.
#' @details DETAILS
#' @examples
#' order_constrain(r=c(5,10),n=c(20,25),order_list=list(1,2))
#' order_constrain(r=c(5,10,13),n=c(20,25,25),order_list=list(1,c(2,3))) ## the strata ordering is 1<2, 1<3 and no ordering between stratum 2 and 3
#' @seealso
#'  \code{\link[quadprog]{solve.QP}}
#' @rdname order_constrain
#' @export
#' @importFrom quadprog solve.QP
order_constrain<-function(r,n,order_list){
  if(length(r)!=length(unlist(order_list))) stop("The number of response and order list must be the same")
  if(length(n)!=length(unlist(order_list))) stop("The number of sample size for strata and order list must be the same")

  p_bar<-r/n
  S<-length(p_bar)
  Dmat<-matrix(0,S,S)
  diag(Dmat)<-n
  A_upper<-diag(S)
  ## partial ordering according to order list
  A_lower <- partial_order(order_list)
  rownames(A_lower) <- NULL
  Amat<-rbind(A_upper,A_lower)
  Amat<-t(Amat)
  b_vec<-c(numeric(S-1+length(p_bar)))

  quadprog::solve.QP(Dmat,dvec=r,Amat,b_vec)$solution}



#' @title Binary outcome estimation for two treatments with two strata.
#' @description Estimation of P_corr and P_amb for selection design between
#' two treatments with two strata using exact binomial method.
#' @param n Total sample size for all strata.
#' @param p_inf A vector of true probability for inferior treatment stratum arm.
#' @param D A vector of two treatment arms differences for each stratum, Default: c(0.15, 0.15)
#' @param d A vector of ambiguous region for each stratum, Default: c(0.05, 0.05)
#' @param prop.strat The sample size proportion for the inferior strata, Default: 0.3
#' @param study Could be either "Constrained" or "Origin" for the two type of study design, Default: NULL
#' @return Return a vector of pcorr and pamb representing $P_{corr}$ and $P_{amb}$
#' @details DETAILS
#' @examples
#' library(constrselect)
#' pickwin_bin_exact(n = 50, p_inf = c(0.25,0.3), D=c(0.15,0.15),d=c(0.05,0.05),prop.strat=0.4,study="Constrained")
#' @rdname pickwin_bin_exact
#' @export
pickwin_bin_exact<- function(n, p_inf,
                             D=c(0.15,0.15),d=c(0.05,0.05),
                             prop.strat=0.3,study=NULL)
{
  n1<- ceiling(prop.strat*n)
  n2<- n-n1
  Nk <- c(n1,n2)
  pp_wrong <- 0
  pp_corr <- 0
  for (i1 in 0:n1){
    for (i2 in 0:n2){
      Rk1 <- c(i1,i2)
      if (study == "Constrained"){
        p_1 <- order_constrain(Rk1,Nk,order_list=list(1,2))

      } else if (study =="Origin"){
        p_1 <- Rk1/Nk

      } else{
        stop("This is an error message")
      }

      for (j1 in 0:n1) {
        for (j2 in 0:n2){
          Rk2 <- c(j1,j2)
          if (study == "Constrained"){
            p_2 <- order_constrain(Rk2,Nk,order_list=list(1,2))

          } else if (study =="Origin"){
            p_2 <- Rk2/Nk

          } else{
            stop("This is an error message")
          }

          ## Wrong cases
          if (all(p_2+d<p_1)){

            pp_wrong <- pp_wrong+prod(dbinom(c(i1,i2), size=c(n1,n2), prob=p_inf))*
              prod(dbinom(c(j1,j2), size=c(n1,n2), prob=p_inf+D))
          }

          ## Correct cases
          if (all(p_2>p_1+d)) {
            pp_corr <- pp_corr+prod(dbinom(c(i1,i2), size=c(n1,n2), prob=p_inf))*
              prod(dbinom(c(j1,j2), size=c(n1,n2), prob=p_inf+D))
          }

        }

      }

    }

  }
  p_amb = 1-pp_wrong-pp_corr
  output <- c(pcorr=pp_corr, pamb=p_amb)

  return(output)
}



#' @title Binary outcome estimation for two treatments with multiple strata.
#' @description Estimation of P_corr and P_amb for selection design between
#' two treatments with multiple strata using simulation method.
#' @param n Total sample size for each treatment arm.
#' @param p_inf A vector of response probabilities for the inferior treatment arm for each stratum.
#' @param D A vector of two treatment arms differences for each stratum, Default: c(0.15, 0.15, 0.15)
#' @param d A vector of ambiguous region for each stratum, Default: c(0.05, 0.05, 0.05)
#' @param prop.strat The sample size proportion for each stratum, Default: c(0.2, 0.3, 0.5)
#' @param study Could be either "Constrained" or "Origin" for the two type of study design, Default: 'Constrained'
#' @param S Number of simulation
#' @param cluster Number of parallel running CPU cores, Default: 6
#' @param order_list A list of strata order allowing for partial ordering grouped in a vector within a list.
#' @param with_seed Random seed for simulation, Default: NULL
#' @return A data frame of binary probability, correct and error decision for each simulated scenario
#' @details DETAILS
#' @examples
#' library(constrselect)
#' out <- pickwin_bin_multiple(n = 50, p_inf = c(0.25,0.28,0.28), D=c(0.15,0.15,0.15),d=c(0.05,0.05,0.05),
#' prop.strat=c(0.3,0.3,0.4),study="Constrained",S = 1000,cluster=2,order_list=list(1,c(2,3)))
#' @rdname pickwin_bin_multiple
#' @export
#' @import doParallel
#' @import doRNG
#' @import foreach
pickwin_bin_multiple <- function(n, p_inf,
                                 D=c(0.15,0.15,0.15), d = c(0.05,0.05,0.05),
                                 prop.strat=c(0.2,0.3,0.5),study="Constrained",
                                 S,cluster = 6,order_list,with_seed=NULL) {
  Nk<- round(prop.strat*n)
  Nk[length(p_inf)] <- n - sum(Nk[1:(length(p_inf)-1)])

  cl <- makeCluster(cluster)
  registerDoParallel(cl)

  if(!is.null(with_seed)){set.seed(with_seed, kind = "L'Ecuyer-CMRG")}

  bin_estimator <- foreach(1:S, .combine = rbind,.packages = c("quadprog"),
                           .export = c("order_constrain","partial_order")) %dorng% {
                             corr <- 0
                             wrong <- 0
                             Rk1 <- mapply(function(x,y) rbinom(1,size=x,y),Nk,p_inf)
                             Rk2 <- mapply(function(x,y) rbinom(1,size=x,y),Nk,p_inf+D)


                             if (study == "Constrained"){
                               p_1 <- order_constrain(Rk1,Nk,order_list)
                               p_2 <- order_constrain(Rk2,Nk,order_list)

                             } else if (study =="Origin"){
                               p_1 <- Rk1/Nk
                               p_2 <- Rk2/Nk

                             } else{
                               stop("This is an error message")
                             }

                             if(all(p_2>p_1+d)){
                               corr <- corr+1
                             }
                             if(all(p_2+d<p_1)){
                               wrong <- wrong+1
                             }

                             result <- c(p_1,p_2,corr,wrong)
                             data.frame(t(result))

                           }
  on.exit(stopCluster(cl))
  colnames(bin_estimator) <- c(paste0("P_A",1:length(p_inf)),
                               paste0("P_B",1:length(p_inf)),"Corr","Wrong")
  return (bin_estimator)

}


#' @title Calculate distinct events
#' @description Calculate the number of distinct events in (0,x], which will be used in the following likelihood calculation
#' @param x Time we are interested in comparing.
#' @param event_time A vector of event times.
#' @param event_ind A vector of event indicator (1, 0)
#' @return A number of distinct events in (0,x]
#' @details DETAILS
#' @examples
#' Mg(x = 6, event_time=c(0.41, 0.8, 2.2, 2.5, 3, 4.5, 6, 6.2, 6.5, 7),event_ind=c(1,1,1,1,1,1,0,0,0,0)) #6
#' @rdname Mg
#' @export
Mg <- function(x,event_time,event_ind){
  ## Number of distinct events in (0,x]
  sum(event_time[!duplicated(event_time)]<=x & event_ind[!duplicated(event_time)]==1)
}

#' @title Calculate number at risk at a specific time
#' @description Calculate number at risk at time x, which will be used in the following likelihood calculation
#' @param x Time we are interested in comparing.
#' @param event_time A vector of event times.
#' @return A number at risk at time x
#' @details DETAILS
#' @examples
#' library(constrselect)
#' Ng(x=6,event_time=c(0.41, 0.8, 2.2, 2.5, 3, 4.5, 6, 6.2, 6.5, 7)) ## 4
#' @rdname Ng
#' @export

Ng <- function(x,event_time){
  ## Number at risk at time t
  sum(event_time>=x)
}




#' @title Function used for constrained likelihood optimization
#' @description Function used for constrained likelihood optimization
#' @param x Time we are interested in comparing.
#' @param qn Parameter we are optimizing which survival probability at time x is exp(qn)
#' @param nrisk A vector of number at risk for each event time
#' @param nevent A vector of 1 representing the number of event
#' @param event_time A vector of event times.
#' @param event_ind A vector of event indicator (1, 0)
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' library(constrselect)
#' Kg(x = 6,qn = -0.1,nrisk = c(15, 14, 13, 12, 11, 10),nevent = c(1, 1, 1, 1, 1, 1),
#' event_time = c(0.47, 0.84, 2.25, 2.5, 2.8, 4.5, 6, 6.16, 6.6, 6.86, 8.17),
#' event_ind = c(1, 1, 1, 1, 1, 1, 0, 0,0,0,0)) #29.99996
#' @rdname Kg
#' @export

Kg <- function(x,qn,nrisk,nevent,event_time,event_ind){
  if(qn>0) stop("qn needs to be smaller than 0")
  mgx <- Mg(x,event_time,event_ind)

  if (mgx ==0 & abs(qn)>0.00001){

    return(-Ng(x,event_time))

  } else {
    ## k=inf if q=0
    if (abs(qn)<0.00001){
      k_hat <- Inf
      return(k_hat)

    } else if (qn == -Inf) {

      k_hat <- nevent[mgx]-nrisk[mgx]
      return(k_hat)
    } else{

      fr <- function(k){
        left <- 0
        for (j in 1:mgx){
          left <- left+log(1-nevent[j]/(nrisk[j]+k))
        }
        return(abs(left-qn))
      }
      k_hat = optimize(fr, c(-50, 50))$minimum

      k_hat <- max(k_hat,-Ng(x,event_time))
      return(k_hat)
    }

  }
}



#' @title Calculate constrained likelihood
#' @description Calculate constrained likelihood for futher get the constrained NPMLE
#' @param x Time we are interested in comparing.
#' @param qn Parameter we are optimizing which survival probability at time x is exp(qn)
#' @param nrisk A vector of number at risk for each event time
#' @param nevent A vector of 1 representing the number of event
#' @param event_time A vector of event times.
#' @param event_ind A vector of event indicator (1, 0)
#' @return Likelihood value
#' @details DETAILS
#' @examples
#' library(constrselect)
#' llkhd(x=6,qn = -0.1,nrisk = c(15, 14, 13, 12, 11, 10),nevent = c(1, 1, 1, 1, 1, 1),
#' event_time = c(0.47, 0.84, 2.25, 2.5, 2.8, 4.5, 6, 6.16, 6.6, 6.86, 8.17),
#' event_ind = c(1, 1, 1, 1, 1, 1, 0, 0,0,0,0))  ## -24.12788
#' @rdname llkhd
#' @export

llkhd <- function(x,qn,nrisk,nevent,event_time,event_ind){

  if(qn>0) stop("qn needs to be smaller than 0")
  mgx <- Mg(x,event_time,event_ind)

  if (mgx!=0){
    event_llkhd <- 0
    k_hat <- Kg(x,qn,nrisk,nevent,event_time,event_ind)

    ## if equal, there will be 0 in denominators which will lead to -inf values
    ## n-Ng=d
    if(k_hat==-Ng(x,event_time)){
      return(-1000)
    } else{
      for (i in 1:mgx){
        event_llkhd <- event_llkhd+(nrisk[i]-nevent[i])*
          log(nrisk[i]+k_hat-nevent[i])-
          nrisk[i]*log(nrisk[i]+k_hat)
      }
      return(event_llkhd)
    }} else{
      return(Ng(x,event_time)*qn)
    }

}



## simulation of survival times
#' @title Simulation of survival times
#' @description Simulate survival times based on poisson arrival and exponential survival.
#' Follow up the patients for addition month after the last patient is accrued.
#' @param nmax Number of patients need to be accrued
#' @param arrival_rate The poisson arrival rate for patients, number of patients accrued each month/year
#' @param event_rate The exponential event rate for patients
#' @param FUP Additional follow up time after the last patient is accrued
#' @return Return a dataframe with time to event/censoring and event status (1,0)
#' @details DETAILS
#' library(constrselect)
#' sim_surv(nmax=12,arrival_rate=4,event_rate=0.08,FUP=6)
#' @rdname sim_surv
#' @export

sim_surv <- function(nmax,arrival_rate,event_rate,FUP){

  wait.t = rexp(nmax,rate = arrival_rate)
  arrival.t = cumsum(wait.t)
  event.t = rexp(nmax,rate=event_rate)
  tobs = arrival.t[nmax]
  tobs = tobs + FUP

  ## tobs being the total follow up time
  n.fail = sum(arrival.t+ event.t <= tobs)
  t.event = rep(0,nmax)
  t.ind = rep(0,nmax)
  for(j in 1:length(t.event)) {
    ## count events in the observed follow up time
    t.event[j] = ifelse(arrival.t[j]+event.t[j]<=tobs,event.t[j],tobs-arrival.t[j])
    t.ind[j] = ifelse(arrival.t[j]+event.t[j]<=tobs,1,0)
  }

  return(data.frame(cbind(time = t.event,ind = t.ind)))

}


#' @title Estimate constrained NPMLE of survival probabilities in stratified selection setting
#' @description Estimate constrained NPMLE of survival probabilities in stratified selection setting using simulation.
#' For each simulation setting, get whether it is correctly or wrongly counted in.
#' @param maxn Number of patients need to be accrued for each treatment arm
#' @param prop The sample size proportion for two strata, Default: c(0.2, 0.3, 0.5)
#' @param surv_inf The exponential survival probability for patients in the inferior treatment arm
#' @param surv_sup The exponential survival probability for patients in the superior treatment arm
#' @param d A vector of ambiguous region for each stratum, Default: c(0.05, 0.05, 0.05)
#' @param arrival_rate The poisson arrival rate for patients, number of patients accrued each month/year
#' @param FUP Additional follow up time after the last patient is accrued
#' @param x Time we are interested in comparing.
#' @param S Number of simulation
#' @param study Could be either "Constrained" or "Origin" for the two type of study design, Default: 'Constrained'
#' @param cluster Number of parallel running CPU cores, Default: 6
#' @param order_list A list of strata order allowing for partial ordering grouped in a vector within a list.
#' @param with_seed Random seed for simulation, Default: NULL
#' @return A data frame of survival probability, correct and error decision for each simulated scenario
#' @details DETAILS
#' @examples
#' library(constrselect)
#' ## Use exponential survival to specify event rate for month 6.
# test <- pickwin_surv_fun(maxn=50,prop=c(0.3,0.3,0.4),surv_inf=c(0.5,0.6,0.6),
# surv_sup=c(0.7,0.8,0.8),d=c(0.05,0.05,0.05), arrival_rate=4,FUP=6,
# x=6,S=100,study = "Constrained",cluster=6,order_list=list(1,c(2,3)),with_seed = 111)
#' @seealso
#'  \code{\link[doParallel]{registerDoParallel}}
#'  \code{\link[foreach]{foreach}}
#'  \code{\link[survival]{survfit}}
#' @rdname pickwin_surv_fun
#' @export
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import doRNG
#' @import survival
#' @import dplyr
pickwin_surv_fun <- function(maxn,prop,surv_inf,
                             surv_sup,d,arrival_rate,FUP,
                             x,S,study = "Constrained",cluster,
                             order_list,with_seed=NULL) {
  #S is simulation times
  n<- round(prop*maxn)
  n[length(prop)] <- maxn - sum(n[1:(length(prop)-1)])

  cl <- makeCluster(cluster)
  doParallel::registerDoParallel(cl)

  if(!is.null(with_seed)){set.seed(with_seed, kind = "L'Ecuyer-CMRG")}

  surv_estimate <- foreach(1:S, .combine = rbind,
                           .export = c("order_constrain","partial_order",
                                       "sim_surv","llkhd","Kg","Ng","Mg"),
                           .packages = c("dplyr","survival")) %dorng% {

                             corr <- 0
                             wrong <- 0
                             trtA <- lapply(mapply(function(x,y) sim_surv(x,arrival_rate,event_rate=y,FUP),n,-log(surv_inf)/x,SIMPLIFY = FALSE),data.frame)
                             trtB <- lapply(mapply(function(x,y) sim_surv(x,arrival_rate,event_rate=y,FUP),n,-log(surv_sup)/x,SIMPLIFY = FALSE),data.frame)
                             names(trtA) <- paste("Strata",1:length(n))
                             names(trtB) <- paste("Strata",1:length(n))

                             trtA <- trtA %>% bind_rows(.id = "column_label") %>%  arrange(column_label,time) %>% rename(`group`=`column_label`)
                             trtB <- trtB %>% bind_rows(.id = "column_label")%>%  arrange(column_label,time) %>% rename(`group`=`column_label`)

                             fitA <- survival::survfit(Surv(time, ind)~group,data=trtA)
                             id_A <- trtA %>% filter(ind ==1)
                             ## Split by strata
                             id_A <- as.numeric(substr(id_A$group,7,8))

                             fitB <- survival::survfit(Surv(time, ind)~group,data=trtB)
                             id_B <- trtB %>% filter(ind ==1)
                             id_B <- as.numeric(substr(id_B$group,7,8))

                             if (study =="Constrained"){


                               ## For treatment A
                               nrisk <- split(summary(fitA, extend = TRUE)$n.risk,id_A)
                               nevent <- split(summary(fitA, extend = TRUE)$n.event,id_A)
                               id_full <- as.numeric(factor(trtA$group))

                               event_time <- split(trtA$time,id_full)
                               event_ind <- split(trtA$ind,id_full)

                               if (length(unique(id_A))==0){
                                 # stop("There is no event for all stratum.")
                                 nrisk <- as.list(n)
                                 nevent <- as.list(rep(0,length(surv_inf)))
                                 names(nrisk) <- as.character(1: length(surv_inf))
                                 names(nevent) <- as.character(1: length(surv_inf))

                               } else if (length(unique(id_A))<length(n)){
                                 missing_id <- setdiff(1:length(n),unique(id_A))
                                 for (i in 1:length(missing_id)){
                                   nrisk <- append(nrisk,list(n[missing_id[i]]),missing_id[i]-1)
                                   names(nrisk)[missing_id[i]] <- as.character(missing_id[i])

                                   nevent <- append(nevent,list(0),missing_id[i]-1)
                                   names(nevent)[missing_id[i]] <- as.character(missing_id[i])

                                 }
                               }

                               fn <- function(q){

                                 sum(mapply(function(a,b,c,d,e) -llkhd(x,a,b,c,d,e), q,nrisk,
                                            nevent,event_time,event_ind))

                               }

                               grr <- function(q){
                                 mapply(function(a,b,c,d,e) Kg(x,a,b,c,d,e), q,nrisk,
                                        nevent,event_time,event_ind)


                               }

                               A_upper<-diag(length(n))*(-1)

                               ## partial ordering according to order list
                               A_lower <- partial_order(order_list)

                               ui <- rbind(A_lower,A_upper)
                               start_surv <- sort(runif(length(n),-0.3,-0.1))

                               ## constrained survival prob
                               q_optim <- constrOptim(start_surv, fn, grad=grr,
                                                      ui = ui, ci = rep(0,nrow(ui)))$par

                               S_A <- exp(q_optim)

                               ## For treatment B
                               nrisk <- split(summary(fitB, extend = TRUE)$n.risk,id_B)
                               nevent <- split(summary(fitB, extend = TRUE)$n.event,id_B)
                               id_full <- as.numeric(factor(trtB$group))

                               event_time <- split(trtB$time,id_full)
                               event_ind <- split(trtB$ind,id_full)

                               if (length(unique(id_B))==0){
                                 # stop("There is no event for all stratum.")
                                 nrisk <- as.list(n)
                                 nevent <- as.list(rep(0,length(surv_inf)))
                                 names(nrisk) <- as.character(1: length(surv_inf))
                                 names(nevent) <- as.character(1: length(surv_inf))

                               } else if  (length(unique(id_B))<length(n)){
                                 missing_id <- setdiff(1:length(n),unique(id_B))
                                 for (i in 1:length(missing_id)){
                                   nrisk <- append(nrisk,list(n[missing_id[i]]),missing_id[i]-1)
                                   names(nrisk)[missing_id[i]] <- as.character(missing_id[i])

                                   nevent <- append(nevent,list(0),missing_id[i]-1)
                                   names(nevent)[missing_id[i]] <- as.character(missing_id[i])

                                 }
                               }

                               ## constrained survival prob
                               q_optim <- constrOptim(start_surv, fn, grad=grr,
                                                      ui = ui, ci = rep(0,nrow(ui)))$par
                               S_B <- exp(q_optim)



                             } else if (study == "Origin"){

                               nrisk <- split(summary(fitA, extend = TRUE)$n.risk,id_A)
                               nevent <- split(summary(fitA, extend = TRUE)$n.event,id_A)

                               event_time <- split(summary(fitA, extend = TRUE)$time,id_A)
                               surv_prob <- split(summary(fitA, extend = TRUE)$surv,id_A)

                               S_A <- NULL
                               for (i in 1:length(surv_inf)){
                                 if (i %in% id_A){

                                   if (max(event_time[[as.character(i)]])<x){
                                     S_A_i=min(surv_prob[[as.character(i)]])
                                   } else{
                                     S_A_i <- summary(fitA,t=x, extend = TRUE)$surv[i]
                                   }

                                   S_A <- c(S_A,S_A_i)
                                 } else{
                                   S_A <- c(S_A,1)

                                 }
                               }


                               nrisk <- split(summary(fitB, extend = TRUE)$n.risk,id_B)
                               nevent <- split(summary(fitB, extend = TRUE)$n.event,id_B)

                               event_time <- split(summary(fitB, extend = TRUE)$time,id_B)
                               surv_prob <- split(summary(fitB, extend = TRUE)$surv,id_B)

                               S_B <- NULL
                               for (i in 1:length(surv_inf)){
                                 if (i %in% id_B){

                                   if (max(event_time[[as.character(i)]])<x){
                                     S_B_i=min(surv_prob[[as.character(i)]])
                                   } else{
                                     S_B_i <- summary(fitB,t=x, extend = TRUE)$surv[i]
                                   }

                                   S_B <- c(S_B,S_B_i)
                                 } else{
                                   S_B <- c(S_B,1)

                                 }
                               }

                             } else{
                               stop("This is an error message.")
                             }

                             # if(length(S_A)!=length(surv_inf)) S_A=c(S_A,rep(1,length(surv_inf)-length(S_A)))
                             # if(length(S_B)!=length(surv_inf)) S_B=c(S_B,rep(1,length(surv_inf)-length(S_B)))

                             if(all(S_A+d<S_B)){
                               corr <- corr+1
                             }
                             if(all(S_A>S_B+d)){
                               wrong <- wrong+1
                             }

                             data.frame(t(S_A),t(S_B),corr,wrong)
                           }
  on.exit(stopCluster(cl))

  colnames(surv_estimate) <- c(paste0("S_A",1:length(surv_inf)),
                               paste0("S_B",1:length(surv_inf)),"Corr","Wrong")
  return (surv_estimate)

}




