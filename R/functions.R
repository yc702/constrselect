#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param order_list PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname partial_order
#' @export
partial_order <- function(order_list){
  ## partial ordering according to order list
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


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param r PARAM_DESCRIPTION
#' @param n PARAM_DESCRIPTION
#' @param order_list PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[quadprog]{solve.QP}}
#' @rdname order_constrain
#' @export
#' @importFrom quadprog solve.QP
order_constrain<-function(r,n,order_list){
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

  quadprog::solve.QP(Dmat,dvec=r,Amat,b_vec)}



#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param n PARAM_DESCRIPTION
#' @param p1 PARAM_DESCRIPTION
#' @param strata_diff PARAM_DESCRIPTION
#' @param D PARAM_DESCRIPTION, Default: c(0.15, 0.15)
#' @param d PARAM_DESCRIPTION, Default: c(0.05, 0.05)
#' @param prop.strat PARAM_DESCRIPTION, Default: 0.3
#' @param study PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname pickwin_bin_exact
#' @export
pickwin_bin_exact<- function(n, p1, strata_diff,
                             D=c(0.15,0.15),d=c(0.05,0.05),
                             prop.strat=0.3,study=NULL)
{
  n1<- ceiling(prop.strat*n)
  n2<- n-n1
  Nk <- c(n1,n2)
  pp_err <- 0
  pp_corr <- 0
  for (i1 in 0:n1){
    for (i2 in 0:n2){
      Rk1 <- c(i1,i2)
      if (study == "Constrained"){
        p_1 <- order_constrain(Rk1,Nk,order_list=list(1,2))$solution

      } else if (study =="Origin"){
        p_1 <- Rk1/Nk

      } else{
        stop("This is an error message")
      }

      for (j1 in 0:n1) {
        for (j2 in 0:n2){
          Rk2 <- c(j1,j2)
          if (study == "Constrained"){
            p_2 <- order_constrain(Rk2,Nk,order_list=list(1,2))$solution

          } else if (study =="Origin"){
            p_2 <- Rk2/Nk

          } else{
            stop("This is an error message")
          }

          ## Error cases
          if (p_2<p_1+d){

            pp_err <- pp_err+dbinom(i1, size=n1, prob=p1)*dbinom(i2, size=n2, prob=p1+strata_diff)*
              dbinom(j1, size=n1, prob=p1+D[1])*dbinom(j2, size=n2, prob=p1+strata_diff+D[2])
          }

          ## Correct cases
          if (p_2>p_1+d) {
            pp_corr <- pp_corr+dbinom(i1, size=n1, prob=p1)*dbinom(i2, size=n2, prob=p1+strata_diff)*
              dbinom(j1, size=n1, prob=p1+D[1])*dbinom(j2, size=n2, prob=p1+strata_diff+D[2])
          }

        }

      }

    }

  }
  p_amb = 1-pp_err-pp_corr
  return(data.frame(pcorr=pp_corr, pamb=p_amb))
}



#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param n PARAM_DESCRIPTION
#' @param pa_list PARAM_DESCRIPTION
#' @param D PARAM_DESCRIPTION, Default: c(0.15, 0.15, 0.15)
#' @param d_diff PARAM_DESCRIPTION, Default: c(0.05, 0.05, 0.05)
#' @param prop.strat PARAM_DESCRIPTION, Default: c(0.2, 0.3, 0.5)
#' @param study PARAM_DESCRIPTION, Default: 'Constrained'
#' @param S PARAM_DESCRIPTION
#' @param cluster PARAM_DESCRIPTION, Default: 6
#' @param order_list PARAM_DESCRIPTION
#' @param with_seed PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname pickwin_bin_multiple
#' @export
#' @import doParallel
#' @import foreach
pickwin_bin_multiple <- function(n, pa_list,
                                 D=c(0.15,0.15,0.15), d_diff = c(0.05,0.05,0.05),
                                 prop.strat=c(0.2,0.3,0.5),study="Constrained",
                                 S,cluster = 6,order_list,with_seed=NULL) {
  Nk<- ceiling(prop.strat*n)
  cl <- makeCluster(cluster)
  registerDoParallel(cl)

  if(!is.null(with_seed)){set.seed(with_seed, kind = "L'Ecuyer-CMRG")}

  bin_estimator <- foreach(1:S, .combine = rbind,.packages = c("quadprog"),
                           .export = c("order_constrain","partial_order")) %dopar% {
                             corr <- 0
                             err <- 0
                             Rk1 <- sapply(mapply(function(x,y) rbinom(x,size=1,y),Nk,pa_list),sum)
                             Rk2 <- sapply(mapply(function(x,y) rbinom(x,size=1,y),Nk,pa_list+D),sum)


                             if (study == "Constrained"){
                               p_1 <- order_constrain(Rk1,Nk,order_list)$solution
                               p_2 <- order_constrain(Rk2,Nk,order_list)$solution

                             } else if (study =="Origin"){
                               p_1 <- Rk1/Nk
                               p_2 <- Rk2/Nk

                             } else{
                               stop("This is an error message")
                             }

                             if(all(p_2>p_1+d_diff)){
                               corr <- corr+1
                             }
                             if(all(p_2<p_1+d_diff)){
                               err <- err+1
                             }

                             result <- c(p_1,p_2,corr,err)
                             data.frame(t(result))

                           }
  return (bin_estimator)

}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param event_time PARAM_DESCRIPTION
#' @param event_ind PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Mg
#' @export
Mg <- function(x,event_time,event_ind){
  ## Number of distinct events in (0,x]
  sum(event_time[!duplicated(event_time)]<=x & event_ind[!duplicated(event_time)]==1)
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param event_time PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Ng
#' @export

Ng <- function(x,event_time){
  ## Number at risk at time t
  sum(event_time>=x)
}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param qn PARAM_DESCRIPTION
#' @param nrisk PARAM_DESCRIPTION
#' @param nevent PARAM_DESCRIPTION
#' @param event_time PARAM_DESCRIPTION
#' @param event_ind PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname Kg
#' @export

Kg <- function(qn,nrisk,nevent,event_time,event_ind){
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
      k_hat = optimize(fr, c(-30, 30))$minimum

      k_hat <- max(k_hat,-Ng(x,event_time))
      return(k_hat)
    }

  }
}

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param qq PARAM_DESCRIPTION
#' @param nrisk PARAM_DESCRIPTION
#' @param nevent PARAM_DESCRIPTION
#' @param event_time PARAM_DESCRIPTION
#' @param event_ind PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname llkhd
#' @export

llkhd <- function(qq,nrisk,nevent,event_time,event_ind){
  mgx <- Mg(x,event_time,event_ind)

  if (mgx!=0){
    event_llkhd <- 0
    k_hat <- Kg(qq,nrisk,nevent,event_time,event_ind)

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
      return(Ng(x,event_time)*qq)
    }

}



## simulation of survival times
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param nmax PARAM_DESCRIPTION
#' @param arrival_rate PARAM_DESCRIPTION
#' @param event_rate PARAM_DESCRIPTION
#' @param FUP PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
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

  return(cbind(time = t.event,ind = t.ind))

}


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param maxn PARAM_DESCRIPTION
#' @param prop PARAM_DESCRIPTION
#' @param event_rate_A PARAM_DESCRIPTION
#' @param trt_diff PARAM_DESCRIPTION
#' @param d_diff PARAM_DESCRIPTION
#' @param arrival_rate PARAM_DESCRIPTION
#' @param FUP PARAM_DESCRIPTION
#' @param x PARAM_DESCRIPTION
#' @param S PARAM_DESCRIPTION
#' @param study PARAM_DESCRIPTION, Default: 'Constrained'
#' @param cluster PARAM_DESCRIPTION
#' @param order_list PARAM_DESCRIPTION
#' @param with_seed PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[doParallel]{registerDoParallel}}
#'  \code{\link[foreach]{foreach}}
#'  \code{\link[survival]{survfit}}
#' @rdname pickwin_surv_fun
#' @export
#' @import doParallel
#' @import foreach
#' @import survival
#' @import dplyr
pickwin_surv_fun <- function(maxn,prop,event_rate_A,
                             trt_diff,d_diff,arrival_rate,FUP,
                             x,S,study = "Constrained",cluster,
                             order_list,with_seed=NULL) {
  #S is simulation times
  n<- ceiling(prop*maxn)

  cl <- makeCluster(cluster)
  doParallel::registerDoParallel(cl)

  if(!is.null(with_seed)){set.seed(with_seed, kind = "L'Ecuyer-CMRG")}

  surv_estimate <- foreach::foreach(1:S, .combine = rbind,
                           .packages = c("dplyr","survival")) %dopar% {

                             corr <- 0
                             err <- 0
                             trtA <- lapply(mapply(function(x,y) sim_surv(x,arrival_rate,event_rate=y,FUP),n,event_rate_A),data.frame)
                             trtB <- lapply(mapply(function(x,y) sim_surv(x,arrival_rate,event_rate=y,FUP),n,event_rate_A+trt_diff),data.frame)
                             names(trtA) <- paste("Strata",1:length(n))
                             names(trtB) <- paste("Strata",1:length(n))

                             trtA <- trtA %>% bind_rows(.id = "column_label") %>%  arrange(column_label,time) %>% rename(`group`=`column_label`)
                             trtB <- trtB %>% bind_rows(.id = "column_label")%>%  arrange(column_label,time) %>% rename(`group`=`column_label`)

                             fitA <- survival::survfit(Surv(time, ind)~group,data=trtA)
                             id_A <- trtA %>% filter(ind ==1)
                             ## Split by strata
                             id_A <- as.numeric(factor(id_A$group))

                             fitB <- survival::survfit(Surv(time, ind)~group,data=trtB)
                             id_B <- trtB %>% filter(ind ==1)
                             id_B <- as.numeric(factor(id_B$group))

                             if (study =="Constrained"){
                               ## For treatment A
                               nrisk <- split(summary(fitA)$n.risk,id_A)
                               nevent <- split(summary(fitA)$n.event,id_A)
                               id_full <- as.numeric(factor(trtA$group))

                               event_time <- split(trtA$time,id_full)
                               event_ind <- split(trtA$ind,id_full)


                               fn <- function(q){

                                 sum(mapply(function(a,b,c,d,e) -llkhd(a,b,c,d,e), q,nrisk,
                                            nevent,event_time,event_ind))

                               }

                               grr <- function(q){
                                 grr_comb <- NULL
                                 for (i in 1:length(q)){
                                   grr_comb <- c(grr_comb,Kg(q[i],nrisk[[i]],nevent[[i]],event_time[[i]],event_ind[[i]]))
                                 }
                                 return(grr_comb)


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
                               id_B <- trtB %>% filter(ind ==1)
                               ## Split by strata
                               id_B <- as.numeric(factor(id_B$group))
                               nrisk <- split(summary(fitB)$n.risk,id_B)
                               nevent <- split(summary(fitB)$n.event,id_B)
                               id_full <- as.numeric(factor(trtB$group))

                               event_time <- split(trtB$time,id_full)
                               event_ind <- split(trtB$ind,id_full)

                               ## constrained survival prob
                               q_optim <- constrOptim(start_surv, fn, grad=grr,
                                                      ui = ui, ci = rep(0,nrow(ui)))$par
                               S_B <- exp(q_optim)



                             } else if (study == "Origin"){

                               nrisk <- split(summary(fit)$n.risk,id_A)
                               nevent <- split(summary(fit)$n.event,id_A)

                               event_time <- split(summary(fit)$time,id_A)
                               surv_prob <- split(summary(fit)$surv,id_A)

                               S_A <- NULL
                               for (i in 1:length(event_time)){

                                 if(length(event_time[[i]])==0){
                                   S_A_i=1
                                 } else{
                                   if (max(event_time[[i]])<x){
                                     S_A_i=min(surv_prob[[i]])
                                   } else{
                                     S_A_i <- summary(fit,t=x)$surv[i]
                                   }
                                 }

                                 S_A <- c(S_A,S_A_i)
                               }


                               nrisk <- split(summary(fit)$n.risk,id_B)
                               nevent <- split(summary(fit)$n.event,id_B)

                               event_time <- split(summary(fit)$time,id_B)
                               surv_prob <- split(summary(fit)$surv,id_B)

                               S_B <- NULL
                               for (i in 1:length(event_time)){

                                 if(length(event_time[[i]])==0){
                                   S_B_i=1
                                 } else{
                                   if (max(event_time[[i]])<x){
                                     S_B_i=min(surv_prob[[i]])
                                   } else{

                                     S_B_i <- summary(fit,t=x)$surv[i]
                                   }
                                 }

                                 S_B <- c(S_B,S_B_i)
                               }


                             } else{
                               stop("This is an error message.")
                             }

                             if(all(S_A>S_B+d_diff)){
                               corr <- corr+1
                            }
                            if(all(S_A<S_B+d_diff)){
                              err <- err+1
                            }


                            data.frame(S_A,S_B,corr,err)
                          }

  return (surv_estimate)

}



