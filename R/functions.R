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
#' @seealso
#'  \code{\link[parallel]{makeCluster}}
#'  \code{\link[doParallel]{registerDoParallel}}
#'  \code{\link[base]{NULL}}, \code{\link[base]{Random}}
#'  \code{\link[foreach]{foreach}}
#' @rdname pickwin_bin_multiple
#' @export
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom base is.null set.seed
#' @importFrom foreach foreach
pickwin_bin_multiple <- function(n, pa_list,
                                 D=c(0.15,0.15,0.15), d_diff = c(0.05,0.05,0.05),
                                 prop.strat=c(0.2,0.3,0.5),study="Constrained",
                                 S,cluster = 6,order_list,with_seed=NULL) {
  Nk<- ceiling(prop.strat*n)
  cl <- parallel::makeCluster(cluster)
  doParallel::registerDoParallel(cl)

  if(!base::is.null(with_seed)){base::set.seed(with_seed, kind = "L'Ecuyer-CMRG")}

  bin_estimator <- foreach::foreach(1:S, .combine = rbind) %dopar% {
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




