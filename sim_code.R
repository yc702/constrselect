library(ggplot2)
library(constrselect)
library(dplyr)
library(Hmisc)
library(foreach)
library(doParallel)
options(digits=10)
set.seed(111, kind = "L'Ecuyer-CMRG")
## binary simulation result

n_list <- c(20,30,40,50,60,70)
d_list <- list(c(0.03,0.03),c(0.07,0.07))
rho = c(0,0.5)

## Without Constraints
output_old <- lapply(n_list, function(y)
  lapply(d_list, function(x) pickwin_bin_exact(n = y, p_inf = c(0.35,0.45),
                                               D=c(0.2,0.2),d=x,
                                               prop.strat=0.4,study="Origin")
))

output_old <- bind_rows(output_old, .id = "column_label")
output_old$n <- rep(n_list,each=2)
output_old$d <- rep(c(0.03,0.07),length(n_list))

## With Constraints
output_new <- lapply(n_list, function(y)
  lapply(d_list, function(x) pickwin_bin_exact(n = y, p_inf = c(0.35,0.45),
                                               D=c(0.2,0.2),d=x,
                                               prop.strat=0.4,study="Constrained")
  ))

output_new <- bind_rows(output_new, .id = "column_label")
output_new$n <- rep(n_list,each=2)
output_new$d <- rep(c(0.03,0.07),length(n_list))


## Create dataset to make plot
output_old <- data.frame(output_old)
output_new <- data.frame(output_new)


output_old$Lambda_0 <- output_old$pcorr
output_old$Lambda_05 <- output_old$pcorr + 0.5*output_old$pamb
output_new$Lambda_0 <- output_new$pcorr
output_new$Lambda_05 <- output_new$pcorr + 0.5*output_new$pamb

## create plot Fig1
plot_data <- data.frame(n = rep(output_old$n,4),
                        Lambda = c(output_old$Lambda_0,output_old$Lambda_05,output_new$Lambda_0,output_new$Lambda_05),
                        rho = c(rep(rho,each=length(n_list)*2),rep(rho,each=length(n_list)*2)),
                        d = c(output_old$d,output_old$d,output_new$d,output_new$d),
                        Group = c(rep("Without Constraints",4*length(n_list)),rep("With Constraints",4*length(n_list))))


plot_data$rho <- factor(plot_data$rho, levels = c("0","0.5"),
                        ordered = TRUE, labels=c(expression(paste(rho,"=0")),expression(paste(rho,"=0.5"))))
plot_data$d <- ifelse(plot_data$d==0.03,"(0.03,0.03)","(0.07,0.07)")


ggplot(plot_data, aes(x=n, y=Lambda, group=interaction(Group,d),color = Group,linetype=d))+
  scale_y_continuous(breaks = seq(0.55, 0.95, by = 0.1))+
  theme_classic() +
  facet_wrap(~rho,
             labeller = label_parsed)+
  scale_x_continuous("Sample size")+
  geom_line()+ylab(expression(paste(lambda)))+
  scale_linetype_discrete(name = expression(paste(theta)))+
  geom_point()+ theme(text = element_text(size = 15))



### simulation for survival
maxn_list <- c(20,30,40,50,60,70)
d_list <- list(c(0.03,0.03),c(0.07,0.07))
rho = c(0,0.5)
prop=0.4

## With Constraints
constr_result <- lapply(maxn_list, function(x)
  lapply(d_list,function(y)
    pickwin_surv_fun(maxn=x,prop=c(0.4,0.6),surv_inf=c(0.55,0.65),
                     surv_sup=c(0.75,0.85),
                     d=y, arrival_rate=4,FUP=6,
                     x=6,S=8000,study = "Constrained",cluster=6,
                     order_list=list(1,2),with_seed = 111))
)

constr_result <- bind_rows(constr_result, .id = "column_label")
constr_result$column_label <- as.numeric(constr_result$column_label)
constr_result<- constr_result %>% select(column_label,Corr,Wrong) %>% group_by(column_label)%>% summarise_all(list(sum))
constr_result$case <- do.call(paste, expand.grid(d_list,maxn_list))
constr_result <- constr_result %>% mutate(corr_pct=Corr/8000, wrong_pct=Wrong/8000)
constr_result$amb_pct = 1-constr_result$corr_pct-constr_result$wrong_pct



## Without Constraints
km_result <- lapply(maxn_list, function(x)
  lapply(d_list,function(y)
    pickwin_surv_fun(maxn=x,prop=c(0.4,0.6),surv_inf=c(0.55,0.65),
                     surv_sup=c(0.75,0.85),
                     d=y, arrival_rate=4,FUP=6,
                     x=6,S=8000,study = "Origin",cluster=6,
                     order_list=list(1,2),with_seed = 111))
)


km_result <- bind_rows(km_result, .id = "column_label")
km_result$column_label <- as.numeric(km_result$column_label)
km_result<- km_result %>% select(column_label,Corr,Wrong) %>% group_by(column_label)%>% summarise_all(list(sum))
km_result$case <- do.call(paste, expand.grid(d_list,maxn_list))
km_result <- km_result %>% mutate(corr_pct=Corr/8000, wrong_pct=Wrong/8000)
km_result$amb_pct = 1-km_result$corr_pct-km_result$wrong_pct


## format the data to create plot
sim_data <- rbind(constr_result,km_result)
sim_data <- sim_data %>% data.frame() %>%
  mutate(Group=c(rep("Contr",nrow(km_result)),rep("KM",nrow(km_result))))
sim_data$n <- rep(rep(maxn_list,each=2),2)
sim_data$d <- rep(rep(c(0.03,0.07),length(maxn_list)),2)


sim_data$Lambda_0 <- sim_data$corr_pct
sim_data$Lambda_05 <- sim_data$corr_pct + 0.5*sim_data$amb_pct

## create plot
plot_data <- data.frame(n = c(sim_data$n,sim_data$n),
                        Lambda = c(sim_data$Lambda_0,sim_data$Lambda_05),
                        rho = c(rep(rho,each=length(sim_data$n))),
                        d = c(sim_data$d,sim_data$d),
                        Group = c(ifelse(sim_data$Group=="Contr","With Constraints","Without Constraints"),
                                  ifelse(sim_data$Group=="Contr","With Constraints","Without Constraints")))


plot_data$rho <- factor(plot_data$rho, levels = c("0","0.5"),
                        ordered = TRUE, labels=c(expression(paste(rho,"=0")),expression(paste(rho,"=0.5"))))
plot_data$d <- ifelse(plot_data$d==0.03,"(0.03,0.03)","(0.07,0.07)")


ggplot(plot_data, aes(x=n, y=Lambda, group=interaction(Group,d),color = Group,linetype=d))+theme_classic() +
  facet_wrap(~rho,
             labeller = label_parsed)+
  scale_x_continuous("Sample size")+
  scale_y_continuous(breaks = c(0.65,0.75,0.85,0.95))+
  geom_line()+ylab(expression(paste(lambda)))+
  scale_linetype_discrete(name = expression(paste(theta)))+
  geom_point()+ theme(text = element_text(size = 15))




## case study:
## binary case
n_list <- c(20,25,30,35,40,45,50)
rho = c(0,0.5)

## Without Constraints
output_old <- lapply(n_list, function(x)
  pickwin_bin_exact(n = x, p_inf=c(0.4,0.5),
                    D=c(0.2,0.2),d=c(0.05,0.05),prop.strat=0.7,study="Origin")
  )

output_old <- bind_rows(output_old, .id = "column_label")
output_old$n <- n_list

## With Constraints
output_new <- lapply(n_list, function(x)
  pickwin_bin_exact(n = x,  p_inf=c(0.4,0.5),
                    D=c(0.2,0.2),d=c(0.05,0.05),prop.strat=0.7,study="Constrained")
)

output_new <- bind_rows(output_new, .id = "column_label")
output_new$n <- n_list


output_old <- data.frame(output_old)
output_new <- data.frame(output_new)


output_old$Lambda_0 <- output_old$pcorr
output_old$Lambda_05 <- output_old$pcorr + 0.5*output_old$pamb
output_new$Lambda_0 <- output_new$pcorr
output_new$Lambda_05 <- output_new$pcorr + 0.5*output_new$pamb

## create plot
plot_data <- data.frame(n = rep(output_old$n,4),
                        Lambda = c(output_old$Lambda_0,output_old$Lambda_05,output_new$Lambda_0,output_new$Lambda_05),
                        rho = c(rep(rho,each=length(n_list)),rep(rho,each=length(n_list))),
                        Group = c(rep("Without Constraints",2*length(n_list)),rep("With Constraints",2*length(n_list))))


plot_data$rho <- factor(plot_data$rho, levels = c("0","0.5"),
                        ordered = TRUE, labels=c(expression(paste(rho,"=0")),expression(paste(rho,"=0.5"))))


## make plot
ggplot(plot_data, aes(x=n, y=Lambda, group=Group,color = Group))+theme_classic() +
  facet_wrap(~rho,
             labeller = label_parsed)+
  scale_x_continuous("Sample size")+
  geom_line()+ylab(expression(paste(lambda)))+
  scale_linetype_discrete(name = expression(paste(theta)))+
  geom_point()+ theme(text = element_text(size = 15)) +
  geom_hline(yintercept=0.8, linetype="dashed",
             color = "grey", size=0.5)



## survival case
## simulation result
maxn_list <- c(20,25,30,35,40,45)
rho = c(0,0.5)

## With Constraints
constr_result <- lapply(maxn_list, function(x)
    pickwin_surv_fun(maxn=x,prop=c(0.3,0.7),surv_inf=c(0.6,0.7),
                     surv_sup=c(0.75,0.85),
                     d=c(0.02,0.02), arrival_rate=8,FUP=2,
                     x=2,S=8000,study = "Constrained",cluster=6,
                     order_list=list(1,2),with_seed = 111))


constr_result <- bind_rows(constr_result, .id = "column_label")
constr_result$column_label <- as.numeric(constr_result$column_label)
constr_result<- constr_result %>% select(column_label,Corr,Wrong) %>% group_by(column_label)%>% summarise_all(list(sum))
constr_result$case <- maxn_list
constr_result <- constr_result %>% mutate(corr_pct=Corr/8000, wrong_pct=Wrong/8000)
constr_result$amb_pct = 1-constr_result$corr_pct-constr_result$wrong_pct

## Without Constraints
km_result <- lapply(maxn_list, function(x)
  pickwin_surv_fun(maxn=x,prop=c(0.3,0.7),surv_inf=c(0.6,0.7),
                   surv_sup=c(0.75,0.85),
                   d=c(0.02,0.02), arrival_rate=8,FUP=2,
                   x=2,S=8000,study = "Origin",cluster=6,
                   order_list=list(1,2),with_seed = 111))


km_result <- bind_rows(km_result, .id = "column_label")
km_result$column_label <- as.numeric(km_result$column_label)
km_result<- km_result %>% select(column_label,Corr,Wrong) %>% group_by(column_label)%>% summarise_all(list(sum))
km_result$case <- maxn_list
km_result <- km_result %>% mutate(corr_pct=Corr/8000, wrong_pct=Wrong/8000)
km_result$amb_pct = 1-km_result$corr_pct-km_result$wrong_pct


## format the data for plot
sim_data <- rbind(constr_result,km_result)
sim_data <- sim_data %>% data.frame() %>%
  mutate(Group=c(rep("Contr",nrow(km_result)),rep("KM",nrow(km_result))))
sim_data$n <- rep(maxn_list,2)

sim_data$Lambda_0 <- sim_data$corr_pct
sim_data$Lambda_05 <- sim_data$corr_pct + 0.5*sim_data$amb_pct


## create plot
plot_data <- data.frame(n = c(sim_data$n,sim_data$n),
                        Lambda = c(sim_data$Lambda_0,sim_data$Lambda_05),
                        rho = c(rep(rho,each=length(sim_data$n))),
                        Group = c(ifelse(sim_data$Group=="Contr","With Constraints","Without Constraints"),
                                  ifelse(sim_data$Group=="Contr","With Constraints","Without Constraints")))


plot_data$rho <- factor(plot_data$rho, levels = c("0","0.5"),
                        ordered = TRUE, labels=c(expression(paste(rho,"=0")),expression(paste(rho,"=0.5"))))


ggplot(plot_data, aes(x=n, y=Lambda, group=Group,color = Group))+theme_classic() +
  facet_wrap(~rho,
             labeller = label_parsed)+
  scale_x_continuous("Sample size")+
  geom_line()+ylab(expression(paste(lambda)))+
  scale_linetype_discrete(name = expression(paste(theta)))+
  geom_point()+ theme(text = element_text(size = 15)) +
  geom_hline(yintercept=0.8, linetype="dashed",
             color = "grey", size=0.5)




### Simulation table
## Binary
p_inf_seq <- list(c(0.05,0.15),c(0.2,0.3),c(0.35,0.45),c(0.5,0.6),c(0.65,0.75))
prop_seq=0.4
D = c(0.2,0.2)
d= c(0.05,0.05)
rho = 0.5
## explore different p_inf seq
cl <- makeCluster(6)
registerDoParallel(cl)

## With Constraints
contr_combo <- function(n, p_inf_seq,
                        D,d,prop_seq,sigma,rho,nstart){
  start.p <- foreach(i=1:length(p_inf_seq),.combine = 'rbind',.packages="constrselect") %dopar% {
      n <- nstart
      lambda <- 0
      while (lambda<sigma){
        n <- n+1
        result <- pickwin_bin_exact(n = n, p_inf = p_inf_seq[[i]],
                                    D=D,d=d,
                                    prop.strat=prop_seq,study="Constrained")
        lambda <- result[1]+result[2]*rho
      }

      c(result[1],result[2],lambda,n,p_inf_seq[[i]]+D[1],prop_seq)


    }
  return(start.p)
  on.exit(stopCluster(cl))
}

n_new1 <- contr_combo (n, p_inf_seq,
                       D,d,prop_seq,sigma=0.8,rho = 0,nstart=25)

## Without Constraints
original_combo <- function(n, p_inf_seq,
                           D,d,prop_seq,sigma,rho,nstart){
  start.p <- foreach(i=1:length(p_inf_seq),.combine = 'rbind',.packages="constrselect") %dopar%  {

      n <- nstart+i
      lambda <- 0
      while (lambda<sigma){
        n <- n+1
        result <- pickwin_bin_exact(n = n, p_inf = p_inf_seq[[i]],
                                    D=D,d=d,
                                    prop.strat=prop_seq,study="Origin")
        lambda <- result[1]+result[2]*rho
      }

      c(result[1],result[2],n,p_inf_seq[[i]]+D[1],prop_seq)


    }
  return(start.p)
  on.exit(stopCluster(cl))
}


n_old1 <- original_combo (n, p_inf_seq, D,d,prop_seq,sigma=0.8,rho = 0,nstart=25)





## survival outcome
## FUP = 4,5,6
sigma=0.8;rho = 0;
FUP = 6
prop_seq <- 0.4

surv_inf_seq <- list(c(0.05,0.15),
                           c(0.2,0.3),
                           c(0.35,0.45),
                           c(0.5,0.6),
                           c(0.65,0.75))


surv_sup_seq <- list(c(0.25,0.35),
                           c(0.4,0.5),
                           c(0.55,0.65),
                           c(0.7,0.8),
                           c(0.85,0.95))


d_diff1 <- 0.05
d_diff2 <- 0.05

nstart = c(25,42,53,43,24)

## constrained survival
output <- NULL
result <- NULL
for(i in 1:length(surv_inf_seq)){
  for (j in 1:length(prop_seq)){
    maxn <- nstart[i]
    lambda <- 0

    while (lambda<sigma){
      maxn <- maxn+1
      corr <- 0
      err <- 0
      n1<- ceiling(prop_seq[j]*maxn)
      n2<- maxn-n1
      surv_inf <- surv_inf_seq[[i]]
      surv_sup <- surv_sup_seq[[i]]

      result <- pickwin_surv_fun(maxn=maxn,prop=c(0.4,0.6),
                                 surv_inf=surv_inf,
                                 surv_sup=surv_sup,
                                 d=c(0.05,0.05), arrival_rate=4,FUP=FUP,
                                 x=6,S=8000,study = "Constrained",cluster=6,
                                 order_list=list(1,2),with_seed = 523)
      amb = (8000-sum(result$Corr)-sum(result$Wrong))/8000

      lambda <- mean(result$Corr)+amb*rho
    }

    output <- rbind(output,c(mean(result$Corr),amb,maxn,lambda,
                             surv_inf_seq[[i]],prop_seq[j],
                             surv_sup_seq[[i]],FUP))


  }
}
output <- data.frame(output)


colnames(output) = c("Pcorr","Pamb","Sample Size","Lambda","Hazard2",
                     "Proportion","Strata Diff","Trt Diff","FUP")




## reviewer's feedback: evaluate bias and variance
## Binary
## Example
n_list <- c(30,20)
p_inf_list <- list(c(0.4, 0.5),c(0.35,0.45))

## constrained
constr_output <- lapply(p_inf_list, function(y)
  lapply(n_list, function(x) pickwin_bin_multiple(n = x, p_inf= y, D=c(0.2,0.2),
                                                  d = c(0.05,0.05),
                                                  prop.strat=c(0.4,0.6),study="Constrained",
                                                  S=8000,cluster = 6,order_list=list(1,2),with_seed = 111))
)

constr_result <- bind_rows(constr_output, .id = "column_label")
constr_result <- constr_result %>% mutate(diff1 = P_B1-P_A1, diff2 = P_B2-P_A2)
constr_result<- constr_result %>% select(-Corr,-Wrong) %>% group_by(column_label)%>% summarise_all(list(mean,var))

## original
orig_output <- lapply(p_inf_list, function(y)
  lapply(n_list, function(x) pickwin_bin_multiple(n = x, p_inf= y, D=c(0.2,0.2),
                                                  d = c(0.05,0.05),
                                                  prop.strat=c(0.4,0.6),study="Origin",
                                                  S=8000,cluster = 6,order_list=list(1,2),with_seed = 111))
)

orig_result <- bind_rows(orig_output, .id = "column_label")
orig_result <- orig_result %>% mutate(diff1 = P_B1-P_A1, diff2 = P_B2-P_A2)
orig_result<- orig_result %>% select(-Corr,-Wrong) %>% group_by(column_label)%>% summarise_all(list(mean,var))

save(constr_result,orig_result,file = "review2_bin.RData")




## Survival
n_list=c(30,20)
surv_inf_list <- list(c(0.55,0.65),
                       c(0.35,0.45))

surv_sup_list <- list(c(0.75,0.85),
                      c(0.55,0.65))



## constrained
constr_output <- lapply(n_list, function(x)
  mapply(function(y,z) pickwin_surv_fun(maxn=x,prop=c(0.4,0.6),surv_inf=y,
                                        surv_sup=z,
                                        d=c(0.05,0.05), arrival_rate=4,FUP=6,
                                        x=6,S=8000,study = "Constrained",cluster=6,
                                        order_list=list(1,2),with_seed = 111),
         surv_inf_list,surv_sup_list,SIMPLIFY = FALSE)
)

constr_result <- bind_rows(constr_output, .id = "column_label")
constr_result <- constr_result %>% mutate(diff1 = S_B1-S_A1, diff2 = S_B2-S_A2)
constr_result<- constr_result %>% select(-Corr,-Wrong) %>% group_by(column_label)%>% summarise_all(list(mean,var))


## original
orig_output <- lapply(n_list, function(x)
  mapply(function(y,z) pickwin_surv_fun(maxn=x,prop=c(0.4,0.6),surv_inf=y,
                                        surv_sup=z,
                                        d=c(0.05,0.05), arrival_rate=4,FUP=6,
                                        x=6,S=8000,study = "Origin",cluster=6,
                                        order_list=list(1,2),with_seed = 111),
         surv_inf_list,surv_sup_list,SIMPLIFY = FALSE)
)
orig_result <- bind_rows(orig_output, .id = "column_label")
orig_result <- orig_result %>% mutate(diff1 = S_B1-S_A1, diff2 = S_B2-S_A2)
orig_result<- orig_result %>% select(-Corr,-Wrong) %>% group_by(column_label)%>% summarise_all(list(mean,var))

save(constr_result,orig_result,file = "review2_surv.RData")





## cases where your monotonicity constraints are violated.
## Binary:
p_inf_list <- list(c(0.35,0.45),c(0.35,0.25))
D_list <- list(c(0.2,0.2),c(0.1,0.2),c(-0.1,0.2),c(-0.2,0.2),c(0.2,-0.1),c(0.2,-0.2))
study_list <- c("Constrained","Origin")

constr_output <- lapply(p_inf_list, function(x)
  lapply(D_list, function(y)
    lapply(study_list, function(z)
      pickwin_bin_exact(n = 30, p_inf = x,
                        D=y,d=c(0.05,0.05),
                        prop.strat=0.4,study=z)
)))
constr_result1 <- bind_rows(constr_output[[1]], .id = "column_label")
constr_result2 <- bind_rows(constr_output[[2]], .id = "column_label")


constr_output_bin <- rbind(constr_result1,constr_result2)
constr_output_bin$case <- do.call(paste, expand.grid(study_list,D_list,p_inf_list))
save(constr_output_bin,file = "review3_bin.RData")


## Survival
study_list <- c("Constrained","Origin")
surv_inf_list <- list(c(0.55,0.65),
                      c(0.55,0.65),
                      c(0.55,0.65),
                      c(0.55,0.65),
                      c(0.55,0.65),
                      c(0.55,0.65),

                      c(0.55,0.45),
                      c(0.55,0.45),
                      c(0.55,0.45),
                      c(0.55,0.45),
                      c(0.55,0.45),
                      c(0.55,0.45))

surv_sup_list <- list(c(0.75,0.85),
                      c(0.65,0.85),
                      c(0.45,0.85),
                      c(0.35,0.85),
                      c(0.75,0.55),
                      c(0.75,0.45),

                      c(0.75,0.65),
                      c(0.65,0.65),
                      c(0.45,0.65),
                      c(0.35,0.65),
                      c(0.75,0.35),
                      c(0.75,0.25)
                      )



name_inf_case <- c(0.65,0.45)
name_sup_case <- c(0.75, 0.65, 0.45, 0.35,0.75,0.25)


## constrained
constr_output_surv <- lapply(study_list, function(x)
  mapply(function(y,z) pickwin_surv_fun(maxn=30,prop=c(0.4,0.6),surv_inf=y,
                                        surv_sup=z,
                                        d=c(0.05,0.05), arrival_rate=4,FUP=6,
                                        x=6,S=8000,study = x,cluster=6,
                                        order_list=list(1,2),with_seed = 111),
         surv_inf_list,surv_sup_list,SIMPLIFY = FALSE)
)

constr_result <- bind_rows(constr_output_surv, .id = "column_label")
constr_result$column_label <- as.numeric(constr_result$column_label)
constr_result<- constr_result %>% select(column_label,Corr,Wrong) %>% group_by(column_label)%>% summarise_all(list(sum))
constr_result$case <- do.call(paste, expand.grid(name_sup_case,name_inf_case,study_list))
constr_result <- constr_result %>% mutate(corr_pct=Corr/8000, wrong_pct=Wrong/8000)
constr_result$amb_pct = 1-constr_result$corr_pct-constr_result$wrong_pct
save(constr_result,file = "review3_surv.RData")


