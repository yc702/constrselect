test_that("Check partial_order",{
  A <- c(-1,1,0)
  B <- c(-1,0,1)
  C <- as.matrix(rbind(A,B))
  dimnames(C)[[1]] <- NULL
  expect_equal(partial_order(list(1,c(2,3))),
               C)
  expect_error(partial_order(c(1,2,3)),"The input or order_list must be a list")
})


test_that("Check order_constrain",{
  expect_equal(order_constrain(r=c(8,10),n=c(20,25),order_list=list(1,2)),
               c(0.4,0.4))
})


test_that("Check pickwin_bin_exact",{
  result <- pickwin_bin_exact(n = 50, p1 = 0.25, strata_diff = 0.05,
                              D=c(0.15,0.15),d=c(0.05,0.05),
                              prop.strat=0.4,study="Constrained")
  expect_vector(result)
  expect_length(result,2)
})

test_that("Check pickwin_bin_multiple",{
  pa_list = c(0.25,0.28,0.28)
  result <- pickwin_bin_multiple(n = 50, pa_list = c(0.25,0.28,0.28),
                                 D=c(0.15,0.15,0.15),d=c(0.05,0.05,0.05),
                                 prop.strat=c(0.3,0.3,0.4),study="Constrained",
                                 S = 100,cluster=2,order_list=list(1,c(2,3)))
  expect_s3_class(result,"data.frame")
  expect_equal(dim(result)[2],2*length(pa_list)+2)
})


test_that("Check Mg",{
  expect_equal(Mg(x = 6, event_time=c(0.41, 0.8, 2.2, 2.5, 3, 4.5, 6, 6.2, 6.5, 7),
                  event_ind=c(1,1,1,1,1,1,0,0,0,0)),6)
})



test_that("Check Ng",{
  expect_equal(Ng(x=6,event_time=c(0.41, 0.8, 2.2, 2.5, 3, 4.5, 6, 6.2, 6.5, 7)),4)
})

test_that("Check Kg",{
  result <- Kg(x = 6,qn = -0.1,nrisk = c(15, 14, 13, 12, 11, 10),nevent = c(1, 1, 1, 1, 1, 1),
               event_time = c(0.47, 0.84, 2.25, 2.5, 2.8, 4.5, 6, 6.16, 6.6, 6.86, 8.17),
               event_ind = c(1, 1, 1, 1, 1, 1, 0, 0,0,0,0))
  expect_equal(round(result),30)
  expect_error(Kg(x = 6,qn = 0.1,nrisk = c(15, 14, 13, 12, 11, 10),nevent = c(1, 1, 1, 1, 1, 1),
                  event_time = c(0.47, 0.84, 2.25, 2.5, 2.8, 4.5, 6, 6.16, 6.6, 6.86, 8.17),
                  event_ind = c(1, 1, 1, 1, 1, 1, 0, 0,0,0,0)),"qn needs to be smaller than 0")
})


test_that("Check Ng",{
  expect_equal(Ng(x=6,event_time=c(0.41, 0.8, 2.2, 2.5, 3, 4.5, 6, 6.2, 6.5, 7)),4)
})


test_that("Check Kg",{
  result <- llkhd(x=6,qn = -0.1,nrisk = c(15, 14, 13, 12, 11, 10),nevent = c(1, 1, 1, 1, 1, 1),
                     event_time = c(0.47, 0.84, 2.25, 2.5, 2.8, 4.5, 6, 6.16, 6.6, 6.86, 8.17),
                     event_ind = c(1, 1, 1, 1, 1, 1, 0, 0,0,0,0))
  expect_error(llkhd(x=6,qn = 0.1,nrisk = c(15, 14, 13, 12, 11, 10),nevent = c(1, 1, 1, 1, 1, 1),
                     event_time = c(0.47, 0.84, 2.25, 2.5, 2.8, 4.5, 6, 6.16, 6.6, 6.86, 8.17),
                     event_ind = c(1, 1, 1, 1, 1, 1, 0, 0,0,0,0)),"qn needs to be smaller than 0")
  expect_equal(round(result),-24)
})


test_that("Check sim_surv",{
  result <- sim_surv(nmax=12,arrival_rate=4,event_rate=0.08,FUP=6)
  expect_type(result,"double")
  expect_equal(dim(result)[2],2)

})


test_that("Check pickwin_surv_fun",{
  event_rate_A=c(0.08,0.05, 0.05)
  result <- pickwin_surv_fun(maxn=50,prop=c(0.3,0.3,0.4),event_rate_A=c(0.08,0.05, 0.05),
                             trt_diff=c(0.1,0.1,0.1),d=c(0.05,0.05,0.05), arrival_rate=3,FUP=6,
                             x=6,S=100,study = "Constrained",cluster=2,order_list=list(1,c(2,3)),with_seed = 111)
  expect_s3_class(result,"data.frame")
  expect_equal(dim(result)[2],2*length(event_rate_A)+2)

})

