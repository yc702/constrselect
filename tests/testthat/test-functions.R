test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


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
  result <- pickwin_bin_multiple(n = 50, pa_list = c(0.25,0.28,0.28),
                                 D=c(0.15,0.15,0.15),d=c(0.05,0.05,0.05),
                                 prop.strat=c(0.3,0.3,0.4),study="Constrained",
                                 S = 1000,cluster=6,order_list=list(1,c(2,3)))
  expect_vector(result)
  expect_length(result,2)
})
