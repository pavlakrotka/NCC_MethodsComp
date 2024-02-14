# CALIBRATION OF TIME MACHINE PRIOR

from_jump_size_to_tau_a_b <- function(e_jump_size = 0.01, max_jump_size = 0.1){
  
  fun_to_solve <- function(tau_b, max_jump_size, e_jump_size){
    pgamma(1 / max_jump_size ^ 2, 1 / (e_jump_size ^ 2) * tau_b, tau_b) - 0.01
    # pgamma(1 / max_jump_size ^ 2, 1 / (e_jump_size ^ 2) * tau_b, 1/tau_b) - 0.01
  }
  
  tau_b <- uniroot(fun_to_solve,
                   c(1e-12, 1e12),
                   max_jump_size = max_jump_size,
                   e_jump_size = e_jump_size)$root
  
  tau_a <- 1 / e_jump_size ^ 2 * tau_b
  
  c(tau_a = tau_a,
    tau_b = tau_b)
}


####
# EXAMPLES
######

max_jump_size = 15
e_jump_size = 10
para_gamma<-from_jump_size_to_tau_a_b (e_jump_size, max_jump_size)
para_gamma[[1]]

x<-rgamma(10000,para_gamma[[1]],para_gamma[[2]])
# plot(x,dist,type="l")
x_std <-sqrt(1/x)
# this plot should be on the scale of e_jump_size
# plot(x_std,dist,type="l")
hist(x_std, xlim=c(0,15),ylim=c(0,80000),main = "max_jump_size = 15
e_jump_size = 10")

######

max_jump_size = 10
e_jump_size = 5
para_gamma<-from_jump_size_to_tau_a_b (e_jump_size, max_jump_size)
para_gamma[[1]]


x<-rgamma(10000,para_gamma[[1]],para_gamma[[2]])
x_std <-sqrt(1/x)
# this plot should be on the scale of e_jump_size
hist(x_std, xlim=c(0,15),ylim=c(0,80000),main = "max_jump_size = 10
e_jump_size = 5")


######

max_jump_size = 1.5
e_jump_size = 1
para_gamma<-from_jump_size_to_tau_a_b (e_jump_size, max_jump_size)
para_gamma[[1]]

x<-rgamma(10000,para_gamma[[1]],para_gamma[[2]])
x_std <-sqrt(1/x)
# this plot should be on the scale of e_jump_size
# hist(x_std, xlim=c(0,15))
hist(x_std, xlim=c(0,15),ylim=c(0,80000),main = "max_jump_size = 1.5
e_jump_size = 1")

######
max_jump_size = 0.15
e_jump_size = 0.01
para_gamma<-from_jump_size_to_tau_a_b (e_jump_size, max_jump_size)

x<-rgamma(10000,para_gamma[[1]],para_gamma[[2]])

x_std <-sqrt(1/x)
# this plot should be on the scale of e_jump_size
# hist(x_std, xlim=c(0,15))
hist(x_std, xlim=c(0,15),ylim=c(0,80000),main = "max_jump_size = 0.15
e_jump_size = 0.01")

######
max_jump_size = 0.015 
e_jump_size = 0.001
para_gamma<-from_jump_size_to_tau_a_b (e_jump_size, max_jump_size)

x<-rgamma(10000,para_gamma[[1]],para_gamma[[2]]) 
x_std <-sqrt(1/x)
# this plot should be on the scale of e_jump_size
# hist(x_std, xlim=c(0,15))
hist(x_std, xlim=c(0,15),ylim=c(0,80000),main = "max_jump_size = 0.015
e_jump_size = 0.001")



######
# TIME MACHINE PAPER
# max_jump_size = 0.015 
# e_jump_size = 0.001
# para_gamma<-from_jump_size_to_tau_a_b (e_jump_size, max_jump_size)

x<-rgamma(10000,0.1,0.01) 
x_std <-sqrt(1/x)
# this plot should be on the scale of e_jump_size
# hist(x_std, xlim=c(0,15))
hist(x_std, xlim=c(0,15),ylim=c(0,80000),main = "a = 0.1
b = 0.01")

## given tau_a and tau_b, what are the expected, 'minimum' and 'maximum' jump in the mean between periods/buckets?
from_tau_a_b_to_jump_size <- function(tau_a, tau_b){
  
  c(e_jump_size = sqrt(1 / (tau_a /  tau_b)),
    min_jump_size = sqrt(1 / qgamma(0.99, tau_a, tau_b)),
    max_jump_size = sqrt(1 / qgamma(0.01, tau_a, tau_b)))
  
}

from_tau_a_b_to_jump_size(0.1,0.01)