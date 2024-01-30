# CALIBRATION OF TIME MACHINE PRIOR

from_jump_size_to_tau_a_b <- function(e_jump_size = 0.01, max_jump_size = 0.1){
  
  fun_to_solve <- function(tau_b, max_jump_size, e_jump_size){
    pgamma(1 / max_jump_size ^ 2, tau_b / (e_jump_size ^ 2) , 1/tau_b) - 0.99
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


max_jump_size = 0.1
e_jump_size = 0.01
para_gamma<-from_jump_size_to_tau_a_b (e_jump_size, max_jump_size)
para_gamma[[1]]

x<-seq(1,1.1*(1/max_jump_size^2),length.out=100)

dist<-dgamma(x,para_gamma[[1]],1/para_gamma[[2]])
#dist<-dgamma(x,0.09639522,963.9522)

plot(x,dist,type="l")
x_std <-sqrt(1/x)
# this plot should be on the scale of e_jump_size
plot(x_std,dist,type="l")


######


max_jump_size = 0.001
e_jump_size = 0.015
para_gamma<-from_jump_size_to_tau_a_b (e_jump_size, max_jump_size)


x<-seq(1,1.1*(1/max_jump_size^2),length.out=100)

dist<-dgamma(x,para_gamma[[1]],1/para_gamma[[2]])
#dist<-dgamma(x,0.09639522,963.9522)

plot(x,dist,type="l")
x_std <-sqrt(1/x)
# this plot should be on the scale of e_jump_size
plot(x_std,dist,type="l")