## Script for the calculation of gamma PDCCPP 
## (see section "Prior Data conflict calibrated power priors" and 
## "Application of PDCCPP in Bayesian SSR")

## Created on 15-02-2018 by T Brakenhoff (t.brakenhoff@gmail.com)
## Used in manuscript: "Bayesian sample size re-estimation using power priors"
## Authors: T.B. Brakenhoff, K.C.B. Roes and S. Nikolakopoulos

## See manuscript for details on the methodology and the implementation
## This script specifically facilitates the calculation of the scaling parameter
## (gamma PDCCPP)

## See provided scripts of the simulation and the analysis of the example for 
## full implementation of the PDCCPP SSR.

### Necessary parameters for calculation (obtained from data collected at interim)

fc.int # desired width of interval set for the prior predictive distribution of M.
       # is equal to 1-c/2

k      # sample size that has been collected per group assuming equal allocation to both groups

a1     # is equal to alpha 0 : a parameter of the prior gamma distribution for the variance

M      # Observed value of M determined from the collected sample size (see eq. 3)



### PROCEDURE FOR DETERMINING SCALING PARAMETER (GAMMA PDCCPP) 
#   depending on observed M and the ppd of M

##Find limits of the prior predictive distribution of M with width set by fc.int (= 1-c/2)
M.ci      <- qf(c(0.5-fc.int/2,0.5+fc.int/2),k*2,2*a1)

##Test if M is within M.ci. If TRUE then it is within M.ci
if(M<M.ci[2] & M>M.ci[1]){
  
  #No changing of original alpha0 so scale parameter (gamma PDCCPP) = 1  
  scal <- 1  
  
  #save upper limit of ci
  ci.u <- M.ci[2]
  
} else {
  
  #Test if M is lower than lower limit
  if(M<M.ci[1]){
    
    #Set scal equal such that second degree of freedom equals 1 
    scal <- 1/a1
    
    #save upper limit of ci
    ci.u <- M.ci[2]
    
  } else { #Where M is higher than upper limit
    
    #Use search procedure to determine the scaling parameter (gamma PDCCPP)
    # PRECISION = 0.001
    scal.m  <-seq(0.001,1,0.001)
    
    #Find corresponding a0 values for each potential scaling value
    a1.m    <- a1*scal.m
    
    #Find upper limit of CI of ppd of M using these scaled alpha0 values
    ci.m.u  <- qf(0.5+(fc.int/2),k*2,2*a1.m)
    
    #Find first occurence where observed M is now included in the intervals 
    # (determined using a scaled alpha0) Take the -1 value because it finds first time M falls out
    ind.ci  <- which(M>ci.m.u)[1]-1
    
    #Take found scaling factor and put into scal object
    scal <- scal.m[ind.ci]
    
    #save upper limit of ci
    ci.u[int.cont-1] <- ci.m.u[ind.ci]
  }
}


## The scale parameter obtained can then be used to downscale the prior information
## if it is in conflict with the variance observed in the collected sample.