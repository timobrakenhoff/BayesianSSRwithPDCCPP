##########################################################
## Frequentist calculation of sample size               ##
##  R version of Listing 2 in Whitehead et al. (2008)   ##
##                                                      ##
## Input:  standard deviation, type 1 error (alpha)     ## 
##         type II error (beta), delta star             ##
## Output: total required sample size for design        ##
##########################################################

freq.white.ss <- function(sd,r=1,alpha=0.025,beta=0.20,d.star=1.5){
  
  for (m in 2:10000){                                # Counter
    n <- (r+1)*m                                     # Sample size
    k <- (((r+1)^2)*sd^2)/(r*d.star^2)               
    rhs <- k*(qt(1-alpha,n-2)+qt(1-beta,n-2))^2      # Right hand side
    
    # Output first time that n is greater than or equal to rhs
    if (n >= rhs) {
      results               <- n
      break
      break
    }}
  results
}


##########################################################
## Sample size caculation for unknown variance  (eq. 7) ##
## R version of  Listing 4 in Whitehead et al. (2008)   ##
##                                                      ##
## Input:  alpha_0, beta_0, eta, zeta, xi, q_0E, q_0C   ##
##         allocation ratio (r)                         ##
## Output: Sample size required to achieve              ##
##         design characteristics                       ##
##########################################################

ss.calc <- function(a0=100, b0=4, eta=0.95, zeta=0.90, 
                    xi=0.80,r=1,q0.e=0,q0.c=0,d.star=0.1) {
  
  for (m in 1:1000000) {                                     # Counter  
    n.calc      <- (r+1)*m                                   # Calculated sample size
    a1.hat      <- a0 + 0.5*n.calc                           # Posterior alpha
    t.eta       <- qt(eta, 2*a1.hat)                         # eta point on t-distribution
    t.zeta      <- qt(zeta, 2*a1.hat)                        # zeta point on t-distribution
    b.xi        <- qbeta(xi, 0.5*n.calc,a0)                  # xi point on beta
    
    n.calc.e    <- (r/(r+1))*n.calc                          # sample size in group e and c   
    n.calc.c    <- (1/(r+1))*n.calc
    D.hat       <- (q0.e+n.calc.e)*(q0.c+n.calc.c)/(q0.e+q0.c+n.calc)  # calculation of D
    
    lhs         <- D.hat*(a1.hat/b0)*(1-b.xi)                # left hand side of equation
    rhs         <- ((t.eta+t.zeta)/d.star)^2                 # right hand side
    
    xi.crit     <- 1- (rhs*b0)/(a1.hat*D.hat)                # Corresponding critical value of xi (see eq. 10) 
    xi.int      <- pbeta(xi.crit,(n.calc/2),a0)              # Xi for this design (should be at least 0.80 for sample size calculated here)
    
    #Output only the first row to satisfy the condition
    if (lhs >= rhs) {
      results               <- matrix(c(m,n.calc.e,n.calc.c,n.calc,t.eta,
                                        t.zeta,rhs,xi.crit,
                                        xi.int,a1.hat,D.hat),nrow=1,ncol=11)
      colnames(results)     <- c("Obs", "n.e", "n.c","n","t.eta","t.zeta",
                                 "rhs","xi.crit","xi.int","a1.hat","D")
      break
      break
    }}
  results
}

###############################################
RES2$g.cond

d.star=0.6
eta=0.95
zeta=0.8
xi=0.9
r=1
prior=10
a0=5
b0=5
alpha=0.05
power=0.8
beta=1-power
q0.e=0
q0.c=0
sd=1


## Initial ss calculation

b0 <- ss.calc(a0,b0,eta,zeta,xi,r,q0.e,q0.c,d.star)[4]

f0 <- freq.white.ss(sd,r,alpha,beta,d.star)

## With interim = 10

b1 <- ss.calc(a0=10,b0=10,eta,zeta,xi,r,q0.e=5,q0.c=5,d.star)[4] + 10

f1 <- freq.white.ss(sd,r,alpha,beta,d.star)

## With interim = 20

b2 <- ss.calc(a0=15,b0=15,eta,zeta,xi,r,q0.e=10,q0.c=10,d.star)[4] + 20

f2 <- freq.white.ss(sd,r,alpha,beta,d.star)

## With interim = 30

b3 <- ss.calc(a0=20,b0=20,eta,zeta,xi,r,q0.e=15,q0.c=15,d.star)[4] + 30

f3 <- freq.white.ss(sd,r,alpha,beta,d.star)

## With interim = 40

b4 <- ss.calc(a0=25,b0=25,eta,zeta,xi,r,q0.e=20,q0.c=20,d.star)[4] + 40

f4 <- freq.white.ss(sd,r,alpha,beta,d.star)

## With interim = 50

b5 <- ss.calc(a0=30,b0=30,eta,zeta,xi,r,q0.e=25,q0.c=25,d.star)[4] + 50

f5 <- freq.white.ss(sd,r,alpha,beta,d.star)


resss <- data.frame(int.size=c(0,10,20,30,40,50),freq=c(f0,f1,f2,f3,f4,f5),
           bayes=c(b0,b1,b2,b3,b4,b5))

xtable(resss)
