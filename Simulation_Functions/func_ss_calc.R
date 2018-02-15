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
