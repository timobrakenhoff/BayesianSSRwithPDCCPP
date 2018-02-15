## NEW APPLY FUNCTION DOING ALL INTERIMS AT SAME TIME
## OUTPUT:
## row with eq5 fulfillment for all percentages and calculated N and re-estimated total N. 

# tm1 <- t.m[1,]
# s.conds <- g.cond[ind.j[j],]
# int.size <- unique(g.cond$int.size)
# fc.int <- 0.50

# tm1      <- t.m[1,]
# int.size <- unique(g.cond$int.size)
# s.conds  <- g.cond[ind.j[j],]
# 

sim.BRA.alt4 <- function(tm1,int.size,s.conds,fc.vec){
  

# Set starting values for parameters
  eta    <- s.conds$eta
  zeta   <- s.conds$zeta
  xi     <- s.conds$xi
  r      <- s.conds$r
  d.star <- s.conds$d.star
  d.pts  <- s.conds$d.pts
  per.sp <- s.conds$per.sp
  tot.sp <- s.conds$tot.sp
  
  q1.c   <- s.conds$q0.c
  q1.e   <- s.conds$q0.e
  u1.e   <- s.conds$u0.e
  u1.c   <- s.conds$u0.c
  D.int  <- 0
  d1.int <- 0
  
# Make matrices for changing pars
  a1     <- matrix(NA,ncol=4,nrow=length(tm1)/2)
  b1     <- matrix(NA,ncol=4,nrow=length(tm1)/2)
  
# Fill first row of matrices
  a1[1,] <- s.conds$a0
  b1[1,] <- s.conds$b0

# Make Matrices for all n.left and tots
  N.left <- matrix(NA,nrow=length(int.size),ncol=2)
  N.tot  <- matrix(NA,nrow=length(int.size),ncol=2)
  
      
## Calculate initial sample size required (N.left is sample size left EXC. interim)
  N.left[1,] <- ss.calc(a1[1,1], b1[1,1], eta, zeta, xi, r, q1.e, q1.c, d.star)[4]
  N.tot[1,]  <- N.left[1,1]
  
# Set pars for scaling factors
  M       <- rep(NA,length(int.size)-1)
  ci.u    <- rep(NA,length(int.size)-1)
  scal    <- rep(NA,length(int.size)-1)
  fc.int  <- rep(NA,length(int.size)-1)
  
# Collection of eq5
  eq5.col <- matrix(NA,ncol=4,nrow=length(tm1)/2)
  suc.con <- matrix(NA,ncol=4,nrow=length(tm1)/2)
  fut.con <- matrix(NA,ncol=4,nrow=length(tm1)/2)
  
  
  # Update parameters per 1 patient per group
  for(k in 1:(length(tm1)/2)){ # K IS THE SAMPLE SIZE COUNT PER GROUP
    
    #Select data from row for control and treatment
    data.c         <- tm1[1:k]
    data.e         <- tm1[(d.pts/2+1):(d.pts/2+k)]
    
    #Calculate means for the interim
    u.c.mean       <- mean(data.c)
    u.e.mean       <- mean(data.e) 
    
    #Update consistent parameters
    H   <- sum((data.e - u.e.mean)^2) +
           sum((data.c - u.c.mean)^2) + 
           (k*q1.c[1]*(u.c.mean-u1.c[1])^2)/(q1.c[1]+k) +
           (k*q1.e[1]*(u.e.mean-u1.e[1])^2)/(q1.e[1]+k)
    
    q1.c[k+1]        <- q1.c[1] + k           
    q1.e[k+1]        <- q1.e[1] + k
    
    u1.c[k+1]        <- u1.c[1]*(q1.c[1]/q1.c[k+1]) + k*u.c.mean/q1.c[k+1]    
    u1.e[k+1]        <- u1.e[1]*(q1.e[1]/q1.e[k+1]) + k*u.e.mean/q1.e[k+1]   
    
    D.int[k+1]       <- (q1.e[k+1]*q1.c[k+1])/(q1.e[k+1]+q1.c[k+1])
    
    d1.int[k+1]      <- u1.e[k+1]-u1.c[k+1]

# Changing pars
    
    #Matrix for a (for all int.sizes)
    a1[k+1,]         <- a1[1,] + k
    
    #matrix for b (for all int.sizes)
    b1[k+1,]         <- b1[1,] + H/2

# Did we reach an interim?
    if((k*2) %in% int.size){
      
      #Which interim size are we at
      int.cont <- which((k*2)==int.size)
      
      #Regular re-estimation (what is left and total) NO ADAP
      N.left[int.cont,1]         <- ss.calc(a1[k+1,1], b1[k+1,1], eta, zeta, xi, r,
                                            q1.e[k+1], q1.c[k+1], d.star)[4]
      
      N.tot[int.cont,1]          <- N.left[int.cont,1]+k*2
      
      
      #Calculate the observed value for M
      M[int.cont-1]    <- (H*a1[1,int.cont])/((k*2)*b1[1,int.cont])
      
      #Calculate fc.int as ratio of int.size/prior
      int.rat <- (k*2)/s.conds$prior
      
      #Find index of vector with probs
      int.rat.int <- which(int.rat==c(0.5,1,2,4))
      
      #Error if int.rat.int = integer(0)
      if(length(int.rat.int)==0){cat("int.rat not one of 4 options!")}
      
      #Determine fc.int based on vector input in function
      fc.int[int.cont-1] <- fc.vec[int.rat.int]
        
      #Find limits of the ppd of M according to width of interval in fc.int
      M.ci      <- qf(c(0.5-(fc.int[int.cont-1]/2),0.5+(fc.int[int.cont-1]/2)),k*2,2*a1[1,int.cont])
      
      #Test if M is within M.ci. If TRUE then it is within M.ci
      if(M[int.cont-1]<M.ci[2] & M[int.cont-1]>M.ci[1]){
      
        #No changing of original alpha0 so scale parameter = 1  
        scal[int.cont-1] <- 1  
        
        #save upper limit of ci
        ci.u[int.cont-1] <- M.ci[2]
        
      } else {
      
        #Test if M is lower than lower limit
        if(M[int.cont-1]<M.ci[1]){
          
          #Set scal equal such that df2=1 
          scal[int.cont-1] <- 1/a1[1,1]
          
          #save upper limit of ci
          ci.u[int.cont-1] <- M.ci[2]
        
          } else {
          
            #DIfferent scale values for search PRECISION = 0.001
            scal.m  <-seq(0.001,1,0.001)
            
            #Find corresponding a0 values
            a1.m    <- a1[1,1]*scal.m
            
            #Find upper limit of CI of ppd
            ci.m.u  <- qf(0.5+(fc.int[int.cont-1]/2),k*2,2*a1.m)
            
            #Find first occurence where M is included in interval (-1 is because it finds first time M falls out)
            ind.ci  <- which(M[int.cont-1]>ci.m.u)[1]-1
            
            #Get corresponding index for the scal
            scal[int.cont-1] <- scal.m[ind.ci]
            
            #save upper limit of ci
            ci.u[int.cont-1] <- ci.m.u[ind.ci]
          }
        }
        
        #update original a0 with scaled parameter for each int.size
        a1[1,int.cont]          <- a1[1,1]*scal[int.cont-1]
        b1[1,int.cont]          <- b1[1,1]*scal[int.cont-1]
        
        #Update parameters using new a0 and b0
        a1[k+1,int.cont]        <- a1[1,int.cont] + k
        b1[k+1,int.cont]        <- b1[1,int.cont] + H/2
        
        #Re-estimated sample size with scale factor
        N.left[int.cont,2] <- ss.calc(a1[k+1,int.cont], b1[k+1,int.cont], eta, zeta, xi, r,
                                      q1.e[k+1], q1.c[k+1], d.star)[4]
        
        N.tot[int.cont,2]  <- N.left[int.cont,2]+k*2
      }
      
    
    #Evaluate criteria eq5 for all settings
    lhs.int         <- (D.int[k+1]*a1[k+1,])/b1[k+1,]
    t.eta.int       <- qt(eta, 2*a1[k+1,])                         
    t.zeta.int      <- qt(zeta, 2*a1[k+1,])                        
    rhs.int         <- ((t.eta.int+t.zeta.int)/d.star)^2
    eq5.col[k+1,]   <- lhs.int>=rhs.int 
    
    #Evaluate suc or fut
    suc.bin         <- pt(d1.int[k+1]*sqrt(D.int[k+1]*a1[k+1,]/b1[k+1,]),2*a1[k+1,])
    fut.bin         <- pt((d.star-d1.int[k+1])*sqrt(D.int[k+1]*a1[k+1,]/b1[k+1,]),2*a1[k+1,])  
    
    suc.con[k+1,]   <- suc.bin>=eta
    fut.con[k+1,]   <- fut.bin>=zeta
    
    #Break if all values of eq5 are true (for all situations)
    if(all(eq5.col[k+1,])==T){break} 
    
}

# Fill out rest of N.left matrix with 0
N.left[is.na(N.left)==T] <- 0   

# Fill out N.tot with interim sizes
N.tot[which(is.na(N.tot[,1])),] <- int.size[which(is.na(N.tot[,1]))]  

# Look up point where eq5 is true for each situation (takes first true value of each column)
eq5.inds <- apply(eq5.col,2,function(x) match(T,x))-1

#make mat for per.left 
eq5.mat  <- cbind(rep(eq5.inds[1],length(int.size)),eq5.inds)*2

# Calculate percentage of re-est ss
per.left <- rounder((eq5.mat-cbind(int.size,int.size))/N.left,inc=0.01,"ceiling")

#Set Inf values to 0
per.left[per.left==-Inf | per.left==Inf] <- 0

#Collect suc or fut results
sucfut   <- c(diag(suc.con[eq5.inds,]),diag(fut.con[eq5.inds,]))

#Put all results in a vector
tot.res <- c(as.vector(N.tot),as.vector(per.left),sucfut,scal,M,ci.u,fc.int)

tot.res
}
