## SIMULATION FUNCTION USING APPLY
## FUNCTION LOOPS OVER ALL SCENARIOS
## Uses func_apply_bayes.R for apply over simulations
## Uses func_apply_plots.R for making graphs 

sim.conditions4 <- function(list.conds,fc.vec){

#Expand grid of conditions
  g.cond <- expand.grid(list.conds)  

#Add variables to grid that need computing
  g.cond$a0       <- g.cond$prior/2
  g.cond$nu.alpha <- g.cond$prior/2
  g.cond$q0.e     <- 0.01 #Uninformative
  g.cond$q0.c     <- 0.01 #Uninformative
  g.cond$b0       <- g.cond$a0*g.cond$sd^2
  g.cond$nu.beta  <- g.cond$nu.alpha*g.cond$sd.data^2
  g.cond$d.star   <- g.cond$DE
  g.cond$d.data   <- g.cond$DE.data
  g.cond$u0.e     <- g.cond$d.data
  
  #set seed
  set.seed(g.cond$sims[1])
  
  #Create results matrix
  sim.res  <- matrix(NA,nrow=nrow(g.cond),ncol=10)
  
  #Save all info
    sim.list <- list()
    tm.list  <- list()
    cum.list <- list()

  #Notification of start sim
  cat("Starting simulation (",nrow(g.cond),"conditions )","\n")
  
  #Indices to get nu.alpha and nu.beta to simulate the data
  sdd.dat <- c(match(unique(g.cond$sd.data),g.cond$sd.data),nrow(g.cond))

#Generate data for the 3 unique sd.data (The simulated data is only different for these 3 conditions For the rest same data can be used)
  for (i in 1:length(unique(g.cond$sd.data))){
   
    #Generate control data points matrix 
    m1     <- matrix(rnorm((g.cond$d.pts[sdd.dat[i]]/2)*g.cond$sims[sdd.dat[i]],g.cond$u0.c[sdd.dat[i]],
                           sqrt(1/(g.cond$nu.alpha[sdd.dat[i]]/g.cond$nu.beta[sdd.dat[i]]))),
                           nrow=g.cond$sims[sdd.dat[i]],byrow=T)
    
    
    #Generate exp data points matrix
    m2     <- matrix(rnorm((g.cond$d.pts[sdd.dat[i]]/2)*g.cond$sims[sdd.dat[i]],g.cond$u0.e[sdd.dat[i]],
                            sqrt(1/(g.cond$nu.alpha[sdd.dat[i]]/g.cond$nu.beta[sdd.dat[i]]))),
                            nrow=g.cond$sims[sdd.dat[i]],byrow=T)
    
    #cbind mat to form total sim dataset
    t.m    <- cbind(m1,m2)  
    
    #Ind for second loop (This is to identify the rows of the g.cond where the prior differs)
    ind.j <- c(sdd.dat[i],sdd.dat[i]+length(unique(g.cond$int.size)))
    
    #apply function for two different priors
    for(j in seq_along(ind.j)){
      
      #Sim counter and scenario counter
      cat(paste("i=",i," and j=",j),"\n", paste("sd.data=",g.cond$sd.data[sdd.dat[i]],"prior=",g.cond$prior[ind.j[j]]),"\n")
      
      #Apply sim.bra.alt4 over all rows of the simulated data
      ap.res <- t(apply(t.m,1,function(tm1) sim.BRA.alt4(tm1,unique(g.cond$int.size),g.cond[ind.j[j],],fc.vec=fc.vec))) 
      
      #Save all results in list
      sim.list[[i*2-2+j]] <- ap.res
      tm.list[[i*2-2+j]]  <- t.m
      
      
      #Take mean and quantiles of normal and adaptation with int size as rows
      m.N.tot    <- t(apply(ap.res[,1:4],2,function(x) c(mean(x),quantile(x,c(0.10,0.90)))))
      m.N.a.tot  <- t(apply(ap.res[,5:8],2,function(x) c(mean(x),quantile(x,c(0.10,0.90)))))
      
      #Take mean of sucfut
      m.sucfut   <- apply(ap.res[,17:24],2,function(x) mean(x))
      
      #Make matrix with nrow as int.size
      mat.sucfut <- matrix(c(rep(m.sucfut[1],5),m.sucfut[2:4],rep(m.sucfut[5],5),m.sucfut[6:8]),nrow=4)
      
      #Store averaged results in sim.res
      sim.res[ind.j[j]:(ind.j[j]+length(unique(g.cond$int.size))-1),] <- cbind(m.N.tot,m.N.a.tot,mat.sucfut)
      
      #Store per.left probabilities in list
      cum.list[[i*2-2+j]] <- ap.res[,9:16]
      
      #Set column names for cum.list matrix
      colnames(cum.list[[i*2-2+j]]) <- c(paste0("int=",unique(g.cond$int.size)),
                                         paste0("A_int=",unique(g.cond$int.size)))
      }
  }
  
  #Give results column names
  colnames(sim.res) <- c("Ntot","Ntot.10","Ntot.90","Atot","Atot.10","Atot.90","Suc","ASuc","Fut","AFut")

  #Give cum.list names
  names(cum.list) <- paste(rep(unique(g.cond$sd.data),each=length(unique(g.cond$prior))),
                           rep(unique(g.cond$prior),length(unique(g.cond$sd.data))),sep="_")
  #Give cum.list colnames
  
## PLOTS 
plot.res <- apply_plots4(sim.res=sim.res,g.cond=g.cond,perleft=cum.list)

## OUTPUT
out.list <- list('sim.res'=sim.res,'g.cond'=g.cond,'plot.res'=plot.res,
                 'sim.list'=sim.list,'tm.list'=tm.list,'cum.list'=cum.list)
out.list
}