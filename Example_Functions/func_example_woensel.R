
ex.conds.sim <- function(list.ex,fc.seq=seq(0.05,0.95,0.1)){

  
    
#Expand grid of conditions
g.cond <- expand.grid(list.ex)  

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


#SET SEED
set.seed(g.cond$seed[1])

#Create results matrix
#sim.res  <- matrix(NA,nrow=nrow(g.cond),ncol=8)

#Collect all sim results
sim.list <- list()
tm.list  <- list()
cum.list <- list()
col.mat  <- NULL

#Indices to get nu.alpha and nu.beta to simulate the data
sdd.dat <- c(match(unique(g.cond$sd.data),g.cond$sd.data),nrow(g.cond))

#Loop over sd.data
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

  #Loop over the different CI intervals
  for (j in seq_along(fc.seq)){
    
    #counter
    cat(paste("i=",i," and j=",j),"\n")
    
    #Apply sim.bra.alt4 over all rows of the simulated data
    ap.res <- t(apply(t.m,1,function(tm1) sim.app.ex(tm1,unique(g.cond$int.size),g.cond[sdd.dat[i],],fc.vec=fc.seq[j]))) 
    
    #Save all results in list
    sim.list[[(i-1)*length(fc.seq)+j]] <- ap.res
    tm.list[[(i-1)*length(fc.seq)+j]]  <- t.m
    
    #Give ap.res columns names
    colnames(ap.res) <- c(paste0("N.tot_",unique(g.cond$int.size)),
                          paste0("A.N.tot_",unique(g.cond$int.size)),
                          paste0("per.left_",unique(g.cond$int.size)),
                          paste0("A.per.left_",unique(g.cond$int.size)),
                          paste0("suc_",unique(g.cond$int.size)),
                          paste0("fut_",unique(g.cond$int.size)),
                          paste0("scal_",unique(g.cond$int.size)[-1]),
                          paste0("M_",unique(g.cond$int.size)[-1]),
                          paste0("ci.u_",unique(g.cond$int.size)[-1]),
                          paste0("fc.int_",unique(g.cond$int.size)[-1]))
                        
    
    #Store per.left probabilities in list
    cum.list[[(i-1)*length(fc.seq)+j]] <- ap.res[,7:12]
    
    #Extract point of re-est ss where xi=0.90
    xi.sim.90 <- matrix(apply(ap.res[,7:12],2,function(x) quantile(x,.90)),ncol=2)
    
    #Results matrix (only int.size prior sd.data and per.left with and without adap)
    res.mat <- cbind(g.cond[sdd.dat[i]:(sdd.dat[i]+2),c(2,3,4)],rep(fc.seq[j],nrow(xi.sim.90)),xi.sim.90)
    
    #Bind al res.mats together in col.mat
    col.mat <- rbind(col.mat,res.mat)
    
    
### BULK ###    
    #NTOT AND NLEFT  
      #Take mean and quantiles of normal and adaptation with int size as rows 
      #m.N.tot    <- t(apply(ap.res[,1:3],2,function(x) c(mean(x),quantile(x,c(0.10,0.90)))))
      #m.N.a.tot  <- t(apply(ap.res[,4:6],2,function(x) c(mean(x),quantile(x,c(0.10,0.90)))))
      
    #SUCFUT
      #Take mean of sucfut
      #m.sucfut   <- apply(ap.res[,13:18],2,function(x) mean(x))
      
      #Make matrix with nrow as int.size
      #mat.sucfut <- matrix(c(rep(m.sucfut[1],5),m.sucfut[2:3],rep(m.sucfut[4],5),m.sucfut[5:6]),nrow=length(unique(g.cond$int.size)))
      #mat.sucfut <- cbind(m.sucfut[1:3],m.sucfut[4:6])
      
    #Store averaged results in sim.res FIX THIS
      #sim.res[ind.j[j]:(ind.j[j]+length(unique(g.cond$int.size))-1),] <- cbind(m.N.tot,m.N.a.tot,mat.sucfut)
###
      
      
    #Names for cum.list
    names(cum.list)[(i-1)*length(fc.seq)+j] <- paste0("sd.data=",unique(g.cond$sd.data)[i],"_fc.int=",fc.seq[j])
  }
}

colnames(col.mat) <- c("int.size","prior","sd.data","f.ci","perleft","Aperleft")

out.list <- list(cum.list=cum.list,col.mat=col.mat)

out.list
}


                                                                                                                    