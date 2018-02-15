## FUNCTION TO MAKE GRAPHS USING OUTPUT FROM func_apply_ss.R
## INPUT is object from previous function

# NEW PLOTTING FUNCTION CREATED ON 29-09-2017 TO SOLVE DEPRACATION OF LABELLER FUNCTION
# THIS FUNCTION CREATES THE PLOTS FOR SUBMISSION TO SMMR

apply_plots4 <- function(sim.res,g.cond,save=T,perleft){
  
out.list <- list()
  
## N LEFT PLOTS (WITHOUT AND WITH PP) (THINK OF WAY TO INCLUDE LINE FOR xi=90)

#Select nleft and scenarios
res     <- sim.res[,c("Ntot","Ntot.10","Ntot.90")]
Ares    <- sim.res[,c("Atot","Atot.10","Atot.90")]
scen    <- g.cond[,2:4]

#Make total df
tot.mat.n <- cbind(scen,res)
Atot.mat.n <- cbind(scen,Ares)

#make into factors
tot.mat.n$prior    <- as.factor(tot.mat.n$prior)
tot.mat.n$sd.data  <- as.factor(tot.mat.n$sd.data)

Atot.mat.n$prior    <- as.factor(Atot.mat.n$prior)
Atot.mat.n$sd.data  <- as.factor(Atot.mat.n$sd.data)

#Labeling of panes (fixed to now include prior label and equal sign)
levels(tot.mat.n$prior)  <- c(expression("n"[0] == 10), expression("n"[0] == 20))
levels(Atot.mat.n$prior) <- c(expression("n"[0] == 10), expression("n"[0] == 20))

#See earlier versions for colour plots code

nleft.plot.GS <- ggplot(data=tot.mat.n, aes(x=int.size, y=Ntot, group=sd.data, 
                                            shape=sd.data,colour=sd.data)) +
                geom_line(aes(linetype=sd.data),size=1.4) +
                geom_point(size=5) +
                geom_ribbon(aes(ymin=Ntot.10,ymax=Ntot.90,linetype=sd.data,colour=sd.data,fill=sd.data),alpha=0.2) +
                facet_grid(~prior,labeller=label_parsed) +
                scale_colour_grey(name=expression(sigma[R]/sigma[0]),start=0.3,end=0.7,guide=guide_legend(keywidth=2.5))+
                scale_fill_grey(name=expression(sigma[R]/sigma[0]),start=0.3,end=0.7,guide=guide_legend(keywidth=2.5))+
                theme_bw() +
                labs(x="Interim timing (total amount of patients)", y="Total required sample size at interim")+ 
                scale_shape_discrete(name=expression(sigma[R]/sigma[0]))+
                scale_linetype_discrete(name=expression(sigma[R]/sigma[0]))+
                theme(text=element_text(size=20,family="Palatino"),legend.title=element_text(face="bold",size=22),
                      legend.title.align=0)

Anleft.plot.GS <- ggplot(data=Atot.mat.n, aes(x=int.size, y=Atot, group=sd.data, 
                                            shape=sd.data,colour=sd.data)) +
                  geom_line(aes(linetype=sd.data),size=1.4) +
                  geom_point(size=5) +
                  geom_ribbon(aes(ymin=Atot.10,ymax=Atot.90,linetype=sd.data,colour=sd.data,fill=sd.data),alpha=0.2) +
                  facet_grid(~ prior,labeller=label_parsed) +
                  scale_colour_grey(name=expression(sigma[R]/sigma[0]),start=0.3,end=0.7,guide=guide_legend(keywidth=2.5))+
                  scale_fill_grey(name=expression(sigma[R]/sigma[0]),start=0.3,end=0.7,guide=guide_legend(keywidth=2.5))+
                  theme_bw() +
                  labs(x="Interim timing (total amount of patients)", y="Total required sample size at interim")+                       
                  scale_shape_discrete(name=expression(sigma[R]/sigma[0]))+
                  scale_linetype_discrete(name=expression(sigma[R]/sigma[0]))+
                  theme(text=element_text(size=20,family="Palatino"),legend.title=element_text(face="bold",size=22),
                        legend.title.align=0)


## XI PLOTS

scen2     <- g.cond[,2:4]
tot.mat   <- cbind(id=seq(1,nrow(scen2)),scen2,sim.res)  


# Select eq5 pers
res2  <- as.vector(abind(lapply(perleft,function(x) t(x)[1:4,]),along=1))
Ares2 <- as.vector(abind(lapply(perleft,function(x) t(x)[5:8,]),along=1))


#Repeat matrix as many times as sims
scen2.rep <- do.call(rbind, replicate(g.cond$sims[1], scen2, simplify=FALSE)) # where m is your matrix

#make plot df
out.df  <- cbind(scen2.rep,res2)
Aout.df <- cbind(scen2.rep,Ares2)

#Make into factors
out.df$prior    <- as.factor(out.df$prior)
out.df$int.size <- as.factor(out.df$int.size)
out.df$sd.data  <- as.factor(out.df$sd.data)

Aout.df$prior    <- as.factor(Aout.df$prior)
Aout.df$int.size <- as.factor(Aout.df$int.size)
Aout.df$sd.data  <- as.factor(Aout.df$sd.data)

#Filter out red line
out.sel  <- out.df[out.df$sd.data!=0.75,]
Aout.sel <- Aout.df[Aout.df$sd.data!=0.75,]

#Make labels for facet grid
levels(out.sel$prior)  <- c(expression("n"[0] == 10), expression("n"[0] == 20))
levels(Aout.sel$prior) <- c(expression("n"[0] == 10), expression("n"[0] == 20))



xi.plot.GS <- ggplot(out.sel, aes(x=res2, group=interaction(sd.data,int.size),
                               colour=sd.data,linetype=int.size))+ 
              stat_ecdf(size=1.2)+
              scale_linetype_manual(values=c("solid", "dotted","dashed","dotdash"))+
              geom_hline(aes(yintercept=0.90),linetype="dashed")+
              geom_vline(aes(xintercept=1.0),linetype="dashed")+
              theme_bw()+
              xlim(0,2) +
              scale_colour_grey(name=expression(sigma[R]/sigma[0]),start=0.35,end=0.65,guide=guide_legend(keywidth=2.5))+
              theme(text=element_text(size=20,family="Palatino"),
                    legend.title=element_text(face="bold",size=22),
                    legend.title.align=0)+
              scale_linetype_discrete(name=expression(paste(N[I])),guide=guide_legend(keywidth=2.5))+
              facet_grid(~ prior,labeller=label_parsed) +            
              labs(x="Ratio of re-estimated sample size", y=expression(xi[emp]))

Axi.plot.GS <- ggplot(Aout.sel, aes(x=Ares2, group=interaction(sd.data,int.size),
                                  colour=sd.data,linetype=int.size))+ 
              stat_ecdf(size=1.2)+
              scale_linetype_manual(values=c("solid", "dotted","dashed","dotdash"))+
              geom_hline(aes(yintercept=0.90),linetype="dashed")+
              geom_vline(aes(xintercept=1.0),linetype="dashed")+
              facet_grid(~ prior,labeller=label_parsed) +
              theme_bw()+
              xlim(0,2) +
              scale_colour_grey(name=expression(sigma[R]/sigma[0]),start=0.35,end=0.65,guide=guide_legend(keywidth=2.5))+
              theme(text=element_text(size=20,family="Palatino"),
                    legend.title=element_text(face="bold",size=22),
                    legend.title.align=0)+
              scale_linetype_discrete(name=expression(paste(N[I])),guide=guide_legend(keywidth=2.5))+
              labs(x="Ratio of re-estimated sample size", y=expression(xi[emp]))


if(save==T){
  ggsave(nleft.plot.GS,filename="plot_nleft_GS.pdf",width=19,height=10.2)
  ggsave(Anleft.plot.GS,filename="plot_Anleft_GS.pdf",width=19,height=10.2)
  ggsave(xi.plot.GS,filename="plot_xi_GS.pdf",width=19,height=10.2)
  ggsave(Axi.plot.GS,filename="plot_Axi_GS.pdf",width=19,height=10.2)}

out.list <- list(nleft.plot,Anleft.plot,xi.plot,Axi.plot)

out.list
}