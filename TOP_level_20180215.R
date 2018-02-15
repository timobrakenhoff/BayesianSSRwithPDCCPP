## Script used for analysis and plots of the simulation and example
## Created on 15-02-2018 by T Brakenhoff (t.brakenhoff@gmail.com)

## Manuscript: "Bayesian sample size re-estimation using power priors"
## Authors: T.B. Brakenhoff, K.C.B. Roes and S. Nikolakopoulos

## See manuscript for details

### This script is a top level function which calls all necessary 
# functions to run the simulation and analyze the example dataset.

#Set your directory

#Packages
require(ggplot2)
require(xtable)
library(extrafont)
require(xtable)
require(microbenchmark)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(Rmisc)
require(abind)


#Makes Times New Roman available in R (Necessary for plots)
windowsFonts(Times=windowsFont("TT Times New Roman"))
windowsFonts(Palatino=windowsFont("Palatino Linotype"))

#Necessary background functions (sourcable)
source("Simulation_Functions\\func_ss_calc.R")

#SIMULATION STUDY
source("Simulation_Functions\\func_apply_ss7.R")
source("Simulation_Functions\\func_apply_bayes_alt11.R") 
source("Plotting_Functions\\func_apply_plots13.R") 

#WOENSEL SIMULATION
source("Example_Functions\\func_example_woensel.R") #Overall function for example
source("Example_Functions\\func_apply_example.R") #Updating mechanism for example
source("Plotting_Functions\\func_example_plot4.R") 

#################################
## SIMULATIONS Conditions list ##
#################################

list.conds <- list(
  
  DE.data  = 0.6,              # Effect size in the simulated data
  int.size = c(0,10,20,40),    #4 different interim sizes (total patients)
  prior    = c(10,20),         #2 different priors (total patients)  
  sd.data  = c(0.75,1,1.5),    #Standard deviation in the data
  sd       = 1,                #SD assumed by design
  alpha    = 0.05,             #Alpha for frequentist ss calculation
  power    = 0.80,             #Power for freq ss calculation
  
  DE       = 0.6,              #Effect size assumed by design
  eta      = 0.95,             #Eta parameter
  zeta     = 0.80,             #Zeta parameter
  xi       = 0.90,             #Value set for xi
  d.star   = 0.6,              #NOT RELEVANT: Is overwritten by DE in simulation
  u0.c     = 0,                #Mean value for the control group (mean for u0.e is = DE)
  r        = 1,                #Not a relevant parameter
  
  interims = 1,                #Amount of interims
  sims     = 100000,           #Amount of iterations in simulation
  seed     = 99,
  per.sp   = 0.01,             #By which percentage steps to collect sample size and assess xi
  tot.sp   = 2,                #Maximum factor of re-estimated sample size that is collected. 
  d.pts    = 800               #Amount of data points to simulate in total (more than enough) 
)

#####################################
## EXAMPLE WOENSEL Conditions list ##
#####################################

list.ex <- list(
  
  DE.data  = 0.746,                #This is actually a DE of 1.5 scaled by an sd of 2.01 --> 1.5/2.01=0.746
  int.size = c(0,20,50),           #3 different interim sizes (total patients)
  prior    = 14,                   #Fixed prior based on woensel (total patients)
  sd.data  = c(0.75,1,1.5,2),      #Standard deviation in simulated data
  sd       = 1,                    #SD assumed by design
  alpha    = 0.05,                 #Alpha for frequentist ss calculation
  power    = 0.80,                 #Power for freq ss calculation
  
  DE       = 0.355,                #Effect size assumed by design based on woensel
  eta      = 0.95,                 #Eta parameter
  zeta     = 0.80,                 #Zeta parameter
  xi       = 0.90,                 #Value set for xi
  d.star   = 0.5,                  #NOT RELEVANT: Is overwritten by DE in simulation
  u0.c     = 0,                    #Mean value for the control group (mean for u0.e is = DE)
  r        = 1,                    #Not a relevant parameter
  
  interims = 1,                    #Amount of interims
  sims     = 100000,               #Amount of iterations in simulation
  seed     = 88,
  per.sp   = 0.01,                 #By which percentage steps to collect sample size and assess xi
  tot.sp   = 2,			               #Maximum factor of re-estimated sample size that is collected. 
  d.pts    = 1000 		             #Amount of data points to simulate in total (more to be safe)   
  )

#Initial sample sizes for woensel
ss.calc(a0=7,b0=28.3,eta=0.95,zeta=0.80,xi=0.90,r=1,q0.e=0.01,q0.c=0.01,d.star=1.5)[4]
#78

#######################
## SIMULATION OUTPUT ##
#######################

## Run simulation function (alternate fc.vec)
RES2 <- sim.conditions4(list.conds,fc.vec=c(0.2,0.4,0.6,0.8))

#Results of above simulation can be found in: 
# Simulations04-03-16_RES2 _ ORIGINAL SIM RESULTS

## To save data
#save.image(file="Simulations02-02-16_RES2.RData")

## Make table with succes and futility probabilities
jj <- RES2$g.cond[,c(2,3,4,5,9,10,11,16,21,22,25,26)]
ll <- cbind(jj,RES2$sim.res)

xtable(ll[,1:11])
xtable(ll[,12:22])

## Plotting
sim.res <- RES2$sim.res
g.cond  <- RES2$g.cond
perleft <- RES2$cum.list

############################
## EXAMPLE WOENSEL OUTPUT ##
############################

#Run woensel simulation
WOE.RES <- ex.conds.sim(list.ex,fc.seq=seq(0.05,0.95,0.1))

## To save data
#save.image(file="Ex_100000_28-04-16.RData")


#####################
## MANUAL PLOTTING ##
#####################

#####
## Simulation
#####

#load data


#set parameters for the plotting function
sim.res <- RES2$sim.res
g.cond  <- RES2$g.cond
perleft <- RES2$cum.list

#Run plotting function
sim.plots <- apply_plots4(sim.res,g.cond,save=T,perleft)


### Put figure 2 and 3 together in ggplot

#take objects from PlotDataSimEPS
out.sel
Aout.sel

#change name of second object so col names match
names(Aout.sel)[4] <- "res2"

#Combine dfs
tot.sel.all <- rbind(out.sel,Aout.sel)

#add col to tot.sel.all
tot.sel.all$PDC <- factor(rep(c("Uncalibrated","Calibrated"),each=nrow(Aout.sel)),
                          levels=c("Uncalibrated","Calibrated"))

#Make labels for facet grid
levels(tot.sel.all$prior)  <- c(expression("n"[0] == 10), expression("n"[0] == 20))
#levels(tot.sel.all$PDC)    <- c("Uncalibrated","Calibrated")

#Make plot using total dataset
tot.sel.plot <- ggplot(tot.sel.all, aes(x=res2, group=interaction(sd.data,int.size),
                     colour=sd.data,linetype=int.size))+ 
                stat_ecdf(size=1.2)+
                scale_linetype_manual(values=c("solid", "dotted","dashed","dotdash"))+
                geom_hline(aes(yintercept=0.90),linetype="dashed")+
                geom_vline(aes(xintercept=1.0),linetype="dashed")+
                facet_grid(PDC ~ prior,labeller=label_parsed) +
                theme_bw()+
                xlim(0,2) +
                scale_colour_grey(name=expression(sigma[R]/sigma[0]),start=0.35,end=0.65,guide=guide_legend(keywidth=2.5))+
                theme(text=element_text(size=20,family="Palatino"),
                      legend.title=element_text(face="bold",size=22),
                      legend.title.align=0)+
                scale_linetype_discrete(name=expression(paste(N[I])),guide=guide_legend(keywidth=2.5))+
                labs(x="Ratio of re-estimated sample size", y=expression(xi[emp]))



#Save plot
ggsave(tot.sel.plot,filename="plot_all_XI.pdf",width=19,height=20.4)

#####
## EXAMPLE
#####

#load data


#set parameters
col.mat <- WOE.RES$col.mat

#Run plotting function
WOE.plots <- plot.example(col.mat,save=F) 

