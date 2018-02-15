## PLOTTING FUNCTION

## FIXED ON 29-09-2017 to update the labelling
## THIS FUNCTION MAKES PLOTS FOR SUBMISSION TO SMMR

# col.mat <- WOE.RES$col.mat


plot.example <- function(col.mat,save=T){

#Make sd.data into factor
col.mat$sd.data <- as.factor(col.mat$sd.data)

#Make f.ci into numeric
col.mat$f.ci <- as.numeric(as.character(col.mat$f.ci))

#Make int.size into factor
col.mat$int.size <- as.factor(col.mat$int.size)

#Labeling of panes (fixed to now include label and equal sign)

levels(col.mat$int.size) <- c(expression("N"[I] == 0), expression("N"[I] == 20),expression("N"[I] == 50))

#SEE OLDER VERSIONS OF THIS FUNCTION FOR COLOUR PLOTS!

### GRAPHS WITH SDD AS LINES ###

sd.line.GS <- 
  
  ggplot(data=col.mat,aes(y=perleft,x=f.ci,group=sd.data,colour=sd.data,shape=sd.data))+
  geom_line(data=col.mat,size=1.2) +
  geom_point(size=5) + 
  facet_grid(.~ int.size,labeller=label_parsed)+ 
  theme_bw() +
  labs(x="Width of predictive interval for M", 
       y=expression(atop(paste("Proportion of total re-estimated sample size"),
                         paste("such that ",xi[emp]," is at least 0.9")))) + 
  scale_colour_discrete(name=expression(sigma[R]/sigma[0]),guide=guide_legend(keywidth=2.5))+
  scale_shape_discrete(name=expression(sigma[R]/sigma[0]),guide=guide_legend(keywidth=2.5))+
  scale_colour_grey(name=expression(sigma[R]/sigma[0]),start=0.2,end=0.8,guide=guide_legend(keywidth=2.5))+
  theme(text=element_text(size=20,family="Palatino"),legend.title=element_text(face="bold",size=22),
        legend.title.align=0)


A.sd.line.GS <- 
  
  ggplot(data=col.mat,aes(y=Aperleft,x=f.ci,group=sd.data,colour=sd.data,shape=sd.data))+
  geom_line(data=col.mat,size=1.2) +
  geom_point(size=5) + 
  facet_grid(~ int.size,labeller=label_parsed)+ 
  theme_bw() +
  labs(x="Width of predictive interval for M", 
       y=expression(atop(paste("Proportion of total re-estimated sample size"),
                         paste("such that ",xi[emp]," is at least 0.9")))) + 
  scale_colour_discrete(name=expression(sigma[R]/sigma[0]),guide=guide_legend(keywidth=2.5))+
  scale_shape_discrete(name=expression(sigma[R]/sigma[0]),guide=guide_legend(keywidth=2.5))+
  scale_colour_grey(name=expression(sigma[R]/sigma[0]),start=0.2,end=0.8,guide=guide_legend(keywidth=2.5))+
  theme(text=element_text(size=20,family="Palatino"),legend.title=element_text(face="bold",size=22),
        legend.title.align=0)

### GRAPHS WITH CI AS LINES ###

# Make selection of CIs to reduce amount of lines   
col.mat.sec <- col.mat[col.mat$f.ci==0.05|
                       col.mat$f.ci==0.35|
                       col.mat$f.ci==0.65|
                       col.mat$f.ci==0.95,]
# Make f.ci into factor
col.mat.sec$f.ci <- as.factor(col.mat.sec$f.ci)  


#This version adds a label to the linegraph  
cols <- c("LINE1"="#000000")

CI.line.GS <- 
              ggplot(data=col.mat.sec,aes(x=sd.data))+
                geom_line(aes(y=perleft,group=1,colour="LINE1"),size=1.2) +
                geom_point(aes(y=perleft,group=1,colour="LINE1"),size=5) + 
                scale_colour_manual(name="PI width",values=cols,guide=guide_legend(keywidth=2.5),
                                    labels=expression(infinity))+
                facet_grid(~ int.size,labeller=label_parsed)+ 
                theme_bw() +
                labs(x=expression(sigma[R]/sigma[0]), 
                     y=expression(atop(paste("Proportion of total re-estimated sample size"),
                                       paste("such that ",xi[emp]," is at least 0.9")))) + 
                theme(text=element_text(size=20,family="Palatino"),legend.title=element_text(face="bold",size=16),
                      legend.text=element_text(size=16),legend.title.align=0)

A.CI.line.GS <- 
          ggplot(data=col.mat.sec,aes(y=Aperleft,x=sd.data,group=f.ci,colour=f.ci,shape=f.ci))+
          geom_line(data=col.mat.sec,size=1.2) +
          geom_point(size=5) + 
          facet_grid(~ int.size,labeller=label_parsed)+ 
          theme_bw() +
          labs(x=expression(sigma[R]/sigma[0]), y=expression(atop(paste("Proportion of total re-estimated sample size"),
                                                      paste("such that ",xi[emp]," is at least 0.9")))) + 
          scale_colour_discrete(name="PI Width",guide=guide_legend(keywidth=2.5))+
          scale_shape_discrete(name="PI Width",guide=guide_legend(keywidth=2.5))+
          scale_colour_grey(name="PI Width",start=0.2,end=0.8,guide=guide_legend(keywidth=2.5))+
          theme(text=element_text(size=20,family="Palatino"),legend.title=element_text(face="bold",size=16),
                legend.title.align=0)


if(save==T){
  ggsave(sd.line.GS,filename="plot_sdline_GS.pdf",width=19,height=10.2)
  ggsave(A.sd.line.GS,filename="plot_Asdline_GS.pdf",width=19,height=10.2)
  ggsave(CI.line.GS,filename="plot_ciline_GS.pdf",width=19,height=10.2)
  ggsave(A.CI.line.GS,filename="plot_Aciline_GS.pdf",width=19,height=10.2)}

out.list <- list(sd.line.GS,A.sd.line.GS,CI.line.GS,A.CI.line.GS)

out.list
}
