boot_analysis_func<-function(
                            tvec,
                            bndupper,
                            bndlower,
                            priormean,
                            n0,
                            Tmax,
                            P,
                            plot.dat.mat.frame,
                            B,
                            pre.rec.cost,
                            fixed.recruit.follow.cost,
                            pp.recruit.cost,
                            pp.followup.cost,
                            trial.end.cost,
                            block)
{tau<-ceiling(95*(Tmax/95)/((638.751*(Tmax/95))/365.25)) # delay
colnames(tvec)<-"tvec"
colnames(bndupper)<-"bndupper"
colnames(bndlower)<-"bndlower"
boundary.dat<-cbind(tvec,bndupper,bndlower)

boundary.dat95<-cbind(t(read.csv("Data/Derived/Tmax 95/tvec.csv",header=F)),
                      t(read.csv("Data/Derived/Tmax 95/bndupper.csv",header=F)),
                      t(read.csv("Data/Derived/Tmax 95/bndlower.csv",header=F)))

colnames(boundary.dat95)  <-c("tvec","bndupper","bndlower")

boundary.dat435<-cbind(t(read.csv("Data/Derived/Tmax 435/tvec.csv",header=F)),
                       t(read.csv("Data/Derived/Tmax 435/bndupper.csv",header=F)),
                       t(read.csv("Data/Derived/Tmax 435/bndlower.csv",header=F)))

colnames(boundary.dat435)  <-c("tvec","bndupper","bndlower")

# Load functions (should be saved in the working directory) --------------------

source("Programs/Functions/Unadjusted_pilot_data_analysis.R")
source("Programs/Functions/CACTUS_HE_model.R")

# Load data --------------------------------------------------------------------

paired.data<-read.csv("Data/Raw/data_set.csv")

  param<-list()
for(i in seq(block*3,270,block*3)){
  param[[i]]<- param_func(paired.data[1:i,])
}


param2<-param[-which(sapply(param, is.null))]
INMB_CSLT_UC<-rep(NA)
INMB_CSLT_AC<-rep(NA)
n_CSLT<-rep(NA)
n_UC<-rep(NA)
n_AC<-rep(NA)
QALY.diff_CSLT_UC<-rep(NA)
cost.diff_CSLT_UC<-rep(NA)
param<-matrix(NA,nrow=length(seq(block*3,270,block*3)),ncol=30)

for(i in 1:length(seq(block*3,270,block*3))){
  library(dplyr)
  X<-BC_HE_model_func(in.age=65,
                   in.prob.death.disease=0.024690, # 3-month probability
                   in.p_relapse_6_9_CSLT=unlist(param2[[i]]$relapse_9_CSLT),
                   in.p_relapse_6_9_UC=unlist(param2[[i]]$relapse_9_UC),
                   in.p_relapse_6_9_AC=unlist(param2[[i]]$relapse_9_AC),
                   # Probability of relapse (9–12 months)
                   in.p_relapse_9_12_CSLT=unlist(param2[[i]]$relapse_12_CSLT),
                   in.p_relapse_9_12_UC=unlist(param2[[i]]$relapse_12_UC),
                   in.p_relapse_9_12_AC=unlist(param2[[i]]$relapse_12_AC),
                   # Probability of relapse (12 months onwards)
                   in.p_relapse_12_CSLT=unlist(param2[[i]]$relapse_12_CSLT),
                   in.p_relapse_12_UC=unlist(param2[[i]]$relapse_12_UC),
                   in.p_relapse_12_AC= unlist(param2[[i]]$relapse_12_AC),
                   # Probability of good response (0–6 months)
                   in.p_good_0_6_CSLT=1-exp(-(-(log(1-unlist(param2[[i]]$response_6_CSLT))/6))*3),
                   in.p_good_0_6_UC=1-exp(-(-(log(1-unlist(param2[[i]]$response_6_UC))/6))*3),
                   in.p_good_0_6_AC=1-exp(-(-(log(1-unlist(param2[[i]]$response_6_AC))/6))*3),
                   # Probability of new good response (6–9 months)
                   in.p_good_6_9_CSLT=unlist(param2[[i]]$response_9_CSLT),
                   in.p_good_6_9_UC=unlist(param2[[i]]$response_9_UC),
                   in.p_good_6_9_AC=unlist(param2[[i]]$response_9_AC),
                   # Probability of new good responseor renewed response in people who responded at 6 months and relapsed at 9 months (9–12 months)
                   in.p_good_12_CSLT=unlist(param2[[i]]$response_12_CSLT),
                   in.p_good_12_UC=unlist(param2[[i]]$response_12_UC),
                   in.p_good_12_AC=unlist(param2[[i]]$response_12_AC),
                   # Costs
                   in.cost_CSLT=unlist(param2[[i]]$cost_CSLT),
                   in.cost_UC=0,
                   in.cost_AC=unlist(param2[[i]]$cost_AC),
                   in.u_aphasia=param2[[i]]$utility0,
                   in.u_response_6=param2[[i]]$utility0+param2[[i]]$u_response_6,
                   in.u_response_9=param2[[i]]$utility0+param2[[i]]$u_response_9,
                   in.u_response_12=param2[[i]]$utility0+param2[[i]]$u_response_12,
                   in.cycle.length=3,
                   in.threshold=20000 ,
                   in.DR=0.035)
  
  param[i,]<-(X$param)
  INMB_CSLT_UC[i]<-X$INB_CSLT_UC
  # INMB_CSLT_AC[i]<-X$INB_CSLT_AC
  QALY.diff_CSLT_UC[i]<-X$QALY.diff_CSLT_UC
  cost.diff_CSLT_UC[i]<-X$cost.diff_CSLT_UC
  
  n_CSLT[i]<-(param2[[i]]$n_CSLT)
  n_UC[i]<-param2[[i]]$n_UC
  # n_AC[i]<-param2[[i]]$n_AC
}

# Storing the data to create a plot of the INB observed during the trial -------

plot.data<-data.frame(cbind(na.omit(unlist(n_CSLT)),
                           # na.omit(unlist(n_AC)),
                           na.omit(unlist(n_UC)),
                           na.omit(param),
                           na.omit(INMB_CSLT_UC),
                           # na.omit(INMB_CSLT_AC),
                           na.omit(QALY.diff_CSLT_UC),
                           na.omit(cost.diff_CSLT_UC),
                           seq(block*2,270/3*2,block*2)))
colnames(plot.data)<-c("n_CSLT",
                      # "n_AC",
                      "n_UC",
                      "in.age",
                      "in.prob.death.disease", # 3-month probability
                      "in.p_relapse_6_9_CSLT",
                      "in.p_relapse_6_9_UC",
                      "in.p_relapse_6_9_AC",
                      # Probability of relapse (9–12 months)
                      "in.p_relapse_9_12_CSLT",
                      "in.p_relapse_9_12_UC",
                      "in.p_relapse_9_12_AC" ,
                      # Probability of relapse (12 months onwards)
                      "in.p_relapse_12_CSLT",
                      "in.p_relapse_12_UC",
                      "in.p_relapse_12_AC",
                      # Probability of good response (0–6 months)
                      "in.p_good_0_6_CSLT",
                      "in.p_good_0_6_UC",
                      "in.p_good_0_6_AC",
                      # Probability of new good response (6–9 months)
                      "in.p_good_6_9_CSLT",
                      "in.p_good_6_9_UC",
                      "in.p_good_6_9_AC",
                      # Probability of new good response,or renewed response in people who responded at 6 months and relapsed at 9 months (9–12 months)
                      "in.p_good_12_CSLT",
                      "in.p_good_12_UC",
                      "in.p_good_12_AC",
                      # Costs
                      "in.cost_CSLT",
                      "in.cost_UC",
                      "in.cost_AC",
                      "in.u_aphasia",
                      "in.u_response_6",
                      "in.u_response_9",
                      "in.u_response_12",
                      "in.cycle.length",
                      "in.threshold",
                      "in.DR",
                      "INMB_CSLT_UC",
                      # "INMB_CSLT_AC",
                      "QALY_CSLT_UC",
                      "Cost_CSLT_UC",
                      "patients")

# Calculating the prior/posterior INB ------------------------------------------

mu<-c(plot.data$INMB_CSLT_UC)
mu_n<-rep(NA,length( seq(block,270/3,block)+1)) # +1 because there is the prior and then the trial allocations
mu_n[1]<-3190.78
m<-c(7,7+seq(block,270/3,block))
t<-c(seq(0,max( seq(block,270/3,block)),block))

for(i in 2:(length( seq(block,270/3,block))+1)){
  mu_n[i]<-((m[1]*mu_n[1])+((t[i])*mu[(i-1)]))/(m[1]+(t[i]))
}

original.post.INMB<-data.frame(mu_n,m,m+round(tau),0) 

colnames(original.post.INMB)<-c("INMB","m","m_delay","label")

# Collecting information about each bootstrap ----------------------------------

original.stop<-rep(NA)
original.cross.inmb<-rep(NA)
original.pipe.inmb<-rep(NA)
original.randomised.N<-rep(NA)
original.trt.select<-rep(NA)

original.MU<-(original.post.INMB %>% select(INMB)) # REMOVED as.vector

for(i in 1:nrow(tvec)){
  # +1 because the MU vector include the prior mean that relates to 0 patients randomised in the trial
  original.stop[i]<-original.MU[(ceiling(tvec[i])-n0-round(tau))+1,1]>bndupper[i]|original.MU[(ceiling(tvec[i])-n0-round(tau))+1,1]<bndlower[i]
  
  # How many randomised when stopping? 
  original.randomised.N[i]<-ifelse(original.stop[i]==TRUE,min(nrow(original.MU)-1,ceiling(tvec[i])-n0),NA)
  # need to subtract 1 from length of MU because this includes the prior mean 
  
  # If yes save the INMB at crossing
  original.cross.inmb[i]<-ifelse(original.stop[i]==TRUE,original.MU[(ceiling(tvec[i]))-n0-round(tau)+1,1],NA)
  # If yes save the INMB at pipeline followed up
  original.pipe.inmb[i]<-ifelse(original.stop[i]==TRUE,original.MU[min(nrow(original.MU),(ceiling(tvec[i]))-n0+1),1],NA)
  
  original.trt.select[i]<-ifelse(original.pipe.inmb[i]>=0,"CSLT","UC") # make decision based on the pipeline INMB
}

original.res<-as.matrix(c(first(na.omit(original.randomised.N)),
                          first(na.omit(original.cross.inmb)),
                          first(na.omit(original.pipe.inmb)),
                          first(na.omit(original.trt.select))
), ncol=4)
rownames(original.res)<-c("randomised.N","cross.inmb","pipe.inmb","trt.select")

### NOTE - this code only works correctly when the block size is a multiple of tau. Checks and changes needed when this isn't the case as
## the pipe line INMB will be between two points in the MU vector. Can use interpolation or rounding up? At the moment the code likely to round up
## as using ceiling but this takes the pipe INMB from more participants 

stop<-rep(NA)
cross.inmb<-rep(NA)
pipe.inmb<-rep(NA)
randomised.N<-rep(NA)
trt.select<-rep(NA)
res<-matrix(rep(NA),ncol=B,nrow=5)

for(j in 1:B){
  MU<-(plot.dat.mat.frame%>%
                  select(label2:label)%>%
                  drop_na()%>%
                  filter(label==j)%>%
                  select(post_INMB)) # REMOVED as.vector
  
  for(i in 1:nrow(tvec)){
    # +1 because the MU vector include the prior mean that relates to 0 patients randomised in the trial
    
    # Will the trial stop early? 
    stop[i]<-MU[(ceiling((tvec[i]-n0-round(tau))/block)+1),1]>bndupper[i]|MU[(ceiling((tvec[i]-n0-round(tau))/block)+1),1]<bndlower[i]
    
    # How many randomised when stopping? 
    randomised.N[i]<-ifelse(stop[i]==TRUE,min((nrow(MU)-1)*block,min(ceiling((tvec[i]-n0)/block)*block,Tmax)),NA) #(nrow(MU)-1)*block is the maximum sample size of the trial
    # need to subtract 1 from length of MU because this includes the prior mean 
    
    # If yes save the INMB at crossing
    cross.inmb[i]<-ifelse(stop[i]==TRUE,MU[(ceiling((tvec[i]-n0-round(tau))/block)+1),1],NA)
    # If yes save the INMB at pipeline followed up
    pipe.inmb[i]<-ifelse(stop[i]==TRUE,MU[(min((nrow(MU))*block,(ceiling((tvec[i]-n0)/block)+1)*block))/block,1],NA)
  }
  res[,j]<-as.matrix(c(first(na.omit(randomised.N)),
                       first(na.omit(cross.inmb)),
                       first(na.omit(pipe.inmb)),
                       first(na.omit(trt.select)),
                       j),ncol=5)
  rownames(res)<-c("randomised.N","cross.inmb","pipe.inmb","trt.select","label")
}


# Value based design results ---------------------------------------------------

formatted.res<-data.frame(t(res))%>%
  mutate(randomised.N=as.numeric(as.character(randomised.N)))%>%
  mutate(cross.inmb=as.numeric(as.character(cross.inmb)))%>%
  mutate(pipe.inmb=as.numeric(as.character(pipe.inmb)))%>%
  mutate(new.cost.sampling=pre.rec.cost+trial.end.cost+randomised.N*pp.recruit.cost*2+randomised.N*pp.followup.cost*2+fixed.recruit.follow.cost)%>%
  mutate(planned.cost.sampling=pre.rec.cost+trial.end.cost+95*pp.recruit.cost*2+95*pp.followup.cost*2+fixed.recruit.follow.cost)%>%
  mutate(change.budget.planned=new.cost.sampling-planned.cost.sampling)%>%
  mutate(P.INMB=P*pipe.inmb)%>%
  mutate(adoption.decision=ifelse(P.INMB-0>0,"CSLT","UC"))%>%
  mutate(exp.value=ifelse(adoption.decision=="UC",-new.cost.sampling-0,P.INMB-new.cost.sampling-0))

## Proportion of trials that were shorter than 95 ------------------------------

formatted.res%>%
  filter(randomised.N<95)%>%
  nrow()/formatted.res%>%
  nrow()

formatted.res%>%
  filter(randomised.N==95)%>%
  nrow()/formatted.res%>%
  nrow()

formatted.res%>%
  filter(randomised.N>95)%>%
  nrow()/formatted.res%>%
  nrow()

length.tab_95<-matrix(c(formatted.res%>%
                       filter(randomised.N<95)%>%
                       nrow()/formatted.res%>%
                       nrow(),
                     
                     formatted.res%>%
                       filter(randomised.N==95)%>%
                       nrow()/formatted.res%>%
                       nrow(),
                     
                     formatted.res%>%
                       filter(randomised.N>95)%>%
                       nrow()/formatted.res%>%
                       nrow()),ncol=3)
colnames(length.tab_95)<-c("shorter","same","longer")

## Proportion of trials that were shorter than 435 -----------------------------

formatted.res%>%
  filter(randomised.N<435)%>%
  nrow()/formatted.res%>%
  nrow()

formatted.res%>%
  filter(randomised.N==435)%>%
  nrow()/formatted.res%>%
  nrow()

formatted.res%>%
  filter(randomised.N>435)%>%
  nrow()/formatted.res%>%
  nrow()

length.tab_435<-matrix(c(formatted.res%>%
                           filter(randomised.N<435)%>%
                           nrow()/formatted.res%>%
                           nrow(),
                         
                         formatted.res%>%
                           filter(randomised.N==435)%>%
                           nrow()/formatted.res%>%
                           nrow(),
                         
                         formatted.res%>%
                           filter(randomised.N>435)%>%
                           nrow()/formatted.res%>%
                           nrow()),ncol=3)
colnames(length.tab_435)<-c("shorter","same","longer")

value.adaptive.tab<-formatted.res%>%
  select(-adoption.decision,-trt.select)%>%
  summarise_all(list(mean = mean, 
                      sd = sd,
                     min = min, 
                      max = max
                      ),na.rm=T)%>% 
  gather(stat, val) %>%
  separate(stat, into = c("var", "stat"), sep = "_") %>%
  spread(stat, val) %>%
  select(var, mean, sd,min, max)

## Summary of the budget change from planned -----------------------------------

value.adaptive.budget<-formatted.res%>%
  summarise(mean.budget.change.planned=mean(change.budget.planned),
            sd.budget.change.planned=sd(change.budget.planned),
            min.budget.change.planned=min(change.budget.planned),
            max.budget.change.planned=max(change.budget.planned),
            perc.total.fixed=mean(change.budget.planned)/mean(planned.cost.sampling)*100)

## Which treatment was favoured ------------------------------------------------

value.adaptive.trt<-formatted.res%>%
  mutate(lower.crossed.CSLT=ifelse(cross.inmb<=0&pipe.inmb>0,1,0))%>%
  mutate(lower.crossed.UC=ifelse(cross.inmb<=0&pipe.inmb<=0,1,0))%>%
  mutate(upper.crossed.CSLT=ifelse(cross.inmb>0&pipe.inmb>0,1,0))%>%
  mutate(upper.crossed.UC=ifelse(cross.inmb>0&pipe.inmb<=0,1,0))%>%
  summarise(perc.lower.CSLT=mean(lower.crossed.CSLT),
            perc.lower.UC=mean(lower.crossed.UC),
            perc.upper.CSLT=mean(upper.crossed.CSLT),
            perc.upper.UC=mean(upper.crossed.UC))


# Fixed design results ---------------------------------------------------------

## Mean INMB -------------------------------------------------------------------

fixed.tab<-plot.dat.mat.frame%>%
  select(label:fixed_m_delay)%>%
  drop_na()%>%
  filter(fixed_m==Tmax+7)%>%
  mutate(fixed.INMB=fixed_INMB)%>%
  mutate(fixed.m=fixed_m-n0)%>%
  mutate(fixed.m.delay=fixed_m_delay)%>%
  select(label,fixed.INMB:fixed.m.delay)%>%
  mutate(cost.sampling=pre.rec.cost+trial.end.cost+Tmax*pp.recruit.cost*2+Tmax*pp.followup.cost*2+fixed.recruit.follow.cost)%>%
  mutate(planned.cost.sampling=pre.rec.cost+trial.end.cost+95*pp.recruit.cost*2+95*pp.followup.cost*2+fixed.recruit.follow.cost)%>%
  mutate(change.budget.planned=cost.sampling-planned.cost.sampling)%>%
  summarise_all(list(mean = mean, 
                     sd = sd,
                     min = min, 
                     max = max))%>% 
  gather(stat, val) %>%
  separate(stat, into = c("var", "stat"), sep = "_") %>%
  spread(stat, val) %>%
  select(var, mean, sd,min, max)

## Budget info -----------------------------------------------------------------

fixed.budget<-plot.dat.mat.frame%>%
  select(label:fixed_m_delay)%>%
  drop_na()%>%
  filter(fixed_m==Tmax+7)%>%
  mutate(fixed.INMB=fixed_INMB)%>%
  mutate(fixed.m=fixed_m-n0)%>%
  mutate(fixed.m.delay=fixed_m_delay)%>%
  select(label,fixed.INMB:fixed.m.delay)%>%
  mutate(cost.sampling=pre.rec.cost+trial.end.cost+Tmax*pp.recruit.cost*2+Tmax*pp.followup.cost*2+fixed.recruit.follow.cost)%>%
  mutate(planned.cost.sampling=pre.rec.cost+trial.end.cost+95*pp.recruit.cost*2+95*pp.followup.cost*2+fixed.recruit.follow.cost)%>%
  mutate(change.budget.planned=cost.sampling-planned.cost.sampling)%>%
  summarise(mean.change.budget.fixed=mean(change.budget.planned),
            sd.change.budget.fixed=sd(change.budget.planned),
            min.change.budgetb.fixed=min(change.budget.planned),
            max.change.budget.fixed=max(change.budget.planned),
            perc.total.fixed=mean(change.budget.planned)/mean(planned.cost.sampling)*100)

## Which treatment was favoured -------------------------------------------------

fixed.trt<-plot.dat.mat.frame%>%
  select(label:fixed_m_delay)%>%
  drop_na()%>%
  filter(fixed_m==Tmax+7)%>%
  mutate(fixed.INMB=fixed_INMB)%>%
  mutate(fixed.m=fixed_m-n0)%>%
  mutate(fixed.m.delay=fixed_m_delay)%>%
  select(label,fixed.INMB:fixed.m.delay)%>%
  mutate(CSLT.picked=ifelse(fixed.INMB>0,1,0))%>%
  mutate(UC.picked=ifelse(fixed.INMB<0,1,0))%>%
  summarise(perc.CSLT=mean(CSLT.picked),
            perc.UC=mean(UC.picked))

## Looking to see whether fixed and adaptive designs make the same decision ----

adaptive.fixed.agree<-plot.dat.mat.frame%>%
  select(label, fixed_m,fixed_INMB)%>%
  drop_na()%>%
  filter(fixed_m==Tmax+7)%>%
  full_join(formatted.res%>%
              select(label,pipe.inmb))%>%
  mutate(agree=ifelse(fixed_INMB>0 &pipe.inmb>0|fixed_INMB<=0&pipe.inmb<=0,1,0))%>%
  summarise(agree.sum=sum(agree)/B)

# Print results ----------------------------------------------------------------

out=list(value.adaptive.tab,
         value.adaptive.budget,
         value.adaptive.trt,
         fixed.tab,
         fixed.budget,
         fixed.trt,
         adaptive.fixed.agree,
         length.tab_95,
         length.tab_435)
print(out)
}
