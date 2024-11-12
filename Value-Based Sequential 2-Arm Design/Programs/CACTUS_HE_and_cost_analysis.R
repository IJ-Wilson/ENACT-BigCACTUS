################################################################################
## ENACT Case Study Analysis                                                  ##
## Laura Flight                                                               ##
## 31Aug2022                                                                  ##
################################################################################

#####################
## Setup workspace ##
#####################

# Load required packages -------------------------------------------------------

rm(list=ls())
require(devtools)

install_version("heemod", version = "0.13.0", repos = "http://cran.us.r-project.org")

pkgs<-c("dplyr","ggplot2","reshape2","tidyr","tibble","xtable","heemod")
inst<-lapply(pkgs,library,character.only=TRUE)

# Set working directory --------------------------------------------------------

setwd("xxxxxx/Value-Based Sequential 2-Arm Design")

# Load functions (should be saved in the working directory) --------------------

source("Programs/Functions/Big_CACTUS_R_model.R")
source("Programs/Functions/BC_model_param_function.R")

################################################################################

###############
## Load data ##
###############

# Based on the imputed Big CACTUS data and organises participants in threes from each arm

DATA<-read.csv("Data/Raw/data_set.csv")

################################################################################

############################
## Primary analysis - ALL ##
############################

BC_HE_model_func(in.age=65,
                 in.prob.death.disease=0.024690, # 3-month probability
                 in.p_relapse_6_9_CSLT=unlist(param_func(DATA)$relapse_9_CSLT),
                 in.p_relapse_6_9_UC=unlist(param_func(DATA)$relapse_9_UC),
                 in.p_relapse_6_9_AC=unlist(param_func(DATA)$relapse_9_AC),
                 # Probability of relapse (9–12 months)
                 in.p_relapse_9_12_CSLT=unlist(param_func(DATA)$relapse_12_CSLT),
                 in.p_relapse_9_12_UC=unlist(param_func(DATA)$relapse_12_UC),
                 in.p_relapse_9_12_AC=unlist(param_func(DATA)$relapse_12_AC),
                 # Probability of relapse (12 months onwards)
                 in.p_relapse_12_CSLT=unlist(param_func(DATA)$relapse_12_CSLT),
                 in.p_relapse_12_UC=unlist(param_func(DATA)$relapse_12_UC),
                 in.p_relapse_12_AC= unlist(param_func(DATA)$relapse_12_AC),
                 # Probability of good response (0–6 months)
                 in.p_good_0_6_CSLT=1-exp(-(-(log(1-unlist(param_func(DATA)$response_6_CSLT))/6))*3),
                 in.p_good_0_6_UC=1-exp(-(-(log(1-unlist(param_func(DATA)$response_6_UC))/6))*3),
                 in.p_good_0_6_AC=1-exp(-(-(log(1-unlist(param_func(DATA)$response_6_AC))/6))*3),
                 # Probability of new good response (6–9 months)
                 in.p_good_6_9_CSLT=unlist(param_func(DATA)$response_9_CSLT),
                 in.p_good_6_9_UC=unlist(param_func(DATA)$response_9_UC),
                 in.p_good_6_9_AC=unlist(param_func(DATA)$response_9_AC),
                 # Probability of new good response or renewed response in people who responded at 6 months and relapsed at 9 months (9–12 months)
                 in.p_good_12_CSLT=unlist(param_func(DATA)$response_12_CSLT),
                 in.p_good_12_UC=unlist(param_func(DATA)$response_12_UC),
                 in.p_good_12_AC=unlist(param_func(DATA)$response_12_AC),
                 # Costs
                 in.cost_CSLT=unlist(param_func(DATA)$cost_CSLT),
                 in.cost_UC=0,
                 in.cost_AC=unlist(param_func(DATA)$cost_AC),
                 in.u_aphasia=param_func(DATA)$utility0,
                 in.u_response_6=param_func(DATA)$utility0+param_func(DATA)$u_response_6,
                 in.u_response_9=param_func(DATA)$utility0+param_func(DATA)$u_response_9,
                 in.u_response_12=param_func(DATA)$utility0+param_func(DATA)$u_response_12,
                 in.cycle.length=3,
                 in.threshold=20000 ,
                 in.DR=0.035)

################################################################################

################################################
## Fitting the model to the data sequentially ##
################################################

# After 6 participants randomised ----------------------------------------------

param<-list()
for(i in seq(6,270,3)){
  param[[i]]<- param_func(DATA[1:i,])
}

param2<-na.omit(param)
INMB_CSLT_UC<-rep(NA)
INMB_CSLT_AC<-rep(NA)
n_CSLT<-rep(NA)
n_UC<-rep(NA)
n_AC<-rep(NA)
QALY.diff_CSLT_UC<-rep(NA)
cost.diff_CSLT_UC<-rep(NA)
param<-matrix(NA,nrow=270,ncol=30)

for(i in seq(6,270,3)){
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
                      # Probability of new good response or renewed response in people who responded at 6 months and relapsed at 9 months (9–12 months)
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
                      in.threshold=20000,
                      in.DR=0.035)
  
  param[i,]<-(X$param)
  INMB_CSLT_UC[i]<-X$INB_CSLT_UC
  INMB_CSLT_AC[i]<-X$INB_CSLT_AC
  QALY.diff_CSLT_UC[i]<-X$QALY.diff_CSLT_UC
  cost.diff_CSLT_UC[i]<-X$cost.diff_CSLT_UC
  
  n_CSLT[i]<-(param2[[i]]$n_CSLT)
  n_UC[i]<-param2[[i]]$n_UC
  n_AC[i]<-param2[[i]]$n_AC
}

# Storing the data to create a plot of the INB observed during the trial -------

plot.dat<-data.frame(cbind(na.omit(unlist(n_CSLT)),
                           na.omit(unlist(n_AC)),
                           na.omit(unlist(n_UC)),
                           na.omit(param),
                           na.omit(INMB_CSLT_UC),
                           na.omit(INMB_CSLT_AC),
                           na.omit(QALY.diff_CSLT_UC),
                           na.omit(cost.diff_CSLT_UC),
                           seq(6,270,3),
                           as.Date(DATA$random_date[seq(6,270,3)]),
                           as.Date(DATA$final_date[seq(6,270,3)])))

colnames(plot.dat)<-c("n_CSLT",
                      "n_AC",
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
                      "in.DR","INMB_CSLT_UC","INMB_CSLT_AC","QALY_CSLT_UC","Cost_CSLT_UC","patients","random_date","final_date")

################################################################################

###################
## Cost analysis ##
###################

# Calculating costs during the trial -------------------------------------------

cost_GRANT_3arm<-read.csv("Data/Raw/BigCACTUS_3arm_costs.csv")%>%
  mutate(date=rep(seq(as.Date("2013-12-01"), by = "month", length.out = 54),3))

# Plot of budget and INMB ------------------------------------------------------

ggplot()+
  geom_line(data=data.frame(cost_GRANT_3arm)%>%
              filter(budget=="NIHR"),aes(x=date,y=(cumulative)))+
  geom_line(data=plot.dat,aes(x=as.Date(final_date,origin="1970-01-01"),y=(INMB_CSLT_UC*600)),color="blue",linetype="dashed")+
  geom_point(data=plot.dat%>%
               filter(patients==270),aes(x=as.Date(final_date,origin="1970-01-01"),y=(INMB_CSLT_UC*600)),color="blue")+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_text(data=plot.dat,label="A",aes(x=as.Date("2014-09-01"),y=1600000),colour=1)+
  geom_vline(xintercept=as.Date("2014-10-14"),linetype="dashed",colour="red")+ #recruitment starts
  # geom_vline(xintercept=as.Date("2015-10-14"),linetype="dashed")+ #first patient followed-up
  geom_text(data=plot.dat,label="B",aes(x=as.Date("2016-07-01"),y=1600000),colour=1)+
  geom_vline(xintercept=as.Date("2016-08-18"),linetype="dashed",colour="red")+ #recruitment ends
  geom_text(data=plot.dat,label="C",aes(x=as.Date("2017-07-01"),y=1600000),colour=1)+
  geom_vline(xintercept=as.Date("2017-08-18"),linetype="dashed")+ #last patient followed-up
  # geom_vline(xintercept=as.Date(plot.dat$final_date[1],origin="1970-01-01"),linetype="dashed",colour="blue")+ #6th person reaches 12 months and first interim
  # geom_vline(xintercept=as.Date("2014-01-01"),colour="black")+ # contract start
  geom_text(data=plot.dat,label="D",aes(x=as.Date("2018-04-15"),y=1600000),colour=1)+
  geom_vline(xintercept=as.Date("2018-06-01"),colour="black")+  # publication date
  theme_bw()+
  theme(axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
        axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position="top",
        axis.title.y = element_text(color = "black"),
        axis.title.y.right = element_text(color = "blue",angle=90))+
  scale_y_continuous(
    # Features of the first axis
    name = "Projected Cumulative Spend (3 arms) £",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./600, name="Cumulative estimate of E[INMB] at 1 year (£)"))+
  xlab("Date")

# Code to create output for MATLAB version of Figure 1.1 in original report ----

# https://eprints.whiterose.ac.uk/180084/7/enactDeliverable3.pdf

fig1.1_data<-full_join(plot.dat%>%
                         mutate(date=as.Date(final_date, origin="1970-01-01"))%>%
                         select(INMB_CSLT_UC,date),
                       data.frame(cost_GRANT_3arm)%>%
                         filter(budget=="NIHR")%>%
                         select(date,cumulative))

write.csv(fig1.1_data%>%
            mutate(year = as.numeric(format(date, format = "%Y")))%>%
            mutate(month = as.numeric(format(date, format = "%m")))%>%
            mutate(day = as.numeric(format(date, format = "%d"))), "Data/Derived/CACTUS_spend_data_matlab.csv")
