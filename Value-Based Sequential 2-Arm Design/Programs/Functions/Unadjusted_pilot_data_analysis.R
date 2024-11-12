################################################################################
## Function that analyses the CACTUS pilot data and calculates parameters     ##
## required for health economic model analysis                                ##
################################################################################

EVSI_analysis_pilot<-function(terminal_analysis, # 
                              terminal_design, # Design of trial
                              terminal_analysis_number, # Sample size at analysis
                              terminal_n_intervention, # Number of participants intervention arm
                              terminal_n_control, # Number of participants control arm
                              terminal_dataset,
                              alpha, # type I error rate 
                              beta, # type II error rate
                              rho_time) # Correlation between time points
{

  level<-c(seq(0.02,0.999,0.005))
  perc<-c(level/2,0.5,rev(1-level/2))
  
  ###################
  # PRIMARY OUTCOME #
  ###################
  
  primary_unadjusted<-terminal_dataset%>%
                        filter(trt==1)%>%
                            summarise(mean=mean(TE.gain.6,na.rm=T))-
    
                                    terminal_dataset%>%
                                        filter(trt==0)%>%
                                          summarise(mean=mean(TE.gain.6,na.rm=T))
  
  primary_sd_pool<-sqrt((terminal_dataset%>%
                          filter(trt==1)%>%
                          summarise(sd_sq=sd(TE.gain.6,na.rm=T)^2)+
                          
                          terminal_dataset%>%
                          filter(trt==0)%>%
                          summarise(sd_sq=sd(TE.gain.6,na.rm=T)^2))/2)
    
  primary_se_pool<-sqrt((1/terminal_n_intervention)+(1/terminal_n_control))*primary_sd_pool
  
  primary_lCI_unadjusted<-primary_unadjusted-qnorm(1-alpha/2)*primary_se_pool
  primary_uCI_unadjusted<-primary_unadjusted+qnorm(1-alpha/2)*primary_se_pool

  #########
  # COSTS #
  #########
  
  cost_unadjusted<-terminal_dataset%>%
    filter(trt==1)%>%
    summarise(mean=mean(Total.C,na.rm=T))-
    
    terminal_dataset%>%
    filter(trt==0)%>%
    summarise(mean=mean(Total.C,na.rm=T))
  
  cost_sd_pool<-sqrt((terminal_dataset%>%
                          filter(trt==1)%>%
                          summarise(sd_sq=sd(Total.C,na.rm=T)^2)+
                          
                          terminal_dataset%>%
                          filter(trt==0)%>%
                          summarise(sd_sq=sd(Total.C,na.rm=T)^2))/2)
  
  cost_se_pool<-sqrt((1/terminal_n_intervention)+(1/terminal_n_control))*cost_sd_pool
  
  cost_lCI_unadjusted<-cost_unadjusted-qnorm(1-alpha/2)*cost_se_pool
  cost_uCI_unadjusted<-cost_unadjusted+qnorm(1-alpha/2)*cost_se_pool
  
  
  
  cost_intervention_unadjusted<-terminal_dataset%>%
                    filter(trt==1)%>%
                    summarise(mean=mean(Total.C,na.rm=T))
  
  cost_control_unadjusted<-terminal_dataset%>%
                      filter(trt==0)%>%
                      summarise(mean=mean(Total.C,na.rm=T))
  
  ########
  # QALY #
  ########
  
  qaly_unadjusted<-terminal_dataset%>%
    filter(trt==1)%>%
    summarise(mean=mean(QALY.6,na.rm=T))-
    
    terminal_dataset%>%
    filter(trt==0)%>%
    summarise(mean=mean(QALY.6,na.rm=T))
  
  qaly_sd_pool<-sqrt((terminal_dataset%>%
                       filter(trt==1)%>%
                       summarise(sd_sq=sd(QALY.6,na.rm=T)^2)+
                       
                       terminal_dataset%>%
                       filter(trt==0)%>%
                       summarise(sd_sq=sd(QALY.6,na.rm=T)^2))/2)
  
  qaly_se_pool<-sqrt((1/terminal_n_intervention)+(1/terminal_n_control))*qaly_sd_pool
  
  qaly_lCI_unadjusted<-qaly_unadjusted-qnorm(1-alpha/2)*qaly_se_pool
  qaly_uCI_unadjusted<-qaly_unadjusted+qnorm(1-alpha/2)*qaly_se_pool
  
  
  
  qaly_intervention_unadjusted<-terminal_dataset%>%
    filter(trt==1)%>%
    summarise(mean=mean(QALY.6,na.rm=T))
  
  qaly_control_unadjusted<-terminal_dataset%>%
    filter(trt==0)%>%
    summarise(mean=mean(QALY.6,na.rm=T))
  
  #######
  # INB #
  #######
  
  inb_unadjusted<-terminal_dataset%>%
    filter(trt==1)%>%
    summarise(mean=mean(NB.6,na.rm=T))-
    
    terminal_dataset%>%
    filter(trt==0)%>%
    summarise(mean=mean(NB.6,na.rm=T))
  
  inb_sd_pool<-sqrt((terminal_dataset%>%
                       filter(trt==1)%>%
                       summarise(sd_sq=sd(NB.6,na.rm=T)^2)+
                       
                       terminal_dataset%>%
                       filter(trt==0)%>%
                       summarise(sd_sq=sd(NB.6,na.rm=T)^2))/2)
  
  inb_se_pool<-sqrt((1/terminal_n_intervention)+(1/terminal_n_control))*inb_sd_pool
  
  inb_lCI_unadjusted<-inb_unadjusted-qnorm(1-alpha/2)*inb_se_pool
  inb_uCI_unadjusted<-inb_unadjusted+qnorm(1-alpha/2)*inb_se_pool

  ####################
  # Model parameters #
  ####################
  
  ## Baseline utility control
  baseline_utility_unadjusted<-terminal_dataset%>%
                filter(trt==1)%>%
                summarise(mean=mean(baseline.U,na.rm=T))-
    
                terminal_dataset%>%
                filter(trt==0)%>%
                summarise(mean=mean(baseline.U,na.rm=T))

  
  baseline_utility_control_unadjusted<-terminal_dataset%>%
                                              filter(trt==0)%>%
                                              summarise(mean=mean(baseline.U,na.rm=T))
  
  baseline_utility_control_sd<-terminal_dataset%>%
                                    filter(trt==0)%>%
                                    summarise(sd=sd(baseline.U,na.rm=T))

  ## Utility of responders
  terminal_dataset$response<-matrix(rep(NA),length(terminal_dataset$TE.gain.6))
  terminal_dataset$response<-ifelse(terminal_dataset$TE.gain.6>=0.17,1,0)
  
  #step where calculate utility of responders?
  utility_responders<-terminal_dataset%>%
                            filter(response==1)%>%
                            summarise(mean=mean(baseline.U,na.rm=T))

  utility_non_responders<-terminal_dataset%>%
                             filter(response==0)%>%
                            summarise(mean=mean(baseline.U,na.rm=T))
  
  utility_increment_unadjusted<-terminal_dataset%>%
                                    filter(response==1)%>%
                                    summarise(mean=mean(E.gain.6,na.rm=T))-
                                            terminal_dataset%>%
                                            filter(response==0)%>%
                                            summarise(mean=mean(E.gain.6,na.rm=T))
    
  ## Probability of a good response at 6 months
  p_good_6_intervention_unadjusted<-  terminal_dataset%>%
                                          filter(trt==1&TE.gain.6>=0.17)%>%
                                            nrow()/terminal_dataset%>%
                                                      filter(trt==1&!is.na(TE.6))%>%
                                                            nrow()
  
  ## Probability of a good response at 9 months
    primary_9_unadjusted<-terminal_dataset%>%
    filter(trt==1)%>%
    summarise(mean=mean(TE.gain.9,na.rm=T))-
    
    terminal_dataset%>%
    filter(trt==0)%>%
    summarise(mean=mean(TE.gain.9,na.rm=T))

  p_good_9_intervention_unadjusted<- terminal_dataset%>%
                                              filter(trt==1&TE.gain.9>=0.17)%>%
                                              nrow()/terminal_dataset%>%
                                              filter(trt==1&!is.na(TE.9))%>%
                                              nrow()

  ## Relapse rate
  relapse_rate_unadjusted<-ifelse(p_good_6_intervention_unadjusted-p_good_9_intervention_unadjusted<0,0,
                                  p_good_6_intervention_unadjusted-p_good_9_intervention_unadjusted)
  month_relapse_rate_unadjusted<--(1/4)*log(1-relapse_rate_unadjusted)
  p_relapse_month_unadjusted<-1-exp(-month_relapse_rate_unadjusted)
 
  ## Resource Costs
  resource_cost_unadjusted<-terminal_dataset%>%
    filter(trt==1)%>%
    summarise(mean=mean(Resource.C,na.rm=T))-
    
    terminal_dataset%>%
    filter(trt==0)%>%
    summarise(mean=mean(Resource.C,na.rm=T))
  
  resource_cost_intervention_unadjusted<-terminal_dataset%>%
    filter(trt==1)%>%
    summarise(mean=mean(Resource.C,na.rm=T))
  
  resource_cost_control_unadjusted<-terminal_dataset%>%
    filter(trt==0)%>%
    summarise(mean=mean(Resource.C,na.rm=T))
  
  resource_cost_intervention_sd_unadjusted<-terminal_dataset%>%
    filter(trt==1)%>%
    summarise(sd=sd(Resource.C,na.rm=T))
  
  resource_cost_control_sd_unadjusted<-terminal_dataset%>%
    filter(trt==0)%>%
    summarise(sd=sd(Resource.C,na.rm=T))
  
  ###########
  # OUTPUTS #
  ###########
  
  unadjusted_list<-list(resource_cost_intervention_unadjusted=resource_cost_intervention_unadjusted,
                        resource_cost_control_unadjusted=resource_cost_control_unadjusted,
                        utility_increment_unadjusted=utility_increment_unadjusted,
                        p_good_6_intervention_unadjusted=p_good_6_intervention_unadjusted,
                        p_relapse_month_unadjusted=p_relapse_month_unadjusted,
                        baseline_utility_control_unadjusted=baseline_utility_control_unadjusted)
  
  data_param<- list(alpha=alpha,
                    beta=beta,
                    baseline_primary_control=terminal_dataset%>%
                                                  filter(trt==0)%>%
                                                  summarise(mean=mean(baseline.TE,na.rm=T)),
                    month_6_primary_control=terminal_dataset%>%
                                                  filter(trt==0)%>%
                                                  summarise(mean=mean(TE.6,na.rm=T)),
                    month_9_primary_control=terminal_dataset%>%
                                                  filter(trt==0)%>%
                                                  summarise(mean=mean(TE.9,na.rm=T)),
                    
                    primary_intervention_sd_pool=sqrt((terminal_dataset%>%
                                         filter(trt==1)%>%
                                         summarise(sd_sq=sd(baseline.TE,na.rm=T)^2)+
                                         
                                         terminal_dataset%>%
                                         filter(trt==1)%>%
                                         summarise(sd_sq=sd(TE.6,na.rm=T)^2)+
                                         
                                         terminal_dataset%>%
                                         filter(trt==1)%>%
                                         summarise(sd_sq=sd(TE.9,na.rm=T)^2))/3),
                      
                    primary_control_sd_pool=sqrt((terminal_dataset%>%
                                          filter(trt==0)%>%
                                          summarise(sd_sq=sd(baseline.TE,na.rm=T)^2)+
                                          
                                          terminal_dataset%>%
                                          filter(trt==0)%>%
                                          summarise(sd_sq=sd(TE.6,na.rm=T)^2)+
                                          
                                          terminal_dataset%>%
                                          filter(trt==0)%>%
                                          summarise(sd_sq=sd(TE.9,na.rm=T)^2))/3),
                    
                    baseline_primary_intervention=terminal_dataset%>%
                                          filter(trt==1)%>%
                                          summarise(mean=mean(baseline.TE,na.rm=T)),
                    month_6_primary_intervention=terminal_dataset%>%
                                          filter(trt==1)%>%
                                          summarise(mean=mean(TE.6,na.rm=T)),
                    month_9_primary_intervention=terminal_dataset%>%
                                          filter(trt==1)%>%
                                          summarise(mean=mean(TE.9,na.rm=T)),
              
                    baseline_utility_control_unadjusted=baseline_utility_control_unadjusted,
                    month_6_utility_control=terminal_dataset%>%
                                          filter(trt==0)%>%
                                          summarise(mean=mean(U.6,na.rm=T)),
                    month_9_utility_control=terminal_dataset%>%
                                          filter(trt==0)%>%
                                          summarise(mean=mean(U.9,na.rm=T)),
                    
                    utility_intervention_sd_pool=sqrt((terminal_dataset%>%
                                          filter(trt==1)%>%
                                          summarise(sd_sq=sd(baseline.U,na.rm=T)^2)+
                                          
                                          terminal_dataset%>%
                                          filter(trt==1)%>%
                                          summarise(sd_sq=sd(U.6,na.rm=T)^2)+
                                          
                                          terminal_dataset%>%
                                          filter(trt==1)%>%
                                          summarise(sd_sq=sd(U.9,na.rm=T)^2))/3),
                    
                    utility_control_sd_pool=sqrt((terminal_dataset%>%
                                           filter(trt==0)%>%
                                           summarise(sd_sq=sd(baseline.U,na.rm=T)^2)+
                                           
                                           terminal_dataset%>%
                                           filter(trt==0)%>%
                                           summarise(sd_sq=sd(U.6,na.rm=T)^2)+
                                           
                                           terminal_dataset%>%
                                           filter(trt==0)%>%
                                           summarise(sd_sq=sd(U.9,na.rm=T)^2))/3),
                    
                    baseline_utility_intervention=terminal_dataset%>%
                                          filter(trt==1)%>%
                                          summarise(mean=mean(baseline.U,na.rm=T)),
                    month_6_utility_intervention=terminal_dataset%>%
                                          filter(trt==1)%>%
                                          summarise(mean=mean(U.6,na.rm=T)),
                    month_9_utility_intervention=terminal_dataset%>%
                                          filter(trt==1)%>%
                                          summarise(mean=mean(U.9,na.rm=T)),
                    
                    resource_cost_control_unadjusted=resource_cost_control_unadjusted,
                    resource_cost_control_sd_unadjusted=resource_cost_control_sd_unadjusted,
                    resource_cost_intervention_sd_unadjusted=resource_cost_intervention_sd_unadjusted,
                    resource_cost_intervention_unadjusted=resource_cost_intervention_unadjusted,
                    
                    rho_time=rho_time # Correlation between time points
             
  )
  
  out=list(unadjusted_list=unadjusted_list,data_param=data_param)
  
}
