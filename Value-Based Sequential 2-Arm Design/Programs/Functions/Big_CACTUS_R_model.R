################################################################################
## Health economic model developed for the analysis of the Big CACTUS trial   ##
## data                                                                       ##
## Laura Flight                                                               ##
## 31Aug22                                                                    ##
################################################################################

# This builds on work from 
####  Flight L. The use of health economics in the design and analysis of adaptive clinical trials (Doctoral dissertation, University of Sheffield).
# and is based on the original analysis conducted in Stata/Excel and reported by 
#### Palmer R, Dimairo M, Latimer N, Cross E, Brady M, Enderby P, Bowen A, Julious S, Harrison M, Alshreef A, Bradley E. Computerised speech and language therapy or attention control added to usual care for people with long-term post-stroke aphasia: the Big CACTUS three-arm RCT. Health Technology Assessment (Winchester, England). 2020 Apr;24(19):1.
#### Latimer NR, Bhadhuri A, Alshreef A, Palmer R, Cross E, Dimairo M, Julious S, Cooper C, Enderby P, Brady MC, Bowen A. Self-managed, computerised word finding therapy as an add-on to usual care for chronic aphasia post-stroke: An economic evaluation. Clinical rehabilitation. 2021 May;35(5):703-17.
# NOTE simplifications have been made to facilitate the replication of this analysis 

BC_HE_model_func<-function(in.age,
                        in.prob.death.disease, # 3-month probability 
                        in.p_relapse_6_9_CSLT,
                        in.p_relapse_6_9_UC,
                        in.p_relapse_6_9_AC,
                        # Probability of relapse (9–12 months)
                        in.p_relapse_9_12_CSLT,
                        in.p_relapse_9_12_UC,
                        in.p_relapse_9_12_AC ,   
                        # Probability of relapse (12 months onwards)
                        in.p_relapse_12_CSLT,
                        in.p_relapse_12_UC,
                        in.p_relapse_12_AC,
                        # Probability of good response (0–6 months)
                        in.p_good_0_6_CSLT,
                        in.p_good_0_6_UC,
                        in.p_good_0_6_AC,
                        # Probability of new good response (6–9 months)
                        in.p_good_6_9_CSLT,
                        in.p_good_6_9_UC,
                        in.p_good_6_9_AC,   
                        # Probability of new good response, or renewed response in people who responded at 6 months and relapsed at 9 months (9–12 months)
                        in.p_good_12_CSLT,
                        in.p_good_12_UC,
                        in.p_good_12_AC, 
                        # Costs
                        in.cost_CSLT, 
                        in.cost_UC,
                        in.cost_AC,
                        in.u_aphasia,
                        in.u_response_6,
                        in.u_response_9, 
                        in.u_response_12, 
                        in.cycle.length,
                        in.threshold, 
                        in.DR){
  set.seed(123)
  #######################
  ## Define Parameters ##
  #######################
  par_mod<-define_parameters(
    cyclelength=in.cycle.length, #cycle length in months
    
    ## DEMOGRAPHICES ##
    # Age at baseline
    age.base=in.age,
    
    # Age at each cycle
    age.cycle=age.base+model_time/(12/cyclelength),
    # model_time is cycle number
    
    
    # Probability of death from disease (literature) per cycle (3 months)
    p_death_disease=in.prob.death.disease,
    
    p_death_disease_cycle=ifelse(trunc(age.cycle)==70,
                                    0.025078,
                                    ifelse(trunc(age.cycle)==71,
                                           0.025519,
                                    ifelse(trunc(age.cycle)==72,
                                           0.026174,
                                    ifelse(trunc(age.cycle)==73,
                                           0.026722,
                                          ifelse(trunc(age.cycle)==74,
                                                 0.027422,
                                          ifelse(trunc(age.cycle)==75,
                                                 0.028042,
                                          ifelse(trunc(age.cycle)==76,
                                                 0.028803,
                                          ifelse(trunc(age.cycle)==77,
                                                 0.029618,
                                          ifelse(trunc(age.cycle)==78,
                                                 0.030631,
                                          ifelse(trunc(age.cycle)==79,
                                                 0.031770,
                                          ifelse(trunc(age.cycle)==80,
                                                 0.033317,
                                          ifelse(trunc(age.cycle)==81,
                                                 0.034717,
                                          ifelse(trunc(age.cycle)==82,
                                                 0.036628,
                                          ifelse(trunc(age.cycle)==83,
                                                 0.038519,
                                          ifelse(trunc(age.cycle)==84,
                                                 0.040706,
                                          ifelse(trunc(age.cycle)==85,
                                                 0.043282,
                                          ifelse(trunc(age.cycle)==86,
                                                 0.045776,
                                          ifelse(trunc(age.cycle)==87,
                                                 0.048749,
                                          ifelse(trunc(age.cycle)==88,
                                                 0.051993,
                                          ifelse(trunc(age.cycle)==89,
                                                 0.055438,
                                          ifelse(trunc(age.cycle)==90,
                                                 0.059212,
                                          ifelse(trunc(age.cycle)==91,
                                                 0.062791,
                                          ifelse(trunc(age.cycle)==92,
                                                 0.067610,
                                          ifelse(trunc(age.cycle)==93,
                                                 0.072767,
                                          ifelse(trunc(age.cycle)==94,
                                                 0.076469,
                                          ifelse(trunc(age.cycle)==95,
                                                 0.081268,
                                          ifelse(trunc(age.cycle)==96,
                                                 0.084970,
                                          ifelse(trunc(age.cycle)==97,
                                                 0.091027,
                                          ifelse(trunc(age.cycle)==98,
                                                 0.097160,
                                          ifelse(trunc(age.cycle)==99,
                                                 0.101883,
                                          ifelse(trunc(age.cycle)==100,
                                                 0.107623,p_death_disease
                                          ))))))))))))))))))))))))))))))
                                          ),
    
    
    
    
    ## RELAPSE ##
    
    # Probability of relapse (6–9 months)
    p_relapse_6_9_CSLT=in.p_relapse_6_9_CSLT,
    p_relapse_6_9_UC=in.p_relapse_6_9_UC,
    p_relapse_6_9_AC=in.p_relapse_6_9_AC,
    
    # Probability of relapse (9–12 months)
    p_relapse_9_12_CSLT=in.p_relapse_9_12_CSLT,
    p_relapse_9_12_UC=in.p_relapse_9_12_UC,
    p_relapse_9_12_AC=in.p_relapse_9_12_AC ,   
   
    # Probability of relapse (12 months onwards)
    p_relapse_12_CSLT=in.p_relapse_12_CSLT,
    p_relapse_12_UC=in.p_relapse_12_UC,
    p_relapse_12_AC=in.p_relapse_12_AC ,
  
    
    p_relapse_CSLT_cycle=ifelse(model_time<=2,0,
                                ifelse(model_time==3,p_relapse_6_9_CSLT,
                                ifelse(model_time==4,p_relapse_9_12_CSLT,p_relapse_12_CSLT))),
    
    p_relapse_UC_cycle=ifelse(model_time<=2,0,
                                ifelse(model_time==3,p_relapse_6_9_UC,
                                ifelse(model_time==4,p_relapse_9_12_UC,p_relapse_12_UC))),
    
    
    p_relapse_AC_cycle=ifelse(model_time<=2,0,
                                ifelse(model_time==3,p_relapse_6_9_AC,
                                ifelse(model_time==4,p_relapse_9_12_AC,p_relapse_12_AC))),
    
    
    ## RESPONSE ##
    
    # Probability of good response (0–6 months)
    p_good_0_6_CSLT=in.p_good_0_6_CSLT,
    p_good_0_6_UC=in.p_good_0_6_UC,
    p_good_0_6_AC=in.p_good_0_6_AC,
    
    # Probability of new good response (6–9 months)
    p_good_6_9_CSLT=in.p_good_6_9_CSLT,
    p_good_6_9_UC=in.p_good_6_9_UC,
    p_good_6_9_AC=in.p_good_6_9_AC,   
    
    # Probability of new good response,or renewed response in people who responded at 6 months and relapsed at 9 months (9–12 months)
    p_good_12_CSLT=in.p_good_12_CSLT,
    p_good_12_UC=in.p_good_12_UC,
    p_good_12_AC=in.p_good_12_AC,  
    
    p_good_CSLT_cycle=ifelse(model_time<=2,p_good_0_6_CSLT,
                                ifelse(model_time==3,p_good_6_9_CSLT,
                                       ifelse(model_time==4,p_good_12_CSLT,0))),
    
    p_good_UC_cycle=ifelse(model_time<=2,p_good_0_6_UC,
                           ifelse(model_time==3,p_good_6_9_UC,
                                  ifelse(model_time==4,p_good_12_UC,0))),
    
    
    p_good_AC_cycle=ifelse(model_time<=2,p_good_0_6_AC,
                           ifelse(model_time==3,p_good_6_9_AC,
                                  ifelse(model_time==4,p_good_12_AC,0))),
                              
    
    
  ## HEALTH ECONOMIC PARAMETERS ##
    # Discount rate for cost
    DR.C=in.DR,
    # As a monthly rate
    R.C=rescale_discount_rate(DR.C,4,1),
    
    # Discount rate for effect
    DR.E=in.DR,
    R.E=rescale_discount_rate(DR.E,4,1),
    
    ## COSTS ##
    # cost of the interventions
    cost_CSLT=in.cost_CSLT, 
    cost_UC=in.cost_UC,
    cost_AC=in.cost_AC,
    
  cost_CSLT_cycle=ifelse(model_time<=2,cost_CSLT/2,0), 
  cost_UC_cycle=ifelse(model_time<=2,cost_UC/2,0) ,
  cost_AC_cycle=ifelse(model_time<=2,cost_AC/2,0) ,
  
    
    ## EFFECTS ##
 
    ## Utility Decrement Multiplier ## 
    u_decrement=
                    ifelse(trunc(age.cycle)==65,
                           1,
                    ifelse(trunc(age.cycle)==66,
                           0.994287186,
                    ifelse(trunc(age.cycle)==67,
                           0.988492050,                
                    ifelse(trunc(age.cycle)==68,
                           0.982614592,
                    ifelse(trunc(age.cycle)==69,
                           0.976654813,
                    ifelse(trunc(age.cycle)==70,
                           0.970612711,
                    ifelse(trunc(age.cycle)==71,
                           0.964488288,
                    ifelse(trunc(age.cycle)==72,
                           0.958281543,
                    ifelse(trunc(age.cycle)==73,
                           0.951992475,      
                    ifelse(trunc(age.cycle)==74,
                           0.945621087,
                    ifelse(trunc(age.cycle)==75,
                           0.939167376,
                    ifelse(trunc(age.cycle)==76,
                           0.932631343,
                    ifelse(trunc(age.cycle)==77,
                           0.926012989,
                    ifelse(trunc(age.cycle)==78,
                           0.919312313,
                    ifelse(trunc(age.cycle)==79,
                           0.912529314 ,
                    ifelse(trunc(age.cycle)==80,
                           0.905663994 ,
                    ifelse(trunc(age.cycle)==81,
                           0.898716353,
                    ifelse(trunc(age.cycle)==82,
                           0.891686389,
                    ifelse(trunc(age.cycle)==83,
                           0.884574103,
                    ifelse(trunc(age.cycle)==84,
                           0.877379496,
                    ifelse(trunc(age.cycle)==85,
                           0.870102567,
                    ifelse(trunc(age.cycle)==86,
                           0.862743316,
                    ifelse(trunc(age.cycle)==87,
                           0.855301743,
                    ifelse(trunc(age.cycle)==88,
                           0.847777848,
                    ifelse(trunc(age.cycle)==89,
                           0.840171631,
                    ifelse(trunc(age.cycle)==90,
                           0.832483093,
                    ifelse(trunc(age.cycle)==91,
                           0.824712232,
                    ifelse(trunc(age.cycle)==92,
                           0.816859050,
                    ifelse(trunc(age.cycle)==93,
                           0.808923546,
                    ifelse(trunc(age.cycle)==94,
                           0.800905720,
                    ifelse(trunc(age.cycle)==95,
                           0.792805572,
                    ifelse(trunc(age.cycle)==96,
                           0.784623103,
                    ifelse(trunc(age.cycle)==97,
                           0.776358311,
                    ifelse(trunc(age.cycle)==98,
                           0.768011198,
                    ifelse(trunc(age.cycle)==99,
                           0.759581763,
                    ifelse(trunc(age.cycle)==100,
                           0.751070006,1
                    )))))))))))))))))))))))))))))))))))),
  
    ## Utility ##
  u_aphasia=in.u_aphasia,
  u_aphasia_cycle=(u_aphasia*u_decrement)*cyclelength/12,
  
  u_response_6=in.u_response_6,
  u_response_9=in.u_response_9, 
  u_response_12=in.u_response_12,
  
  u_response_cycle=ifelse(model_time<=2,(u_response_6*u_decrement)*cyclelength/12,
                          ifelse(model_time==3,(u_response_9*u_decrement)*cyclelength/12,(u_response_12*u_decrement)*cyclelength/12)))
  
 
  
  #######################
  ## Transition Matrix ##
  #######################
  
  mat_CSLT<-define_transition(
    state_names=c("APHASIA",
                  "RESPONSE",
                  "DEAD"),
   
    1-(p_good_CSLT_cycle*(1-p_death_disease_cycle)+p_death_disease_cycle),p_good_CSLT_cycle*(1-p_death_disease_cycle),p_death_disease_cycle,
   p_relapse_CSLT_cycle*(1-p_death_disease_cycle),1-(p_death_disease_cycle+p_relapse_CSLT_cycle*(1-p_death_disease_cycle)), p_death_disease_cycle,
    0,0,1
   )
   
  
  
  mat_AC<-define_transition(
    state_names=c("APHASIA",
                  "RESPONSE",
                  "DEAD"),
    
    1-(p_good_AC_cycle*(1-p_death_disease_cycle)+p_death_disease_cycle),p_good_AC_cycle*(1-p_death_disease_cycle),p_death_disease_cycle,
    p_relapse_AC_cycle*(1-p_death_disease_cycle),1-(p_death_disease_cycle+p_relapse_AC_cycle*(1-p_death_disease_cycle)),p_death_disease_cycle,
    0,0,1
  )
  
  
  
  mat_UC<-define_transition(
    state_names=c("APHASIA",
                  "RESPONSE",
                  "DEAD"),
    
    1-(p_good_UC_cycle*(1-p_death_disease_cycle)+p_death_disease_cycle),p_good_UC_cycle*(1-p_death_disease_cycle),p_death_disease_cycle,
    p_relapse_UC_cycle*(1-p_death_disease_cycle),1-(p_death_disease_cycle+p_relapse_UC_cycle*(1-p_death_disease_cycle)),p_death_disease_cycle,
    0,0,1
  )
  
  
  
  ##################
  ## STATE VALUES ##
  ##################
  # Define costs and effects in each of the states
  state_DEAD<-define_state(
    cost=dispatch_strategy(
      CSLT=cost_CSLT_cycle,
      AC=cost_AC_cycle,
      UC=cost_UC_cycle),
    cost_total=discount(cost,R.C,first=FALSE),
    utility=dispatch_strategy(
      CSLT=0,
      AC=0,
      UC=0),
    utility_total=discount(utility,R.E,first=FALSE))
  
  state_APHASIA<-define_state(
    cost=dispatch_strategy(
      CSLT=cost_CSLT_cycle,
      AC=cost_AC_cycle,
      UC=cost_UC_cycle),
    # resource cost in the DEAD state
    cost_total=discount(cost,R.C,first=FALSE),
    # total discounted cost in the DEAD state
    utility=dispatch_strategy(
      CSLT=u_aphasia_cycle,
      AC=u_aphasia_cycle,
      UC=u_aphasia_cycle),
    # effect in the intervention arm for the DEAD state
    utility_total=discount(utility,R.E,first=FALSE)
    # total discounted effect in the DEAD state in each arm
  )
  
  state_RESPONSE<-define_state(
    cost=dispatch_strategy(
      CSLT=cost_CSLT_cycle,
      AC=cost_AC_cycle,
      UC=cost_UC_cycle),
    # resource cost in the DEAD state
    cost_total=discount(cost,R.C,first=FALSE),
    # total discounted cost in the DEAD state
    utility=dispatch_strategy(
      CSLT=u_response_cycle,
      AC=u_response_cycle,
      UC=u_response_cycle),
    # effect in the intervention arm for the DEAD state
    utility_total=discount(utility,R.E,first=FALSE)
    # total discounted effect in the DEAD state in each arm
  )
  
 
  
  ################
  ## STRATEGIES ##
  ################
  strat_CSLT<-define_strategy(
    transition=mat_CSLT,
    
    APHASIA=state_APHASIA,
    RESPONSE=state_RESPONSE,
    DEAD=state_DEAD)
  
  strat_UC<-define_strategy(
    transition = mat_UC,
    
    APHASIA=state_APHASIA,
    RESPONSE=state_RESPONSE,
    DEAD=state_DEAD
    )
  
  strat_AC<-define_strategy(
    transition = mat_AC,
        APHASIA=state_APHASIA,
    RESPONSE=state_RESPONSE,
    DEAD=state_DEAD)
  
  
  ###############
  ## Run Model ##
  ###############
  res_mod<-run_model(
    parameters=par_mod,
    CSLT=strat_CSLT,
    AC=strat_AC,
    UC=strat_UC,
    cycles=142, # run for 50 years
    cost=cost_total,
    effect=utility_total,
    method="end",
    init=c(1000,0,0)
    # life-table correction or beginning or end (no half cycle option)
  )
  

QQ<-res_mod$eval_strategy_list$UC
 
  
  total.cost_CSLT<-sum(res_mod$eval_strategy_list$CSLT$values$cost_total)
  total.QALY_CSLT<-sum(res_mod$eval_strategy_list$CSLT$values$utility_total)
  total.cost_AC<-sum(res_mod$eval_strategy_list$AC$values$cost_total)
  total.QALY_AC<-sum(res_mod$eval_strategy_list$AC$values$utility_total)
  total.cost_UC<-sum(res_mod$eval_strategy_list$UC$values$cost_total)
  total.QALY_UC<-sum(res_mod$eval_strategy_list$UC$values$utility_total)
  
  
  
  ICER_CSLT_AC<-(total.cost_CSLT-total.cost_AC)/(total.QALY_CSLT-total.QALY_AC)
  ICER_CSLT_UC<-(total.cost_CSLT-total.cost_UC)/(total.QALY_CSLT-total.QALY_UC)
  
  n.individ_CSLT<-res_mod$eval_strategy_list$CSLT$n_indiv
  n.individ_AC<-res_mod$eval_strategy_list$AC$n_indiv
  n.individ_UC<-res_mod$eval_strategy_list$UC$n_indiv
  # 
  cost.diff_CSLT_AC<-total.cost_CSLT/n.individ_CSLT-total.cost_AC/n.individ_AC
  QALY.diff_CSLT_AC<-total.QALY_CSLT/n.individ_CSLT-total.QALY_AC/n.individ_AC
  cost.diff_CSLT_UC<-total.cost_CSLT/n.individ_CSLT-total.cost_UC/n.individ_UC
  QALY.diff_CSLT_UC<-total.QALY_CSLT/n.individ_CSLT-total.QALY_UC/n.individ_UC
  
  
  QALY_CSLT<-total.QALY_CSLT/n.individ_CSLT
  QALY_AC<-total.QALY_AC/n.individ_AC
  QALY_UC<-total.QALY_UC/n.individ_UC
  
  Cost_CSLT<-total.cost_CSLT/n.individ_CSLT
  Cost_AC<-total.cost_AC/n.individ_AC
  Cost_UC<-total.cost_UC/n.individ_UC
  
  NB_CSLT<-(in.threshold*total.QALY_CSLT-total.cost_CSLT)/n.individ_CSLT
  NB_AC<-(in.threshold*total.QALY_AC-total.cost_AC)/n.individ_AC
  NB_UC<-(in.threshold*total.QALY_UC-total.cost_UC)/n.individ_UC
  
  INB_CSLT_AC<-(in.threshold*(QALY.diff_CSLT_AC)-
                cost.diff_CSLT_AC  )
  INB_CSLT_UC<-(in.threshold*(QALY.diff_CSLT_UC)-
                  cost.diff_CSLT_UC  )
  param<-c(in.age,
           in.prob.death.disease, # 3-month probability 
           in.p_relapse_6_9_CSLT,
           in.p_relapse_6_9_UC,
           in.p_relapse_6_9_AC,
           # Probability of relapse (9–12 months)
           in.p_relapse_9_12_CSLT,
           in.p_relapse_9_12_UC,
           in.p_relapse_9_12_AC ,   
           # Probability of relapse (12 months onwards)
           in.p_relapse_12_CSLT,
           in.p_relapse_12_UC,
           in.p_relapse_12_AC,
           # Probability of good response (0–6 months)
           in.p_good_0_6_CSLT,
           in.p_good_0_6_UC,
           in.p_good_0_6_AC,
           # Probability of new good response (6–9 months)
           in.p_good_6_9_CSLT,
           in.p_good_6_9_UC,
           in.p_good_6_9_AC,   
           # Probability of new good response, or renewed response in people who responded at 6 months and relapsed at 9 months (9–12 months)
           in.p_good_12_CSLT,
           in.p_good_12_UC,
           in.p_good_12_AC, 
           # Costs
           in.cost_CSLT, 
           in.cost_UC,
           in.cost_AC,
           in.u_aphasia,
           in.u_response_6,
           in.u_response_9, 
           in.u_response_12, 
           in.cycle.length,
           in.threshold, 
           in.DR)
  out=list(NB_CSLT=NB_CSLT,NB_AC=NB_AC,NB_UC=NB_UC,
           INB_CSLT_AC=INB_CSLT_AC,INB_CSLT_UC=INB_CSLT_UC,
           ICER_CSLT_AC=ICER_CSLT_AC,ICER_CSLT_UC=ICER_CSLT_UC,
           res_mod=res_mod,
           cost.diff_CSLT_AC=cost.diff_CSLT_AC,cost.diff_CSLT_UC=cost.diff_CSLT_UC,
           QALY.diff_CSLT_AC=QALY.diff_CSLT_AC,QALY.diff_CSLT_UC=QALY.diff_CSLT_UC,
           total.cost_CSLT=total.cost_CSLT,
           total.QALY_CSLT=total.QALY_CSLT,
           total.cost_AC=total.cost_AC,
           total.QALY_AC=total.QALY_AC,
           total.cost_UC=total.cost_UC,
           total.QALY_UC=total.QALY_UC,
           param=param)
           
  out  
}
