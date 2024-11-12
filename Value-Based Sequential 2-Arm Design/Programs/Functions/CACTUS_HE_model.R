# Function based on the original CACTUS pilot Excel based health economic model

# There have been some simplifications to this application for ease of conduct
# Full details and comparison given in Chapter 5 Flight L. The use of health economics in the design and analysis of adaptive clinical trials (Doctoral dissertation, University of Sheffield).
# Original analysis reported by Latimer NR, Dixon S, Palmer R. Cost-utility of self-managed computer therapy for people with aphasia. International journal of technology assessment in health care. 2013 Oct;29(4):402-9. 

HE_model_func<-function(in_age,
                       in_prob_death_disease,
                       in_prob_relapse_month_int,
                       in_prob_good_response_int, 
                       in_cost_intervention,
                       in_cost_cont, # cost not including intervention cost
                       in_cost_int_resp, # cost not including intervention cost
                       in_cost_int_non_resp,
                       in_dis_no_resp_cont,
                       in_utility_change_int,
                       threshold)
{
  
  set.seed(123)
  
  #######################
  ## Define Parameters ##
  #######################
  
  par_mod<-define_parameters(
    cyclelength=12, # cycle length is months
    
    # Demographics -------------------------------------------------------------
    
    # Age at baseline
    age_base=in_age,
    # Age at each cycle
    age_cycle=age_base+model_time/cyclelength, # model_time is cycle number
    # Gender
    sex_male="MLE",
    # Probability of death from disease (literature) per month
    prob_death_disease=in_prob_death_disease,
    prob_death_disease_cycle=ifelse(trunc(age_cycle-0.07)==73,
                                    0.0113740263,
                                  ifelse(trunc(age_cycle-0.07)==74,
                                         0.0116037087125002,
                                  ifelse(trunc(age_cycle-0.07)==75,
                                         0.0118961575630027,
                                  ifelse(trunc(age_cycle-0.07)==76,
                                         0.012222392402447,
                                  ifelse(trunc(age_cycle-0.07)==77,
                                         0.0125523652497712,
                                  ifelse(trunc(age_cycle-0.07)==78,
                                         0.0129242959435807,
                                  ifelse(trunc(age_cycle-0.07)==79,
                                         0.013427307526383,
                                  ifelse(trunc(age_cycle-0.07)==80,
                                         0.013945178719201,
                                  ifelse(trunc(age_cycle-0.07)==81,
                                         0.0144983961217108,
                                  ifelse(trunc(age_cycle-0.07)==82,
                                         0.0150780074615502,
                                  ifelse(trunc(age_cycle-0.07)==83,
                                         0.0157420385460144,
                                  ifelse(trunc(age_cycle-0.07)==84,
                                         0.0165333900569085,
                                  ifelse(trunc(age_cycle-0.07)==85,
                                         0.0173376782101331,
                                  ifelse(trunc(age_cycle-0.07)==86,
                                         0.0182089366562969,
                                  ifelse(trunc(age_cycle-0.07)==87,
                                         0.0190907135893689,
                                  ifelse(trunc(age_cycle-0.07)==88,
                                         0.0196915610947849,
                                  ifelse(trunc(age_cycle-0.07)==89,
                                         0.0206586949271182,
                                  ifelse(trunc(age_cycle-0.07)==90,
                                         0.0215852581532521,
                                  ifelse(trunc(age_cycle-0.07)==91,
                                         0.0232433895309773,
                                  ifelse(trunc(age_cycle-0.07)==92,
                                         0.02497065466537,
                                  ifelse(trunc(age_cycle-0.07)==93,
                                         0.0263318908753054,
                                  ifelse(trunc(age_cycle-0.07)==94,
                                         0.0276624039008831,
                                  ifelse(trunc(age_cycle-0.07)==95,
                                         0.029774315827079,
                                  ifelse(trunc(age_cycle-0.07)==96,
                                         0.0310891942245149,
                                  ifelse(trunc(age_cycle-0.07)==97,
                                         0.0327137783402275,
                                  ifelse(trunc(age_cycle-0.07)==98,
                                         0.034242655440335,
                                  ifelse(trunc(age_cycle-0.07)==99,
                                         0.0355802603399481,
                                  ifelse(trunc(age_cycle-0.07)==100,
                                         0.0372704690687986
                                         ,prob_death_disease
                                         )))))))))))))))))))))))))))
                                  ),
                                              
    # Relapse ------------------------------------------------------------------
    
    # Probability of relapse per month
    # added in after month 6
    prob_relapse_month=in_prob_relapse_month_int,
    
    prob_relapse_month_cycle=ifelse(model_time<=5,
                                    0,
                                    prob_relapse_month),
    
    # Probability of death from disease (literature) per month
    prob_death_disease=in_prob_death_disease,
    
    # Response -----------------------------------------------------------------
    
    # Probability of a good response intervention arm applied at the first month
    prob_good_response_int=in_prob_good_response_int,
    prob_good_int=ifelse(model_time==1,
                         prob_good_response_int,
                         0),
    
    # Probability of a good response control arm applied at the first month
    prob_good_response_cont=0,
    prob_good_cont=ifelse(model_time==1,
                          prob_good_response_cont,
                          0),
    
    # Health economic parameters -----------------------------------------------
    
    # Discount rate for cost
    DR_C=0.035,
    # As a monthly rate
    R_C=rescale_discount_rate(DR_C,cyclelength,1),
    # Discount rate for effect
    DR_E=0.035,
    R_E=rescale_discount_rate(DR_E,cyclelength,1),
    
    # Cost ---------------------------------------------------------------------
    
    # Cost of the intervention
    cost_intervention=in_cost_intervention, 
    # Total cost in the control arm
    cost_cont=in_cost_cont,
    # Resource costs control arm per month  
    cost_resource_cont=cost_cont, 
    # Resource costs intervention arm per arm
    cost_resource_int_resp=in_cost_int_resp, 
    cost_resource_int_non_resp=in_cost_int_non_resp, 
    # Intervention group - no response
    # Intervention cost is added during the first cycle only
    total_cost_intervention_no_cycle=ifelse(model_time==1,
                                   cost_resource_int_non_resp+cost_intervention,
                                   ifelse(model_time>=2&model_time<=5,
                                          cost_resource_int_non_resp,
                                          cost_resource_cont)),
    # Intervention group - good response
    # Intervention cost is added during the first cycle only
    total_cost_intervention_good_cycle=ifelse(model_time==1,
                                     cost_resource_int_resp+cost_intervention,
                                     ifelse(model_time>=2&model_time<=5,
                                            cost_resource_int_resp,
                                            cost_resource_cont)),
    # Control group - good response
    # Control cost is added during the first cycle only
    total_cost_control_good_cycle=ifelse(model_time==1,
                                      cost_resource_cont,
                                      cost_resource_cont),
    # Control group - no response
    # Control cost is added during the first cycle only
    total_cost_control_no_cycle=ifelse(model_time==1,
                                    cost_resource_cont,
                                    cost_resource_cont),
    # Intervention group - dead
    total_cost_intervention_dead_cycle=ifelse(model_time==1,
                                     cost_intervention,
                                     0),
    
    # Effects ------------------------------------------------------------------
    
    # Utility Decrements Multiplier
    utility_decrement=ifelse(trunc(age_cycle-0.07)==68,
                         1,
                  ifelse(trunc(age_cycle-0.07)==69,
                         0.993936937133514,
                  ifelse(trunc(age_cycle-0.07)==70,
                         0.987790125775351,
                  ifelse(trunc(age_cycle-0.07)==71,
                         0.981559565925514,
                  ifelse(trunc(age_cycle-0.07)==72,
                         0.975245257584,
                  ifelse(trunc(age_cycle-0.07)==73,
                         0.96884720075081,      
                  ifelse(trunc(age_cycle-0.07)==74,
                         0.962365395425945,
                  ifelse(trunc(age_cycle-0.07)==75,
                         0.955799841609404,
                  ifelse(trunc(age_cycle-0.07)==76,
                         0.949150539301187,
                  ifelse(trunc(age_cycle-0.07)==77,
                         0.942417488501294,
                  ifelse(trunc(age_cycle-0.07)==78,
                         0.935600689209726,
                  ifelse(trunc(age_cycle-0.07)==79,
                         0.928700141426481,
                  ifelse(trunc(age_cycle-0.07)==80,
                         0.921715845151561,
                  ifelse(trunc(age_cycle-0.07)==81,
                         0.914647800384966,
                  ifelse(trunc(age_cycle-0.07)==82,
                         0.907496007126694,
                  ifelse(trunc(age_cycle-0.07)==83,
                         0.900260465376747,
                  ifelse(trunc(age_cycle-0.07)==84,
                         0.892941175135123,
                  ifelse(trunc(age_cycle-0.07)==85,
                         0.885538136401824,
                  ifelse(trunc(age_cycle-0.07)==86,
                         0.878051349176849,
                  ifelse(trunc(age_cycle-0.07)==87,
                         0.870480813460199,
                  ifelse(trunc(age_cycle-0.07)==88,
                         0.862826529251872,
                  ifelse(trunc(age_cycle-0.07)==89,
                         0.85508849655187,
                  ifelse(trunc(age_cycle-0.07)==90,
                         0.847266715360192,
                  ifelse(trunc(age_cycle-0.07)==91,
                         0.839361185676839,
                  ifelse(trunc(age_cycle-07)==92,
                         0.831371907501809,
                  ifelse(trunc(age_cycle-0.07)==93,
                         0.823298880835104,
                  ifelse(trunc(age_cycle-0.07)==94,
                         0.815142105676723,
                  ifelse(trunc(age_cycle-0.07)==95,
                         0.806901582026666,
                  ifelse(trunc(age_cycle-0.07)==96,
                         0.798577309884933,
                  ifelse(trunc(age_cycle-0.07)==97,
                         0.790169289251524,
                  ifelse(trunc(age_cycle-0.07)==98,
                         0.78167752012644,
                  ifelse(trunc(age_cycle-0.07)==99,
                         0.77310200250968,
                  ifelse(trunc(age_cycle-0.07)==100,
                         0.764442736401244
                  ,1
                  ))))))))))))))))))))))))))))))))),
    
    # Effect in the control arm when no response
    dis_no_resp_cont=in_dis_no_resp_cont,
    utility_no_resp_cont= (1-dis_no_resp_cont),
    utility_no_resp_cont_cycle=  (utility_no_resp_cont)*utility_decrement/cyclelength,
    
    # Effect in the intervention arm when no response
    dis_no_resp_int=in_dis_no_resp_cont,
    utility_no_resp_int=(1-dis_no_resp_int), #Value taken from pilot
    utility_no_resp_int_cycle=(utility_no_resp_int)*utility_decrement/cyclelength,
    
    # Change from baseline in effect for the intervention arm
    utility_change_int=in_utility_change_int,
    
    # Effect in intervention arm when good response (see "Model Experimental Cell J8")
    utility_good_resp_int1=(utility_no_resp_int+utility_change_int), 
    
    # If utility in the good response state is less than unwell then use unwelll
    utility_good_resp_int=(ifelse(utility_good_resp_int1>utility_no_resp_int,
                                 utility_good_resp_int1,
                                 utility_no_resp_int)),

    # Add the change from baseline in effect for intervention arm from the trial
    utility_good_resp_int_cycle=(utility_good_resp_int[[1]])*utility_decrement/cyclelength,
    
    # Effect in the control arm when good response
    utility_good_resp_cont= utility_no_resp_cont,
    utility_good_resp_cont_cycle=  (utility_good_resp_cont)*utility_decrement/cyclelength
      )
  
  #######################
  ## Transition Matrix ##
  #######################
  
  # Intervention arm -----------------------------------------------------------
  
  # min(prob_death_disease_cycle,1)
  
    mat_intervention<-define_transition(
      state_names=c("UNWELL","WELL","DEAD"),
      
      1-(prob_good_int+prob_death_disease_cycle),
      
      prob_good_int, prob_death_disease_cycle,
      
      prob_relapse_month_cycle*(1-prob_death_disease_cycle),
      
      1-( prob_relapse_month_cycle*(1-prob_death_disease_cycle)+prob_death_disease_cycle),
      
      prob_death_disease_cycle,
      
      0,0,1)
    # state_names=c("UNWELL","WELL","DEAD"),
    # 
    #     1-(min((1-0.0372704690687986),prob_good_int)+prob_death_disease_cycle),
    # 
    # min((1-0.0372704690687986),prob_good_int), prob_death_disease_cycle,
    # 
    # prob_relapse_month_cycle*(1-prob_death_disease_cycle),
    # 
    # 1-( prob_relapse_month_cycle*(1-prob_death_disease_cycle)+prob_death_disease_cycle),
    # 
    # prob_death_disease_cycle,
    # 
    # 0,0,1)
  
  # Control arm ----------------------------------------------------------------
    
  mat_control<-define_transition(
    
    state_names=c("UNWELL","WELL","DEAD"),
    
    1-(prob_good_cont+prob_death_disease_cycle),
    
    prob_good_cont, prob_death_disease_cycle,
    
    0, 1-prob_death_disease_cycle,prob_death_disease_cycle,
    
    0,0,1)
  
  ##################
  ## State Values ##
  ##################
  
  # Define costs and effects in each of the states
  state_DEAD<-define_state(
    cost=dispatch_strategy(
      control=0,
      intervention=total_cost_intervention_dead_cycle), # resource cost in the DEAD state
    
    cost_total=discount(cost,DR_C,first=FALSE,period=12), # total discounted cost in the DEAD state
   
    utility=dispatch_strategy(
      control=0, # effect in the control arm for the DEAD state
      intervention=0), # effect in the intervention arm for the DEAD state
    
    utility_total=discount(utility,DR_E,first=FALSE,period=12) # total discounted effect in the DEAD state in each arm
    
  )
  
  state_WELL<-define_state(
    cost=dispatch_strategy(
      control=total_cost_control_good_cycle, # total cost in the WELL state for control arm
     
      intervention=total_cost_intervention_good_cycle), # total cost in the WELL state for the intervention arm
   
    cost_total=discount(cost,DR_C,first=FALSE,period=12),
    utility=dispatch_strategy(
      control=utility_good_resp_cont_cycle, # effect in the WELL state for control arm
      
      intervention= utility_good_resp_int_cycle), # effect in the WELL state for the intervention arm
    
    utility_total=discount(utility,DR_E,first=FALSE,period=12))
  
  state_UNWELL<-define_state(
    cost=dispatch_strategy(
      control=total_cost_control_no_cycle, # total cost in the UNWELL state for control arm
     
      intervention=total_cost_intervention_no_cycle), # total cost in the UNWELL state for intervention arm
   
    cost_total=discount(cost,DR_C,first=FALSE,period=12),
    utility=dispatch_strategy(
      control=utility_no_resp_cont_cycle, # effect in the UNWELL state for the control arm
      
      intervention=utility_no_resp_int_cycle), # effect in the UNWELL state for the intervention arm
   
    utility_total=discount(utility,DR_E,first=FALSE,period=12))
  
  ################
  ## Strategies ##
  ################
  
  strat_intervention<-define_strategy(
    transition=mat_intervention,
    UNWELL=state_UNWELL,
    WELL=state_WELL,
    DEAD=state_DEAD)
  
  strat_control<-define_strategy(
    transition = mat_control,
    DEAD=state_DEAD,
    WELL=state_WELL,
    UNWELL=state_UNWELL)
  
  ###############
  ## Run Model ##
  ###############
  
  res_mod<-run_model(
    parameters=par_mod,
    intervention=strat_intervention,
    control=strat_control,
    cycles=385, # run for 50 years
    cost=cost_total,
    effect=utility_total,
    method="end" # life-table correction or beginning or end (no half cycle option)
   
  )
  
  total_cost_intervention<-sum(res_mod$eval_strategy_list$intervention$values$cost_total)
  total_qaly_interventionervention<-sum(res_mod$eval_strategy_list$intervention$values$utility_total)
  total_cost_control<-sum(res_mod$eval_strategy_list$control$values$cost_total)
  total_qaly_controlrol<-sum(res_mod$eval_strategy_list$control$values$utility_total)

  model_ICER<-(total_cost_intervention-total_cost_control)/(total_qaly_interventionervention-total_qaly_controlrol)
  
  n_intervention<-res_mod$eval_strategy_list$intervention$n_indiv
  n_control<-res_mod$eval_strategy_list$control$n_indiv
   
  cost_difference<-total_cost_intervention/n_intervention-total_cost_control/n_control
  qaly_difference<-total_qaly_interventionervention/n_intervention-total_qaly_controlrol/n_control
  
  qaly_intervention<-total_qaly_interventionervention/n_intervention
  qaly_control<-total_qaly_controlrol/n_control
  cost_intervention<-total_cost_intervention/n_intervention
  cost_control<-total_cost_control/n_control
  
  nb_intervention<-(threshold*total_qaly_interventionervention-total_cost_intervention)/n_intervention
  nb_control<-(threshold*total_qaly_controlrol-total_cost_control)/n_control
  
  model_INB<-(threshold*(qaly_difference)-
                cost_difference  )
  
  out=list(nb_intervention=nb_intervention,nb_control=nb_control,model_INB=model_INB,model_ICER=model_ICER,
           #res_mod=res_mod,#creating warning when applying to boostraps so remove for now (LF: 30Aug22)
           cost_difference=cost_difference,
           qaly_difference=qaly_difference,qaly_intervention=qaly_intervention,qaly_control=qaly_control,cost_intervention=cost_intervention,cost_control=cost_control)
  out  
  }
