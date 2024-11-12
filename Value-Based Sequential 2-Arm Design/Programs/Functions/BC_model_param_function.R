################################################################################
## Function to calculate health economic model parameters required            ##
## Laura Flight                                                               ## 
## 31Aug22                                                                    ##
################################################################################

param_func<-function(data){
  library(dplyr)
  
  # CSLT total cost for model --------------------------------------------------
  cost_CSLT<-data%>%
    summarize(mean_SLT=mean(hours_therapist,na.rm=T),
              mean_SLT_training=mean(hours_therapist_training,na.rm=T),
              mean_SLTA=mean(hours_slta_ppt,na.rm=T),
              mean_SLTA_training=mean(hours_slta_training,na.rm=T),
              mean_cost_travel=mean(cost_travel,na.rm=T))%>%
    mutate(cost_SLT=mean_SLT*45)%>%
    mutate(cost_SLT_training=mean_SLT_training*45)%>%
    mutate(cost_SLTA=mean_SLTA*25)%>%
    mutate(cost_SLTA_training=mean_SLTA_training*25)%>%
    mutate(cost_therapist=cost_SLT+cost_SLT_training+cost_SLTA+cost_SLTA_training)%>%
    mutate(cost_travel=mean_cost_travel)%>%
    mutate(cost_therapist_training=10.33)%>%
    mutate(cost_hardware=46.64845+4.783505+109.84)%>%
    mutate(cost_CSLT=cost_therapist+cost_travel+cost_hardware+cost_therapist_training)%>%
    select(cost_CSLT)
  
  # AC total cost for model ----------------------------------------------------
  cost_AC<-data%>%
    summarize(mean_book=mean(book_total,na.rm=T),mean_hours=mean(hours_puzzle,na.rm=T))%>%
    mutate(cost_book=mean_book*2.5)%>%
    mutate(cost_therapist=mean_hours*34.2)%>%
    mutate(AC_totalcost=cost_book+cost_therapist)%>%
    select(AC_totalcost)  
  
  ########################################
  # Calculating transition probabilities #
  ########################################
  
  # Calculating improvements in primary outcomes -------------------------------
  transition.prob<-data%>%
    mutate(changewords6 = (perc_wordfindscore6 - perc_wordfindscore0))%>%
    mutate(changetoms6 = score6 - score0)%>%
    mutate( changewords9 = (perc_wordfindscore9 - perc_wordfindscore0))%>%
    mutate( changetoms9 = score9 - score0)%>%
    mutate( changewords12 = (perc_wordfindscore12 - perc_wordfindscore0))%>%
    mutate( changetoms12 = score12 - score0)%>%
    # *6 month transition probabilities* 
    mutate(response_6=ifelse(changewords6>=10 | changetoms6>=0.5,1,0))%>%
    # *9 month new response transition probabilities* 
    mutate(response_9=ifelse((changewords9>=10 | changetoms9>=0.5),1,0))%>%
    # *12 month new response transition probabilities*
    mutate(response_12=ifelse((changewords12>=10 | changetoms12>=0.5),1,0))%>%
    group_by(random_group)%>%
    summarize(mean_response_6=mean(response_6,na.rm=T),
              mean_response_9=mean(response_9[response_6==0],na.rm=T),
              mean_relapse_9=1-mean(response_9[response_6==1],na.rm=T),
              mean_response_12=mean(response_12[response_9==0],na.rm=T),
              mean_relapse_12=1-mean(response_12[response_9==1],na.rm=T))%>% 
    replace(is.na(.),0)
  
  ##################################
  # Calculating utility parameters #
  ##################################
  
  utilities<-data%>%
    mutate(changewords6 = (perc_wordfindscore6 - perc_wordfindscore0))%>%
    mutate(changetoms6 = score6 - score0)%>%
    mutate( changewords9 = (perc_wordfindscore9 - perc_wordfindscore0))%>%
    mutate( changetoms9 = score9 - score0)%>%
    mutate( changewords12 = (perc_wordfindscore12 - perc_wordfindscore0))%>%
    mutate( changetoms12 = score12 - score0)%>%
    # *6 month transition probabilities* 
    mutate(response_6=ifelse(changewords6>=10 | changetoms6>=0.5,1,0))%>%
    # *9 month new response transition probabilities* 
    mutate(response_9=ifelse((changewords9>=10 | changetoms9>=0.5),1,0))%>%
    mutate(response_9_new=ifelse((changewords9>=10 | changetoms9>=0.5)&response_6==0,1,0))%>%
    # *12 month new response transition probabilities*
    mutate(response_12=ifelse((changewords12>=10 | changetoms12>=0.5),1,0))%>%
    mutate(response_12_new=ifelse((changewords12>=10 | changetoms12>=0.5)&response_9==0,1,0))%>%
    mutate(change_eq5d_6 = EQ5D3L6 - EQ5D3L0)%>%
    mutate(change_eq5d_9 = EQ5D3L9 - EQ5D3L0)%>%
    mutate(change_eq5d_12 = EQ5D3L12 - EQ5D3L0)%>%
    mutate(change_eq5d_6_good = ifelse(response_6==1,change_eq5d_6,NA))%>%
    mutate(change_eq5d_6_no = ifelse(response_6==0,change_eq5d_6,NA))%>%
    mutate(change_eq5d_9 = EQ5D3L9 - EQ5D3L0)%>%
    mutate(change_eq5d_9_good = ifelse(response_9==1,change_eq5d_9,NA))%>%
    mutate(change_eq5d_9_no = ifelse(response_9==0,change_eq5d_9,NA))%>%
    mutate(change_eq5d_12 = EQ5D3L12 - EQ5D3L0)%>%
    mutate(change_eq5d_12_good = ifelse(response_12==1,change_eq5d_12,NA))%>%
    mutate(change_eq5d_12_no = ifelse(response_12==0,change_eq5d_12,NA))%>%
    summarize(mean_EQ5D3L0=mean(EQ5D3L0,na.rm=T),# mean utility at baseline for all participants
              mean_eq5d_6_good=mean(change_eq5d_6_good,na.rm=T),# mean utility at 6 months for responders
              mean_eq5d_6_no=mean(change_eq5d_6_no,na.rm=T),# mean utility at 6 months for non-responders
              mean_eq5d_9_good=mean(change_eq5d_9_good,na.rm=T),# mean utility at 9 months for responders
              mean_eq5d_9_no=mean(change_eq5d_9_no,na.rm=T),# mean utility at 9 months for n0n-responders
              mean_eq5d_12_good=mean(change_eq5d_12_good,na.rm=T), # mean utility at 12 months for responders
              mean_eq5d_12_no=mean(change_eq5d_12_no,na.rm=T))%>%# mean utility at 12 months for non-responders
    replace(is.na(.),0) %>%
    mutate(utility0=mean_EQ5D3L0)%>%
    mutate(u_response_6=mean_eq5d_6_good-mean_eq5d_6_no)%>%
    mutate(u_response_9=mean_eq5d_9_good-mean_eq5d_9_no)%>%
    mutate(u_response_12=mean_eq5d_12_good-mean_eq5d_12_no)%>%
    select(utility0,u_response_6,u_response_9,u_response_12)
  
  n_CSLT<- data%>%
    filter(random_group==2)%>%
    summarize(n_CSLT=n())
  n_UC<- data%>%
    filter(random_group==1)%>%
    summarize(n_UC=n())
  n_AC<- data%>%
    filter(random_group==3)%>%
    summarize(n_AC=n())
  
  out=list(n_CSLT=n_CSLT[1],
           n_UC=n_UC[1],
           n_AC=n_AC[1],
           cost_AC=cost_AC,
           cost_CSLT=cost_CSLT,
           utility0=utilities$utility0,
           u_response_6=utilities$u_response_6,
           u_response_9=utilities$u_response_9,
           u_response_12=utilities$u_response_12,
           response_6_CSLT=transition.prob%>%
             filter(random_group==2)%>%
             select(mean_response_6),
           response_9_CSLT=transition.prob%>%
             filter(random_group==2)%>%
             select(mean_response_9),
           response_12_CSLT=transition.prob%>%
             filter(random_group==2)%>%
             select(mean_response_12),
           relapse_9_CSLT=transition.prob%>%
             filter(random_group==2)%>%
             select(mean_relapse_9),
           relapse_12_CSLT=transition.prob%>%
             filter(random_group==2)%>%
             select(mean_relapse_12),
           response_6_UC=transition.prob%>%
             filter(random_group==1)%>%
             select(mean_response_6),
           response_9_UC=transition.prob%>%
             filter(random_group==1)%>%
             select(mean_response_9),
           response_12_UC=transition.prob%>%
             filter(random_group==1)%>%
             select(mean_response_12),
           relapse_9_UC=transition.prob%>%
             filter(random_group==1)%>%
             select(mean_relapse_9),
           relapse_12_UC=transition.prob%>%
             filter(random_group==1)%>%
             select(mean_relapse_12),
           response_6_AC=transition.prob%>%
             filter(random_group==3)%>%
             select(mean_response_6),
           response_9_AC=transition.prob%>%
             filter(random_group==3)%>%
             select(mean_response_9),
           response_12_AC=transition.prob%>%
             filter(random_group==3)%>%
             select(mean_response_12),
           relapse_9_AC=transition.prob%>%
             filter(random_group==3)%>%
             select(mean_relapse_9),
           relapse_12_AC=transition.prob%>%
             filter(random_group==3)%>%
             select(mean_relapse_12))
  out
}