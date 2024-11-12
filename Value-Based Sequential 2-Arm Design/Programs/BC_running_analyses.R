################################################################################
## Annotated code to analyse the bootstrap samples                            ##
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

pkgs<-c("dplyr","ggplot2","reshape2","tidyr","tibble","xtable","heemod","foreach","doParallel","scales")
inst<-lapply(pkgs,library,character.only=TRUE)

# Set working directory --------------------------------------------------------

setwd("xxxxxx/Value-Based Sequential 2-Arm Design")

# Load functions (should be saved in the working directory) --------------------

source("Programs/Functions/BC_model_param_function.R")
source("Programs/Functions/Big_CACTUS_R_model.R")
source("Programs/Functions/BC_bootstrap_analysis_function.R")

################################################################################

# Load data --------------------------------------------------------------------

DATA <- read.csv("Data/Raw/data_set.csv")

################################################################################

# Merge the bootstrap files ----------------------------------------------------

## 95 allocations --------------------------------------------------------------

boot_95 <- rbind(read.csv("Data/Derived/Bootstraps/block5_CSLTvsUC_Tmax_95_attempt2_B_2500a_0.csv"),
               read.csv("Data/Derived/Bootstraps/block5_CSLTvsUC_Tmax_95_attempt2_B_2500a_2500.csv") %>%
        mutate(label=label+2500))

output_95 <- boot_analysis_func(tvec=t(read.csv("Data/Derived/Tmax 95/tvec.csv",header=F)),
                                bndupper=t(read.csv("Data/Derived/Tmax 95/bndupper.csv",header=F)),
                                bndlower=t(read.csv("Data/Derived/Tmax 95/bndlower.csv",header=F)),
                                priormean=3190.78,
                                
                                n0=7,
                                Tmax=95,
                                plot.dat.mat.frame=boot_95,
                                B=5000,
                                P=215378,
                                pre.rec.cost=0,
                                fixed.recruit.follow.cost=0,
                                pp.recruit.cost=1647,
                                pp.followup.cost=706,
                                trial.end.cost=0,
                                
                                block=5)

## 435 allocations -------------------------------------------------------------

boot_435 <- rbind(read.csv("Data/Derived/Bootstraps/block5_CSLTvsUC_Tmax_435_attempt2_B_2500a_0.csv"),
                  read.csv("Data/Derived/Bootstraps/block5_CSLTvsUC_Tmax_435_attempt2_B_2500a_2500.csv") %>%
                    mutate(label=label+2500))

output_435 <- boot_analysis_func(tvec=t(read.csv("Data/Derived/Tmax 435/tvec.csv",header=F)),
                                 bndupper=t(read.csv("Data/Derived/Tmax 435/bndupper.csv",header=F)),
                                 bndlower=t(read.csv("Data/Derived/Tmax 435/bndlower.csv",header=F)),
                                 priormean=3190.78,
                                 
                                 n0=7,
                                 Tmax=435,
                                 plot.dat.mat.frame=boot_435,
                                 B=5000,
                                 P=215378,
                                 pre.rec.cost=0,
                                 fixed.recruit.follow.cost=0,
                                 pp.recruit.cost=1647,
                                 pp.followup.cost=706,
                                 trial.end.cost=0,
                                 
                                 block=5)
