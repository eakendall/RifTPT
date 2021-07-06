library(tidyverse)
library(readxl)
library(abind)
library(patchwork)
library(reshape2)
library(gridExtra)
library(RColorBrewer)
library(lhs)
library(purrr)
library(sensitivity)
library(betafunctions)
library(cowplot)
library(ggpubr)

tag <- "20210603"

## Setup and utility functions ####

# Read and format external param file, 
read.params <- function()
{
  params <- read_excel(paste0("RIF TPT params ",tag,".xlsx"), 
                       col_types = c("text","text","text",
                                      "numeric","numeric","numeric",
                                      "text","text","text","text","text","text"))
  params[,2:3][is.na(params[,2:3])] <- ""
  return(params)
}


# converts OR for incidence versus placebo into a RR of progression with 
# TPT versus placebo if adherent
tptEfficacyConversion <- function(OR, percentReinfection, 
                                  percentNonadherent, cumulativeIncidence=0.05)
{
  RR <- OR/
    (1 - cumulativeIncidence + cumulativeIncidence*OR)
  
  RR_progression <- 1 - ( 1 -(RR-percentReinfection)/(1-percentReinfection) )/
                              (1-percentNonadherent)
  return(RR_progression)
}

# take a probability, convert to odds, apply an odds ratio, and convert back to probability
probwithOR <- function(prob, OR)
{
  odds <- prob/(1-prob)
  newodds <- odds*OR
  newprob <- newodds/(1+newodds)
  return(newprob)
}

# abind along 3rd dimension without needing additional arguments (for use in reduce)
abshort <- function(x,y) {return(abind(x,y,along=3))}

# Will use the same samplemat to run each of multiple scenarios over the sample probabilistic sample of other params.
# So sample the other parameters first, then use the same sample for each scenario.
sample.paramvals <- function(n=10000, set.seed=12345, samplemat)
  # Shift indicates change in the support, to start at low (gamma)
{
  params <- read.params()
  paramvals <- as.list(params$estimate);
  names(paramvals) <- paste0(params$paramname, params$setting, params$subpopulation)
  
  if(missing(samplemat)) 
    samplemat <- randomLHS(n=n, k=nrow(params))
  
  if(length(setdiff(c("gamma","beta","none","shifted gamma"), params$shape))>0) 
    stop(paste0("haven't defined handling of shape(s) ",
         setdiff(c("gamma","beta","none","shifted gamma"), params$shape)))
  
  sampledparams <- do.call("rbind", replicate(n, unlist(paramvals), simplify = FALSE))
  
  for (r in 1:length(paramvals))
  {
    if (params$sd[r] == 0 | params$shape[r]=='none') sampledparams[,r] <- params$estimate[r] else
      {
      if (params$shape[r]=="gamma") sampledparams[,r] <- 
        qgamma(p=samplemat[,r], 
               shape = params$estimate[r]^2/params$sd[r]^2, #mean^2/var
               rate = params$estimate[r]/params$sd[r]^2 ) #mean/var
      if (params$shape[r]=="shifted gamma") sampledparams[,r] <- 
        params$low[r] + 
          qgamma(p=samplemat[,r], 
                shape = (params$estimate[r]-params$low[r])^2/params$sd[r]^2, 
                rate = (params$estimate[r]-params$low[r])/params$sd[r]^2 ) 
      if (params$shape[r]=="beta") sampledparams[,r] <-   
          qbeta(p=samplemat[,r], 
                shape1 = AMS(mean=params$estimate[r], var=params$sd[r]^2), 
                shape2 = BMS(mean=params$estimate[r], var=params$sd[r]^2))   
      }
    }
  return(list("params"=sampledparams, "samplemat"=samplemat))
}



###### main model functions #####  
# Regimen: 4R, 6H, or none
# Load: latent (LTBI, future progression) or subclinical (current undetected active TB)
# pendingresistance: a vector, the distribution of resistance among those starting regimen (DS, HR, RR, MDR). 
  # May use future or subclinical from above, or outputs from a prior run of this function.
TPToutcome <- function(regimen="none", load="latent", pendingresistance=c(1,0,0,0), paramvals)
{
  # subclinical_progression_OR is OR of progression of subclinical TB with TPT relative to placebo
  # and in the populations we're modeling, the risk of progression (of latnet to active, or subclinical to clinical) with placebo is 1.
  latent_progression_prob_4r <- 
    tptEfficacyConversion(OR = paramvals$tptOR_4r, 
                                percentReinfection = paramvals$reinfectionfraction, 
                                percentNonadherent = paramvals$nonadherencefraction)
  latent_progression_prob_6h <- latent_progression_prob_4r + 
    paramvals$tptprog_6h_vs_4r * (1-latent_progression_prob_4r)
  subclinical_progression_prob_4r <- latent_progression_prob_4r + 
    paramvals$rxprog_subclinical_butnot_latent_4r*(1-latent_progression_prob_4r)
  subclinical_progression_prob_6h <- subclinical_progression_prob_4r + 
    paramvals$tptprog_6h_vs_4r * (1-subclinical_progression_prob_4r)
  
  becomes_a_case <- NA
  case_res <- c(NA, NA, NA, NA); names(case_res) <- c("DS", "HR", "RR", "MDR")
  if (regimen=="none")  
    if (load %in% c("latent", "subclinical")) 
    {   becomes_a_case <- pendingresistance # without TPT, all will eventually develop TB...
        case_res <- becomes_a_case  # ...with no new acquired resistance
    }
  if (regimen=="6H")
  {
    if (load=="latent")
    {
      # of those initially susceptible, eff (after TPT efficacy conversion) progress, of whom aRR progress with resistance. 
      becomes_a_case <- pendingresistance * c( latent_progression_prob_6h,
                                              1, 
                                              latent_progression_prob_6h,
                                              1)
      case_res <- becomes_a_case %*% # assuming HR higher risk than RR by some factor
                  t(array(c(1-probwithOR(prob=paramvals$aRR_4R_ltbi_if_progresses, 
                                         OR=paramvals$aHR_vs_aRR_OR),
                            probwithOR(prob=paramvals$aRR_4R_ltbi_if_progresses, 
                                       OR=paramvals$aHR_vs_aRR_OR),
                            0,
                            0,
                          0, 1, 0, 0,
                          0, 0, 1-probwithOR(prob=paramvals$aRR_4R_ltbi_if_progresses, 
                                             OR=paramvals$aHR_vs_aRR_OR), 
                          probwithOR(prob=paramvals$aRR_4R_ltbi_if_progresses, 
                                     OR=paramvals$aHR_vs_aRR_OR),
                          0, 0, 0, 1), dim=c(4,4))) 
    }
    if (load=="subclinical")
    {
      becomes_a_case <- pendingresistance * c(subclinical_progression_prob_6h, 
                                              1, 
                                              subclinical_progression_prob_6h, 
                                              1) # of those who would have progressed, which do (by resistance profile)
      case_res <- (becomes_a_case ) %*%  # and what is the resulting resistance?
        t(array(c(1-probwithOR(prob=paramvals$aRR_fail_subclinical, 
                               OR=paramvals$aHR_vs_aRR_OR), 
                      probwithOR(prob=paramvals$aRR_fail_subclinical, 
                                 OR=paramvals$aHR_vs_aRR_OR), 
                      0, 
                      0, 
                  0, 1, 0, 0,
                  0, 0, 1-probwithOR(prob=paramvals$aRR_fail_subclinical, 
                                     OR=paramvals$aHR_vs_aRR_OR), 
                        probwithOR(prob=paramvals$aRR_fail_subclinical, 
                                   OR=paramvals$aHR_vs_aRR_OR),
                  0, 0, 0, 1), dim=c(4,4))) 
    }
  }
  if (regimen=="4R")
  {
    if (load=="latent")
    {
      becomes_a_case <- pendingresistance * c(latent_progression_prob_4r,
                                              latent_progression_prob_4r,
                                              1,
                                              1)
      case_res <- becomes_a_case %*% 
        t(array(c(1 - paramvals$aRR_4R_ltbi, 0, paramvals$aRR_4R_ltbi, 0,
                  0, 1 - paramvals$aRR_4R_ltbi, 0, paramvals$aRR_4R_ltbi, 
                  0, 0, 1, 0,
                  0, 0, 0, 1), dim=c(4,4))) 
    }
    if (load=="subclinical")
    {
      becomes_a_case <- pendingresistance * c(subclinical_progression_prob_4r, 
                                              subclinical_progression_prob_4r, 
                                              1, 
                                              1) # of those who would have progressed with resistance, which do (by resistance profile)
      case_res <- (becomes_a_case) %*%  
        t(array(c(1-paramvals$aRR_fail_subclinical, 0, paramvals$aRR_fail_subclinical, 0, 
                  0, 1-paramvals$aRR_fail_subclinical, 0, paramvals$aRR_fail_subclinical,
                  0, 0, 1, 0,
                  0, 0, 0, 1), dim=c(4,4))) 
    }
  }
  return(case_res)
}


# also track resistance after a round of treatment: for a vector of cases and resistance, 
 # get vector of fail/retreatment cases and their new resistance after treatment
# will apply to subclinical diagnoses too
treatonce <- function(res, paramvals)
{
  return((res * c(paramvals$rxfail_ds, paramvals$rxfail_hr, 
                  paramvals$rxfail_rr, paramvals$rxfail_rr)) %*% 
    t(array(
      c( 1-paramvals$aRR_rxfailds - paramvals$aHR_rxfailds - paramvals$aMDR_rxfailds, 
         paramvals$aHR_rxfailds, 
         paramvals$aRR_rxfailds, 
         paramvals$aMDR_rxfailds,
       0, 1-paramvals$aMDR_rxfailhr, 0, paramvals$aMDR_rxfailhr,
       0, 0, 1-paramvals$aMDR_rxfailrr, paramvals$aMDR_rxfailrr,
       0, 0, 0, 1), dim=c(4,4))
    ))
}


## Run model and process outputs ####

# 1. sample array of param sets once for all scenarios,
sampled <- sample.paramvals(n = 5000)

# 2. choose cohort and coverage scenario (will run this below)
set.scenario <- function(setting = "_pak", sens=0.9, uptake=0.8)
{
  scenario <- list()
  scenario$setting <- setting
  scenario$screensens <- sens #0.9 #sensitivity for subclinical (asymptomatic, bacteriologically positive) TB, of the screening method used		
   # Will vary by analysis, rather than probabilistically. 0 for symptoms alone. 
  scenario$tptuptake <- uptake #0.8
  if (scenario$setting=="_kzn") scenario$subpopulations <- 
    c("_lt100", "_100to200", "_200to350", "_gt350")
  if (scenario$setting=="_pak") scenario$subpopulations <- 
    c("_lt5yo", "_5to15yo", "_gt15yo")
  return(scenario)
}

# 3. within each scenario, loop over all subpopulations and states, generating outcomes for each 

runmodel <- function(scenario, paramvals, paramnumber=1)
{
  
  results <- list()
  
  interventions <- c("none","4R","6H")
  
  paramvals <- as.list(paramvals)
  
  setting <- scenario$setting
  for (subpopulation in scenario$subpopulations)
  {
    popres <- c(1-(paramvals[[paste0("hrprev",setting)]] + 
                     paramvals[[paste0("rrprev",setting)]] + 
                     paramvals[[paste0("mdrprev",setting)]]), 
                paramvals[[paste0("hrprev",setting)]], 
                paramvals[[paste0("rrprev",setting)]], 
                paramvals[[paste0("mdrprev",setting)]]) 
    future <- paramvals[[paste0("cumulativeprogressionincidence",setting,subpopulation)]] * 
            paramvals[["cumulativeprogression_multiplier"]] * 
                popres*paramvals[[paste0("subpopsize",setting,subpopulation)]]
    subclinical <- future * paramvals[[paste0("subclin_per_progressor",setting,subpopulation)]]*
                    paramvals[["subclin_multiplier"]]
    # whether these are detected and treated, 
     # or eligible to get TPT, will depend on screensens. This is the number of asymptomatic who may or may not be identified.
    
    # outcomes after 100% TPT 
    # (outcomes among latent will be new incident [clinical] cases; 
    # outcomes among subclinical will be progressions to clinical disease;
    # total TB incidence will be clinical cases, including new incident TB and subclinical progressions despite TPT)
    latent_initial <- sapply(X=interventions, function(x) 
      TPToutcome(regimen=x, load = "latent", pendingresistance = future, paramvals=paramvals))
    subclinical_initial <- sapply(X=interventions, function(x) 
      TPToutcome(regimen=x, load = "subclinical", 
                pendingresistance = subclinical, paramvals=paramvals))

    # but really, for screensens>0, those subclinicals who get detected will be treated immediately, so the outcome of the TPT step is the outcome of treatonce 
    # - in other words, subclinicals who got detected would have definitely develop clinical TB if not tested, but may or may not be cured by initial treatment.
    # Incidence, then, includes LTBI activations (not prevented), subclinical progressions (not detected or prevented), and subclinical treatment failures (when detected) 
    # and the rest will have outcomes as per TPT. Latent TB outcomes won't be affected.
    # Screening comes with TPT, so for "none", subclinical_initial_withscreening should be the same as subclinical_initial.
     subclinical_initial_withscreening <- t(scenario$screensens*
                                              (colnames(subclinical_initial)!="none") * # screensens only gets offered to those considered for a TPT regimen (not in the no-TPT scenraio)
      do.call("rbind", replicate(length(interventions), treatonce(subclinical, paramvals), 
                                 simplify=FALSE))) + # and their outcome is that of immediate treatment.
            t(t(subclinical_initial) * (1-scenario$screensens*
                                          (colnames(subclinical_initial)!="none"))) # those who aren't detected have the same outcome they'd have with no screening
                                                                       
    # and really, for tptuptake<1, outcomes will be a mix of screening+TPT, versus no screening or TPT (for both latent and subclinical). 
    # so incidence includes (subclinical not screened, sublinical screened and not detected or prevented, subclinical screened and detected and not cured, 
    # ... latent not screened or treated, and latent screened and treated but not prevented))
    latent_initial_reality <- scenario$tptuptake*latent_initial + 
      (1-scenario$tptuptake)*latent_initial[,"none"]
    #** FOR COMPARISONS WITH "NO TPT", WE WANT "NONE" to result in neither screening nor TPT -- SO... 
    ##  ...need to make sure that xx_initial[,'none'] doesn't include screening. I.e., even if high tptuptake,...
    ##  ..."none" should include no uptake. i.e. subclinical_initial_withscreening[,'none'] shouldn't include screening.
    subclinical_initial_reality <- scenario$tptuptake*subclinical_initial_withscreening + 
      (1-scenario$tptuptake)*subclinical_initial[,"none"]
    
    # everyone then moves through a round of treatment. 
    latent_postreatment <- apply(latent_initial_reality, 2, treatonce, paramvals=paramvals)
    subclinical_postreatment <- apply(subclinical_initial_reality, 2, treatonce, paramvals=paramvals)

    initial <- abind(subclinical_initial_reality, latent_initial_reality, along = 3)
    posttreatment <- abind(subclinical_postreatment, latent_postreatment, along = 3)
    outcomes <- abind(initial, posttreatment, along = 4)
    dimnames(outcomes) <- list("resistance" = c("DS", "HR", "RR", "MDR"), 
                               "regimen"=interventions,  
                               "initialTB"=c("subclinical","latent"), 
                               "casetiming"=c("initial","posttreatment"))
    
    results[[subpopulation]] <- outcomes # could then abind this over another dimension with lapply to get rid of lists 
  }
    
  # 4. collate population and analyze for key outputs:
  results$all <- Reduce("+", results)
  
  all <- melt(results$all) %>% pivot_wider(names_from = resistance)
  all %>% group_by(regimen, initialTB, casetiming) %>% mutate(rrmdr = sum())
  
  all <- all %>% mutate(TB = DS + HR + RR + MDR)
  all <- all %>% mutate(allRR = RR + MDR)
  all <- all %>% mutate(allHR = HR + MDR)
  all$paramnumber <- paramnumber
  
  return(all)
}

# # 5.  run for all param sets and collate outcomes
# scenario <- set.scenario()
# a <- apply(sampled$params,1,function(x) runmodel(scenario = scenario, paramvals = x))
# 
# 
# # 6. compare regimens within scenario


# 7. Do the same for additional scnearios (varying tptuptake, screensens, setting).
# I’ll plan to make 4R versus 6H the primary comparison, 
#  and to show the same comparison under a few different screening scenarios: 
#  4R vs 6H with symptoms only (screensens=0, tptuptake=100%),
#  4R vs 6H with higher sensitivity screening but complete uptake (screensens=0.8, tptuptake=100%), 
#  4R vs 6H with higher sensitivity screening and a reduction in uptake (screensens=0.8, tptuptake=80%). 
# For screening, assume a sensitivity slightly worse than what CXR has for symptomatic undiagnosed cx+ TB in prevalence surveys (84% in annex 2 of https://apps.who.int/iris/bitstream/handle/10665/252424/9789241511506-eng.pdf?sequence=1), or what Xpert has in sympatomic patients (90%).
# And for impact of screening on uptake, look at dropoff in Pakistan from contacts doing clinical evaluation (and not rec'd for Rx) to those starting TPT: 3523/(13103-262)


### ANALYSIS for manuscript #### 


#### Part 1: Empiric data on subclinical TB ####

# First results table: primary data and subclinical TB prevalence ### 
datatab <- as.data.frame(read_excel(paste0("RIF TPT params ",tag,".xlsx"), sheet = 2, col_names = T, col_types = c("text",rep("numeric",8))))
datatab <- datatab[rowSums(is.na(datatab)) != ncol(datatab),]
rownames(datatab) <- datatab$obs; datatab <- datatab[,-c(1,2)]
(displaytab <- datatab[1,])
colnames(displaytab) <- c("CD4 <100", "CD4 100-200", "CD4 200-350", "CD4 >350", "Contact age <5y", "Contact age 5-14y", "Contact age >=15y")
displaytab[2,] <- paste0(datatab[2,]," (",round(datatab[2,]/datatab[1,]*100,1),"%)") 
rownames(displaytab)[1:2] <- c("Total evaluated","Symptomatic TB at baseline") # includes some Xpert- in KZN and some lcinical diagnoses in Pakistan
displaytab[3,] <- paste0(colSums(datatab[3:5,], na.rm=T)," (",round(colSums(datatab[3:5,], na.rm=T)/datatab[1,]*100, 1),"%)") 
rownames(displaytab)[3] <- c("Asymptomatic TB diagnosed through thorough baseline evaluation (Pak) or 3-month follow-up (KZN)") # taken to represent subclinical TB
displaytab[4,] <- datatab[6,]
displaytab[6,] <- datatab[8,]
displaytab[5,] <- paste0(datatab[7,]," (",round(datatab[7,]/datatab[6,]*100,1),"%)")
displaytab[7,] <- paste0(datatab[9,]," (",round(datatab[9,]/datatab[8,]*100,1),"%)")
rownames(displaytab)[4:7] <- c("Followed to 6 months", "TB diagnoses, 3 to 6 months", "Followed to 12 months", "TB diagnoses, 6 to 12 months")
displaytab[8,] <- paste0(round(100*datatab["cumulativeprogressionincidence",],1),"%")
rownames(displaytab)[8] <- "Estimate of cumulative future (>90d) incidence from infections present at enrollment" 
# "*extrapolated from month 4-12 incidence in KZN contacts, and based on published household contact cohorts
displaytab[9,] <- round(datatab["ratio",],1)
rownames(displaytab)[9] <- "Estimated subclinical prevalent cases per future LTBI progression" 
displaytab
write.csv(displaytab, file="Table2_cohortsubclinical.csv")

sum(datatab[3,1:4])/sum(datatab[1,1:4])
sum(datatab[3,1:4])/(sum(datatab[2:3,1:4]))

sum(datatab[4:5,5:7])/sum(datatab[1,5:7])
sum(datatab[4:5,5:7])/(sum(datatab[2,5:7]) + sum(datatab[4:5,5:7]))

# Supplemental table: Comparison of 3 month incidence to cx+ prevalence
cultured <- c(52,81,100,192)
cxpos <- c(18,9,7,9)
notcultured <- c(327,361,685,1200)
tb3notcultured <- c(106,68,70,83)
cxtab <- as.data.frame(t(cultured))
colnames(cxtab) <- c("CD4 <100", "CD4 100-200", "CD4 200-350", "CD4 >350")
cxtab[2,] <- paste0(cxpos," (",round(cxpos/cultured*100),"%)")
cxtab[3,] <- notcultured
cxtab[4,] <- paste0(tb3notcultured," (",round(tb3notcultured/notcultured*100),"%)")
rownames(cxtab) <- c("Included in culture substudy, N", "Culture positive at enrollment, N (%)", 
                     "Not included in culture substudy, N", "Diagnosed with TB within 3 months, not in substudy, N (%)")
cxtab
write.csv(cxtab, file = "tableS1_culture.csv")



#### Part 2: Modeled outcomes  #####

# run <- runmodel(scenario, paramvals)
# runlong <- melt(run, id.vars = c('regimen','initialTB','casetiming','paramnumber'))

getresults <- function(setting = '_kzn', #or '_pak'
                       uptake = 1, # proportion who get TPT (and any recommended screening)
                       sens = 0.9, # sensitivity of screening if it's recommended, or proportion of subclinical TB taht gets detected by screening. 0 means no screening. 
                       sampled)
{
  scenario <- set.scenario(setting = setting, uptake = uptake, sens = sens); 
  
  all <- do.call(rbind, lapply(1:nrow(sampled$params), 
                               FUN = function(x) runmodel(scenario = scenario,
                                                          paramvals = sampled$params[x,],
                                                          paramnumber=x) ))
  alllong <- melt(all, id.vars = c('regimen','initialTB','casetiming','paramnumber'))
  
  
  diffs <- alllong %>% filter(variable %in% c("TB","allRR","HR")) %>% 
    pivot_wider(names_from = regimen, values_from=value) %>%
    mutate(diff4R = `4R` - none, diff6H = `6H` - none) %>% 
    select(-none, -`4R`,-`6H`) %>% 
    pivot_longer(cols = contains('diff'), names_to = 'regimen' ) %>% 
    pivot_wider(names_from=c(initialTB, casetiming), values_from = c(value)) %>% 
    mutate(posttreatment = latent_posttreatment + subclinical_posttreatment) %>% 
    select(-latent_posttreatment, -subclinical_posttreatment) %>% 
    mutate(variable2 = fct_relevel(variable, c("TB","allRR","HR")),
           variable3 = fct_recode(variable2, "Any TB"='TB', "RR/MDR TB"='allRR', "Isoniazid monoresistant TB"='HR')) %>%
    mutate(total_TB = subclinical_initial + latent_initial + posttreatment) %>% 
    pivot_longer(cols=c("subclinical_initial":"latent_initial", "posttreatment", "total_TB"))
  

  return(list("alllong" = alllong, "diffs"=diffs))
  
}

# P90 <- getresults("_pak",uptake=0.9, sens=0.9,  sampled=sampled)
P80 <- getresults("_pak",uptake=0.8, sens=0.9,  sampled=sampled)
K80 <- getresults("_kzn",uptake=0.8, sens=0.9,  sampled=sampled)
# K75 <- getresults("_kzn",uptake=0.75, sens=0.9,  sampled=sampled)
P50 <- getresults("_pak",uptake=0.5, sens=0.9,  sampled=sampled)
K50 <- getresults("_kzn",uptake=0.5, sens=0.9,  sampled=sampled)
P100 <- getresults("_pak",uptake=1, sens=0.9,  sampled=sampled)
K100 <- getresults("_kzn",uptake=1, sens=0.9,  sampled=sampled)
P0 <- getresults("_pak",uptake=1, sens=0,  sampled=sampled) # no screening
K0 <- getresults("_kzn",uptake=1, sens=0,  sampled=sampled) # no screening

save(P80, K80, P50, K50, P100, K100, P0, K0, 
     file="TPTsims.Rdata" )

#expected case numbers, no intervention
baseplotdata <- rbind(
  P80$alllong %>% mutate(Cohort="Household contacts"),
  K80$alllong %>% mutate(Cohort="People with HIV")) %>% 
    filter(regimen=="none", variable %in% c("TB","allRR","HR")) %>% 
  pivot_wider(names_from=c(initialTB, casetiming), values_from = c(value)) %>% 
  mutate(posttreatment = latent_posttreatment + subclinical_posttreatment) %>% 
  select(-latent_posttreatment, -subclinical_posttreatment) %>% 
  mutate(total_TB = subclinical_initial + latent_initial + posttreatment) %>% 
  pivot_longer(cols=c("subclinical_initial":"latent_initial", "posttreatment", "total_TB")) %>%
    mutate(variable2 = fct_drop(variable),
           variable3 = fct_relevel(variable2, c("TB","allRR","HR")),
           variable4 = fct_recode(variable3, "Any TB"='TB', "RR/MDR TB"='allRR', "Isoniazid monoresistant TB"='HR'),
           name2 = fct_relevel(name, c("latent_initial","subclinical_initial",
                                              "posttreatment","total_TB")),
            name3 = fct_recode(name2, "From latent" = "latent_initial",
                                     "From subclinical"="subclinical_initial", 
                                     "After future treatment"="posttreatment",
                                     "All averted cases combined" = "total_TB"))
  
tbbaseplot <- ggplot(baseplotdata %>% 
                       filter(name3 !="All averted cases combined",
                              variable4=="Any TB"), 
  aes(x=name3, y=value, col=name3)) +
  geom_boxplot() + 
  ylab("Future incidence of symptomatic TB, if no intervention, per cohort of 1000") + 
  xlab("") +
  theme_bw() + 
  theme(legend.text = element_text(size = 8),
        legend.margin =margin(0.1,0.1,0.1,2, unit='cm'),
        axis.title = element_text(size=9)) +
  facet_grid(variable4 ~ Cohort,  scales="free_y", space="free_y") + 
  theme(legend.position = 'bottom', legend.title=element_blank(), 
        axis.text.x=element_blank(), axis.title.x=element_blank()) + 
  scale_colour_manual(values=hcl(seq(15,300,length.out=4), 
                                 c = 100 , l = 65)[1:3])
resbaseplot <- ggplot(baseplotdata %>% 
                       filter(name3 !="All averted cases combined",
                              variable4!="Any TB"), 
                     aes(x=name3, y=value, col=name3)) +
  geom_boxplot() + 
  ylab("Future incidence of drug-resistant TB, if no intervention, per cohort of 1000") + 
  xlab("") +
  theme_bw() + 
  theme(legend.text = element_text(size = 8),
        legend.margin =margin(0.1,0.1,0.1,2, unit='cm'),
        axis.title = element_text(size=9)) +
  facet_grid(variable4 ~ Cohort,  scales="free_y", space="free_y") + 
  theme(legend.position = 'bottom', legend.title=element_blank(), 
        axis.text.x=element_blank(), axis.title.x=element_blank()) + 
  scale_colour_manual(values=hcl(seq(15,300,length.out=4), 
                                 c = 100 , l = 65)[1:3])

(Fig2 <- ggarrange( tbbaseplot, resbaseplot, 
             common.legend = T, legend='bottom')
)

pdf(file="TPT model Fig 2 no intervention.pdf", width = 8, height = 6)
Fig2
dev.off()


baseline <- P80$alllong %>% filter(regimen=="none")
Kbaseline <- K80$alllong %>% filter(regimen=="none")

# number of cases
baseline %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  summarize(all=sum(value)) %>% summarize(median(all), quantile(all,.025), quantile(all, 0.975))
# number from latent etc
baseline %>% filter(variable=="TB", initialTB=="latent") %>% group_by(paramnumber) %>%
  summarize(all=sum(value)) %>% summarize(median(all), quantile(all,.025), quantile(all, 0.975))
baseline %>% filter(variable=="TB", initialTB=="subclinical") %>% group_by(paramnumber) %>%
  summarize(all=sum(value)) %>% summarize(median(all), quantile(all,.025), quantile(all, 0.975))
# proportoin from latent 
baseline %>% filter(variable=="TB") %>% group_by(paramnumber) %>%
  summarize(all=weighted.mean(initialTB=="latent" & casetiming=="initial", value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))
baseline %>% filter(variable=="TB") %>% group_by(paramnumber) %>%
  summarize(all=weighted.mean(initialTB=="subclinical" & casetiming=="initial", value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))
baseline %>% filter(variable=="TB") %>% group_by(paramnumber) %>%
  summarize(all=weighted.mean(casetiming=="posttreatment", value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))

Kbaseline %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  summarize(all=sum(value)) %>% summarize(median(all), quantile(all,.025), quantile(all, 0.975))
Kbaseline %>% filter(variable=="TB") %>% group_by(paramnumber) %>%
  summarize(all=weighted.mean(initialTB=="latent" & casetiming=="initial", value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))
Kbaseline %>% filter(variable=="TB") %>% group_by(paramnumber) %>%
  summarize(all=weighted.mean(initialTB=="subclinical" & casetiming=="initial", value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))
Kbaseline %>% filter(variable=="TB") %>% group_by(paramnumber) %>%
  summarize(all=weighted.mean(casetiming=="posttreatment", value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))
# resistance
baseline %>% filter(variable=="allRR") %>% group_by(paramnumber) %>%
  summarize(all=weighted.mean(casetiming=="posttreatment", value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))
baseline %>% filter(variable=="allRR") %>% group_by(paramnumber) %>%
  summarize(all=sum(value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))
Kbaseline %>% filter(variable=="allRR") %>% group_by(paramnumber) %>%
  summarize(all=weighted.mean(casetiming=="posttreatment", value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))
Kbaseline %>% filter(variable=="allRR") %>% group_by(paramnumber) %>%
  summarize(all=sum(value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))
baseline %>% filter(variable=="HR") %>% group_by(paramnumber) %>%
  summarize(all=weighted.mean(casetiming=="posttreatment", value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))
baseline %>% filter(variable=="HR") %>% group_by(paramnumber) %>%
  summarize(all=sum(value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))
Kbaseline %>% filter(variable=="HR") %>% group_by(paramnumber) %>%
  summarize(all=weighted.mean(casetiming=="posttreatment", value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))
Kbaseline %>% filter(variable=="HR") %>% group_by(paramnumber) %>%
  summarize(all=sum(value)) %>% 
  summarize(median(all), quantile(all,.025), quantile(all, 0.975))


# interventions in single setting (put other in supplement by changing first two rows)
setting <- "_pak"
withscreening <- P80$diffs
noscreening <- P0$diffs

# setting <- "_kzn"
# withscreening <- K80$diffs
# noscreening <- K0$diffs

plotdata <-
  rbind(withscreening %>% mutate(screen=TRUE), 
               noscreening %>% mutate(screen=FALSE)) %>% 
  mutate(name2 = fct_relevel(name, c("latent_initial","subclinical_initial",
                                     "posttreatment","total_TB")),
         name3 = fct_recode(name2, "From latent" = "latent_initial",
                            "From subclinical"="subclinical_initial", 
                            "After future treatment"="posttreatment",
                            "All averted cases combined" = "total_TB"),
         name4 = fct_recode(paste0(name3, screen),
                            "From latent" = "From latentFALSE",
                            "From subclinical"="From subclinicalFALSE", 
                            "After future treatment"="After future treatmentFALSE",
                            "All TB combined" = "All averted cases combinedFALSE",
                            " "="From latentTRUE",
                            "  "="From subclinicalTRUE", 
                            "   "="After future treatmentTRUE",
                            "    "="All averted cases combinedTRUE"),
         name5 = fct_relevel(name4,
                             c("From latent", " ",
                               "From subclinical", "  ",
                               "After future treatment", "   ",
                               "All TB combined", "    ")),
         regimen2 = fct_recode(regimen, "4R versus no TPT" ="diff4R",
                               "6H versus no TPT" ="diff6H"))


gres <- ggplot(plotdata %>% filter(variable3 != "Any TB"),
              aes(x=name3, y=value, col=name5)) + 
  geom_hline(yintercept = 0, col='gray') +
  geom_boxplot(position = position_dodge(width = 0.3), outlier.size=0.4, lwd=0.5) + 
  ylab("Change in drug-resistant TB cases") + 
  theme_bw() + 
  theme(legend.text = element_text(size = 8),
        legend.margin =margin(0.1,0.1,0.1,2, unit='cm'),
        axis.title = element_text(size=9)) +
  facet_grid(variable3 ~ regimen2,  scales="free_y", space="free_y") + 
  theme(legend.position = 'bottom', legend.title=element_blank(), 
        axis.text.x=element_blank(), axis.title.x=element_blank()) + 
  scale_colour_manual(values=hcl(rep(seq(15,300,length.out=4),each=2), 
                                 c = 100 , l = 65+ rep(c(-15,15), times=4))) 

gtb <- ggplot(plotdata %>% filter(variable3 == "Any TB"),
               aes(x=name3, y=value, col=name5)) + 
  geom_hline(yintercept = 0, col='gray') +
  geom_boxplot(position = position_dodge(width = 0.3), outlier.size=0.4, lwd=0.5) + 
  ylab("Change in overall TB cases") + 
  theme_bw() + 
  theme(legend.text = element_text(size = 8),
        legend.margin =margin(0.1,0.1,0.1,2, unit='cm'),
        axis.title = element_text(size=9)) +
  facet_grid(variable3 ~ regimen2,  scales="free_y", space="free_y") + 
  theme(legend.position = 'bottom', legend.title=element_blank(), 
        axis.text.x=element_blank(), axis.title.x=element_blank()) + 
  scale_colour_manual(values=hcl(rep(seq(15,300,length.out=4),each=2), 
                                 c = 100 , l = 65+ rep(c(-15,15), times=4))) 

(diffFig <- ggdraw( 
  ggarrange( gtb + labs(tag="A"), gres + labs(tag="B"), 
  common.legend = T, legend='bottom')
  ) +
  draw_label("No screening\n \nScreening + 20% reduced uptake", 
             x = 0.26, y = 0.021, 
             hjust=1, vjust=0, size = 8, fontface = 2)
)

if(setting=="_pak")
{
  pdf(file="TPT model Fig 3 Pakistan diffs.pdf", width = 9, height=6)
  diffFig
  dev.off()
} else
{
  pdf(file="TPT model Fig Sx PWH diffs.pdf", width = 9, height=6)
  diffFig
  dev.off()
}

# numeric outcomes , all TB prevented:
# with 4r:
P0$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
                      summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
                                quantile(number, probs = c(0.5,0.025,0.975)))
# from latent versus subclinical
P0$alllong %>% filter(variable=="TB", initialTB=='latent') %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
P0$alllong %>% filter(variable=="TB", initialTB=='subclinical') %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))

K0$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
# with 6H:
P0$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`6H`))/sum(none), number=sum(none)-sum(`6H`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
K0$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`6H`))/sum(none), number=sum(none)-sum(`6H`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))

# isoniazid monoresistant cases prevented by 4R:
P0$alllong %>% filter(variable=="HR") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
K0$alllong %>% filter(variable=="HR") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
# RR/MDR added
P0$alllong %>% filter(variable=="allRR") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
P0$alllong %>% filter(variable %in% c("RR","allRR")) %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = c("initialTB", "casetiming")) %>% 
  mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment ) %>% 
  select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment )) %>% 
  pivot_longer(all) %>% 
  pivot_wider(names_from = "regimen") %>%
  mutate(diff = `4R` - none) %>% select(-c(none,`4R`,`6H`)) %>% 
  pivot_wider(names_from = "variable", values_from="diff") %>%
  summarize("proportionRifMono" = ifelse(RR>0 & allRR>RR, RR/allRR, NA)) %>% 
  summarize(
    quantile(proportionRifMono, probs = c(0.5,0.025,0.975), na.rm = T))
K0$alllong %>% filter(variable=="allRR") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))

#resistance outcomes with 6H:
# isoniazid monoresistant cases added :
P0$alllong %>% filter(variable=="HR") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`6H`))/sum(none), number=sum(none)-sum(`6H`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
K0$alllong %>% filter(variable=="HR") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`6H`))/sum(none), number=sum(none)-sum(`6H`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
# RR/MDR averted
P0$alllong %>% filter(variable=="allRR") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`6H`))/sum(none), number=sum(none)-sum(`6H`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
K0$alllong %>% filter(variable=="allRR") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`6H`))/sum(none), number=sum(none)-sum(`6H`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))


# 4R versus 6H
# RR added per TB, HR, and MDR averted
P0$alllong %>% filter(variable %in% c("TB","HR", "allRR","MDR")) %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = c("initialTB", "casetiming")) %>% 
  mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment ) %>% 
  select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment )) %>% 
  pivot_longer(all) %>% 
  pivot_wider(names_from = "regimen") %>%
  mutate(diff = `4R` - `6H`) %>% select(-c(none,`4R`,`6H`)) %>% 
  pivot_wider(names_from = "variable", values_from="diff") %>%
  summarize("HmonoperRmono" = ifelse(allRR>0 & HR<0, -HR/allRR, NA),
            "TBperRmono" = ifelse(allRR>0 & TB<0, -TB/allRR, NA),
              "allRRpos" = (allRR>0),
              "HRneg" = (HR<0),
              "TBneg" = (TB<0)) %>% 
  summarize(
    quantile(HmonoperRmono, probs = c(0.5,0.025,0.975), na.rm = T),
          quantile(TBperRmono, probs = c(0.5,0.025,0.975), na.rm = T),
    mean(allRRpos), mean(HRneg), mean(TBneg)
    )

K0$alllong %>% filter(variable %in% c("TB","HR", "allRR","MDR")) %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = c("initialTB", "casetiming")) %>% 
  mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment ) %>% 
  select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment )) %>% 
  pivot_longer(all) %>% 
  pivot_wider(names_from = "regimen") %>%
  mutate(diff = `4R` - `6H`) %>% select(-c(none,`4R`,`6H`)) %>% 
  pivot_wider(names_from = "variable", values_from="diff") %>%
  summarize("HmonoperRmono" = ifelse(allRR>0 & HR<0, -HR/allRR, NA),
            "TBperRmono" = ifelse(allRR>0 & TB<0, -TB/allRR, NA),
            "allRRpos" = (allRR>0),
            "HRneg" = (HR<0),
            "TBneg" = (TB<0)) %>% 
  summarize(
    quantile(HmonoperRmono, probs = c(0.5,0.025,0.975), na.rm = T),
    quantile(TBperRmono, probs = c(0.5,0.025,0.975), na.rm = T),
    mean(allRRpos), mean(HRneg), mean(TBneg))


# Hmono added per TB pr3vented, 6H
P0$alllong %>% filter(variable %in% c("TB","allHR","MDR")) %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = c("initialTB", "casetiming")) %>% 
  mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment ) %>% 
  select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment )) %>% 
  pivot_longer(all) %>% 
  pivot_wider(names_from = "regimen") %>%
  mutate(diff = `6H` - `none`) %>% select(-c(none,`4R`,`6H`)) %>% 
  pivot_wider(names_from = "variable", values_from="diff") %>%
  summarize("TBperHR" = ifelse(allHR>0 & TB<0, -TB/allHR, NA)) %>% 
  summarize(
    quantile(TBperHR, probs = c(0.5,0.025,0.975), na.rm = T))

K0$alllong %>% filter(variable %in% c("TB","allHR","MDR")) %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = c("initialTB", "casetiming")) %>% 
  mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment ) %>% 
  select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment )) %>% 
  pivot_longer(all) %>% 
  pivot_wider(names_from = "regimen") %>%
  mutate(diff = `6H` - `none`) %>% select(-c(none,`4R`,`6H`)) %>% 
  pivot_wider(names_from = "variable", values_from="diff") %>%
  summarize("TBperHR" = ifelse(allHR>0 & TB<0, -TB/allHR, NA)) %>% 
  summarize(
    quantile(TBperHR, probs = c(0.5,0.025,0.975), na.rm = T))


########
# Table: summarize outcome measures incl ratios.
# format like figure
# columns: 4R sxs, 6H sxs, 4R cxr, 6H cxr, ?? also 4R w/ CXR + reduced uptake 80% -> 60%? 
# rows: INHmono created/prevented, Rmono created/.., MDR created, RR/MDR created, total DR created,
#        TB incidence, TB prevented (vs none), TB prevented (vs 6H),
#        TB prevented per RR created (vs 6H), TB prevented per RR created (vs 6H)

Pakistan <- TRUE

if(Pakistan)
{
  P80t <- P80$alllong %>% mutate(Scenario="Subclinical TB screening, 20% reduced uptake")
  P0t <- P0$alllong %>% mutate(Scenario="Symptom screening only")
  P100t <- P100$alllong %>% mutate(Scenario="Subclinical TB screening, same uptake")
  Pt <- rbind(rbind(P0t, P100t, P80t) %>% filter(regimen!='none'),
              P0$alllong %>% filter(regimen=='none') %>% mutate(Scenario="Baseline without TPT"))
} else
{
  P80t <- K80$alllong %>% mutate(Scenario="Subclinical TB screening, 20% reduced uptake")
  P0t <- K0$alllong %>% mutate(Scenario="Symptom screening only")
  P100t <- K100$alllong %>% mutate(Scenario="Subclinical TB screening, same uptake")
  Pt <- rbind(rbind(P0t, P100t, P80t) %>% filter(regimen!='none'),
              K0$alllong %>% filter(regimen=='none') %>% mutate(Scenario="Baseline without TPT"))
}


tempnames <- t(Pt %>% filter(variable=='TB') %>% 
                 pivot_wider(names_from = c(initialTB, casetiming)) %>% 
                 mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment) %>% 
                 select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment)) %>%
                 
                 group_by(regimen, Scenario) %>%
                 summarise(mean(all)))

colorder <- c(1,4,3,2,7,6,5)

resulttable <- 
  rbind(
    
    # TB incidence
    as.numeric(t(
      Pt %>% filter(variable=='TB') %>% 
        pivot_wider(names_from = c(initialTB, casetiming)) %>% 
        mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment) %>% 
        select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment)) %>%
        group_by(regimen, Scenario) %>%
        summarise(median(all)))[3,]),
    
    # TB prevented, vs none
    c(NA,as.numeric(t(
      Pt %>% filter(variable=='TB') %>% 
        pivot_wider(names_from = c(initialTB, casetiming)) %>% 
        mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment) %>% 
        select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment)) %>%
        pivot_wider(names_from = c(regimen, Scenario), values_from=c(all)) %>%
        mutate(across(.cols=contains('4R')|contains('6H'), .fns = ~ `none_Baseline without TPT` - .x)) %>%
        pivot_longer(cols = contains('4R')|contains('6H')|contains('none'),
                     names_to = c("regimen","Scenario"), 
                     names_pattern = "(.*)_(.*)",
                     values_to = "averted") %>% 
        mutate(regimen=factor(regimen, levels=c("none","4R","6H"))) %>%
        group_by(regimen, Scenario) %>%
        summarise(median(averted)))[3,-1])),
    
    
    #Hmono prevented, vs none
    c(NA,as.numeric(t(
      Pt %>% filter(variable=='HR') %>% 
        pivot_wider(names_from = c(initialTB, casetiming)) %>% 
        mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment) %>% 
        select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment)) %>%
        pivot_wider(names_from = c(regimen, Scenario), values_from=c(all)) %>%
        mutate(across(.cols=contains('4R')|contains('6H'), .fns = ~ `none_Baseline without TPT` - .x)) %>%
        pivot_longer(cols = contains('4R')|contains('6H')|contains('none'),
                     names_to = c("regimen","Scenario"), 
                     names_pattern = "(.*)_(.*)",
                     values_to = "averted") %>% 
        filter(regimen != 'none') %>% 
        mutate(regimen=factor(regimen, levels=c("none","4R","6H"))) %>%
        group_by(regimen, Scenario) %>%
        summarise(median(averted)))[3,])),
    
    # R mono revented, vs none
    c(NA,as.numeric(t(
      Pt %>% filter(variable=='RR') %>% 
        pivot_wider(names_from = c(initialTB, casetiming)) %>% 
        mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment) %>% 
        select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment)) %>%
        pivot_wider(names_from = c(regimen, Scenario), values_from=c(all)) %>%
        mutate(across(.cols=contains('4R')|contains('6H'), .fns = ~ `none_Baseline without TPT` - .x)) %>%
        pivot_longer(cols = contains('4R')|contains('6H')|contains('none'),
                     names_to = c("regimen","Scenario"), 
                     names_pattern = "(.*)_(.*)",
                     values_to = "averted") %>% 
        filter(regimen != 'none') %>% 
        mutate(regimen=factor(regimen, levels=c("none","4R","6H"))) %>%
        group_by(regimen, Scenario) %>%
        summarise(median(averted)))[3,])),
    
    # MDR averted vs none
    c(NA,as.numeric(t(
      Pt %>% filter(variable=='MDR') %>% 
        pivot_wider(names_from = c(initialTB, casetiming)) %>% 
        mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment) %>% 
        select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment)) %>%
        pivot_wider(names_from = c(regimen, Scenario), values_from=c(all)) %>%
        mutate(across(.cols=contains('4R')|contains('6H'), .fns = ~ `none_Baseline without TPT` - .x)) %>%
        pivot_longer(cols = contains('4R')|contains('6H')|contains('none'),
                     names_to = c("regimen","Scenario"), 
                     names_pattern = "(.*)_(.*)",
                     values_to = "averted") %>% 
        filter(regimen != 'none') %>% 
        mutate(regimen=factor(regimen, levels=c("none","4R","6H"))) %>%
        group_by(regimen, Scenario) %>%
        summarise(median(averted)))[3,])),
    
    
    # any DR averted vs none
    c(NA,as.numeric(t(
      Pt %>% filter(variable=='HR'|variable=="allRR") %>% 
        pivot_wider(names_from = c(initialTB, casetiming)) %>% 
        mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment) %>% 
        select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment)) %>%
        pivot_wider(names_from = c(regimen, Scenario), values_from=c(all)) %>%
        mutate(across(.cols=contains('4R')|contains('6H'), .fns = ~ `none_Baseline without TPT` - .x)) %>%
        pivot_longer(cols = contains('4R')|contains('6H')|contains('none'),
                     names_to = c("regimen","Scenario"), 
                     names_pattern = "(.*)_(.*)",
                     values_to = "averted") %>% 
        filter(regimen != 'none') %>% 
        mutate(regimen=factor(regimen, levels=c("none","4R","6H"))) %>%
        group_by(regimen, Scenario) %>%
        summarise(median(averted)))[3,])),                 
    
    
    
    
    # # TB prevented, vs 6H
    # c(NA,as.numeric(t(
    #   Pt %>% filter(variable=='TB', regimen %in% c("4R","6H")) %>% 
    #     pivot_wider(names_from = c(initialTB, casetiming)) %>% 
    #     mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment) %>% 
    #     select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment)) %>%
    #     pivot_wider(names_from = c(regimen), values_from=c(all)) %>%
    #     mutate(difference = `6H` - `4R`) %>%
    #     group_by(Scenario) %>%
    #     summarise(median(difference)))[2,]),NA,NA,NA),
    
    # TB prevented per RR created (vs none)
    c(NA,
      as.numeric(t(
        Pt %>% filter(variable %in% c('TB','allRR'), regimen!='6H') %>% 
          pivot_wider(names_from = c(initialTB, casetiming)) %>% 
          mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment) %>% 
          select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment)) %>%
          pivot_wider(names_from = c(regimen, Scenario), values_from=c(all)) %>%
          mutate(across(.cols=contains('4R'), .fns = ~ `none_Baseline without TPT` - .x)) %>%
          pivot_longer(cols = contains('4R')|contains('none'),
                       names_to = c("regimen","Scenario"), 
                       names_pattern = "(.*)_(.*)",
                       values_to = "averted") %>% 
          filter(regimen=='4R') %>% 
          pivot_wider(names_from=variable, values_from=averted) %>%
          group_by(regimen, Scenario) %>%
          mutate(RRperTB = ifelse(allRR<0,-TB/allRR,NA)) %>% 
          summarise(ifelse(mean(is.na(RRperTB))<0.1, median(RRperTB, na.rm=T), -99), mean(is.na(RRperTB))))[3,]),-99,-99,-99),
    # no RR added 4% of the time for Pak symptoms only. others >50% none added. 
    
    
    # H mono prevented per RR created (vs no TPT)
    c(NA,
      as.numeric(t(
        Pt %>% filter(variable %in% c('HR','allRR'), regimen!='6H') %>% 
          pivot_wider(names_from = c(initialTB, casetiming)) %>% 
          mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment) %>% 
          select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment)) %>%
          pivot_wider(names_from = c(regimen, Scenario), values_from=c(all)) %>%
          mutate(across(.cols=contains('4R'), .fns = ~ `none_Baseline without TPT` - .x)) %>%
          pivot_longer(cols = contains('4R')|contains('none'),
                       names_to = c("regimen","Scenario"), 
                       names_pattern = "(.*)_(.*)",
                       values_to = "averted") %>% 
          filter(regimen=='4R') %>% 
          pivot_wider(names_from=variable, values_from=averted) %>%
          group_by(regimen, Scenario) %>%
          mutate(RRperHR = ifelse(allRR<0,-HR/allRR,NA)) %>% 
          summarise(ifelse(mean(is.na(RRperHR))<0.2, median(RRperHR, na.rm=T), -99), mean(is.na(RRperHR))))[3,]),-99,-99,-99),
    
    # TB prevented per RR created (vs 6H)
    c(NA,
      as.numeric(t(
        Pt %>% filter(variable %in% c('TB','allRR'), regimen!='none') %>% 
          pivot_wider(names_from = c(initialTB, casetiming)) %>% 
          mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment) %>% 
          select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment)) %>%
          pivot_wider(names_from = c(regimen), values_from=c(all)) %>%
          mutate(difference = `6H` - `4R`) %>%
          select(-c(`4R`,`6H`)) %>%
          pivot_wider(names_from=variable, values_from=difference) %>%
          mutate(RRperTB = ifelse(allRR<0,-TB/allRR,NA)) %>% 
          group_by(Scenario) %>%
          summarise(ifelse(mean(is.na(RRperTB))<0.2, median(RRperTB, na.rm=T), -99), mean(is.na(RRperTB))))[2,]),NA,NA,NA),
    # no RR added 2% of the time for Pak **. others 11% with none added. 
    
    
    # H mono prevented per RR created (vs 6H)
    c(NA,
      as.numeric(t(
        Pt %>% filter(variable %in% c('HR','allRR'), regimen!='none') %>% 
          pivot_wider(names_from = c(initialTB, casetiming)) %>% 
          mutate(all = subclinical_initial + latent_initial + subclinical_posttreatment + latent_posttreatment) %>% 
          select(-c(subclinical_initial, latent_initial, subclinical_posttreatment, latent_posttreatment)) %>%
          pivot_wider(names_from = c(regimen), values_from=c(all)) %>%
          mutate(difference = `6H` - `4R`) %>%
          select(-c(`4R`,`6H`)) %>%
          pivot_wider(names_from=variable, values_from=difference) %>%
          mutate(RRperHR = ifelse(allRR<0,-HR/allRR,NA)) %>% 
          group_by(Scenario) %>%
          summarise(ifelse(mean(is.na(RRperHR))<0.2, median(RRperHR, na.rm=T), -99), mean(is.na(RRperHR))))[2,]),NA,NA,NA)
    # no RR added 2% of the time for Pak **. others 11% with none added. 
  )[,colorder]


rownames(resulttable) <- c(
  "Incident TB  cases (absolute)",
  "Symptomatic TB cases averted (vs no TPT)",
  "INH monoresistant cases added (versus no TPT)", 
  "RIF monoresistant cases added (versus no TPT)", 
  "MDR cases added (versus no TPT)", 
  "Total DR cases added (versus no TPT)", 
  "Symptomatic TB cases averted per RR added (vs no TPT)",
  "INH monoresistance prevented per RR added (vs no TPT)",
  "Symptomatic TB cases averted per RR added (vs 6H)",
  "INH monoresistance prevented per RR added (vs 6H)"
)

colnames(resulttable) <-
  paste0(ifelse(tempnames[1,]=='none','',paste0(tempnames[1,],', ')),tempnames[2,])[colorder]

roundedresulttable <- t(as.data.frame(round(signif(t(resulttable),digits = 2),2)) %>% mutate_all(as.character))
roundedresulttable[roundedresulttable=="-99"]<-"†"

colnames(roundedresulttable) <- colnames(resulttable)

roundedresulttable

if(Pakistan) write.csv(roundedresulttable, file=paste0("Table3_modeloutcomes","_pak",".csv")) else
  write.csv(roundedresulttable, file=paste0("Table3_modeloutcomes","_kzn",".csv"))


##### Part 3: symptom-only versus more sensitive screening, still uptake = 1-noncompletion ###########
P50$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
P0$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
P80$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
# P90$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
#   pivot_wider(names_from = "regimen") %>%
#   summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
#   summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
#             quantile(number, probs = c(0.5,0.025,0.975)))

P0$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`6H`))/sum(none), number=sum(none)-sum(`6H`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
P80$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`6H`))/sum(none), number=sum(none)-sum(`6H`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
# P90$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
#   pivot_wider(names_from = "regimen") %>%
#   summarize(proportion=(sum(none)-sum(`6H`))/sum(none), number=sum(none)-sum(`6H`)) %>% 
#   summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
#             quantile(number, probs = c(0.5,0.025,0.975)))


K0$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
K80$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`4R`))/sum(none), number=sum(none)-sum(`4R`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))

K0$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`6H`))/sum(none), number=sum(none)-sum(`6H`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
K80$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
  pivot_wider(names_from = "regimen") %>%
  summarize(proportion=(sum(none)-sum(`6H`))/sum(none), number=sum(none)-sum(`6H`)) %>% 
  summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
            quantile(number, probs = c(0.5,0.025,0.975)))
# K75$alllong %>% filter(variable=="TB") %>% group_by(paramnumber) %>% 
#   pivot_wider(names_from = "regimen") %>%
#   summarize(proportion=(sum(none)-sum(`6H`))/sum(none), number=sum(none)-sum(`6H`)) %>% 
#   summarize(quantile(proportion, probs = c(0.5,0.025,0.975)),
#             quantile(number, probs = c(0.5,0.025,0.975)))
# 
# 
# 
P50$alllong %>% filter(regimen=="4R") %>% 
  filter(variable=="TB") %>% group_by(paramnumber) %>% 
  summarize(all=sum(value)) %>% summarize(median(all), quantile(all,.025), quantile(all, 0.975))
P0$alllong %>% filter(regimen=="4R") %>% 
  filter(variable=="TB") %>% group_by(paramnumber) %>% 
  summarize(all=sum(value)) %>% summarize(median(all), quantile(all,.025), quantile(all, 0.975))
P0$alllong %>% filter(regimen=="none") %>% 
  filter(variable=="TB") %>% group_by(paramnumber) %>% 
  summarize(all=sum(value)) %>% summarize(median(all), quantile(all,.025), quantile(all, 0.975))
P80$alllong %>% filter(regimen=="4R") %>% 
  filter(variable=="TB") %>% group_by(paramnumber) %>% 
  summarize(all=sum(value)) %>% summarize(median(all), quantile(all,.025), quantile(all, 0.975))



# TB per RR or HR, summarized
P0$diffs %>% filter(variable %in% c('TB','allRR'), regimen=='diff4R', name=='total_TB') %>% 
  select(-c(variable2, variable3, name, regimen)) %>% 
  pivot_wider(names_from=variable) %>%
  mutate(RRperTB = ifelse(allRR>0,-TB/allRR,NA)) %>% 
  summarise(quantile(RRperTB, c(0.025,0.5,0.975), na.rm=T))
K0$diffs %>% filter(variable %in% c('TB','allRR'), regimen=='diff4R', name=='total_TB') %>% 
  select(-c(variable2, variable3, name, regimen)) %>% 
  pivot_wider(names_from=variable) %>%
  mutate(RRperTB = ifelse(allRR>0,-TB/allRR,NA)) %>% 
  summarise(quantile(RRperTB, c(0.025,0.5,0.975), na.rm=T))
P0$diffs %>% filter(variable %in% c('TB','HR'), regimen=='diff6H', name=='total_TB') %>% 
  select(-c(variable2, variable3, name, regimen)) %>% 
  pivot_wider(names_from=variable) %>%
  mutate(HRperTB = ifelse(HR>0,-TB/HR,NA)) %>% 
  summarise(quantile(HRperTB, c(0.025,0.5,0.975), na.rm=T))
K0$diffs %>% filter(variable %in% c('TB','HR'), regimen=='diff6H', name=='total_TB') %>% 
  select(-c(variable2, variable3, name, regimen)) %>% 
  pivot_wider(names_from=variable) %>%
  mutate(HRperTB = ifelse(HR>0,-TB/HR,NA)) %>% 
  summarise(quantile(HRperTB, c(0.025,0.5,0.975), na.rm=T))


#### Effects of reduced uptake #########

# For various tptuptake, generate estimates of TB inc and RR-TB inc, 
# with 4R and screensens=0.9

# want to consider uptakes from 20% to 100% (relative to the 100% of sx only screening)

uptakes <- seq(0.9,0.2,by=-0.1)

setting <- "_pak"
allout <- P100$diffs %>% mutate(reduction=0)

setting <- "_kzn"
allout <- K100$diffs %>% mutate(reduction=0)

for(uptake in uptakes)
{
  sout <- getresults(setting = setting,uptake=uptake, sens=0.9,  sampled=sampled)$diffs %>% 
    mutate(reduction=1-uptake)
  allout <- rbind(allout,sout)
}

save(uptake, file=paste0("TPT varying uptake",setting,".Rdata"))
# load(file=paste0("TPT varying uptake",setting,".Rdata"), verbose = T)


# Figure, varying uptake
fig4data <-
  allout %>% 
  filter(name=="total_TB") %>%
  mutate(regimen2 = fct_recode(regimen, "4R versus no TPT" ="diff4R",
                               "6H versus no TPT" ="diff6H"),
         variable2 = fct_recode(variable, "Any TB" = "TB",
                               "Isoniazid monoresistant TB" = "HR",
                               "MDR/RR TB" = "allRR")) 

if (setting=="_pak") tempout <- P0$diffs else tempout <- K0$diffs
noout <- tempout %>% 
  filter(name=="total_TB") %>%
  mutate(regimen2 = fct_recode(regimen, "4R versus no TPT" ="diff4R",
                               "6H versus no TPT" ="diff6H"),
         variable2 = fct_recode(variable, "Any TB" = "TB",
                                "Isoniazid monoresistant TB" = "HR",
                                "MDR/RR TB" = "allRR"),
         reduction = NA)

nout_summary <- noout %>% group_by(regimen2, variable2) %>% 
  summarize(med=median(value))

grangetb <- ggplot(fig4data %>% filter(variable =="TB"),
               aes(x=factor(paste0(100*reduction,"%")), y=value)) + 
  geom_hline(yintercept = 0, col='gray') +
  geom_boxplot(data=noout %>% filter(variable =="TB"),
               outlier.size=0.4, lwd=0.5,
               # col="black", fill="lightgray",
               col=hcl(rep(seq(15,300,length.out=4),each=2),
                       c = 100 , l = 65+ rep(c(-15,15), times=4))[7],
               fill = 'grey94',
                 # hcl(rep(seq(15,300,length.out=4),each=2),
                 #   c = 100 , l = 65+ rep(c(-10,10), times=4), alpha=0.2)[7],
               aes(x="40%"), width=9) + 
  geom_boxplot(outlier.size=0.4, lwd=0.5,
               col=hcl(rep(seq(15,300,length.out=4),each=2), 
                       c = 100 , l = 65+ rep(c(-15,15), times=4))[8]) + 
  ylab("Change in TB cases") + 
  xlab("Reduction in uptake with subclinical screening requirement") +
  theme_bw() + 
  facet_grid(variable2 ~ regimen2,  scales="free_y", space="free_y") +
  geom_text(data = nout_summary %>% filter(variable2=="Any TB"), aes(y = med, x="60%"), 
             label = "Symptom screening only", 
            vjust=0, nudge_y = 0.3, hjust=0.5, size=3, 
            col='black')

ylim <- ifelse(setting=="_pak",4,7)

grangeres <- ggplot(fig4data %>% filter(variable %in% c("allRR","HR")),
                    aes(x=factor(paste0(100*reduction,"%")), y=value)) + 
  geom_hline(yintercept = 0, col='gray') +
  geom_boxplot(data=noout %>% filter(variable  %in% c("allRR","HR")),
               outlier.size=0.4, lwd=0.5,
               col=hcl(rep(seq(15,300,length.out=4),each=2),
                       c = 100 , l = 65+ rep(c(-15,15), times=4))[7],
               fill = 'grey94',
               aes(x="40%"), width=9) + 
  geom_boxplot(outlier.size=0.4, lwd=0.5,
               col=hcl(rep(seq(15,300,length.out=4),each=2), 
                       c = 100 , l = 65+ rep(c(-10,10), times=4))[8]) + 
  ylab("Change in drug-resistant TB cases") + 
  xlab("Reduction in uptake with subclinical screening requirement") +
  theme_bw() + 
  facet_grid(variable2 ~ regimen2) +
  ylim(-ylim,ylim) 


pdf(file = paste0('uptakesensis',setting,'.pdf'), width = 8, height = 9)
ggarrange( grangetb + labs(tag="A"), grangeres + labs(tag="B"), ncol = 1,
           common.legend = T, legend='bottom')
dev.off()



### Sensitivity analyses ####

# choose params to vary
params <- read.params()
targetsetting <- "_pak"
targetsetting <- "_kzn"
varyparams <- params$sensis=='y' & params$setting %in% c(targetsetting,"")
varynames <- paste0(unlist(params[varyparams,1]), unlist(params[varyparams,2]))
names(varynames) <- unlist(params[varyparams,'nicename'])

outcomenames <- c(paste0("Change in ",c("all TB", #1
                                        "RR/MDR TB"),",\nsymptom-only screening,\n4R vs no TPT"), #2
                  "Change in RR/MDR TB,\nsubclinical TB screening,\n4R vs no TPT",
                  "TB prevented per\nRR/MDR added,\nsymptom-only screening,\n4R vs no TPT",
                  "TB prevented per\nRR/MDR added,\nsymptom-only screening,\n4R vs 6H")

nicenames <- rev(names(varynames))

# identify outcomes of interest
sensis <- function(paramname, nicename, diffs, diffs2)
{
  hiindices <- which(sampled$params[,paramname] >= quantile(sampled$params[,paramname],0.9))
  loindices <- which(sampled$params[,paramname] <= quantile(sampled$params[,paramname],0.1))
  hidiffs <- diffs %>% filter(paramnumber %in% hiindices) %>% filter(variable %in% c('TB','allRR'), name=="total_TB") %>%
    select(-c(variable2, variable3))
  lodiffs <- diffs %>% filter(paramnumber %in% loindices) %>% filter(variable %in% c('TB','allRR'), name=="total_TB") %>%
    select(-c(variable2, variable3))
  hidiffs2 <- diffs2 %>% filter(paramnumber %in% hiindices) %>% filter(variable %in% c('TB','allRR'), name=="total_TB") %>%
    select(-c(variable2, variable3))
  lodiffs2 <- diffs2 %>% filter(paramnumber %in% loindices) %>% filter(variable %in% c('TB','allRR'), name=="total_TB") %>%
    select(-c(variable2, variable3))
  outcomes <- 
  rbind(
    # TB and RR prevented
    hidiffs %>% 
        filter(regimen=='diff4R') %>% select(variable, value) %>% 
      mutate(level="high", param=paramname,
             outcome = ifelse(variable=="TB",outcomenames[1],outcomenames[2])),
    lodiffs %>% 
      filter(regimen=='diff4R') %>% select(variable, value) %>% 
      mutate(level="low", param=paramname,
             outcome = ifelse(variable=="TB",outcomenames[1],outcomenames[2])),
             
    
    hidiffs2 %>% 
      filter(regimen=='diff4R', variable=="allRR") %>% select(variable, value) %>% 
      mutate(level="high", param=paramname,
             outcome = outcomenames[3]),
    lodiffs2 %>% 
      filter(regimen=='diff4R', variable=="allRR") %>% select(variable, value) %>% 
      mutate(level="low", param=paramname, 
             outcome = outcomenames[3]),
    
    
    hidiffs  %>%
      pivot_wider(names_from = c(variable), values_from=c(value)) %>%
      filter(regimen=='diff4R') %>% 
      mutate(RRperTB = ifelse(allRR>0,-TB/allRR,NA)) %>% 
      select(-c(TB,allRR)) %>% 
      pivot_longer(cols = RRperTB, names_to="variable", values_to="value") %>%
      select(variable, value) %>% 
      mutate(level="high", param=paramname,
             outcome = outcomenames[4]),
    lodiffs  %>%
      pivot_wider(names_from = c(variable), values_from=c(value)) %>%
      filter(regimen=='diff4R') %>% 
      mutate(RRperTB = ifelse(allRR>0,-TB/allRR,NA)) %>% 
      select(-c(TB,allRR)) %>% 
      pivot_longer(cols = RRperTB, names_to="variable", values_to="value") %>%
      select(variable, value) %>% 
      mutate(level="low", param=paramname,
             outcome = outcomenames[4]),
    
    hidiffs  %>%
      pivot_wider(names_from = c(regimen), values_from=c(value)) %>%
      mutate(value = diff4R - diff6H) %>% select(-c(diff4R, diff6H)) %>% 
      pivot_wider(names_from = c(variable), values_from=c(value)) %>%
      mutate(RRperTBvs6H =  case_when(TB<0 & allRR>0 ~ -TB/allRR,
                                      TB>0 & allRR>0 ~ 0)) %>%  
      select(-c(TB,allRR)) %>% 
      pivot_longer(cols = RRperTBvs6H, names_to="variable", values_to="value") %>%
      select(variable, value) %>% 
      mutate(level="high", param=paramname,
             outcome = outcomenames[5]),

    lodiffs  %>%
      pivot_wider(names_from = c(regimen), values_from=c(value)) %>%
      mutate(value = diff4R - diff6H) %>% select(-c(diff4R, diff6H)) %>% 
      pivot_wider(names_from = c(variable), values_from=c(value)) %>%
      mutate(RRperTBvs6H =  case_when(TB<0 & allRR>0 ~ -TB/allRR,
                                      TB>0 & allRR>0 ~ 0)) %>% 
      select(-c(TB,allRR)) %>% 
      pivot_longer(cols = RRperTBvs6H, names_to="variable", values_to="value") %>%
      select(variable, value) %>% 
      mutate(level="low", param=paramname,
             outcome = outcomenames[5])
  )
  outcomes$nicename = nicename
  return(outcomes)
}

if(targetsetting=="_kzn") diffs <- K0$diffs else diffs <- P0$diffs
if(targetsetting=="_kzn") diffs2 <- K80$diffs else diffs2 <- P80$diffs

alloutcomes <- do.call(rbind, lapply(X = 1:length(varynames), 
                                     FUN = function(x) sensis(paramname=varynames[x], 
                                                                nicename=names(varynames)[x],
                                                                diffs=diffs, diffs2=diffs2)))

# # now plot comparing high and low deiles
# ggplot(data = alloutcomes %>% 
#          filter(value<500,
#            variable!="RRperTBvs6H" | value<20) %>% 
#          mutate(outcome = fct_relevel(outcome, outcomenames)) ,
#        aes(x=nicename, col=level, y=value)) + 
#   geom_boxplot(outlier.colour = NA) + 
#   facet_grid(cols = vars(outcome), scales="free_x") + 
#     coord_flip()  

# or just compare the medians of these deciles:
plotdata <- alloutcomes %>% 
  mutate(outcome2 = fct_relevel(outcome, outcomenames),
         nicename2 = fct_relevel(nicename, nicenames)) %>%
  group_by(level, nicename2, outcome2) %>% 
  summarise(median = median(value, na.rm=T),
            q25 = quantile(value, probs = 0.25, na.rm=T),
            q75 = quantile(value, probs = 0.75, na.rm=T)) %>%
  pivot_wider(names_from = level, values_from = c(median, q25, q75))
# %>%
#   mutate(high2 = case_when(median_high<500 ~ median_high,
#                            median_high >=500 ~ 500),
#          low2 = case_when(median_low<500 ~ median_low,
#                           median_low >=500 ~ 500),
#          makearrow = case_when(median_high>500 ~ 1))
         


(gsensis <- ggplot(data = plotdata) +
  facet_grid(cols = vars(outcome2), scales="free_x") +
  geom_segment(aes(x=nicename2, xend=nicename2, y=median_low, yend=median_high),
               size=3, col='gray50') +
  geom_segment(aes(x=nicename2, xend=nicename2, y=q25_low, yend=q75_low), size=0.3,
               col="red3", arrow = arrow(angle = 90, ends = 'both', length = unit(0.1,"cm"))) +
  geom_segment(aes(x=nicename2, xend=nicename2, y=q25_high, yend=q75_high), size=0.3,
                 col="blue3", arrow = arrow(angle = 90, ends = 'both', length = unit(0.1,"cm"))) +
  # geom_point(aes(x=nicename, y=low2), col="blue") +
  # geom_point(aes(x=nicename, y=high), col="red") +
  # geom_segment(data = plotdata %>% filter(median_high>500),
  #              aes(x=nicename2, xend=nicename2, y=high2, yend=high2),
  #               arrow = arrow(angle = 45, ends = "last")) +
  xlab("") + ylab("IQRs for highest (blue) versus lowest (red) deciles of parameter, and difference in deciles' medians (gray bar)") +
  coord_flip() + 
    theme_bw() +
    ggtitle(ifelse(targetsetting=='_pak', "Household contact cohort", "PWH Cohort"))
)

pdf(paste0("TPT Fig Sx sensis",targetsetting,".pdf"), width = 13, height = 7)
gsensis
dev.off()


