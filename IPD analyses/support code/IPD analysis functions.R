# helper code for IPD analysis

###################################################
# one-stage meta-analysis

# ma.type, "RE" for random or "FE" for fixed effect MA

one_stage_MA <- function(outcome,trial,arm,measure="RR",ma.type="RE",baseline=NULL){
  
  require(lme4)
  
  if (measure=="OR"){
    if (ma.type=="FE")
      ma <- glm(outcome ~ factor(trial) + factor(arm),family=binomial)
    if (ma.type=="RE")
      ma <- glmer(outcome ~ factor(trial) + factor(arm) + (arm-1|trial),family=binomial)
    if (ma.type=="RE_both")
      ma <- glmer(outcome ~ factor(arm) + (1+arm|trial),family=binomial)
    
  }
  
  if (measure=="RR"){
    if (ma.type=="FE")
      ma <- glm(outcome ~ factor(trial) + factor(arm),family=binomial(link="log"))
    if (ma.type=="RE")
      ma <- glmer(outcome ~ factor(trial) + factor(arm) + (arm-1|trial),family=binomial(link="log"))
    if (ma.type=="RE_both")
      ma <- glmer(outcome ~ factor(arm) + (1+arm|trial),family=binomial(link="log"))
  }
  
  if (measure=="MD"){
    if (ma.type=="FE")
      ma <- lm(outcome ~ factor(trial) + factor(arm))
    if (ma.type=="RE")
      ma <- lmer(outcome ~ factor(trial) + factor(arm) + (arm-1|trial))
    if (ma.type=="RE_both")
      ma <- lmer(outcome ~ factor(arm) + (1+arm|trial))
  }
  
  if (measure=="anova"){
    if (ma.type=="FE")
      ma <- lm(outcome ~ factor(trial) + baseline + factor(arm))
    if (ma.type=="RE")
      ma <- lmer(outcome ~ factor(trial) + baseline + factor(arm) + (arm-1|trial))
    if (ma.type=="RE_both")
      ma <- lmer(outcome ~ baseline + factor(arm) + (1+arm|trial))
  }
  
  
  return(summary(ma))
  
}

# for continuous outcomes

one_stage_cont_MA <- function(baseline,final,study,arm,method="ancova",scale="MD",ma.type="RE"){
  
  require(lme4)
  change <- final - baseline
  
  # convert data
  if (scale=="MD"){
    data.ma <- data.frame(study,arm,baseline,final,change)
  }
  
  if (scale=="SMD"){
    # scaling by sd (testing)
    studysize <- tapply(study,study,length)
    Nt <- length(studysize)
    sd.base <- tapply(baseline,factor(study),sd,na.rm=T)
    scale.factor <- rep(sd.base,studysize)
    data.ma <- data.frame(study,arm,baseline=baseline/scale.factor,final=final/scale.factor,change=change/scale.factor)
  }
  
  data.ma <- subset(data.ma,!is.na(baseline) & !is.na(final))
  
  # analysis models
  if (method=="change"){
    if (ma.type=="RE")
      ma <- lmer(change ~ factor(study) + factor(arm) + (arm-1|study), data=data.ma)
    if (ma.type=="RE_both")
      ma <- lmer(change ~ factor(arm) + (1+arm|study), data=data.ma)
  }
  
  if (method=="final"){
    if (ma.type=="RE")
      ma <- lmer(final ~ factor(study)+ factor(arm) + (arm-1|study), data=data.ma)
    if (ma.type=="RE_both")
      ma <- lmer(final ~ factor(arm) + (1+arm|study), data=data.ma)
  }
  
  if (method=="ancova"){
    if (ma.type=="RE")
      ma <- lmer(final ~ factor(study) + baseline + factor(arm) + (arm+baseline-1|study), data=data.ma)
    if (ma.type=="RE_both")
      ma <- lmer(final ~ baseline + factor(arm) + (1+baseline+arm|study), data=data.ma)
  }
  
  return(ma)
  
}



###############################################################
# one-stage meta-analysis with treat-covariate interaction

# ma.type, "RE"/"RE_both" for random or "FE" for fixed effect MA

one_stage_MA_inter <- function(outcome,trial,arm,covar,measure="RR",ma.type="RE",centre=T,baseline=NULL){
  
  require(lme4)
  
  if(centre==T)
    covar <- covar - mean(covar,na.rm=T)
  
  if (measure=="OR"){
    if (ma.type=="FE")
      ma <- glm(outcome ~ factor(trial) + factor(arm)*covar,family=binomial)
    if (ma.type=="RE")
      ma <- glmer(outcome ~ factor(trial) + factor(arm)*covar + (arm-1|trial),family=binomial)
    if (ma.type=="RE_both")
      ma <- glmer(outcome ~ factor(arm)*covar + (1+arm|trial),family=binomial)
  }
  
  if (measure=="RR"){
    if (ma.type=="FE")
      ma <- glm(outcome ~ factor(trial) + factor(arm)*covar,family=binomial(link="log"))
    if (ma.type=="RE")
      ma <- glmer(outcome ~ factor(trial) + factor(arm)*covar + (arm-1|trial),family=binomial(link="log"))
    if (ma.type=="RE_both")
      ma <- glmer(outcome ~ factor(arm)*covar + (1+arm|trial),family=binomial(link="log"))
  }
  
  if (measure=="MD"){
    if (ma.type=="FE")
      ma <- lm(outcome ~ factor(trial) + factor(arm)*covar)
    if (ma.type=="RE")
      ma <- lmer(outcome ~ factor(trial) + factor(arm)*covar + (arm-1|trial))
    if (ma.type=="RE_both")
      ma <- lmer(outcome ~ factor(arm)*covar + (1+arm|trial))
  }
  
  if (measure=="anova"){
    if (ma.type=="FE")
      ma <- lm(outcome ~ factor(trial) + baseline + factor(arm)*covar)
    if (ma.type=="RE")
      ma <- lmer(outcome ~ factor(trial) + baseline + factor(arm)*covar + (arm-1|trial))
    if (ma.type=="RE_both")
      ma <- lmer(outcome ~ baseline + factor(arm)*covar + (1+arm|trial))
  }
  
  return(summary(ma))
  
}


# repeated measures (no covar)

rep_measures_MA <- function(outcome,treat,trial,time,IDVar,
                             model="linear",measure="OR",baseline=NULL){
  
  require(lme4)

  if (measure=="OR"){
    if (model=="linear")
      res <- glmer(outcome ~  treat*time + (1|IDVar/trial) + (1+treat|trial), family=binomial)
    if (model!="linear")
      res <- glmer(outcome ~  treat*factor(time) -treat + (1|IDVar/trial) + (1+treat|trial), family=binomial)
  }
  
  if (measure=="MD"){
    if (model=="linear")
      res <- lmer(outcome ~  treat*time + (1|IDVar/trial) + (1+treat|trial))
    if (model!="linear")
      res <- lmer(outcome ~  treat*factor(time) -treat + (1|IDVar/trial) + (1+treat|trial))
  }
  
  if (measure=="anova"){
    if (model=="linear")
      res <- lmer(outcome ~ baseline + treat*time + (1|IDVar/trial) + (1+treat|trial))
    if (model!="linear")
      res <- lmer(outcome ~ baseline + treat*factor(time) -treat + (1|IDVar/trial) + (1+treat|trial))
  }
  
  return(summary(res))
  
}

##########################################
# convert to odds ratios

cont_convert <- function(result){
  
  coeffs <- result$coefficients
  out <- data.frame(Estimate=coeffs[,1],
                    CIlow=coeffs[,1]-1.96*coeffs[,2],
                    CIhigh=coeffs[,1]+1.96*coeffs[,2])
  out <- round(out,3)
  return(out)
}

log_convert <- function(result){
  
  coeffs <- result$coefficients
  out <- data.frame(Estimate=exp(coeffs[,1]),
                    CIlow=exp(coeffs[,1]-1.96*coeffs[,2]),
                    CIhigh=exp(coeffs[,1]+1.96*coeffs[,2]))
  out <- round(out,3)
  return(out)
}
