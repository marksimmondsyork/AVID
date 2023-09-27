# Two-stage MAoI

two_stage_MAoI <- function(outcome,trial,arm,covar,measure="RR",centre=T,filename=NULL){
  
  # in one trial
  trial_model <- function(outcome,arm,covar,measure,centre){
    if(centre==T)
      covar <- covar - mean(covar,na.rm=T)
    
    if (measure=="OR"){
      trial.res <- glm(outcome ~ factor(arm)*covar,family=binomial)
    }
    if (measure=="RR"){
      trial.res <- glm(outcome ~ factor(arm)*covar,family=binomial(link="log"))
    }
    if (measure=="MD"){
      trial.res <- lm(outcome ~ factor(trial) + factor(arm)*covar)
    }
    out <- data.frame(summary(trial.res)$coefficients) %>%
      rownames_to_column(var="param")
    return(out)
  }
  
  # select trials with data
  data.all <- data.frame(outcome,trial,arm,covar)
  
  trials.in <- data.all %>%
    group_by(trial) %>%
    summarise(n.covar=length(unique(covar))) %>%
    filter(n.covar>1)
  
  # run for all trials
  res.all <- data.all %>%
    filter(trial %in% trials.in$trial) %>%
    group_by(trial) %>%
    nest() %>%
    mutate(results=map(data,~trial_model(.$outcome,.$arm,.$covar,measure,centre))) %>%
    select(-data) %>%
    unnest(cols=results)
  
  # MAoI
  data.ma <- res.all %>%
    filter(grepl("\\:",param))
  
  MAoI <- data.ma %$% metagen(Estimate,`Std..Error`,sm="OR",studlab=data.ma$trial)
  
  # forest plot
  if(!is.null(filename)){
    filename <- paste("figures/",filename," ","MAoI forest plot.png",sep="")
    png(filename,width=800,height=500)
    forest(MAoI,xlab=measure,leftcols="studlab",rightcols=c("effect", "ci"))
    dev.off()
  }
  
  return(MAoI)
  
}
