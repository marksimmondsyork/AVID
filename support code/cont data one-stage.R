# Continuous data MA using LMER

Cont_data_MA <- function(Study,Arm,N,Mean,SD){
  
  require(tidyverse)
  require(lme4)
  
  data.IPD <- data.frame(study=rep(Study,N), arm=rep(Arm,N),
                         mean=rep(Mean,N), sd=rep(SD,N))
  
  # generate pseudo-IPD
  set.seed(4533776)
  data.IPD$ytemp <- rnorm(nrow(data.IPD), 0, 1)#
  
  data.IPD.summ <- data.IPD %>%
    group_by(study,arm) %>%
    summarise(ytemp.m=mean(ytemp),
              ytemp.sd=sd(ytemp))
  
  data.IPD.2 <- left_join(data.IPD,data.IPD.summ) %>%
    mutate(y=mean + (ytemp-ytemp.m)*(sd/ytemp.sd),
           arm2=1000*study + arm)
  
  # model with lmer
  lmer.res <- lmer(y ~ arm + (1+arm|study), data=data.IPD.2)
  
  lmer.coeffs <- summary(lmer.res)$coefficients
  lmer.hetgty <- as.numeric(attr(summary(lmer.res)$varcor$study,"stddev")[2])^2
  lmer.out <- data.frame(lmer.est=lmer.coeffs[2,1],
                         lmer.sd=lmer.coeffs[2,2],
                         lmer.hetgty)

  out <- lmer.out
  return(out)
  
  }
