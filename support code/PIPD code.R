# Code for pseudo-IPD creation

# generate Pseudo-IPD for BCVA and other cont. outcomes

BCVA_PIPD <- function(Trial,Treat,Time,N,Mean,SD){
  
  data.IPD <- data.frame(trial=rep(Trial,N), treat=rep(Treat,N), time=rep(Time,N),
                         mean=rep(Mean,N), sd=rep(SD,N))
  
  # generate pseudo-IPD
  set.seed(4533776)
  data.IPD$ytemp <- rnorm(nrow(data.IPD), 0, 1)
  
  data.IPD.summ <- data.IPD %>%
    group_by(trial,treat,time) %>%
    summarise(ytemp.m=mean(ytemp),
              ytemp.sd=sd(ytemp))
  
  data.IPD.2 <- left_join(data.IPD,data.IPD.summ) %>%
    mutate(y=mean + (ytemp-ytemp.m)*(sd/ytemp.sd))
  
  return(data.IPD.2)
}

# generate pseudo-IPD from binary outcomes

Binary_PIPD <- function(Trial,Treat,Time,N, Outcome){
  
  data.IPD <- data.frame(trial=rep(Trial,N), treat=rep(Treat,N), time=rep(Time,N))
  
outcome_pipd <- function(N,Outcome){
  if(is.na(Outcome))
    out <- rep(NA,N)
  if(!is.na(Outcome))
    out <- rep(c(1,0),c(Outcome,N-Outcome))
  return(out)
}
  data.IPD$outcome <- map2(N, Outcome, outcome_pipd) %>%
    do.call(c,.)
  return(data.IPD)
}
