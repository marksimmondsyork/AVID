# convert IPD to summary data for MA

IPD_to_summary <- function(outcome,trial,arm,type){
  
  data <- tibble(trial,arm,outcome)
  npats <- data %>%
    group_by(trial,arm) %>%
    summarise(N=length(outcome[!is.na(outcome)])) %>%
    pivot_wider(id_cols=trial,names_from=arm,values_from=N,names_prefix="N")
  
  if(type=="binary"){
    nevents <- data %>%
      group_by(trial,arm) %>%
      summarise(E=sum(outcome,na.rm=T)) %>%
      pivot_wider(id_cols=trial,names_from=arm,values_from=E,names_prefix="E")
    agg.data <- left_join(npats,nevents)
  }
  
  if(type=="continuous"){
    meanval <- data %>%
      group_by(trial,arm) %>%
      summarise(mean=mean(outcome,na.rm=T),
                sd=sd(outcome,na.rm=T)) %>%
      pivot_wider(id_cols=trial,names_from=arm,values_from=mean:sd)
    agg.data <- left_join(npats,meanval)
  }
  return(agg.data)
}