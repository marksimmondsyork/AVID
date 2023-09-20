# code to run and save NMAs after set-up


run_NMA <- function(NMA.model,filename,lower=T){
  
  png(paste(filename," - network.png"),width=600,height=400)
  print(plot(NMA.model))
  dev.off()
  
  NMA.re <- nma(NMA.model,
                trt_effects="random",
                prior_intercept=normal(scale=10),
                prior_trt=normal(scale=10),
                prior_het=half_normal(scale=5))
  
  NMA.res <- relative_effects(NMA.re,all_contrasts=TRUE)$summary
  NMA.ranks <- posterior_ranks(NMA.re,lower_better=lower)$summary
  NMA.probs <- posterior_rank_probs(NMA.re,lower_better=lower)$summary
  
  # all results
  sink(paste(filename," - all results.txt"))
  print(NMA.re)
  sink()
  
  # results as tables
  write_csv(NMA.res, paste(filename," - relative effects.csv"))
  write_csv(NMA.ranks,paste(filename," - probability rankings A.csv"))
  write_csv(NMA.probs,paste(filename," - probability rankings B.csv"))
  
  # plots
  tiff(paste(filename," - relative effects forest.tiff"),width=1200,height=800,res=200,compression="lzw")
  print(plot(relative_effects(NMA.re,all_contrasts=TRUE), ref_line=0, .width=c(0,0.95)))
  dev.off()
  
  png(paste(filename," - probability ranking.png"),width=600,height=400)
  print(plot(posterior_rank_probs(NMA.re,lower_better=lower)))
  dev.off()
  
  return(NMA.re)
}
