# code to run threshold analysis from NMA


threshold_run <- function(NMA,data,trt.names,filename,adjusted=FALSE,ref.data=NULL){
  
  if (adjusted==FALSE){
    
    # extract results
    nma.res <- relative_effects(NMA)
    nma.means <- nma.res$summary$mean
    nma.cov <- cov(as.matrix(nma.res$sims[,1,]))
    nma.all <- as.data.frame(relative_effects(NMA,all_contrasts=TRUE)$summary)
    
  }
  
  if (adjusted==TRUE){
    
    nma.res <- relative_effects(NMA,newdata=ref.data)
    nma.means <- nma.res$summary$mean
    nma.cov <- cov(as.matrix(nma.res$sims[,1,]))
    nma.all <- as.data.frame(relative_effects(NMA,newdata=ref.data,all_contrasts=TRUE)$summary)
    
  }
  
  # contrast matrix using external function
  lik.mat <- network_contrast(data$Trialn,data$Trtn)
  
  # # prior covariance
  lik.cov <- recon_vcov(nma.cov,
                        prior.prec=0.001,
                        X=lik.mat)
  
  # threshold analysis
  nma.thresh <- nma_thresh(mean.dk=nma.means,
                           lhood=lik.cov,
                           post=nma.cov,
                           nmatype="fixed",
                           X=lik.mat,
                           opt.max=FALSE,
                           trt.code=trt.names)
  
  c.set <- get_contrasts(lik.mat) # from external code
  
  # plot
  plotdat.nma <- data.frame(lab=paste0(trt.names[c.set$d.b], " vs. ", trt.names[c.set$d.a]),
                            contr.mean=nma.all[c.set$d.i,"mean"],
                            CI2.5 = nma.all[c.set$d.i,"2.5%"],
                            CI97.5 = nma.all[c.set$d.i,"97.5%"])
  
  png(filename,width=1000,height=400)
  thresh_forest(nma.thresh, contr.mean, CI2.5, CI97.5, 
                label = lab, data = plotdat.nma,
                label.title = "Contrast", 
                y.title= "logMAR",
                xlab = "MD", 
                CI.title = "95% Credible Interval",
                xlim = c(-1, 1),
                refline = 0, digits = 2, cutoff = 30,
                add.columns.hjust = 0)
  dev.off()
  
}

