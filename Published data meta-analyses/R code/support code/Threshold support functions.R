###############################
# construct network likelihood matrix 
# for threshold analysis
###############################

# input:
# study, vector of study codes
# trt, vector of treatment arms (one row per arm)
# treatments must be numerical 1...t
# base.trt, baseline (reference) treatment

network_contrast <- function(study,trt,base.trt=1){
  
  data <- data.frame(study,trt)
  
  # comparisons in study
  comp_study <- function(x,data){
    datai <- filter(data,study==x)
    comps <- expand_grid(base=datai$trt,comp=datai$trt) %>%
      filter(comp>base)
    return(comps)
  }
  
  # all comparisons
  all.comps <- map_dfr(unique(data$study), ~comp_study(.,data)) %>%
    arrange(base,comp)
  inc.comps <- distinct(all.comps,base,comp)
  
  # construct matrix
  n.comps <- length(inc.comps$base)
  n.trt <- length(unique(trt))
  lik.matrix <- matrix(0,n.comps,n.trt)
  
  for (i in 1:n.comps){
    lik.matrix[i,inc.comps$base[i]] <- -1
    lik.matrix[i,inc.comps$comp[i]] <- 1
  }
  
  lik.matrix <- lik.matrix[,-base.trt]
  return(lik.matrix)
}

# get contrasts from the matrix produced above for plotting

get_contrasts <- function(net.matrix){
  
  K <- ncol(net.matrix) + 1
  all.cols <- expand.grid(b=1:K,a=1:K)
  all.cols <- filter(all.cols, b>a)
  all.cols <- subset(all.cols, b>a)
  
  d.a <- d.b <- d.i <- vector(length=nrow(net.matrix))
  for (i in 1:nrow(net.matrix)){
    d.a[i] <- ifelse(any(net.matrix[i, ] == -1), which(net.matrix[i, ] == -1), 0) + 1
    d.b[i] <- ifelse(any(net.matrix[i, ] == 1), which(net.matrix[i, ] == 1), 0) + 1
    for (j in 1:nrow(all.cols)){
      if(d.a[i]==all.cols$a[j] & d.b[i]==all.cols$b[j])
        d.i[i] <- j
    }
  }
  
  # Transform from d_ab style contrast references into d[i] style
  # d.i <- nmathresh::d_ab2i(d.a, d.b, K = ncol(net.matrix)+1)
  
  out <- data.frame(d.a,d.b,d.i)
  return(out)
}
