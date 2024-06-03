library(parallel)
library(igraph)
library(matrixcalc)
library(MASS)
install.packages("/Users/satabdisaha/Downloads/QUIC", 
                 repos = NULL, 
                 type = "source")
# library(QUIC)
library(Hmisc)
library(robustbase)
library(GENIE3)
library(fitdistrplus)
library(ZIM)

library(ggplot2); theme_set(theme_bw())
library(tidyr)
library(dplyr)

library(zoo) 
library(VGAM) 
library(corpcor) 
library(Matrix) 

#files.sources = list.files("./pcorSimulator", full.names = TRUE)
#sapply(files.sources, source)

source("pcorSimulator.R")
source("controlsPcorSimulator.R")
source("graphStructure.R")
source("lowerTriInd.R")
source("funs-get-res.R")
source("get_mix_parameters.R")

rdir = "/Users/satabdisaha/Documents/Documents/Summer_2019/scSGL/scSGL/simulation_scLink"
dir.create(rdir, recursive = TRUE)



### ------------------------------------------------------------------
### simulate data
### ------------------------------------------------------------------

nc = 5000 # number of cells
B = 30
ncores = 1
thre_no = 10

simulate_hub = function(p, nc = 100, nclust, rho = 0.3,
                        pattern = "powerLaw", 
                        seed = 1, pdf = TRUE, path, ...){
  set.seed(seed)
  sim_cor = pcorSimulatorSimple(nobs = 50, nclusters=nclust, 
                                nnodesxcluster=rep(p,nclust), 
                                pattern=pattern, ...)
  
  theta = sim_cor$omega
  sum(abs(theta) < 1e-5)
  
  # pars = readRDS("~/Box/scSimilarity/code/suppfuns/mixpara.rds")
  # pars = pars[complete.cases(pars), ]
  # pars = pars[, c("rate", "mu","sigma")]
  # pars = pars[order(pars[,"mu"],decreasing = TRUE)[1:1000], ]
  # pars_selected = pars[sample(1:1000, p*nclust), ]
  # sigma = pars_selected[, "sigma"]
  sigma = pbmc_300f_CD4T_count_sd[1:(p*nclust)] 
  
  cov = solve(theta)
  dem = sqrt(diag(cov))
  cov = sweep(cov, MARGIN = 2, dem, FUN = "/")
  cov = sweep(cov, MARGIN = 1, dem, FUN = "/")
  cov = sweep(cov, MARGIN = 2, sigma, FUN = "*")
  cov = sweep(cov, MARGIN = 1, sigma, FUN = "*")
  theta = solve(cov)
  theta[abs(theta) < 1e-5] = 0
  if(!is.positive.definite(theta)) stop("precision matrix is not PD!")  
  
  adj_true = ((theta)!=0)*1
  
  if(pdf == TRUE){
    g = graph_from_adjacency_matrix(adj_true, mode = "undirected", diag = FALSE, weighted = TRUE)
    cols = rep("white", length(V(g)))
    V(g)$color = cols
    ecols = rep("#5C88DA", length(E(g)))
    # ecols[E(g)$weight < 3] = "#CC0C00" # false negative
    E(g)$color = ecols
    ly = layout_nicely(g)
    pdf(path, width = 6, height = 6)
    # plot(g, vertex.label.dist=0, vertex.label.color="black", vertex.label=NA,
    #      vertex.cex = 2, edge.width = 2, layout = ly)
    plot(g, vertex.label=NA, vertex.size = 7, edge.width = 3, layout = ly)
    dev.off()
  }
  
  ### simulate count matrices for all time points -----------------------
  mu = pbmc_300f_CD4T_count_mean[1:(p*nclust)]
  cnt = mvrnorm(n = nc, mu = mu, Sigma = cov)
  cnt[cnt < 0] = 0
  cnt_zero = cnt
  droprate = exp(-rho * cnt^2)
  dropind = rbinom(prod(dim(droprate)), size = 1, prob = as.vector(droprate))
  dropvals = sample(c(0,log10(2)), sum(dropind == 1), prob = c(0.7,0.3), replace = TRUE)
  cnt_zero[dropind == 1] = dropvals
  # cnt_zero = cnt
  # droprate = exp(-rho * mu^2)
  # for(gid in 1:(p*nclust)){
  #   ntp = ceiling(nc*droprate[gid])
  #   cnt_zero[sample(1:nc, ntp),gid] =  sample(c(0,log10(2)), ntp, prob = c(0.7,0.3), replace = TRUE)
  # }
  return(list(cov_true = cov, theta_true = theta, 
              adj_true = adj_true, count = cnt, countd = cnt_zero))
}

theta2pcor = function(theta){
  denom = sqrt(diag(theta))
  pcor = sweep(theta, MARGIN = 1, denom, FUN = "/")
  pcor = sweep(pcor, MARGIN = 2, denom, FUN = "/")
  return(-pcor)
}


nr = 10
l1 = seq(0.5, 0.001, length=nr)
p = 10
nclustvals = c(10, 20, 25)
rhovals = c(0.07, 0.10, 0.13, 0.16)


### simulate data
sim_list = lapply(nclustvals, function(nclust){
  print(paste("nc", nclust))
  tp_res_list = lapply(rhovals, function(rho){
    print(paste("rho", rho))
    rep = mclapply(1:B, function(b){
      if(b %% 10 == 0) print(paste("rep", b)); gc()
      ### simulate data ---------------------------------
      sim_data = simulate_hub(p = p, nc = nc, nclust = nclust, rho = rho, seed = b,
                              pattern = "powerLaw", low.strength = 0.4, sup.strength = 0.7,
                              pdf = FALSE)
      # pdfpath = paste0(rdir, "graph-nc", nclust, "-b-", b, ".pdf")
      # if(rho == 0.1){
      #   sim_data = simulate_hub(p = p, nc = nc, nclust = nclust, rho = rho, seed = b,
      #                           pattern = "powerLaw", low.strength = 0.4, sup.strength = 0.7,
      #                           pdf = TRUE, path = pdfpath)
      # }else{
      # }
      
      # countf = sim_data$count
      countd = sim_data$countd
      
      # save countd
      saveda = t(countd)
      colnames(saveda) = paste0("cell", 1:ncol(saveda))
      write.table(saveda, file = paste0(rdir, "graph-nc", nclust, "-rho-", rho, "-b-", b, ".txt"),
                  quote = FALSE)
      
      return(sim_data)
    }, mc.cores = ncores)
    saveRDS(rep, file = paste0(rdir, "rep-nc", nclust, "-rho", rho, ".rds"))
    gc()
    return(rep)
  })
  return(0)
})



### ------------------------------------------------------------------
### Network Inference
### ------------------------------------------------------------------

res_list = lapply(nclustvals, function(nclust){
  print(paste("nc", nclust))
  tp_res_list = lapply(rhovals, function(rho){
    print(paste("rho", rho))
    simdata = readRDS(paste0(rdir, "rep-nc", nclust, "-rho", rho, ".rds"))
    rep = mclapply(1:B, function(b){
      # if(b %% 10 == 0)
      print(paste("rep", b)); gc()
      ### simulate data ---------------------------------
      
      sim_data = simdata[[b]]
      countf = sim_data$count
      countd = sim_data$countd
      adj_true = sim_data$adj_true
      theta_true = sim_data$theta_true
      # pcor_true = theta2pcor(theta_true)
      epercent = mean(adj_true[upper.tri(adj_true)])
      zpercent = mean(countd == 0) - mean(countf == 0)
      
      ### estimate covariance matrices
      cor_pearson = cor(countd)
      x.q=apply(countd,2,sd)
      cov_pearson=diag(x.q)%*%cor_pearson%*%diag(x.q)
      
      nobs = nrow(countd)
      cor_pearson.c = est_pearson.c(countd, cor.p = cor_pearson, thre_no = 10)
      
      ### estimate precision matrices
      aucscore = t(sapply(1:length(methods), function(k){
        method = methods[k]
        # print(method)
        if(method == "scLink"){
          weights_p = 1-abs(cor_pearson.c)
          cov_pearson.c=diag(x.q)%*%cor_pearson.c%*%diag(x.q)
          cov_pearson.c = easy.psd(cov_pearson.c,method="perturb")
          res_seq = get_res_wGL_quic(cov_pearson.c, weights_p, nobs, l1, adj_true)
        }
        if(method == "glasso"){res_seq = get_res_GL(cov_pearson, nobs, l1, adj_true)}
        if(method == "pearson"){
          res_seq = get_res_cors(countd, adj_true, method = method,
                                 cmat = cor_pearson, thre_no = thre_no)
        }
        if(method == "spearman"){
          res_seq = get_res_cors(countd, adj_true, method = method, thre_no = thre_no)
        }
        if(method == "PIDC"){
          wpath = paste0(rdir, "PIDC/weight-p", p, "-rho-", rho, "-b-", b, ".txt")
          res_seq = get_res_pidc(wpath, adj_true)
        }
        if(method == "GENIE3"){
          res_seq = get_res_genie3(countd, adj_true)
        }
        
        
        # theta_list = res_seq$theta_list
        accuracy = res_seq$accuracy
        
        colnames(accuracy) = c("nE", "bic", "precision", "recall", "FP", "lambda1")
        accuracy = as.data.frame(accuracy)
        accuracy = accuracy[complete.cases(accuracy), , drop = FALSE]
        if(nrow(accuracy) < 2) return(c(NA, NA))
        
        acc = accuracy
        tp = acc[, c("precision","recall")]
        tp = tp[order(tp[,1]), ]
        tp = rbind(c(tp[1,1],1), tp)
        tp = rbind(tp, c(1,0))
        AUPR = sum(diff(tp[,1])*rollmean(tp[,2],2))
        tp = acc[, c("FP","recall")]
        tp = tp[order(tp[,1]), ]
        tp = rbind(c(0,0), tp)
        tp = rbind(tp, c(1,1))
        AUROC = sum(diff(tp[,1])*rollmean(tp[,2],2))
        return(c(AUPR, AUROC))
      }))
      
      colnames(aucscore) = c("AUPR", "AUROC")
      aucscore = as.data.frame(aucscore)
      aucscore$method = methods
      aucscore$b = b
      aucscore$rho = rho
      aucscore$p = p
      aucscore$zpercent = zpercent
      aucscore$epercent = epercent
      return(aucscore)
    }, mc.cores = ncores)
    rep = Reduce(rbind, rep)
    saveRDS(rep, file = paste0(rdir, "est-nc", nclust, "-rho", rho, ".rds"))
    gc()
    return(rep)
  })
  tp_res = Reduce(rbind, tp_res_list)
  return(tp_res)
})
