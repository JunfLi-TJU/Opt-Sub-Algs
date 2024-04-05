setwd("D:/experiment/JMLR/JMLR2024/code")
rm(list = ls())

d_index <- 7

dpath          <- file.path("D:/experiment/online learning dataset/regression/")   

Dataset       <- c("elevators_all","bank_all", "Year_test","ailerons_all","calhousing","N-cpusmall",
                   "N-parkinsons","N-slice_all")               

savepath1      <- paste0("D:/experiment/JMLR/JMLR2024/Result/",
                         paste0("AOMD-OGS-s-",Dataset[d_index],".txt"))

traindatapath  <- file.path(dpath, paste0(Dataset[d_index], ".train"))

traindatamatrix <- as.matrix(read.table(traindatapath))                       
trdata     <- traindatamatrix[ ,-1]
ylabel     <- traindatamatrix[ ,1] 

length_tr  <- nrow(trdata)                                               
feature_tr <- ncol(trdata)              

##############################################################################

# -4 -3 -2 -1 0 1 2 3 4
sigma     <- 8
B0        <- round(sqrt(feature_tr^2+4*feature_tr*length_tr)/2-feature_tr/2)
U         <- 4
coe       <- 0.5

reptimes  <- 10
runtime   <- c(rep(0, reptimes))
errorrate <- c(rep(0, reptimes))
All_bud   <- c(rep(0, reptimes))
sumtdelta <- c(rep(0, reptimes))
bar_t     <- c(rep(length_tr, reptimes))

for( re in 1:reptimes)
{
  order      <- sample(1:length_tr,length_tr,replace = F)   #dis
  #  order      <- c(1:length_tr)
  k          <- 0
  error      <- 0
  alpha      <- 5/length_tr^(1/2)
  alpha_t    <- 0
  sum_delta  <- 0
  sum_t_delta<- 0.04*(U+1)^2
  eta_t      <- U*coe
  Norm       <- 0

  svmat      <- matrix(0,nrow = feature_tr,ncol=1)
  
  Inver_K    <- matrix(0,nrow = 1,ncol=1)           # The inverse kernel matrix
  Gram       <- matrix(0,nrow = 1,ncol=1)           # The inverse kernel matrix
  beta_ast   <- array(0,1)                          # The optimal parameter d
  delta      <- 0                                   # The difference of f''-f'
  kt         <- array(0,1)
  indx       <- array(0,1)
  
  svpara     <- array(0,1)
  t1         <- proc.time()  #proc.time()
  
  ### the first instance
  error      <- (ylabel[order[1]])^2
  svmat[,1]  <- trdata[order[1], ]
  k          <- 1
  Inver_K[1,1] <- 1
  Gram[1,1]    <- 1 
  sum_t_delta  <- sum_t_delta+4*(ylabel[order[1]])^2
  eta_t        <- coe*U/sqrt(sum_t_delta)
  svpara[1]    <- -eta_t*2*(0-ylabel[order[1]])
  Norm  <- abs(svpara[1])
  i     <- 1
  
  ### from the second instance
  while(bar_t[re]>=length_tr && i< length_tr)
  {
    i=i+1
    diff  <- svmat- trdata[order[i], ]
    tem   <- colSums(diff*diff)
    kt    <- exp(tem/(-2*(sigma)^2))
    fx    <- crossprod(svpara[1:k],kt)[1,1]
    
    error <- error + (fx - ylabel[order[i]])^2
    #### compute alpha_t
    beta_ast   <- Inver_K%*%kt
    alpha_t    <- 1-crossprod(beta_ast,kt)[1,1]
    if(alpha_t<0)
      alpha_t  <- 0
    sq_alpha_t <- sqrt(alpha_t)
    if(sq_alpha_t <= alpha)
    {
      tem         <- Gram%*%beta_ast
      tem1        <- crossprod(tem,beta_ast)[1,1]
      sum_t_delta <- sum_t_delta+4*(fx - ylabel[order[i]])^2*tem1
      eta_t       <- coe*U/sqrt(sum_t_delta)
      tem3        <- crossprod(tem,svpara)[1,1]
      g_t         <- eta_t*2*(fx-ylabel[order[i]])
      svpara      <- svpara - g_t*as.vector(beta_ast)
      Norm        <- sqrt(Norm^2-2*g_t*tem3+g_t^2*tem1)
      if(Norm >U)
      {
        svpara <- svpara*U/Norm
        Norm <- U
      }
    }else{
      k           <- k+1
      svmat       <- cbind(svmat,trdata[order[i],])
      sum_t_delta <- sum_t_delta+4*(fx - ylabel[order[i]])^2
      eta_t       <- coe*U/sqrt(sum_t_delta)
      g_t         <- eta_t*2*(fx-ylabel[order[i]])
      svpara[k]   <- -g_t
      
      #        update the inverse kernel matrix 
      tem_d       <- beta_ast
      tem_d[k]    <- -1
      incre       <- tem_d %*% t(tem_d)/alpha_t
      incre[1:(k-1),1:(k-1)] <- incre[1:(k-1),1:(k-1)]+Inver_K
      
      Inver_K     <- incre
      Gram        <- cbind(Gram,kt)
      Gram        <- rbind(Gram,c(kt,1))
      Norm        <- sqrt(Norm^2-2*g_t*fx+g_t^2)
      if(Norm >U)
      {
        svpara <- svpara*U/Norm
        Norm   <- U
      }
    }
    if(k==B0)
      bar_t[re] <- i
  }
  if(bar_t[re]<length_tr)
  {
    for(t in (i+1):length_tr)
    {
      diff  <- svmat- trdata[order[t], ]
      kt    <- exp(colSums(diff*diff)/(-2*(sigma)^2))
      fx    <- crossprod(svpara[1:k],kt)[1,1]
      error <- error + (fx - ylabel[order[t]])^2
      #### compute alpha_t
      k     <- k+1
      svmat        <- cbind(svmat,trdata[order[t],])
      sum_t_delta  <- sum_t_delta + 4*(fx - ylabel[order[t]])^2
      eta_t        <- coe*U/sqrt(sum_t_delta)
      g_t          <- eta_t*2*(fx-ylabel[order[t]])
      svpara[k]    <- -g_t
      Norm         <- sqrt(Norm^2-2*g_t*fx+g_t^2)
      if(Norm > U)
      {
        svpara  <- svpara*U/Norm
        Norm    <- U
      }
    }
  }
  t2 <- proc.time()
  runtime[re]   <- (t2 - t1)[3]
  errorrate[re] <- error/length_tr
  All_bud[re]   <- k
}

save_result <- list(
  note     = c("the next term are:alg_name--dataname--sam_num--sigma--sv_num--run_time--err_num--tot_run_time--ave_run_time--ave_err_rate--sd_time--sd_err"),
  alg_name = c("AOMD-OGS-s-"),
  dataname = paste0(Dataset[d_index], ".train"),
  ker_para = sigma,
  coe      = coe,
  sv_num   = sum(All_bud)/re,
  run_time = as.character(runtime),
  err_num = errorrate,
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time      <- sd(runtime),
  sd_err       <-sd(errorrate)
)

write.table(save_result,file=savepath1,row.names =TRUE, col.names =FALSE, quote = T) 

sprintf("the candidate kernel parameter are :")
sprintf("%.5f", sigma)
sprintf("the number of sample is %d", length_tr)
sprintf("the number of support vectors is %d", round(sum(All_bud)/re))
sprintf("total training time is %.4f in dataset", sum(runtime))
sprintf("average training time is %.5f in dataset", sum(runtime)/reptimes)
sprintf("the average MSE is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of MSE is %.5f in dataset", sd(errorrate))
