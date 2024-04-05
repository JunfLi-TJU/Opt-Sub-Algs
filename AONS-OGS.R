setwd("D:/experiment/JMLR/JMLR2024/code")
rm(list = ls())

d_index <- 7

dpath          <- file.path("D:/experiment/online learning dataset/regression/")   

Dataset        <- c("elevators_all","bank_all", "Year_test","ailerons_all","calhousing","N-cpusmall",
                   "N-parkinsons","N-slice_all")                 

savepath1      <- paste0("D:/experiment/JMLR/JMLR2024/Result/",
                         paste0("AONS-OGS-",Dataset[d_index],".txt"))

traindatapath  <- file.path(dpath, paste0(Dataset[d_index], ".train"))

traindatamatrix <- as.matrix(read.table(traindatapath))                       
trdata     <- traindatamatrix[ ,-1]
ylabel     <- traindatamatrix[ ,1] 

length_tr  <- nrow(trdata)                                               
feature_tr <- ncol(trdata)              

##############################################################################

# -4 -3 -2 -1 0 1 2 3 4
sigma     <- 2^3
C         <- 1

reptimes  <- 10
runtime   <- c(rep(0, reptimes))
errorrate <- c(rep(0, reptimes))
All_bud   <- c(rep(0, reptimes))

for( re in 1:reptimes)
{
  order      <- sample(1:length_tr,length_tr,replace = F)   #dis
  #  order      <- c(1:length_tr)
  k          <- 0
  error      <- 0
  alpha      <- 5/length_tr^(1/2)
  alpha_t    <- 0
  eta_t      <- 1/(4*C^2+4)
  mu         <- 1

  svmat      <- matrix(0,nrow = feature_tr,ncol=1)
  
  Inver_K    <- matrix(0,nrow = 1,ncol=1)           # The inverse kernel matrix
  Gram       <- matrix(0,nrow = 1,ncol=1)           # The inverse kernel matrix
  copy_u     <- Gram
  copy_D     <- Gram
  beta_ast   <- array(0,1)                          # The optimal parameter d
  kt         <- array(0,1)
  
  t1         <- proc.time()  #proc.time()
  
  ### the first instance
  error      <- (ylabel[order[1]])^2
  svmat[,1]  <- trdata[order[1], ]
  k          <- 1
  Inver_K[1,1] <- 1
  Gram[1,1]    <- 1 
  
  diff    <- svmat- trdata[order[2], ]
  tem     <- colSums(diff*diff)
  kt      <- exp(tem/(-2*(sigma)^2))
  A       <- mu
  inver_A <- 1/mu
  wt      <- 0
  T       <- svd(Gram)
  U       <- T$u
  Sig     <- T$d  ## vector
  
  D       <- Sig^{-0.5}
  tem     <- t(U)%*%kt
  
  copy_u  <- U
  copy_D  <- D
  
  ### from the second instance
  for(t in 2:(length_tr-1))
  {

    #### compute alpha_t
    beta_ast   <- Inver_K%*%kt
    alpha_t    <- 1-crossprod(beta_ast,kt)[1,1]
    if(alpha_t<0)
      alpha_t  <- 0
    if(alpha_t <= alpha^2)
    {
      ph_t    <- D*tem 
      hatyt   <- crossprod(wt,ph_t)[1,1]
      
      gt      <- 2*(hatyt-ylabel[order[t]])*ph_t
      t_gt    <- t(gt)
      A       <- A + eta_t*gt%*%t_gt
      error   <- error + (hatyt-ylabel[order[t]])^2
      
      ########### solving A^-1
      A_gt    <- inver_A%*%gt
      tem1    <- 1+eta_t*(t_gt%*%A_gt)[1,1]
      tem2    <- eta_t*A_gt%*%t_gt%*%inver_A
      inver_A <- inver_A - tem2/tem1
      
      ########### projection 
      vt      <- wt - A_gt

      diff    <- svmat- trdata[order[t+1], ]
      tem     <- colSums(diff*diff)
      kt      <- exp(tem/(-2*(sigma)^2))
      tem     <- t(U)%*%kt
      ph_t_1  <- D*tem
      
      tildeyt <- crossprod(vt,ph_t_1)[1,1]
      z2      <- inver_A%*%ph_t_1
      z1      <- (t(ph_t_1)%*%z2)[1,1]
      if(z1==0)
        z1 <- 0.00001
      wt      <- vt - sign(tildeyt)*max(abs(tildeyt)-C,0)/z1*z2
    }else{
      k           <- k+1
      Gram        <- cbind(Gram,kt)
      Gram        <- rbind(Gram,c(kt,1))
      T           <- svd(Gram)
      U           <- T$u
      Sig         <- T$d  ## type:vector
      t_U         <- t(U)
      
      D        <- Sig^{-0.5}
      Gram_k_1 <- Gram[,1:(k-1)]%*%copy_u
      Q        <- diag(D) %*% t_U %*% Gram_k_1 %*% diag(copy_D)
      A        <- (A-mu*diag(k-1))
      A        <- mu*diag(k)+Q %*% A %*% t(Q)
      wt       <- Q %*% wt
      
      tem      <- t_U %*% Gram[,k]
      ph_t     <- D * tem 
      hatyt    <- crossprod(wt,ph_t)[1,1]
      gt       <- 2*(hatyt-ylabel[order[t]])*ph_t
      A        <- A + eta_t*gt%*%t(gt)
      error    <- error + (hatyt-ylabel[order[t]])^2
      
      ########### solving A^-1
      inver_A <- solve(A)
      
      ########### projection 
      vt      <- wt - inver_A%*%gt
      
      svmat   <- cbind(svmat,trdata[order[t],])
      
      diff    <- svmat- trdata[order[t+1], ]
      tem     <- colSums(diff*diff)
      kt      <- exp(tem/(-2*(sigma)^2))
      tem     <- t_U%*%kt
      ph_t_1  <- D*tem
      
      tildeyt <- crossprod(vt,ph_t_1)[1,1]
      z2      <- inver_A%*%ph_t_1
      z1      <- (t(ph_t_1)%*%z2)[1,1]
      if(z1==0)
        z1 <- 0.00001
      wt      <- vt - sign(tildeyt)*max(abs(tildeyt)-C,0)/z1*z2
      
      copy_u   <- U
      copy_D   <- D
      
      #        update the inverse kernel matrix 
      tem_d       <- beta_ast
      tem_d[k]    <- -1
      incre       <- tem_d %*% t(tem_d)/alpha_t
      incre[1:(k-1),1:(k-1)] <- incre[1:(k-1),1:(k-1)]+Inver_K
      Inver_K     <- incre  
    }
  }
  t       <- length_tr
  hatyt   <- crossprod(wt,ph_t_1)[1,1]
  error   <- error + (hatyt-ylabel[order[t]])^2
  
  t2 <- proc.time()
  runtime[re]   <- (t2 - t1)[3]
  errorrate[re] <- error/length_tr
  All_bud[re]   <- k
}

save_result <- list(
  note     = c("the next term are:alg_name--dataname--sam_num--sigma--sv_num--run_time--err_num--tot_run_time--ave_run_time--ave_err_rate--sd_time--sd_err"),
  alg_name = c("AONS-OGS-"),
  dataname = paste0(Dataset[d_index], ".train"),
  ker_para = sigma,
  mu      = mu,
  sv_num   = sum(All_bud)/re,
  run_time = as.character(runtime),
  err_num = errorrate,
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time      <- sd(runtime),
  sd_err       <-sd(errorrate)
)

write.table(save_result,file=savepath1,row.names =TRUE, col.names =FALSE, quote = F) 

sprintf("the candidate kernel parameter are :")
sprintf("%.5f", sigma)
sprintf("the number of sample is %d", length_tr)
sprintf("the number of support vectors is %d", round(sum(All_bud)/re))
sprintf("total training time is %.4f in dataset", sum(runtime))
sprintf("average training time is %.5f in dataset", sum(runtime)/reptimes)
sprintf("the average MSE is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of MSE is %.5f in dataset", sd(errorrate))
