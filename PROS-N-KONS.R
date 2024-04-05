setwd("D:/experiment/JMLR/JMLR2024/code")

rm(list = ls())

d_index <- 8

dpath          <- file.path("D:/experiment/online learning dataset/regression/") 

Dataset       <- c("elevators_all","bank_all", "Year_test","ailerons_all","calhousing","N-cpusmall",
                   "N-parkinsons","N-slice_all")

savepath1      <- paste0("D:/experiment/JMLR/JMLR2024/Result/",
                         paste0("PROS-sq-",Dataset[d_index],".txt"))

traindatapath  <- file.path(dpath, paste0(Dataset[d_index], ".train"))

traindatamatrix <- as.matrix(read.table(traindatapath))                       
trdata     <- traindatamatrix[ ,-1]
ylabel     <- traindatamatrix[ ,1] 

length_tr  <- nrow(trdata)                                               
feature_tr <- ncol(trdata)              

baseline   <- crossprod(ylabel,ylabel)[1,1]/length_tr

##############################################################################

# -4 -3 -2 -1 0 1 2 3 4
sigma     <- 2^6
reptimes  <- 10

gamma     <- 5
alpha     <- 1
beta      <- 1
varepsilon   <- 0.5
C         <- 1
sigma_    <- 1/(8*C^2)

runtime   <- c(rep(0, reptimes))
errorrate <- c(rep(0, reptimes))
B         <- c(rep(0, reptimes))

for(re in 1:reptimes)
{
  order      <- sample(1:length_tr,length_tr,replace = F)   #dis
  #  order       <- c(1:length_tr)
  
  svmat      <- matrix(0,nrow = feature_tr,ncol=1)
  svpara     <- array(0,1)
  t1         <- proc.time()  #proc.time()
  
  ### the first instance
  error      <- (ylabel[order[1]])^2
  Gram       <- matrix(1,nrow = 1,ncol=1)
  A          <- matrix(alpha,nrow = 1,ncol=1)
  inver_A    <- matrix(1/alpha,nrow = 1,ncol=1)
  S          <- array(1,1)
  k          <- 1
  
  svmat[,1]  <- trdata[order[1], ]
  zt_1       <- 1
  Gramt_1    <- Gram
  inver_K    <- (1+gamma)^{-1}
  
  ### from the second instance
  for(t in 2:length_tr)
  {
    if(zt_1 == 1)
    {
      T <- svd(Gramt_1)
      U <- T$u
      D <- T$d  ##vector
      
      D <- D^{-0.5}
      
      num     <- ncol(Gramt_1)
      A       <- alpha*diag(num)
      inver_A <- 1/alpha*diag(num)
      wt      <- c(rep(0, num))
      
      diff    <- svmat - trdata[order[t], ]
      tem     <- colSums(diff*diff)
      kt      <- exp(tem/(-2*(sigma^2)))
      tem     <- t(U)%*%kt
      ph_t    <- D*tem
    }else{
      diff    <- svmat - trdata[order[t], ]
      tem     <- colSums(diff*diff)
      kt      <- exp(tem/(-2*(sigma^2)))
      tem     <- t(U)%*%kt
      ph_t    <- D*tem
      
      z       <- crossprod(vt,ph_t)[1,1]
      tem2    <- inver_A%*%ph_t
      tem1    <- (t(ph_t)%*%tem2)[1,1]
      if(tem1==0)
        tem1 <- 0.00001
      wt      <- vt - sign(z)*max(abs(z)-C,0)/tem1*tem2
    }
    hatyt   <- crossprod(wt,ph_t)[1,1]
    gt      <- 2*(hatyt-ylabel[order[t]])*ph_t
    A       <- A + sigma_/2*gt%*%t(gt)
    error   <- error + (hatyt-ylabel[order[t]])^2
    
    ########### solving A^-1
    A_gt    <- inver_A%*%gt
    tem3    <- 1+sigma_/2*(t(gt)%*%A_gt)[1,1]
    tem4    <- sigma_/2*A_gt%*%t(gt)%*%inver_A
    inver_A <- inver_A - tem4/tem3
    
    vt      <- wt - A_gt
    
    ############ kernelized online row sampling    
    k         <- k+1
    copy_S    <- S
    Gram      <- cbind(Gram,kt)
    Gram      <- rbind(Gram,c(kt,1)) 
    S[k]      <- 1
    S_mat     <- diag(S)
    kt[k]     <- 1
    tem3      <- kt*S     #    tem3  <- diag(S)%*%kt
    tem1      <- t(tem3)  #    tem1  <- t(kt)%*%diag(S)
    tem_inver_K <- inver_K
    
    kt_ <- tem3[1:(k-1)]
    AB <- inver_K%*%kt_
    E  <- (gamma+1-t(kt_)%*%AB)^{-1}
    E  <- E[1,1]
    E2 <- -E*AB
    E3 <- t(E2)
    E1 <- inver_K + E*AB%*%t(AB)
    inver_K  <- cbind(E1,E2)
    inver_K  <- rbind(inver_K,c(E3,E))
    kk_inver_K <-inver_K
    tem4  <- 1-tem1%*%inver_K%*%tem3
    pt    <- min(beta*(1+varepsilon)*tem4/gamma,1)
    bt    <- rbinom(1,1,pt)

    if(bt==1)
    {
      S[k]    <- 1/sqrt(pt)
      svmat   <- cbind(svmat,trdata[order[t], ])
      Gramt_1 <- Gram
      zt_1    <- 1
      
      kt_ <- tem3[1:(k-1)]/sqrt(pt)
      AB  <- tem_inver_K%*%kt_
      E   <- (gamma+1/pt-t(kt_)%*%AB)^{-1}
      E  <- E[1,1]
      E2 <- -E*AB
      E3 <- t(E2)
      E1 <- tem_inver_K + E*AB%*%t(AB)
      inver_K  <- cbind(E1,E2)
      inver_K  <- rbind(inver_K,c(E3,E))
    }else{
      k    <- k-1
      Gram <- Gramt_1
      S    <- copy_S
      zt_1 <- 0
      inver_K <- tem_inver_K 
    }
  }
  t2 <- proc.time()
  B[re]         <- k
  runtime[re]   <- (t2 - t1)[3]
  errorrate[re] <- error/length_tr
}

save_result <- list(
  note     = c("the next term are:alg_name--dataname--sam_num--sigma--sv_num--run_time--err_num--tot_run_time--ave_run_time--ave_err_rate--sd_time--sd_err"),
  alg_name = c("PROS-r-"),
  dataname = paste0(Dataset[d_index], ".train"),
  sam_num  = length_tr,
  ker_para = sigma,
  gamma    = gamma,
  alpha    = alpha,
  beta     = beta,
  varepsilon = varepsilon,
  sv_num   = mean(B),
  run_time = as.character(runtime),
  err_num = errorrate,
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time  <- sd(runtime),
  sd_err    <-sd(errorrate)
)

write.table(save_result,file=savepath1,row.names =TRUE, col.names =FALSE, quote = F) 

sprintf("the candidate kernel parameter are :")
sprintf("%.5f", sigma)
sprintf("the number of sample is %d", length_tr)
sprintf("the number of support vectors is %.4f", mean(B))
sprintf("total training time is %.4f in dataset", sum(runtime))
sprintf("average training time is %.5f in dataset", sum(runtime)/reptimes)
sprintf("the MSE is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of MSE is %.5f in dataset", sd(errorrate))
