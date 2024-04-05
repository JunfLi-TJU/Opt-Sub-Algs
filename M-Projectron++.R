setwd("D:/experiment/JMLR/JMLR2024/code")
rm(list = ls())


d_index <- 10


dpath         <- file.path("D:/experiment/online learning dataset/multi-class") 

Dataset       <- c("usps_all", "protein_all", "mnist_test","cifar10","letter",
                   "poker","shuttle_all","segmentation","waveform","covtype")

savepath1      <- paste0("D:/experiment/JMLR/JMLR2024/Result/",
                         paste0("M-Projectron++-",Dataset[d_index],".txt"))

traindatapath  <- file.path(dpath, paste0(Dataset[d_index], ".train"))

traindatamatrix <- as.matrix(read.table(traindatapath))                       
trdata     <- traindatamatrix[ ,-1]
ylabel     <- traindatamatrix[ ,1] 

length_tr  <- nrow(trdata)                                               
feature_tr <- ncol(trdata)              
if(min(ylabel)<= 0)
  ylabel <- -min(ylabel)+ylabel+1
##############################################################################

# -4 -3 -2 -1 0 1 2 3 4
sigma     <- 2^(11)
eta       <- 500/length_tr^(1)
M         <- max(ylabel)

reptimes  <- 10
runtime   <- c(rep(0, reptimes))
errorrate <- c(rep(0, reptimes))
All_bud   <- c(rep(0, reptimes))

for( re in 1:reptimes)
{
  order      <- sample(1:length_tr,length_tr,replace = F)   #dis
  k          <- c(rep(1, M))
  sum        <- c(rep(0, M))
  error      <- 0
  
  svmat      <- matrix(0,nrow = feature_tr,ncol=1)
  
  Inver_K    <- matrix(1,nrow = 1,ncol=1)           # The inverse kernel matrix
  Gram       <- matrix(1,nrow = 1,ncol=1)           # The inverse kernel matrix
  d_ast      <- array(0,1)                          # The optimal parameter d
  delta      <- 0                                   # The difference of f''-f'
  kt         <- array(0,1)
  
###################################################################################
  
  kt_list <- list(
    kt1   = array(0,1),
    kt2   = array(0,1),
    kt3   = array(0,1),
    kt4   = array(0,1),
    kt5   = array(0,1),
    kt6   = array(0,1),
    kt7   = array(0,1),
    kt8   = array(0,1),
    kt9   = array(0,1),
    kt10   = array(0,1),
    kt11   = array(0,1),
    kt12   = array(0,1),
    kt13   = array(0,1),
    kt14   = array(0,1),
    kt15   = array(0,1),
    kt16   = array(0,1),
    kt17   = array(0,1),
    kt18   = array(0,1),
    kt19   = array(0,1),
    kt20   = array(0,1),
    kt21   = array(0,1),
    kt22   = array(0,1),
    kt23   = array(0,1),
    kt24   = array(0,1),
    kt25   = array(0,1),
    kt26   = array(0,1)
  )
  
  sv_coe_list <- list(
    svpara1   = array(-1,1),
    svpara2   = array(-1,1),
    svpara3   = array(-1,1),
    svpara4   = array(-1,1),
    svpara5   = array(-1,1),
    svpara6   = array(-1,1),
    svpara7   = array(-1,1),
    svpara8   = array(-1,1),
    svpara9   = array(-1,1),
    svpara10   = array(-1,1),
    svpara11   = array(-1,1),
    svpara12   = array(-1,1),
    svpara13   = array(-1,1),
    svpara14   = array(-1,1),
    svpara15   = array(-1,1),
    svpara16   = array(-1,1),
    svpara17   = array(-1,1),
    svpara18   = array(-1,1),
    svpara19   = array(-1,1),
    svpara20   = array(-1,1),
    svpara21   = array(-1,1),
    svpara22   = array(-1,1),
    svpara23   = array(-1,1),
    svpara24   = array(-1,1),
    svpara25   = array(-1,1),
    svpara26   = array(-1,1)
  )
  
  Gram_list <- list(
    Gram1  = matrix(1,nrow = 1,ncol=1),
    Gram2  = matrix(1,nrow = 1,ncol=1),
    Gram3  = matrix(1,nrow = 1,ncol=1),
    Gram4  = matrix(1,nrow = 1,ncol=1),
    Gram5  = matrix(1,nrow = 1,ncol=1),
    Gram6  = matrix(1,nrow = 1,ncol=1),
    Gram7  = matrix(1,nrow = 1,ncol=1),
    Gram8  = matrix(1,nrow = 1,ncol=1),
    Gram9  = matrix(1,nrow = 1,ncol=1),
    Gram10  = matrix(1,nrow = 1,ncol=1),
    Gram11  = matrix(1,nrow = 1,ncol=1),
    Gram12  = matrix(1,nrow = 1,ncol=1),
    Gram13  = matrix(1,nrow = 1,ncol=1),
    Gram14  = matrix(1,nrow = 1,ncol=1),
    Gram15  = matrix(1,nrow = 1,ncol=1),
    Gram16  = matrix(1,nrow = 1,ncol=1),
    Gram17  = matrix(1,nrow = 1,ncol=1),
    Gram18  = matrix(1,nrow = 1,ncol=1),
    Gram19  = matrix(1,nrow = 1,ncol=1),
    Gram20  = matrix(1,nrow = 1,ncol=1),
    Gram21  = matrix(1,nrow = 1,ncol=1),
    Gram22  = matrix(1,nrow = 1,ncol=1),
    Gram23  = matrix(1,nrow = 1,ncol=1),
    Gram24  = matrix(1,nrow = 1,ncol=1),
    Gram25  = matrix(1,nrow = 1,ncol=1),
    Gram26  = matrix(1,nrow = 1,ncol=1)
  )
  
  InvGram_list <- list(
    InvGram1  = matrix(1,nrow = 1,ncol=1),
    InvGram2  = matrix(1,nrow = 1,ncol=1),
    InvGram3  = matrix(1,nrow = 1,ncol=1),
    InvGram4  = matrix(1,nrow = 1,ncol=1),
    InvGram5  = matrix(1,nrow = 1,ncol=1),
    InvGram6  = matrix(1,nrow = 1,ncol=1),
    InvGram7  = matrix(1,nrow = 1,ncol=1),
    InvGram8  = matrix(1,nrow = 1,ncol=1),
    InvGram9  = matrix(1,nrow = 1,ncol=1),
    InvGram10  = matrix(1,nrow = 1,ncol=1),
    InvGram11  = matrix(1,nrow = 1,ncol=1),
    InvGram12  = matrix(1,nrow = 1,ncol=1),
    InvGram13  = matrix(1,nrow = 1,ncol=1),
    InvGram14  = matrix(1,nrow = 1,ncol=1),
    InvGram15  = matrix(1,nrow = 1,ncol=1),
    InvGram16  = matrix(1,nrow = 1,ncol=1),
    InvGram17  = matrix(1,nrow = 1,ncol=1),
    InvGram18  = matrix(1,nrow = 1,ncol=1),
    InvGram19  = matrix(1,nrow = 1,ncol=1),
    InvGram20  = matrix(1,nrow = 1,ncol=1),
    InvGram21  = matrix(1,nrow = 1,ncol=1),
    InvGram22  = matrix(1,nrow = 1,ncol=1),
    InvGram23  = matrix(1,nrow = 1,ncol=1),
    InvGram24  = matrix(1,nrow = 1,ncol=1),
    InvGram25  = matrix(1,nrow = 1,ncol=1),
    InvGram26  = matrix(1,nrow = 1,ncol=1)
  )
  
  sv_max_list <- list(
    svmat1  = matrix(0,nrow = feature_tr,ncol=1),
    svmat2  = matrix(0,nrow = feature_tr,ncol=1),
    svmat3  = matrix(0,nrow = feature_tr,ncol=1),
    svmat4  = matrix(0,nrow = feature_tr,ncol=1),
    svmat5  = matrix(0,nrow = feature_tr,ncol=1),
    svmat6  = matrix(0,nrow = feature_tr,ncol=1),
    svmat7  = matrix(0,nrow = feature_tr,ncol=1),
    svmat8  = matrix(0,nrow = feature_tr,ncol=1),
    svmat9  = matrix(0,nrow = feature_tr,ncol=1),
    svmat10  = matrix(0,nrow = feature_tr,ncol=1),
    svmat11  = matrix(0,nrow = feature_tr,ncol=1),
    svmat12  = matrix(0,nrow = feature_tr,ncol=1),
    svmat13  = matrix(0,nrow = feature_tr,ncol=1),
    svmat14  = matrix(0,nrow = feature_tr,ncol=1),
    svmat15  = matrix(0,nrow = feature_tr,ncol=1),
    svmat16  = matrix(0,nrow = feature_tr,ncol=1),
    svmat17  = matrix(0,nrow = feature_tr,ncol=1),
    svmat18  = matrix(0,nrow = feature_tr,ncol=1),
    svmat19  = matrix(0,nrow = feature_tr,ncol=1),
    svmat20  = matrix(0,nrow = feature_tr,ncol=1),
    svmat21  = matrix(0,nrow = feature_tr,ncol=1),
    svmat22  = matrix(0,nrow = feature_tr,ncol=1),
    svmat23  = matrix(0,nrow = feature_tr,ncol=1),
    svmat24  = matrix(0,nrow = feature_tr,ncol=1),
    svmat25  = matrix(0,nrow = feature_tr,ncol=1),
    svmat26  = matrix(0,nrow = feature_tr,ncol=1)
  ) 
  
  t1         <- proc.time()  #proc.time()
  
  ### the first instance
  
  error      <- 1
  for(h in 1:M)
  {
    svmat      <- sv_max_list[[h]]
    svmat[,1]  <- trdata[order[1], ]
    sv_max_list[[h]]  <- svmat
  }
  y_t          <- ylabel[order[1]]
  sv_coe_list[[y_t]][1]    <- 1
  
###################################################################################
  t1           <- proc.time()  #proc.time()

  ### from the second instance
  for (t in 2:length_tr)
  {
    for(h in 1:M)
    {
      svmat  <- sv_max_list[[h]]
      svpara <- sv_coe_list[[h]]
      diff   <- svmat - trdata[order[t],]
      kt     <- exp(colSums(diff*diff)/(-2*sigma^2))
      sum[h] <- crossprod(svpara[1:k[h]],kt)[1,1]
      kt_list[[h]] <- kt
    }
    
    hat_y <- which(sum==max(sum))[1]
    y_t   <- ylabel[order[t]]

    tem_sum      <- sum
    tem_sum[y_t] <- min(sum)-0.1
    i_t          <- which(tem_sum==max(tem_sum))[1]
    if(hat_y != y_t)
    {
      error = error+1
      for(h in c(y_t,i_t))
      {
        svmat   <- sv_max_list[[h]]
        svpara  <- sv_coe_list[[h]]
        Gram    <- Gram_list[[h]]
        Inver_K <- InvGram_list[[h]]
        kt      <- kt_list[[h]]
        
        #### compute delta
        d_ast   <-Inver_K%*%kt
        delta   <- 1-crossprod(d_ast,kt)[1,1]
        if(abs(delta)<1e-8)
          delta <- 0
        if(delta <= eta)
        {
          if(h == y_t)
            svpara <- svpara + as.vector(d_ast)
          else
            svpara <- svpara - as.vector(d_ast)
          sv_coe_list[[h]]  <- svpara
        }else{
          k[h]        <- k[h]+1
          svmat       <- cbind(svmat,trdata[order[t],])
          
          if(h == y_t)
            svpara[k[h]]   <- 1
          else
            svpara[k[h]]   <- -1
          
          #update the inverse kernel matrix 
          tem_d       <- d_ast
          tem_d[k[h]] <- -1
          incre       <- tem_d %*% t(tem_d)/delta
          incre[1:(k[h]-1),1:(k[h]-1)] <- incre[1:(k[h]-1),1:(k[h]-1)]+Inver_K
          Inver_K     <- incre
          Gram        <- cbind(Gram,kt)
          Gram        <- rbind(Gram,c(kt,1)) 
          
          sv_max_list[[h]]  <- svmat
          sv_coe_list[[h]]  <- svpara
          Gram_list[[h]]    <- Gram
          InvGram_list[[h]] <- Inver_K
        }
      }
    }
    if(hat_y == y_t && (sum[y_t]-sum[i_t])<1)
    {
      for(h in c(y_t,i_t))
      {
        svpara  <- sv_coe_list[[h]]
        Inver_K <- InvGram_list[[h]]
        kt      <- kt_list[[h]]
        
        #### compute delta
        d_ast  <-Inver_K%*%kt
        tem_t  <- crossprod(d_ast,kt)[1,1]
        delta  <- 1-tem_t
        if(abs(delta)<1e-8)
          delta <- 0
        if((sum[y_t]-sum[i_t]) < (1-sqrt(delta/eta)))
        {
          tau_t1  <- (1-(sum[y_t]-sum[i_t]))/tem_t
          tau_t2  <- 2*(1-(sum[y_t]-sum[i_t])-sqrt(delta/eta))/tem_t
          tau_t   <- min(tau_t1,tau_t2,1)
          if(h == y_t)
            svpara <- svpara + tau_t*as.vector(d_ast)
          else
            svpara <- svpara - tau_t*as.vector(d_ast)
          sv_coe_list[[h]]  <- svpara
        }
      }
    }
  }
  t2 <- proc.time()
  runtime[re]   <- (t2 - t1)[3]
  errorrate[re] <- error/length_tr
  All_bud[re]   <- sum(k)
}

save_result <- list(
  note     = c("the next term are:alg_name--dataname--sam_num--sigma--sv_num--run_time--err_num--tot_run_time--ave_run_time--ave_err_rate--sd_time--sd_err"),
  alg_name = c("M-Projectron++-"),
  dataname = paste0(Dataset[d_index], ".train"),
  eta  = eta,
  ker_para = sigma,
  sv_num   = sum(All_bud)/re,
  run_time = as.character(runtime),
  err_num = errorrate,
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time  <- sd(runtime),
  sd_err    <-sd(errorrate)
)

write.table(save_result,file=savepath1,row.names =TRUE, col.names =FALSE, quote = T) 

sprintf("the candidate kernel parameter are :")
sprintf("%.5f", sigma)
sprintf("the number of sample is %d", length_tr)
sprintf("the number of support vectors is %d", round(sum(All_bud)/re))
sprintf("total training time is %.4f in dataset", sum(runtime))
sprintf("average training time is %.5f in dataset", sum(runtime)/reptimes)
sprintf("the average error rate is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of error is %.5f in dataset", sd(errorrate))


