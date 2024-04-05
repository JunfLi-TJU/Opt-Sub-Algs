setwd("D:/experiment/JMLR/JMLR2024/code")
rm(list = ls())

d_index <- 10

dpath           <- file.path("D:/experiment/online learning dataset/multi-class") 

Dataset         <- c("usps_all", "protein_all", "mnist_test","cifar10","letter","poker","shuttle_all",
                   "segmentation","waveform","covtype")

savepath1       <- paste0("D:/experiment/JMLR/JMLR2024/Result/",paste0("MAVP-OGS-",Dataset[d_index],".txt"))

traindatapath   <- file.path(dpath, paste0(Dataset[d_index], ".train"))

traindatamatrix <- as.matrix(read.table(traindatapath))                       
trdata          <- traindatamatrix[ ,-1]
ylabel          <- traindatamatrix[ ,1]                                        
if(min(ylabel)<= 0)
  ylabel <- -min(ylabel)+ylabel+1

length_tr       <- nrow(trdata)                                               
feature_tr      <- ncol(trdata)              

##############################################################################

sigma     <- 2^11
U         <- 10
alpha     <- 500/length_tr^(1)
#alpha     <- 0.01
M         <- max(ylabel)

epsilon   <- 0.9
coe       <- 0.1
reptimes  <- 10   
runtime   <- c(rep(0, reptimes))
errorrate <- c(rep(0, reptimes))
All_bud   <- c(rep(0, reptimes))

for( re in 1:reptimes)
{
  order      <- sample(1:length_tr,length_tr,replace = F)   #dis
  lambda_t   <- coe*sqrt(M)*U/sqrt(M*U^2+1)
  k          <- 1
  sum        <- c(rep(0, M))
  norm_f     <- c(rep(lambda_t, M))
  cum_delta  <- 1
  
  sv_coe_list <- list(
    svpara1   = array(-lambda_t,1),
    svpara2   = array(-lambda_t,1),
    svpara3   = array(-lambda_t,1),
    svpara4   = array(-lambda_t,1),
    svpara5   = array(-lambda_t,1),
    svpara6   = array(-lambda_t,1),
    svpara7   = array(-lambda_t,1),
    svpara8   = array(-lambda_t,1),
    svpara9   = array(-lambda_t,1),
    svpara10   = array(-lambda_t,1),
    svpara11   = array(-lambda_t,1),
    svpara12   = array(-lambda_t,1),
    svpara13   = array(-lambda_t,1),
    svpara14   = array(-lambda_t,1),
    svpara15   = array(-lambda_t,1),
    svpara16   = array(-lambda_t,1),
    svpara17   = array(-lambda_t,1),
    svpara18   = array(-lambda_t,1),
    svpara19   = array(-lambda_t,1),
    svpara20   = array(-lambda_t,1),
    svpara21   = array(-lambda_t,1),
    svpara22   = array(-lambda_t,1),
    svpara23   = array(-lambda_t,1),
    svpara24   = array(-lambda_t,1),
    svpara25   = array(-lambda_t,1),
    svpara26   = array(-lambda_t,1)
  )
  
  svmat       <- matrix(0,nrow = feature_tr,ncol=1) 
  Gram        <- matrix(1,nrow = 1,ncol=1)           # The inverse kernel matrix
  Inver_K     <- matrix(1,nrow = 1,ncol=1)           # The inverse kernel matrix
  beta_ast    <- array(0,1)                          # The optimal parameter d
  
  t1          <- proc.time()  #proc.time()
  
  ### the first instance
  error      <- 1
  svmat[,1]  <- trdata[order[1], ]
  y_t        <- ylabel[order[1]]
  sv_coe_list[[y_t]][1]    <- lambda_t
  
  ### from the second instance
  for(t in 2:length_tr)
  {
    diff   <- svmat - trdata[order[t],]
    kt     <- exp(colSums(diff*diff)/(-2*sigma^2))
    for(h in 1:M)
    {
      svpara <- sv_coe_list[[h]]
      sum[h] <- crossprod(svpara[1:k],kt)[1,1]
    }
    
    hat_y <- which(sum==max(sum))[1]
    y_t   <- ylabel[order[t]]
#    if(hat_y != y_t)
#    {
#      error = error+1
#    }
    tem_sum      <- sum
    tem_sum[y_t] <- min(sum)-0.1
    i_t          <- which(tem_sum==max(tem_sum))[1]
    
    if((sum[y_t]-sum[i_t])<1-epsilon)
    {
      #### compute alpha_t
      beta_ast   <- Inver_K%*%kt
      alpha_t    <- 1-crossprod(beta_ast,kt)[1,1]
      if(alpha_t<0)
        alpha_t  <- 0
      if(alpha_t <= alpha)
      {
        tem         <- Gram %*% beta_ast
        tem1        <- crossprod(tem,beta_ast)[1,1]
        if(hat_y != y_t)
        {
          error = error+1
          cum_delta    <- cum_delta + tem1
          lambda_t     <- coe*sqrt(M)*U/sqrt(M*U^2+cum_delta)
        }
        ########################
        h <- y_t 
        svpara         <- sv_coe_list[[h]]
        tem3           <- crossprod(tem,svpara)[1,1]
        svpara         <- svpara + lambda_t*as.vector(beta_ast)
        norm_f[h]      <- sqrt(norm_f[h]^2+lambda_t^2*tem1+2*lambda_t*tem3)
        if(norm_f[h]>U)
        {
          svpara    <- svpara*U/norm_f[h]
          norm_f[h] <- U
        }
        sv_coe_list[[h]]  <- svpara 
        
        h <- i_t
        svpara        <- sv_coe_list[[h]]
        tem3          <- crossprod(tem,svpara)[1,1]
        svpara        <- svpara - lambda_t*as.vector(beta_ast) 
        norm_f[h]     <- sqrt(norm_f[h]^2+lambda_t^2*tem1-2*lambda_t*tem3)
        if(norm_f[h]>U)
        {
          svpara      <- svpara*U/norm_f[h]
          norm_f[h]   <- U
        }
        sv_coe_list[[h]]  <- svpara
      }else
      {
        k           <- k + 1
        svmat       <- cbind(svmat,trdata[order[t],])
        
        Gram        <- cbind(Gram,as.vector(kt))
        Gram        <- rbind(Gram,c(kt,1))
        
        tem_d       <- beta_ast
        tem_d[k]    <- -1
        incre       <- tem_d %*% t(tem_d)/alpha_t
        incre[1:(k-1),1:(k-1)] <- incre[1:(k-1),1:(k-1)]+Inver_K
        Inver_K     <- incre

        if(hat_y != y_t)
        {
          error = error+1
          cum_delta    <- cum_delta + 1
          lambda_t     <- coe*sqrt(M)*U/sqrt(M*U^2+cum_delta)
        }        
        ########################
        h            <- y_t 
        svpara       <- sv_coe_list[[h]]
        svpara[k]    <- lambda_t
        norm_f[h]    <- sqrt(norm_f[h]^2+lambda_t^2+2*lambda_t*sum[h])
        if(norm_f[h]>U)
        {
          svpara     <- svpara*U/norm_f[h]
          norm_f[h]  <- U
        }
        sv_coe_list[[h]]  <- svpara 
        
        h            <- i_t
        svpara       <- sv_coe_list[[h]]
        svpara[k]    <- -lambda_t 
        norm_f[h]    <- sqrt(norm_f[h]^2+lambda_t^2-2*lambda_t*sum[h])
        if(norm_f[h]>U)
        {
          svpara     <- svpara*U/norm_f[h]
          norm_f[h]  <- U
        }
        sv_coe_list[[h]]    <- svpara 
        
        subset              <- setdiff(c(1:M),c(y_t,i_t))
        for( h in subset)
        {
          svpara            <- sv_coe_list[[h]]
          svpara[k]         <- 0
          sv_coe_list[[h]]  <- svpara
        }
      }
    }
  }
  t2            <- proc.time()
  runtime[re]   <- (t2 - t1)[3]
  errorrate[re] <- error/length_tr
  All_bud[re]   <- k
}

save_result <- list(
  note     = c("the next term are:alg_name--dataname--sam_num--sigma--sv_num--run_time--err_num--tot_run_time--ave_run_time--ave_err_rate--sd_time--sd_err"),
  alg_name = c("MAVP-OGS"),
  dataname = paste0(Dataset[d_index], ".train"),
  U        = U,
  Coe      = coe,
  varepsilon   = epsilon,
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
