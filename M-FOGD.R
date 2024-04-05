setwd("D:/experiment/JMLR/JMLR2024/code")
rm(list = ls())
library(MASS)

dpath         <- file.path("D:/experiment/online learning dataset/multi-class")  

d_index       <- 5

Dataset       <- c("usps_all", "protein_all", "mnist_test","cifar10","letter","poker","shuttle_all",
                   "segmentation","waveform","covtype")

savepath      <- paste0("D:/experiment/JMLR/JMLR2024/Result/",
                        paste0("M-FOGD-all-",Dataset[d_index],".txt"))

traindatapath    <- file.path(dpath, paste0(Dataset[d_index], ".train"))                
traindatamatrix  <- as.matrix(read.table(traindatapath))
trdata           <- traindatamatrix[ ,-1]
ylabel           <- traindatamatrix[ ,1]
if(min(ylabel)<= 0)
  ylabel <- -min(ylabel)+ylabel+1

length_tr        <- nrow(trdata)    
feature_tr       <- ncol(trdata)  

p_setting <-list(
  gamma = 2^3,
  d     = feature_tr,
  D     = 8896,
  zx    = array(0,dim=8896*2),
  eta   = 100/sqrt(length_tr)
) 

reptimes <- 10
M        <- max(ylabel)

runtime   <- c(rep(0, reptimes))
errorrate <- c(rep(0, reptimes))


for(re in 1:reptimes)
{
  
  order <- sample(1:length_tr,length_tr,replace = F)   #dis
  error <- 0
  sigma <- (1/p_setting$gamma)^2 *diag(p_setting$d)
  u     <- mvrnorm(p_setting$D,rep(0,p_setting$d),sigma)   # w--->D*d
  w     <- matrix(0,nrow=M,ncol= 2*p_setting$D)
  t1    <- proc.time()                                     #proc.time()
  
  for (i in 1:length_tr)
  {
#    p_setting$eta <- 2/sqrt(i)
    tem = u%*%trdata[order[i],]
    coszx <- cos(tem)/sqrt(p_setting$D)
    sinzx <- sin(tem)/sqrt(p_setting$D)
    p_setting$zx <- c(coszx,sinzx)
    
    sum   <- (w%*%p_setting$zx)[,1]
    hat_y <- which(sum==max(sum))[1]

    y_t   <- ylabel[order[i]]
    
    if(hat_y != y_t)
    {
      error = error+1
    }
    tem_sum      <- sum
    tem_sum[y_t] <- min(sum)-0.1
    i_t          <- which(tem_sum==max(tem_sum))[1]

    epsilon <- 0
    if((sum[y_t]-sum[i_t])<1-epsilon)
    {
      w[y_t,] <- w[y_t,] + p_setting$eta*p_setting$zx
      w[i_t,] <- w[i_t,] - p_setting$eta*p_setting$zx
    }
  }
  t2 <- proc.time()
  runtime[re] <- (t2 - t1)[3]
  errorrate[re] <- error/length_tr
}

save_result <- list(
  note     = c(" the next term are:alg_name--dataname--ker_para--sam_num--RFF--run_time--tot_run_time--ave_run_time--err_num--all_err_rate--ave_err_rate--sd_time--sd_err"),
  alg_name = c("M-FOGD"),
  dataname = paste0(Dataset[d_index], ".train"),
  ker_para = p_setting$gamma,
  sam_num  = length_tr,
  eta      = p_setting$eta,
  rff_num  = p_setting$D,
  run_time = as.character(runtime),
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  err_num  = errorrate,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time  <- sd(runtime),
  sd_err    <-sd(errorrate)
)
write.table(save_result,file=savepath,row.names =TRUE, col.names =FALSE, quote = T)

sprintf("the kernel parameter is %f", p_setting$gamma)
sprintf("the number of sample is %d", length_tr)
sprintf("the number of random features are %d", p_setting$D)
sprintf("total running time is %.1f in dataset", sum(runtime))
sprintf("average running time is %.1f in dataset", sum(runtime)/reptimes)
sprintf("the average AMR is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of AMR is %.5f in dataset", sd(errorrate))
###############################################################################################
