rm(list=ls())

###############################################################################
##FUNCTION:donor.deciles()
##*****************************************************************************
##PURPOSE: CALCULATES THE PERCENTAGE OF FUNDS RAISED FROM DONORS IN SPECIFIED
##         QUANTILES
##ARGUMENTS: ip,cm.mat
##OUTPUT: MATRIX THAT SUMMARIZES FUNDRAISING FROM QUANTILES FOR EACH CANDIDATE
###############################################################################
donor.deciles <- function(ip,cm.mat,nbin =10,ncross=12,challs=FALSE){
  
  
  ##SET CHALLENGER IPS TO NA
  ip[grepl('chall_',colnames(cm.mat))] <- NA
  
  set.seed(1234)
  sip <- sample(1:length(ip),length(ip))
  vvect <- round(quantile(1:length(ip),pr=seq(0,1,1/(ncross+1))))
  
  XX <- ceiling(cm.mat/1e10)
  keep <- (rowSums(XX[,!is.na(ip)]) > 0)
  cm.mat <- cm.mat[keep,]
  XX <- ceiling(cm.mat/1e10)
  ps <- Sys.time()
  c.deciles <- matrix(NA,length(ip),nbin)
  c.ips <- rep(NA,length(ip))
  for(xx in 2:(ncross+2)){
    cat(xx-1,'...')
    tuse <- sip[(vvect[xx-1] +1):vvect[xx]]
    ip.use <- ip
    ip.use[tuse] <- NA
    tuse2 <- (gsub('chall_','',colnames(cm.all)) %in% gsub('chall_','',colnames(cm.all)[tuse]))
    ip.use[tuse2] <- NA
    ##
    MM <- XX[,!is.na(ip.use)]
    keep <- (rowSums(MM) > 0)
    CC <- cm.mat[keep,]
    contrib.ips <- get.contrib.means(CC,weights=NULL,ip.use)
    new.cand.means <- get.cand.means(CC,contrib.ips)
    ord <- order(contrib.ips)
    qq1 <- quantile(contrib.ips,pr=unique(seq(0,1,1/nbin)))
    if(sum(duplicated(qq1))>0){
      jit <- (runif(length(contrib.ips))-.5)*(sd(contrib.ips)) * .001
      qq1 <- quantile(contrib.ips+jit,pr=unique(seq(0,1,1/nbin))) 
      ##qq1[which(duplicated(qq1))] <- (qq1[which(duplicated(qq1))-1] + qq1[which(duplicated(qq1))+1])/2
    }
    oo <- cut(contrib.ips,
              breaks=qq1,
              labels = FALSE,
              include.lowest=TRUE)
    cand.quantiles <- NULL
    csums <- colSums(CC)
    for(x in 1:max(oo)){
      s1 <- colSums(CC[oo==x,]) / csums
      cand.quantiles <- cbind(cand.quantiles,s1)
    }
    c.deciles[tuse,] <- cand.quantiles[tuse,]
    c.ips[tuse] <- new.cand.means[tuse]
  }
  
  if(challs){
    tuse.ch <- grep('chall_',colnames(cm.all))
    c.ips[tuse.ch] <- NA
    set.seed(1234)
    sip.chall <- sample(1:length((tuse.ch)),length((tuse.ch)))
    vvect.chall <- round(quantile(1:length((tuse.ch)),pr=seq(0,1,1/(ncross+1))))
    for(xx in 2:(ncross+2)){
      tuse <- sip.chall[(vvect.chall[xx-1] +1):vvect.chall[xx]]
      tuse2 <- (gsub('chall_','',colnames(cm.all)) %in%
                  gsub('chall_','',colnames(cm.all)[tuse.ch][tuse]))
      ip.use <- ip
      ip.use[tuse.ch[tuse]] <- NA
      ip.use[tuse2] <- NA
      
      MM <- XX[,!is.na(ip.use)]
      keep <- (rowSums(MM) > 0)
      CC <- cm.mat[keep,]
      contrib.ips <- get.contrib.means(CC,weights=NULL,ip.use)
      new.cand.means <- get.cand.means(CC,contrib.ips)
      
      
      ord <- order(contrib.ips)
      qq1 <- quantile(contrib.ips,pr=unique(seq(0,1,1/nbin)))
      qq1[which(duplicated(qq1))] <- (qq1[which(duplicated(qq1))-1] + qq1[which(duplicated(qq1))+1])/2
      oo <- cut(contrib.ips,
                breaks=qq1,
                labels = FALSE,include.lowest=TRUE)
      cand.quantiles <- NULL
      csums <- colSums(CC)
      for(x in 1:max(oo)){
        s1 <- colSums(CC[oo==x,]) / csums
        cand.quantiles <- cbind(cand.quantiles,s1)
      }
      c.deciles[tuse.ch[tuse],] <- cand.quantiles[tuse.ch[tuse],]
      c.ips[tuse.ch[tuse]] <- new.cand.means[tuse.ch[tuse]]
    }
  }
  ps1 <- Sys.time()
  ps1-ps
  cat('\n')
  return(list(c.deciles=c.deciles,c.ips=c.ips))
}

##CALCULATE CONTRIBUTOR MEANS
get.contrib.means <- function(cm,cand.ips,contrib.ips=NULL,weights=NULL,upper.limit=2000,
                              get.sds=FALSE,cores=1){
  numconts <- nrow(cm)
  use.weights <- !is.null(weights)
  use <- !is.na(as.numeric(cand.ips))
  t1 <- cm[,use] %*% (as.numeric(cand.ips[use]))
  t2 <- rowSums(cm[,use])
  contrib.ips <- as.numeric(t1/t2)
  return(contrib.ips)
}

##CALCULATE CANDIDATE MEANS
get.cand.means <- function(cm,contrib.ips,weights=NULL,dynamic.cands=FALSE,
                           upper.limit=30,cores=1,get.sds=F){
  numconts <- ncol(cm)
  if(numconts %% upper.limit == 0){upper.limit <- upper.limit + 2}
  if(numconts %% upper.limit == 1){upper.limit <- upper.limit + 1}
  if(is.null(weights)){weights=rep(1,length(contrib.ips))}
  if(length(weights)>1){weights <- weights[!is.na(contrib.ips)]}
  cm <- cm[!is.na(contrib.ips),]
  contrib.ips <- contrib.ips[!is.na(contrib.ips)]
  count <- 1
  
  if(get.sds){
    if(cores > 1){
      print(cores)
      require('snow')
      require('doSNOW',quietly=TRUE,warn.conflicts=FALSE)
      c3 <- makeCluster(cores, type = "SOCK")
      registerDoSNOW(c3)
    }
    calc.c.means <- function(x){
      cmean <- sum(x*contrib.ips*weights,na.rm=T)/sum(x*weights,na.rm=T)
      out <- sqrt(sum(x*((cmean-contrib.ips)^2),na.rm=T)/sum(x,na.rm=T))
      out
    }
    cind <- c(seq(upper.limit,numconts,upper.limit))
    cand.ips <- foreach(i = cind,.combine='c',
                        .packages=c('Matrix')) %dopar% {
                          out <- apply(cm[,(i-(upper.limit-1)):i],2,calc.c.means)
                        }
    if(ncol(cm) - max(cind) > 0){
      if(ncol(cm) - max(cind) > 1){
        cips2 <- apply(cm[,(max(cind)+1):ncol(cm)],2,calc.c.means)
      }else{
        cips2 <- calc.means(cm[,ncol(cm)])
      }
      cand.ips <- c(cand.ips,cips2)
    }
    if(cores > 1){stopCluster(c3)}
  }else{
    cand.ips <- rep(NA,numconts)
    t1 <- t(cm) %*% (as.numeric(contrib.ips))
    t2 <- colSums(cm)
    cand.ips <- as.numeric(t1/t2)
  }
  return(cand.ips)
}

list2env(readRDS('bonica_data.RDS'), globalenv())

# From Bonica:
fake_y <- rep(NA, nrow(x.test))
y.in <- y.train[
  match(colnames(cm.all),
        c(rownames(x.train), rownames(x.test)))
]
# The "deciles" dataset
cdec <- donor.deciles(y.in, cm.all, ncross=10) 
colnames(cdec$c.deciles) <- paste('cdec_decile_',1:ncol(cdec$c.deciles),sep='')
cdec.in <- cbind(cdec_ips=cdec$c.ips,cdec$c.deciles)
mc.train <- match(rownames(x.train),colnames(cm.all))
mc.test <- match(rownames(x.test),colnames(cm.all))

cdec.train <- cdec.in[mc.train,]
cdec.test <- cdec.in[mc.test,]

fit_lm <- lm(y.train ~ as.matrix(x.train[,1:28]))
library(glmnet)
fit_LASSO <- cv.glmnet(x = x.train, y = y.train)

pred2 <- as.vector(as.matrix(cbind(1,x.test[,1:28])) %*% coef(fit_lm))
pred3 <- as.vector(predict(fit_LASSO, newx = x.test))

# Prepare your data as follows. 
# Give your models short informative names, e.g. "model_LASSO".
# You may hand in up to 10 models
# You *must* beat a naive linear model *and* LASSO
suppressPackageStartupMessages(library(caret))

# model 1: ranger
library(SuperLearner)
library(tidyverse)

x.df <- as.data.frame(as.matrix(x.train))
x.sub <- x.df[,1:28] 
x.sub <- x.sub %>% 
  rename('poly1'='poly(cf,4)1',
         'poly2'='poly(cf,4)2',
         'poly3'='poly(cf,4)3',
         'poly4'='poly(cf,4)4')

ranger <- SuperLearner(Y = y.train,
                         X = x.sub %>% select(poly1, poly2,
                                              tiactive_1980, pvs,
                                              gender_indicator,
                                              sen.cand, pres.cand,
                                              rep_indicator, dem_indicator),
                         SL.library = 'SL.ranger')
saveRDS(ranger, 'ranger.RDS')

# prediction
x.df.test <- as.data.frame(as.matrix(x.test))
x.sub.test <- x.df.test[,1:28] 
x.sub.test <- x.sub.test %>% 
  rename('poly1'='poly(cf,4)1',
         'poly2'='poly(cf,4)2',
         'poly3'='poly(cf,4)3',
         'poly4'='poly(cf,4)4')

ranger_pred <- predict(ranger, newdata=x.sub.test)

# outcome
output <- data.frame(cand_id = rownames(x.sub.test), 
                     model1 = mean(y.train),
                     model2 = ranger_pred,
                     stringsAsFactors = F)

saveRDS(output, 'model_ranger_Orly & Tong.RDS')







fit_glmnet <- SuperLearner(Y = y.train,
                       X = x.sub %>% select(poly1, poly2,poly3, poly4, 
                                            tiactive_1980, pvs,
                                            gender_indicator,
                                            sen.cand, pres.cand,
                                            rep_indicator, dem_indicator),
                       SL.library = 'SL.glmnet')

saveRDS(ranger, 'ranger.RDS')

# model 2: lmer

library(lme4)

fit_lmer <- lmer(y.train ~ (1 | x.sub$gender_indicator) +
                   (1 | x.sub$rep_indicator))

# model 3: random forest

library(caret)
alldata <- cbind(x.sub, y.train)
modelrf <- train(
  y.train ~ .,
  data = alldata,
  method = 'rf'
)

