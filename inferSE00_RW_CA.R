# Key functions for inferring self-esteem by evidence accummulation and mapping
# beliefs to responses using a sigmoid function
# New model - as of discussion on 3rd September 2018
#
# Alexis (An Yee) Low & Michael Moutoussis, 2019
# remove(list=ls()); # drastic clearup!
try(remove(whoami,codeDir,baseDir,seDirs,bestsofar10)) # conservative clearup

# ---------------   Cluster related menu / argument options ------------------
runOnCluster = 1;   # 0 for single off cluster run, 1 for cluster run,
# 2 for cluster-like run (with Rscript and ptN argument, 
# e.g.$ Rscript inferSE04_9p.R 30) but off cluster.
Debug = 0;  if (Debug & !(runOnCluster)){ participant = 25};                 

if (runOnCluster ==0){
  # Then 'participant' must already be in the workspace! 
  print(paste('participant = ',participant)) # participant <- 1;  
} else {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args)!=1) {
    stop("Wrong arguments supplied - you must supply just one, for pt number")
  } else {
    participant <- as.numeric(args[1])
  }
}
## ------------------------------------------------------------------------- 

# ---- Now find out from where this source code and so set paths -----
# Set the base and code directories consistently for different
# users of the main code, LikeMe.R etc.
sewd <- getwd()
# find out where we are
# Depending on what the path contains, decide who is the user. Student PLS EDIT YOUR ENTRY:
if (grepl('Alexis', sewd)) {  whoami <- 'Student' }
if (grepl('/home/hopper', sewd)) {  whoami <- 'WillLinux' }
# if (grepl('C:/Users/mmoutou', sewd)) {  whoami <- 'SpectreMM' }
if (grepl('michael', sewd)) {   whoami <- 'LinuxMM' }
if (grepl('geert-jan', sewd)) {  whoami <- 'GeertJanMac' }
# if (grepl('C_mmoutou', sewd)) {   whoami <- 'SpectreVM' }
if (grepl('D:/mmoutou', sewd)) {  whoami <- 'cedricD' }
if (grepl('sf_',sewd)){ 
  setwd('/media/sf_CHERhomeYdrive/'); 
  print('NB: We are in kikidi vm: switching to CHERhomeYdrive');
}
if (grepl('sf_CHERhomeYdrive', sewd)) {  whoami <- 'kikidi' }
if ((runOnCluster ==1)){   whoami <- 'MMhal9000';  }

# Adjust the base directdory accordingly.  Student PLS EDIT YOUR ENTRY :
switch(
  whoami,
  Student  = {
    baseDir <-
      "X:/OneDrive - University College London/Summer Project - Alexis An Yee Low/"
  },
  WillLinux = {
    baseDir <- "/home/hopper/Dropbox/SelfEvalMEV"
  },
  SpectreMM = {
    baseDir <-
      "C:/Users/mmoutou/OneDrive - University College London/SharePoint/Low, An Yee/Summer Project - Alexis An Yee Low/"
  },
  SpectreVM = {
    baseDir <-
      "/media/michael/C_mmoutou/OneDrive/SharePoint/Low, An Yee/Summer Project - Alexis An Yee Low/"
  },
  LinuxMM = {
    baseDir <- "~/Dropbox/BASOR/BASOR_output/AYLwork/"    # "/home/michael/gitwork/LikeMe/"
  },
  GeertJanMac = {
    baseDir <- "/Users/geert-janwill/Dropbox/GJW_LikeMe/"
  },
  kikidi = {
    baseDir <- "/media/sf_CHERhomeYdrive/R/R4hal/"
  },
  cedricD = {
    baseDir <-
      "D:/mmoutou/Dropbox/BASOR/BASOR_output/AYLwork/"
  },
  MMhal9000 = {
    baseDir <- "~/R/R4hal/"
  }
)

outDir <- baseDir;
codeDir <- paste(baseDir, "likeme-Socio3/", sep = '')
if (whoami == 'cedricD'){ 
  codeDir <- "D:/mmoutou/Dropbox/FIL_aux/R_scripts/";
  outDir <- "D:/mmoutou/Dropbox/BASOR/BASOR_output/AYLwork/" };
if (whoami == 'kikidi'){ 
  codeDir <- '/media/sf_CHERhomeYdrive/R/R4hal/scripts/';
  outDir <- '/media/sf_CHERhomeYdrive/R/R4hal/currentWork/socio4/' };
if (whoami == 'LinuxMM'){ 
  codeDir <- '~/Dropbox/FIL_aux/R_scripts/';
  outDir <- baseDir };
if (whoami == 'MMhal9000'){
  codeDir <- '~/R/R4hal/scripts/';
  outDir <- '~/R/R4hal/currentWork/socio4/' };


seDirs = list(baseDir=baseDir,codeDir=codeDir,outDir=outDir,whoami=whoami)
# --------------------------------------------------------------------------

# LikeMe.R has most of the functions that for this project. Whoever it
# may change the working directories we use here, so then we restored them.
print(getwd())
print(codeDir)
source(paste(codeDir, 'LikeMe.R', sep = ''))
source(paste(codeDir, 'gen_ut.R', sep = ''))
baseDir=seDirs$baseDir; codeDir=seDirs$codeDir; outDir=seDirs$outDir; whoami=seDirs$whoami;


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SLPsocioCA_t <- function(parMat, datAr, onlySLP=0, fixAcc=NA, check=1){
  # sum-log-posterior function for self-esteem - edited by William Hopper - Competence-Acceptance Model
  # to analyse data by Geert-Jan Will et al. 
  #
  # parMat has rows with the parameters for each pt.,
  # par. row: SEb, acc0max, acc0min, etapos, gam, wexp, wrpe, Tau, sig
  # datAr has a page for each pt, and Ntr rows.
  #
  # For testing & example use see roughLikeMe.R, wherein 
  # can use parMat = par0a; datAr = D; onlySLP=0; fixAcc=gpAcc; check = 1; 
  # and at the end remove(parMat,datAr,onlySLP,fixAcc,check,hapPol)
  
  M <- 4; # Number of rater groups.
  knmax <- 200; # 10;  # max. depth of kernel summation. A: ?
  gamNorm = 0;  # 0 ;= don't normalize the decay kernel A: ?
  
  if (check){
    datDim <- dim(datAr);
    if (length(datDim)<3){   # if 2D object was provided, convert to 3D.
      hd <- colnames(datAr);
      datAr <- array(as.matrix(datAr),c(datDim,1));
      colnames(datAr) <- hd;
    }
    if (is.null(dim(parMat))) {
      parMat <- array(parMat,c(1,length(parMat)));
    }
  }
  # To check and add accHalf as last entry of parMat, so that 50% response becomes:
  # predictApp <- 1/(1+exp((parMat[ptN,'accHalf']-hapPol[trN,P_aInd[gpI],ptN])/
  #                                    (parMat[ptN,'accHalf']*parMat[ptN,'Tau'])));
  colnames(parMat) <- c('SEb','acc0max','acc0min','eta','gam','wexp','wrpe','wa','Tau','SSB','sig');
  Nptot <- dim(datAr)[3];
  
  Ntrtot = dim(datAr)[1];  # this will sadly have to be fiddled with later, usually to
  
  # Cater somehow for the fact that scanned task contains interruptions.
  bakD <- datAr;     # backup copy to restore later ...
  stimTr <- array(NA,c(Ntrtot,Nptot));  # will store positions of stimulus-endowed trials.
  
  # Array to hold expectations for each group of raters, RPE,
  # as well as the probability / prob. dens. for the responses emitted,
  # and generated acceptance prediction and SE. So each completed row
  # refers to the values at the END OF THAT TRIAL, i.e. includes the
  # RPE experienced on the basis of the EV JUST BEFORE the trial but the
  # updated EVs AFTER updating with the observation.
  # Note has one more row than trials as starts from 'priors'.
  hapPol <- array(dim=c(Ntrtot+1,M+7,Nptot)); 
  dimnames(hapPol)[[2]] <- c(paste('exp',1:M,sep=''), 'RPE_A','RPE_C',
                             'obsP','SEPD','accP','genSE','expSE');
  P_aInd <- 1:M    # Indices of group expectations in hapPol 
  ArpeInd <- M+1;   # index for  RPE_Acceptance, 
  CrpeInd <- ArpeInd+1 # index for RPE_Competence
  obsPI <- CrpeInd+1;       SEPDI  <- obsPI+1 ;
  accPI <- SEPDI+1;       genSEI <- accPI+1;   expSEI <- genSEI+1;
  
  for (ptN in 1:Nptot){
    
    # Cater somehow for the fact that scanned task contains interruptions
    # including at the very start. Awkward - to start with, just exclude:
    stimTr[,ptN] <- !is.na(datAr[,'pred',ptN]); #altered GJW #A: if there is a prediction indicated, a trial happened, so put 'true'. result - column of TRUEs
    Ntrtot <- sum(stimTr[,ptN]);    # update to exclude non-stimulus trials. Was dim(datAr)[1] above.
    datAr[1:Ntrtot,,ptN] <- datAr[stimTr[,ptN],,ptN];  # change datAr to only include stimulus trials
    if (Ntrtot<dim(datAr)[1]) {
      datAr[(Ntrtot+1):dim(datAr)[1],,ptN] <- NA;      # set all non-stimulus trials to NA
    }
    
    # Initial beliefs about ratings by others (groups):
    if (!is.na(sum(fixAcc))){ #AL fixAcc is something provided beforehand... but what?
      parMat[ptN,'eta'] <- 0.0; #A: set learning rate to 0
      # hapPol[1,P_aInd,ptN] <- fixAcc[P_aInd];   # this will hopefully fail if rubbish fixAcc is given.
      # Scaled a la Geert-Jan: 
      hapPol[1,P_aInd,ptN] <- 2*fixAcc[P_aInd]-1; #A: set initial expectations for each group  
      # this will hopefully fail if rubbish fixAcc is given.
    } else {
      accP0 <-  ((M-1):0)*(parMat[ptN,'acc0max']-parMat[ptN,'acc0min'])/(M-1) + parMat[ptN,'acc0min'];
      hapPol[1,P_aInd,ptN] <- accP0;  # Apply same to SSB below. #A: now expected acceptance prob has been calculated.
      # no longer scaled between 1,-1 - now 1,0
    }
    
    # Stuff to apply to all trials:
    gam <- as.numeric(parMat[ptN,'gam']); #A: recall gam is forgetting factor
    # Poss. include a normalizing constant to gam, equal to 1/sum(kernel vector), or not ...#A: what
    # z = 1; #
    if (gamNorm) {
      z = 1/sum(c(1, cumprod(rep(gam, (sumN-1)))));
    } else {
      z = 1;
    }
    # ... to be included in the weights:  
    We = 0.0001 #We = z*parMat[ptN,'wexp'];    
    Wr = z*parMat[ptN,'wrpe']; #A: weight on expectations, weight on PE
    Wdiff <- We-Wr;    Wa = z*parMat[ptN,'wa']; #A: weight on acceptance?
    # Now for acceptance prediction response function parameters:
    Tau = parMat[ptN,'Tau'];      #   this and SSB below scaled to be directly
    SSB = parMat[ptN,'SSB'];    #   comparable with the values and returns ...
    #   no longer scaled between 1,-1 - now between 1,0
    
    eta <- parMat[ptN, 'eta'];
    # etapos <- parMat[ptN, 'etapos'];     # learning rate for positive SPE (50% higher than average eta)
    # etaneg <- parMat[ptN, 'etaneg'];     # learning rate for negative SPE (50% lower than average eta)
    SEb <- as.numeric(parMat[ptN,'SEb']);  # initialise SEb 
    SEs <- 0;                              # initialise SEs #A: SEs is state self esteem
    
    # loop over trials
    for (trN in 1:Ntrtot){
      hapPol[trN+1,,ptN] <- hapPol[trN,,ptN];   # Initialise-most will remain same.
      
      gpI <- datAr[trN,1,ptN];  # which group rated this pt
      
      P_a <- hapPol[trN,P_aInd[gpI],ptN];    # probability of approval
      
      # Prob. of 'accept' response emitted (this is BEFORE rating seen),
      # if present. NB we will use beliefs after last trial, i.e. hapPol[trN,...
      predictApp <- invlogit((SSB+P_a-1)/Tau); #A: why include the -1
      
      hapPol[trN+1,accPI,ptN] <- predictApp;  # Store      
      
      # Now calc. obsP, the probability under the current model that the
      # participant emits the acceptance prediction in datAr (the experimentally observed one).
      if (datAr[trN,2,ptN] == 0){  # rbinom gives 0 for disapproval, change 0 to -1 for consistency
        datAr[trN,2,ptN] <- -1
      }
      
      predAcc <- datAr[trN,2,ptN]; # rating that actual participant predicted.
      if (!is.na(predAcc)){    
        if (predAcc > 0.5){  # i.e. if it's 1 
          hapPol[trN+1,obsPI,ptN] <- predictApp;
        } else {
          hapPol[trN+1,obsPI,ptN] <- 1-predictApp;
        }
      } # End if valid prediction predAcc
      
      # RPE_Acceptance for this trial. Note that the EV is from the end
      # of the trial before ( hapPol[trN,... ) but the outcome is from this one,
      #  datAr[trN, ... :
      nofb = datAr[trN,'nofb',ptN];  # 0 if feedback given, 1 otherwise
      # so if no feedback given the outcome is said to be 0 acc. to GJW / RR model.
      if (nofb == 1){
        Outcome <- 0;    # if no feedback outcome is 0 -- COME BACK TO THIS re: MM's e-mail
      } else {
        Outcome <- diracd(1,datAr[trN,'obs',ptN]); # outcome = 1 (approval) or 0 (disapproval)
      }
      hapPol[trN+1,ArpeInd,ptN] <-  Outcome - P_a;
      
      # Competence
      # Expected Value of the Action taken
      if (!is.na(predAcc)){ 
        if (predAcc > 0.5){ # i.e. if it is 1 (approval)
          Q <- P_a
        } else {            # i.e. disapproval
          Q <- 1-P_a
        }
      }
      # Competence Prediction Error
      Competence <- diracd(predAcc,datAr[trN,'obs',ptN]);    # dirac delta(action, outcome)
      hapPol[trN+1,CrpeInd,ptN] <- Competence - Q; # RPE_competence
      
      # Now for the SE sum:
      # j loop removed and only state update tracked
      #SE <- SE + gam^(j-1)*
      # (  We* hapPol[trN-j+1,P_aInd[datAr[trN-j+1,'gp',ptN]],ptN] +
      #     Wr* hapPol[trN-j+2,ArpeInd,ptN] )
      
      #SEs <- gam*SEs + Wr*(Outcome*Wa + (1-Wa)*Competence) +
      #Wdiff*(Wa*P_a + (1-Wa)*Q);
      
      SEs <- gam*SEs + Wr*(Wa*hapPol[trN+1,ArpeInd,ptN] + (1-Wa)*hapPol[trN+1,CrpeInd,ptN]) + We*(Wa*P_a + (1-Wa)*Q);
      
      # SEs <-  gam*SEs + Wr*Approval + 
      #           Wdiff*hapPol[trN,P_aInd[datAr[trN,'gp',ptN]],ptN];
      
      SE <- SEb + SEs;
      
      # Store this SE as the current expected SE:
      hapPol[trN+1,expSEI,ptN] <- SE;
      
      # rpeSign <- sign(hapPol[[trN+1,rpeInd,ptN]])    # pos or neg RPE for eta selection
      
      # if loop below zapped as in the first instance NA group trials excluded already ...
      # if (!is.na(datAr[trN,'gp',ptN])){  # if valid group stimulus was presented
      
      
      # Determine density at the actual response:
      if (!is.na( datAr[trN, 'SE',ptN] )) {
        hapPol[trN+1,SEPDI,ptN] <- dnorm( datAr[trN, 'SE',ptN], SE, parMat[ptN,'sig'] );
      } else {
        hapPol[trN+1,SEPDI,ptN] <- NA;
      }
      
      # generate pseudodatum:
      hapPol[trN+1,genSEI,ptN] <- rnorm(1, SE, parMat[ptN,'sig'] );
      
      # Don't forget to update the EV for the gp we've just been rated by:
      # if (rpeSign > 0){       # if the RPE is +ve, use first learning rate
      
      
      hapPol[trN+1,P_aInd[gpI],ptN] <- (1-eta)*P_a + eta*Outcome;
      
      #  } else {                                         # else, if RPE is -ve, use second learning rate
      #     hapPol[trN+1,P_aInd[gpI],ptN] <- hapPol[trN+1,P_aInd[gpI],ptN] + 
      #      etaneg*hapPol[trN+1,rpeInd,ptN] ;} 
      
      #} else {  # if no group stimulus was presented, fill in with NAs. Note
      #          # that the group expectations are not touched - they just stay as they are.
      #    hapPol[trN+1,c(expSEI,SEPDI,genSEI),ptN] = c(NA,NA,NA) ; 
      # }
      
    } # end of loop over trials
    
    
  } # end loop over pts
  
  # Create objects for output
  SLP1 <- sum(log(na.omit(as.vector(hapPol[,obsPI,]))));
  SLP2 <- sum(log(na.omit(as.vector(hapPol[,SEPDI,]))));
  
  if (onlySLP ==1){          # sLP pertaining to all observations
    return( SLP1+SLP2 );
  } else {
    if (onlySLP ==2) {  # sLP pertaining to SE only
      return(SLP2) ; 
    } else {
      SLPetc <- list();
      SLPetc[[1]] <- SLP1;
      SLPetc[[2]] <- SLP2;
      
      datAr = bakD;  Ntrtot=dim(datAr)[1]; # restore in original form.
      # Next to combine both exp. data, beliefs, policies etc. 
      # Therefore has extra gp, pred, obs, SE, nofb and genPred cols :
      colN <- dim(hapPol)[2]+6;
      DatBelPol <- array(NA,c(Ntrtot+1,colN,Nptot));
      dimnames(DatBelPol)[[2]] <- c(colnames(datAr), colnames(hapPol), 'genPred')
      for (ptN in 1:Nptot) {
        offs = 1*is.na(datAr[1,'gp',ptN]);
        DatBelPol[(2-offs):(Ntrtot+1-offs),1:5,ptN] <- datAr[,,ptN];
        stTr = stimTr[,ptN];
        if (offs == 1){ 
          stTr[1] = TRUE; stTr=c(stTr,FALSE);
        } else {
          stTr=c(TRUE, stTr);  
        }
        DatBelPol[stTr,(5+1):(5+dim(hapPol)[2]),ptN] <- hapPol[1:sum(stTr),,ptN];
        for (trN in 1:Ntrtot) {
          if (!is.na(DatBelPol[trN+1,'accP',ptN])) {
            DatBelPol[trN+1,'genPred',ptN] <- rbinom(1,1,DatBelPol[trN+1,'accP',ptN]);
          }
        }
      }
      SLPetc[[3]] <- DatBelPol;
      # 4th element to have generated data only :
      SLPetc[[4]] <- array(NA,c(dim(DatBelPol)[1], 5, dim(DatBelPol)[3]));
      SLPetc[[4]][,,] <- DatBelPol[, c('gp','genPred','obs','genSE','nofb'),] ;
      colnames(SLPetc[[4]]) <- c('gp','pred','obs','SE','nofb'); # Just like real expt. data ...
      SLPetc[[5]] <- parMat;
      colnames(SLPetc[[5]]) <- c('SEb','acc0max','acc0min','eta','gam','wexp','wrpe','wa','Tau','SSB','sig');
      names(SLPetc) <-c('predSLnP','SESLnP','DatBelPol','genD','ptPar');
      
      
      return(SLPetc);
    }
  }
  
} # end of SLPsocioCA - 'Happiness-like model' (Robb's) (Edited by William Hopper)




# end of SLPsocioCA_t - 'Happiness-like model' (Robb's)


msLP0tr <- function(trParM, datAr, gamPri=NA, fixAcc=c(0.85,0.70,0.30,0.15), predRat=0, check=1){
  # trParM only to have atanh(2*SEb-1), log(gam,  wexp,  wrpe, Tresp,  sig)
  if (check){
    if (is.null(dim(trParM))){
      trParM <- matrix(trParM,nrow=1,byrow=TRUE)
    }
  }
  parM <- matrix(NA,nrow=dim(trParM)[1],ncol=3+dim(trParM)[2]);   # long form ...
  
  # fill in long form :
  if (dim(trParM)[2] == 7){  
    parM[,c(1,5:10)]  <- tr2natLP0(trParM);  # keep NAs for acc0max, acc0min, eta
  } else {
    if (dim(trParM)[2] == 11) { 
      parM  <- tr2natLP0(trParM);
    } else {
      print(trParM);
      stop('tr2natLP0 not ready for this transformed parameter matrix, dim(trParM) ');
    }
  }
  
  
  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  mSLPrior <- 0;
  if (length(gamPri)>1){  # legit prior must have 12 elements or so!
    for (ptN in 1:dim(trParM)[1]) {
      # in the line below na.omit bec. usually we don't have acc0max, acc0min, eta
      mSLPrior <- mSLPrior - sum(dgammaMS(na.omit(parM[ptN,]),   
                                          gamPri[1,],gamPri[2,], log=TRUE)); 
    }
  } 
  
  if (mSLPrior == Inf){  # If we are in an a priori prohibited parameter region
    # do not attempt to calculated the likelihood - it will be nonsense anyway.
    return(Inf); 
  } else {
    #if (predRat==1) {  # i.e. if we are to consider pts predictions of ratings
      return ( mSLPrior - SLPsocioCA_t( parM, datAr, onlySLP=1, fixAcc, check=1) );
    #} else {  # this is the default, for predRat=0 : only return sLP pertaining to SE, not ratings.
    #  return ( mSLPrior - SLPsocioCA_t( parM, datAr, onlySLP=2, fixAcc, check=1) );
    #}
  }
  
  # was:  mSLP <- - SLPsocioCA_t( parM, datAr, onlySLP=1, fixAcc, check=1)
  #  return(mSLP);
  
} # end of msLP0tr



# - - - - - - transform from native to fitting space for happinesque - - - - -
nat2trLP0 <- function(parM,check=1){
  if (check){
    if (is.null(dim(parM))){
      parM <- matrix(parM,nrow=1,byrow=TRUE)
    }
  }  
  if (dim(parM)[2] == 11){  # This means acc0max, acc0min, eta given.
    #  c('SEb','acc0max','acc0min','eta','gam','wexp','wrpe','wa','Tau','SSB','sig')
    # stop('tr2natLP0 not ready for non-ready-made general group acceptance probabilities');
    trParM <- matrix(NA,nrow=dim(parM)[1],ncol=dim(parM)[2]); 
    trParM[,c(1:4,8,10)] <- as.matrix(atanh(2*parM[,c(1:4,8,10)]-1));
    trParM[,c(5:7,9,11)] <- as.matrix(log(parM[,c(5:7,9,11)])) ;
    #trParM[,9] <- as.matrix(atanh(2*parM[,9]-1));
  } else {
    if (dim(parM)[2] == 7) {
      #  c('SEb', 'gam', 'wexp', 'wrpe', 'Tresp','Hresp','sig')
      trParM <- matrix(NA,nrow=dim(parM)[1],ncol=dim(parM)[2]); 
      trParM[,c(1,6)] <- as.matrix(atanh(2*parM[,c(1,6)]-1));
      trParM[,c(2:5,7)] <- as.matrix(log(parM[,c(2:5,7)])) ;
    } else {
      print(parM);
      stop('tr2natLP0 not ready for this param. vector length provided'); 
    }
  }
  return(trParM);
}
# - - - - - - transform from fitting space to native for happinesque - - - - -
tr2natLP0 <- function(trParM,check=1){
  if (check){
    if (is.null(dim(trParM))){
      trParM <- matrix(trParM,nrow=1,byrow=TRUE)
    }
  }
  if (dim(trParM)[2] == 11){  # This means acc0max, acc0min, eta given.
    parM <- matrix(NA,nrow=dim(trParM)[1],ncol=dim(trParM)[2]);
    colnames(parM) <- c('SEb','acc0max','acc0min','eta','gam','wexp','wrpe','wa','Tau','SSB','sig');
    parM[,c(1:4,8,10)] <- as.matrix(0.5*(tanh(trParM[,c(1:4,8,10)])+1));
    parM[,c(5:7,9,11)] <- as.matrix(exp(trParM[,c(5:7,9,11)])) ;
    
  } else {
    if (dim(trParM)[2] == 7) {
      parM <- matrix(NA,nrow=dim(trParM)[1],ncol=dim(trParM)[2]);
      colnames(parM) <- c('SEb','gam','wexp','wrpe','Tresp','Hresp','sig');
      parM[,c(1,6)] <- as.matrix(0.5*(tanh(trParM[,c(1,6)])+1));
      parM[,c(2:5,7)] <- as.matrix(exp(trParM[,c(2:5,7)])) ;
    } else {
      stop('tr2natLP0 not ready for this length of trParM');
    }
    
  }
  
  return(parM);
}


# ################################ start of nlm fit #################################

try({load(paste(baseDir,"loadfornlm_RW_set1.RData",sep='')); #change this depending on which set
  bestsofar10 <- importCSV(paste(baseDir,"soc00fit9par01to61_1.csv",sep=''))}) #change this #this contains bestsofar9, tryPmatrix, datArW61, and priors
try({load(paste(baseDir,"/currentWork/loadfornlm_RW_set1.RData",sep='')); 
  bestsofar10 <- importCSV(paste(baseDir,"/currentWork/soc04fit9par01to61_1.csv",sep=''))}) # for kikidi / hal ?

###^change this later

RWModelFit <- function(pts,Par0=NULL,nlmprintlev=0) {
  
  #priors <- priors[,2:10]
  
  #if (is.null(Par0)){ Par0 <- priors; }
  
  ml1fit <- list();
  ml1res <- matrix(NA,nrow=dim(datArW61)[3],ncol=16);
  dimnames(ml1res)[[2]] <- c('SEb','acc0max','acc0min','eta', 'gam','wexp','wrpe','wa','Tresp','Hresp','sig', 'predSLnP', 'SESLnP','SEcor','predProb','BIC');
  #put back n0 above if 10p, also ncol = 15
  
  for (ptN in pts){ #1:dim(datArW61)[3] ){
    
    
    D <- array(NA,c(dim(datArW61)[1:2],1));  # Create & clear the working array
    D[,,1] <- datArW61[,,ptN];
    dimnames(D)[[2]] <- c('gp','pred','obs','SE','nofb');
    dimnames(D)[[3]] <- ptN;
    ml1fit[[ptN]] <- list();
    mPD <- Inf;
    
    #tryPmatrixwbest<-matrix(NA,nrow=129,ncol=10)
    #tryPmatrixwbest[1:128,] <- tryPmatrix
    #tryPmatrixwbest[129,] <-bestsofar[ptN,]
    #allsets <- matrix(NA,nrow=129,ncol=12)
    
    #from 10p versions, now leave out n0 :
    #tryPmatrix9 <- tryPmatrix[,2:10] 
    #bestsofar9 <-  bestsofar[,2:10]
    
    tryPmatrixwbest<-matrix(NA,nrow=129,ncol=11)
    tryPmatrixwbest[1:128,] <- tryPmatrix
    tryPmatrixwbest[129,] <-as.numeric(bestsofar11[ptN,])
    allsets <- matrix(NA,nrow=129,ncol=13) #stopped here
    
    attempts <- c(129,1:128); if (Debug){ attempts <- c(129)} # , 73) } #attempts <- c(1,25,50,75,100,125) for testing
    for (set in attempts) { 
      tryP = nat2trLP0( tryPmatrixwbest[set,] )
      iniLen=length(tryP);
      
      ## this part included if further randomised attempts are wanted - remove the double #s and add attempts as an argument
      ##for (attempt in 1:attempts){  # 2-10 for testing; try (10*iniLen) for real { #put attempts back into function as an argument if you want
      ##if (attempt == 1){
      ##    iniTrPar <- nat2trLP4(tryP);
      ##} else {
      ##    iniTrPar <- nat2trLP4(tryP) * runif(10,0.75,1.25);
      ##}
      ## # atI <- attempt %% iniLen; if (atI == 0){atI <- iniLen;};
      ## #if (attempt > 11) {
      ## #    iniTrPar <- as.numeric(iniTrParM[atI,] * runif(10,0.8,1.2));
      ## #}
      
      print(paste('ptN:',ptN,';  fit attempt:', set));  print(paste('Init. Cond:', paste(round(tr2natLP0(tryP),3),collapse=',')));
      try( fitAttempt <- nlm(msLP0tr, tryP, D, print.level=nlmprintlev, iterlim=500) #Par0?
           
      ); # Par0, print.level=2, iterlim=500) 
      
      
      
      if (vecTRUE(length(fitAttempt$estimate)>1)){
        if ( vecTRUE(fitAttempt$minimum < mPD) || !(vecTRUE(length(ml1fit[[ptN]][[1]])>1)) ){
          mPD <- fitAttempt$minimum;
          ml1fit[[ptN]][[1]] <- fitAttempt;
        }
        #now to save estp and summed log likelihoods for all trials
        estpOfAttempt <- (tr2natLP0(fitAttempt$estimate)) ;
        #allsets[set,1:10] <- estpOfAttempt
        #attemptSLP <- SLPsocio4(estpOfAttempt, D);
        #allsets[set,11]   <- attemptSLP$SESLnP
        #allsets[set,12]   <- attemptSLP$predSLnP
        
        allsets[set,1:11] <- estpOfAttempt
        allsets[set,12]   <- SLPsocioCA_t(estpOfAttempt, D)$SESLnP
        allsets[set,13]   <- SLPsocioCA_t(estpOfAttempt, D)$predSLnP
        ml1fit[[ptN]][[3]] <- allsets 
      }
      
      
      ## } 
    } # End exploration of initial conditions
    
    est10p <- (tr2natLP0(ml1fit[[ptN]][[1]]$estimate)) ;
    ml1fit[[ptN]][[2]] <- SLPsocioCA_t(est10p, D);
    names(ml1fit[[ptN]]) <- c('NLM','SLP','alltrials')
    
    #est9p <- (tr2natLP4(ml1fit[[ptN]][[1]]$estimate)) ;
    #ml1fit[[ptN]][[2]] <- SLPsocio4(est9p, D);
    #names(ml1fit[[ptN]]) <- c('NLM','SLP','alltrials')
    
    # output array storage
    #ml1res[ptN,1:10] <- tr2natLP4(ml1fit[[ptN]][[1]]$estimate);
    #ml1res[ptN,11]   <- ml1fit[[ptN]][[2]][[1]];
    #ml1res[ptN,12]   <- ml1fit[[ptN]][[2]][[2]];
    co <- cor(na.omit(ml1fit[[ptN]][[2]][[3]][,,1][,c('SE','expSE')])) #SE correlation, 2sf
    ml1res[ptN,14]   <- round(co[1,2],digits=4);    
    #ml1res[ptN,14]   <- exp(ml1res[ptN,11]/192) 
    
    ml1res[ptN,1:11] <- tr2natLP0(ml1fit[[ptN]][[1]]$estimate);
    ml1res[ptN,12]   <- ml1fit[[ptN]][[2]][[1]];
    ml1res[ptN,13]   <- ml1fit[[ptN]][[2]][[2]];
    #v <- D[,'SE',1]; v <- 1+v; v<- v/v;
    #expSE <- ml1fit[[ptN]][[2]][[3]][,'expSE',1]*c(NA,v);  expSE <- expSE[-1]; #previously from line 669
    #ml1res[ptN,12]   <- round(cor(na.omit(data.frame(D[,'SE',1], expSE)))[1,2],2) #SE correlation, 2sf
    ml1res[ptN,15]   <- exp(ml1res[ptN,12]/192)  #cross-entropy 
    
    #now the individual BIC
    LnforBIC = ml1fit[[ptN]][[2]]$predSLnP + ml1fit[[ptN]][[2]]$SESLnP 
    nforBIC = length(na.omit(D[,'pred',1]))+ length(na.omit(D[,'SE',1])) 
    kforBIC = 10 #number of params - for 10p, 10
    ml1res[ptN,16]   <- log(nforBIC)*kforBIC - 2*LnforBIC #BIC formula #for 10p, 15 not 14
    
    # filename stem useful for saving stuff:
    if (ptN < 10){  fname=paste(outDir,"soc00fitPt0",ptN,sep='')
    } else {        fname=paste(outDir,"soc00fitPt",ptN,sep='')    };
    
    #now to save images - here they are saved in outDir
    mypath <- file.path(paste(fname, "_1.png", sep = "")) 
    png(file=mypath, width = 912, height = 742, units = "px")
    
    # Prepare for graphs with real & randomly generated data for visual inspection:
    d <- ml1fit[[ptN]][[2]][[3]][,,1] ;  d <- na.omit(TDSE(d));
    c <- round(cor(d[,c('sAp','TDSE')])[1,2],2);
    plot(na.omit(d[,c('sAp','TDSE')]), main=paste('pt',ptN,'  cor=',c));
    # Plotting of expected, measured and generated SE :
    c <- round(ml1res[ptN,14],3)
    d2plot <- (ml1fit[[ptN]][[2]][[3]][,c('SE','expSE','genSE','obsP','pred','obs'),1])
    plot(d2plot[,'SE'],t='p',col='green4',pch=19,lwd=5,
         main=paste('MAP fit, RW model (SLPsocioCA_t):   pt',ptN,';  cor=',c,'\n[green: measured;   blue:fitted,  pink: generated from fit]'),
         xlab='trial number',
         ylab='Self Evaluation'); 
    lines(d2plot[,'expSE'],t='l',col='blue',lwd=3);
    lines(d2plot[,'genSE'],t='l',col='pink3');
    dev.off()  # End of plot 1
    
    # Plot log p for predictions
    mypath <- file.path(paste(fname, "_2.png", sep = "")) 
    png(file=mypath, width = 912, height = 400, units = "px"); 
    plot(log(d2plot[,'obsP']),ylim=c(log(0.005),0),col='cyan3',pch=19, t='b',main=paste('pt:',ptN,'  Cyan: ln(pred.lik.); red=pred; pink=approval')); 
    lines(d2plot[,'pred']-4,col='red3',t='b'); lines(d2plot[,'obs']-4.1,col='pink3',t='p',pch=19);
    dev.off()
    # And their histogram:
    mypath <- file.path(paste(fname, "_3.png", sep = "")) 
    png(file=mypath, width = 912, height = 500, units = "px"); 
    try(hist(log(d2plot[,'obsP']),30,main=paste('pt:',ptN,'  histogram of ln(P(prediction))'),col='gray'));
    dev.off()
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  MAIN SAVING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print('ABOUT TO SAVE:',quote = FALSE)
    print(fname,quote = FALSE)
    y = data.frame(t(as.matrix(ml1res[ptN,],1,length(ml1res[ptN,]))));
    y$ptN = ptN;
    exportCSV(paste(fname,'res.csv',sep=''),y);
    exportCSV(paste(fname,'evo.csv',sep=''),ml1fit[[ptN]][[2]][[3]][,,1]);
    save.image(paste(fname,'all.RData',sep=''));
    
    ml1fit <<- ml1fit #so that these objects have global assignment
    ml1res <<- ml1res #otherwise, only exist within the function
    
  }
}    # end of function definition RWModelFit


RWModelFit(participant,NULL,Debug) 
write.csv(ml1res,"ml1res.csv")
saveRDS(ml1fit, "ml1fit.rds")


#end of nlm fitting script


