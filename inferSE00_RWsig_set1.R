# Key functions for inferring self-esteem by evidence accummulation and mapping
# inferSE00_RWsig_set1 is for doing the original RW model for eLife set complete with wexp on the cluster.
# ON THE CEDRIC VM, SETWD TO setwd("/media/sf_CHERhomeYdrive/R/R4hal") BEFORE TESTING THIS!!
#
# Alexis (An Yee) Low & Michael Moutoussis, 2019
# remove(list=ls()); # drastic clearup!
try(remove(whoami,codeDir,baseDir,seDirs,bestsofar10,bestsofar9,bestsofar)) # conservative clearup

# ---------------   Cluster related menu / argument options ------------------
runOnCluster = 0;   # 0 for off cluster run, 1 for cluster run,
# 2 for cluster-like run (with Rscript and ptN argument, 
# e.g.$ Rscript inferSE00_matrix.R 30) but off cluster.
Debug =0 ;   
# First, determine the size of grid of initial conditions:
# If run serially:
if (runOnCluster == 0){
  participant = 1:61;                  # All pts
  gridAtt = round((1:10)*12.34);   # selection of e.g. 10 grid attempts, max is 1:128 below .... 
} else {  # if run on cluster
  gridAtt = 1:128;
}

# Next, determine wich pt(s) to run:
if (Debug & !(runOnCluster)){ participant = 50; gridAtt = c()};    #A: changing this when testing - remember to change back             
# If running serially:
if (runOnCluster == 0){
  # Then 'participant' must already be in the workspace! 
  print(paste('participant = ',participant)) # participant <- 1;  
  
} else  # if running on cluster or cluster-like
{
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
  outDir <- '/media/sf_CHERhomeYdrive/R/R4hal/currentWork/socioRWsig_set1/' };
if (whoami == 'LinuxMM'){ 
  codeDir <- '~/Dropbox/FIL_aux/R_scripts/';
  outDir <- baseDir };
if (whoami == 'MMhal9000'){
  codeDir <- '~/R/R4hal/scripts/';
  outDir <- '~/R/R4hal/currentWork/socioRWsig_set1/' };


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


SLPsocio0_sig <- function( parMat, datAr, onlySLP=0, fixAcc=NA, check=1){ 
  # parMat has rows with the parameters for each pt.,
  # par. row: sensi, acc0max, acc0min, eta, gam, sesh, Tresp, Hresp, sig e.g. parMat = c(1,0.81,0.47,0.04,0.66,1,0.12,0.42,0.08) or parMat = c(1,0.83,0.59,0.019,0.67,1,0.28,0.26,0.07)
  # datAr has a page for each pt, and Ntr rows.
  #
  # For testing & example use see roughLikeMe.R, wherein 
  # can use parMat = par0a; datAr = D; onlySLP=0; fixAcc=gpAcc; check = 1; 
  # and at the end remove(parMat,datAr,onlySLP,fixAcc,check,hapPol)
  
  M <- 4; # Number of rater groups.
  knmax <- 200; # 10;  # max. depth of kernel summation.
  gamNorm = 0;  # 0 ;= don't normalize the decay kernel
  
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
  # ratingP <- 1/(1+exp((parMat[ptN,'accHalf']-hapPol[trN,expInd[gpI],ptN])/
  #                                    (parMat[ptN,'accHalf']*parMat[ptN,'Tresp'])));
  colnames(parMat) <- c('sensi','acc0max','acc0min','eta','gam','sesh','Tresp','Hresp','sig'); #colnames(parMat) <- c('SEb','acc0max','acc0min','eta','gam','wexp','wrpe','Tresp','Hresp','sig');
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
  hapPol <- array(dim=c(Ntrtot+1,M+6,Nptot)); 
  dimnames(hapPol)[[2]] <- c(paste('exp',1:M,sep=''), 'RPE',
                             'obsP','SEPD','accP','genSE','modeSE');
  expInd <- 1:M    # Indices of group expectations in hapPol 
  rpeInd <- M+1;   # index for  RPE, 
  obsPI <- rpeInd+1;       SEPDI  <- obsPI+1 ;
  accPI <- SEPDI+1;       genSEI <- accPI+1;   modeSEI <- genSEI+1;
  
  for (ptN in 1:Nptot){
    
    # Cater somehow for the fact that scanned task contains interruptions
    # including at the very start. Awkward - to start with, just exclude:
    stimTr[,ptN] <- !is.na(datAr[,'pred',ptN]); #altered GJW
    Ntrtot <- sum(stimTr[,ptN]);    # update to exclude non-stimulus trials. Was dim(datAr)[1] above.
    datAr[1:Ntrtot,,ptN] <- datAr[stimTr[,ptN],,ptN];  
    if (Ntrtot<dim(datAr)[1]) {
      datAr[(Ntrtot+1):dim(datAr)[1],,ptN] <- NA; 
    }
    
    # Initial beliefs about ratings by others (groups):
    if (!is.na(sum(fixAcc))){
      parMat[ptN,'eta'] <- 0.0;
      # hapPol[1,expInd,ptN] <- fixAcc[expInd];   # this will hopefully fail if rubbish fixAcc is given.
      # Scaled a la Geert-Jan: 
      hapPol[1,expInd,ptN] <- 2*fixAcc[expInd]-1;   # this will hopefully fail if rubbish fixAcc is given.
    } else {
      accP0 <-  ((M-1):0)*(parMat[ptN,'acc0max']-parMat[ptN,'acc0min'])/(M-1) + parMat[ptN,'acc0min'];
      hapPol[1,expInd,ptN] <- 2*accP0 - 1;  # Apply same to Hresp below.
    }
    # Stuff to apply to all trials:
    gam <- as.numeric(parMat[ptN,'gam']);
    # Poss. include a normalizing constant to gam, equal to 1/sum(kernel vector), or not ...
    # z = 1; #
    if (gamNorm) {
      z = 1/sum(c(1, cumprod(rep(gam, (sumN-1)))));
    } else {
      z = 1;
    }
    # ... to be included in the weights:  
    We = 0.00001;          Wr = 1 #Wr = z*parMat[ptN,'wrpe']; #A: for 9p RW model just set We = 0.00001
    
    # Now for acceptance prediction response function parameters:
    Tresp = 2*parMat[ptN,'Tresp'];    # this and Hresp below scaled to be directly
    Hresp = 2*parMat[ptN,'Hresp']-1;  #   comparable with the values and returns ...
    
    # loop over trials
    for (trN in 1:Ntrtot){
      hapPol[trN+1,,ptN] <- hapPol[trN,,ptN];   # Initialise-most will remain same.
      
      gpI <- datAr[trN,1,ptN];  # which group rated this pt
      
      # Prob. of 'accept' response emitted (this is BEFORE rating seen),
      # if present. NB we will use beliefs after last trial, i.e. hapPol[trN,...
      ratingP <- 1/(1+exp((-Hresp-hapPol[trN,expInd[gpI],ptN])/Tresp));
      hapPol[trN+1,accPI,ptN] <- ratingP;  # Store      
      
      # Now calc. obsP, the probability under the current model that the
      # participant emits the acceptance prediction in datAr (the experimentally observed one).
      predAcc <- datAr[trN,2,ptN]; # rating that actual participant predicted.
      if (!is.na(predAcc)){    
        if (predAcc > 0.5){  # i.e. if it's 1 
          hapPol[trN+1,obsPI,ptN] <- ratingP;
        } else {
          hapPol[trN+1,obsPI,ptN] <- 1-ratingP;
        }
      } # End if valid prediction predAcc
      
      # Sum happinesSE eqn
      # First find out how many timepoints to sum over. 
      sumN <- knmax; if ( trN < knmax ){ sumN <- trN; };
      
      # RPE for this trial. Note that the EV is from the end
      # of the trial before ( hapPol[trN,... ) but the outcome is from this one,
      #  datAr[trN, ... :
      nofb = datAr[trN,'nofb',ptN];  # 0 if feedback given, 1 otherwise
      # so if no feedback given the outcome is said to be 0 acc. to GJW / RR model.
      EV = hapPol[trN,expInd[datAr[trN,'gp',ptN]],ptN]; 
      hapPol[trN+1,rpeInd,ptN] <-  (1-nofb)*datAr[trN,'obs',ptN]  - EV ;
      # Now for the SE sum:
      sensi <- parMat[ptN,1]
      sesh <- parMat[ptN,6]
      SE <-  0 #we start from SE=0, accumulate the sum, and then pass it through the sigmoid after the 'while (j>0)' loop
      j = sumN;
      while (j>0) {  # done with a while loop to skip it if sumN already 0, 
        # which may be the case if we want to skip trials by using remainder
        # to detect block starts for fmri version ... :-(
        SE <- SE + gam^(j-1)*
          (  We* hapPol[trN-j+1,expInd[datAr[trN-j+1,'gp',ptN]],ptN] +
               Wr* hapPol[trN-j+2,rpeInd,ptN] ) #note that here, in sigmoid version, Wr = 1, not a param
        
        #Now adding sigmoid response function
        
        
        j=j-1;
      }
      # Store the current modal SE:
      hapPol[trN+1,modeSEI,ptN] <- state2seen(SE,sensi,sesh)   # The accummulated sum is passed through the sigmoid response
      
      # if loop below zapped as in the first instance NA group trials excluded already ...
      # if (!is.na(datAr[trN,'gp',ptN])){  # if valid group stimulus was presented
      #hapPol[trN+1,expSEI,ptN] <- SE;
      
      # Determine density at the actual response:
      SEdat <- datAr[trN, 'SE',ptN]
      if (is.na(SEdat)) {
        hapPol[trN+1,SEPDI,ptN] <- NA;
        
        
      } else {
        SEspt <- seen2state(SEdat,sensi,sesh)
        
        SESptdens <- dnorm( SEspt, SE, parMat[ptN,'sig'] ) #check where/how SEs is generated
        hapPol[trN+1,SEPDI,ptN] <- slopeseen2state(SEdat,sensi)*SESptdens
      }
      
      
      # generate pseudodatum:
      genSEs <- rnorm(1, SE, parMat[ptN,'sig'] );
      hapPol[trN+1,genSEI,ptN] <- state2seen(genSEs,sensi,sesh) 
      
      # Don't forget to update the EV for the gp we've just been rated by:
      hapPol[trN+1,expInd[gpI],ptN] <- hapPol[trN+1,expInd[gpI],ptN] + 
        parMat[ptN,'eta']* hapPol[trN+1,rpeInd,ptN] ; 
      
      #} else {  # if no group stimulus was presented, fill in with NAs. Note
      #          # that the group expectations are not touched - they just stay as they are.
      #    hapPol[trN+1,c(expSEI,SEPDI,genSEI),ptN] = c(NA,NA,NA) ; 
      #}
      
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
      colnames(SLPetc[[5]]) <- c('sensi','acc0max','acc0min','eta', 'gam','sesh','Tresp','Hresp','sig');
      names(SLPetc) <-c('predSLnP','SESLnP','DatBelPol','genD','ptPar');
      
      
      return(SLPetc);
    }
  }
  
} # end of SLPsocio0 - 'Happiness-like model' (Robb's)


msLP0tr <- function(trParM, datAr, Pri=NA, Check=0){
  # trParM only to have atanh(2*SEb-1), log(gam,  wexp,  wrpe, Tresp,  sig)
  # if (Check){
    if (is.null(dim(trParM))){
      trParM <- matrix(trParM,nrow=1,byrow=TRUE)
    }
  # }
  parM <- matrix(NA,nrow=dim(trParM)[1],ncol=3+dim(trParM)[2]);   # long form ...
  
  # fill in long form :
  if (dim(trParM)[2] == 7){  
    parM[,c(1,5:10)]  <- tr2natLP0(trParM);  # keep NAs for acc0max, acc0min, eta
  } else {
    if (dim(trParM)[2] == 9) { 
      parM  <- tr2natLP0(trParM);
    } else {
      print(trParM);
      stop('tr2natLP0 not ready for this transformed parameter matrix, dim(trParM) ');
    }
  }
  
  
  # Cacl. the log prior for MAP purposes etc, all calc'd in short form:
  mSLPrior <- 0;
  if (length(Pri)>1){  # legit prior must have 12 elements or so!
    for (ptN in 1:dim(trParM)[1]) {
      # in the line below na.omit bec. usually we don't have acc0max, acc0min, eta

      mSLPrior <- mSLPrior - sum(dbetaMS(parM[ptN,2:3], Pri[1,2:3],Pri[2,2:3], log=TRUE)); 
      
      mSLPrior <- mSLPrior - sum(dnorm(parM[ptN,c(1,4:9)], Pri[1,c(1,4:9)],Pri[2,c(1,4:9)], log=TRUE)); 
    }
  } 
  if ((mSLPrior == Inf) || is.na(mSLPrior)){  # If we are in an a priori prohibited parameter region
    # do not attempt to calculated the likelihood - it will be nonsense anyway.
    return(Inf); 
    
  }  
  else {
  
    SLP <- SLPsocio0_sig( parM, datAr, check=Check)
    predSLnP <- SLP[[1]]
    SESLnP <- SLP[[2]]
    value <- mSLPrior - SLP[[1]] - SLP[[2]] 

  
     # else if (value == 0) {    return(Inf);  }
     #   else {
     
     if (is.na(SLP$DatBelPol[,'obsP',][2]) & 
              is.na(SLP$DatBelPol[,'obsP',][3])) {
       return(Inf); }
       
      else {
        return (value);
        }
     }
     
     #so the lower the better

  
  # was:  mSLP <- - SLPsocio0( parM, datAr, onlySLP=1, fixAcc, check=1)
  #  return(mSLP);
  
} # end of msLP0tr


nat2trLP0 <- function(parM,check=1){
  if (check){
    if (is.null(dim(parM))){
      parM <- matrix(parM,nrow=1,byrow=TRUE)
    }
  }  
  if (dim(parM)[2] == 9){  # This means acc0max, acc0min, eta given.
    #  c('sensi','acc0max','acc0min','eta','gam','sesh','Tresp','Hresp','sig')
    # stop('tr2natLP0 not ready for non-ready-made general group acceptance probabilities');
    trParM <- matrix(NA,nrow=dim(parM)[1],ncol=dim(parM)[2]); 
    trParM[,c(2:4,8)] <- as.matrix(atanh(2*parM[,c(2:4,8)]-1));
    trParM[,c(1,5,7,9)] <- as.matrix(log(parM[,c(1,5,7,9)])) ;
    trParM[,6] <- as.matrix(parM[,6]) ;
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
  if (dim(trParM)[2] == 9 ){  # This means acc0max, acc0min, eta given.
    parM <- matrix(NA,nrow=dim(trParM)[1],ncol=dim(trParM)[2]);
    colnames(parM) <- c('sensi','acc0max', 'acc0min', 'eta', 'gam','sesh','Tresp','Hresp','sig');
    parM[,c(2:4,8)] <- as.matrix(0.5*(tanh(trParM[,c(2:4,8)])+1));
    parM[,c(1,5,7,9)] <- as.matrix(exp(trParM[,c(1,5,7,9)])) ;
    parM[,6] <- as.matrix(trParM[,6]) ;
    
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
parhd =  c('sensi','acc0max','acc0min','eta','gam','sesh','Tresp','Hresp','sig');
npar = length(parhd)
try({load(paste(baseDir,"loadfornlm_RW_sigmoid_set1.RData",sep='')); #change this depending on which set
  bestsofar <- importCSV(paste(baseDir,"ml0res1lr.csv",sep=''))})  
#this contains bestsofar9, tryPmatrix, datArW61, and priors
try({load(paste(baseDir,"/currentWork/loadfornlm_RW_sigmoid_set1.RData",sep='')); 
  bestsofar <- importCSV(paste(baseDir,"currentWork/ml0res1lr.csv",sep=''))}) # for kikidi / hal ?
# ml0res1CNNP.csv for basic CNNP fit with 1 lr
bestsofar <- bestsofar[,parhd];

RWModelFit <- function(pts,Par0=NULL,nlmprintlev=0,grAtt=c()) {
  
  priors <- priors_RWsig
  
  if (is.null(Par0)){ Par0 <- priors; }
  
  ml1fit <- list();
  ml1res <- matrix(NA,nrow=dim(datArW61)[3],ncol=npar+5);
  dimnames(ml1res)[[2]] <- c(parhd, 'predSLnP', 'SESLnP','SEcor','predProb','BIC');
  
  for (ptN in pts){ #1:dim(datArW61)[3] ){
    
    
    D <- array(NA,c(dim(datArW61)[1:2],1));  # Create & clear the working array
    D[,,1] <- datArW61[,,ptN];
    dimnames(D)[[2]] <- c('gp','pred','obs','SE','nofb');
    dimnames(D)[[3]] <- ptN;
    ml1fit[[ptN]] <- list();
    mPD <- Inf;
    
    tryPmatrixwbest<-matrix(NA,nrow=1+length(grAtt),ncol=npar)
    try({ tryPmatrixwbest[1:length(grAtt),] <- tryPmatrix[grAtt,] }) #for the first 128 rows
    tryPmatrixwbest[length(grAtt)+1,] <-as.numeric(bestsofar[ptN,])
    allsets <- matrix(NA,nrow=length(grAtt)+1,ncol=12) 
    
    attempts <- c(length(grAtt)+1,grAtt); 
    attempts <- c(length(grAtt)+1,grAtt)
    if (Debug){ attempts <- c(length(grAtt)+1)} 
    ;# , 73) } #attempts <- c(1,25,50,75,100,125) for testing
    for (set in length(attempts):1) { 
      # correct for the possibility of values rounded to too small:
      p2try <- tryPmatrixwbest[set,]; 
      smallind <- (abs(p2try) < 1e-7);  # We don't care for such small intial values ... 
      p2try[smallind] <- sign(p2try[smallind]+1e-20)*1e-6;  # The 'sign' stuff is probably redundant - defensive coding.
      
      tryP = nat2trLP0( p2try );       
      # It is now convenient to check for overflow approximations:
      bigind <- (abs(tryP) > 100);
      tryP[bigind] <- sign(tryP[bigind])*16;  # Note very small 'large value' as these get exponentiated in a few param. transforms.
      
      iniLen=length(tryP);
      
      print(paste('ptN:',ptN,';  fit attempt:', set));  print(paste('Init. Cond:', paste(round(tr2natLP0(tryP),3),collapse=',')));
      try( fitAttempt <- nlm(msLP0tr, tryP, D, Par0, print.level=nlmprintlev, iterlim=500) #Par0?
           
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
        
        allsets[set,1:9] <- estpOfAttempt
        apply <- SLPsocio0_sig(estpOfAttempt, D)
        allsets[set,10]   <- apply$SESLnP
        allsets[set,11]   <- apply$predSLnP
        ml1fit[[ptN]][[3]] <- allsets 
      }
      
      
      ## } 
    } # End exploration of initial conditions
    
    est9p <- (tr2natLP0(ml1fit[[ptN]][[1]]$estimate)) ;
    ml1fit[[ptN]][[2]] <- SLPsocio0_sig(est9p, D);
    names(ml1fit[[ptN]]) <- c('NLM','SLP','alltrials')
    
    # output array storage
    #ml1res[ptN,1:10] <- tr2natLP4(ml1fit[[ptN]][[1]]$estimate);
    #ml1res[ptN,11]   <- ml1fit[[ptN]][[2]][[1]];
    #ml1res[ptN,12]   <- ml1fit[[ptN]][[2]][[2]];
    co <- cor(na.omit(ml1fit[[ptN]][[2]][[3]][,,1][,c('SE','modeSE')])) #SE correlation, 2sf
    ml1res[ptN,12]   <- round(co[1,2],digits=4);    
    #ml1res[ptN,14]   <- exp(ml1res[ptN,11]/192) 
    
    ml1res[ptN,1:9] <- tr2natLP0(ml1fit[[ptN]][[1]]$estimate);
    ml1res[ptN,10]   <- ml1fit[[ptN]][[2]][[1]];
    ml1res[ptN,11]   <- ml1fit[[ptN]][[2]][[2]];
    #v <- D[,'SE',1]; v <- 1+v; v<- v/v;
    #expSE <- ml1fit[[ptN]][[2]][[3]][,'expSE',1]*c(NA,v);  expSE <- expSE[-1]; #previously from line 669
    #ml1res[ptN,12]   <- round(cor(na.omit(data.frame(D[,'SE',1], expSE)))[1,2],2) #SE correlation, 2sf
    ml1res[ptN,13]   <- exp(ml1res[ptN,10]/192)  #cross-entropy 
    
    #now the individual BIC
    LnforBIC = ml1fit[[ptN]][[2]]$predSLnP + ml1fit[[ptN]][[2]]$SESLnP 
    nforBIC = length(na.omit(D[,'pred',1]))+ length(na.omit(D[,'SE',1])) 
    kforBIC = length(parhd) 
    ml1res[ptN,14]   <- log(nforBIC)*kforBIC - 2*LnforBIC #BIC formula #for 10p, 15 not 14
    
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
    c <- round(ml1res[ptN,12],3)
    d2plot <- (ml1fit[[ptN]][[2]][[3]][,c('SE','modeSE','genSE','obsP','pred','obs'),1])
    plot(d2plot[,'SE'],t='p',col='green4',pch=19,lwd=5,
         main=paste('MAP fit, RW model (SLPsocio0_sig):   pt',ptN,';  cor=',c,'\n[green: measured;   blue:fitted,  pink: generated from fit]'),
         xlab='trial number',
         ylab='Self Evaluation'); 
    lines(d2plot[,'modeSE'],t='l',col='blue',lwd=3);
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


RWModelFit(participant,NULL,Debug,gridAtt) 
#RWModelFit(26,NULL,Debug,1:5)
write.csv(ml1res,"ml1res.csv")
saveRDS(ml1fit, "ml1fit.rds")


#end of nlm fitting script


