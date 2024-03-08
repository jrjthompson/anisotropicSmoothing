## Load packages
require(np)
require(parallel)

options(np.messages=FALSE)

## Data generation function
data.generate <- function(n,data.type,sigma){
  if (data.type=="uniformJump"){# 3 uniforms
    X <- c(seq(0,1,3/n),seq(1+3/n,2,3/n),seq(2+3/n,3,3/(n-1)))
    oracle <- c(rep(1,length(seq(0,1,3/n))),rep(7,length(seq(1+3/n,2,3/n))),rep(3,length(seq(2+3/n,3,3/n))))
    Y <- oracle+rnorm(length(oracle),sd=sigma)
  } else if (data.type=="continuous") { #continuous function
    X <- seq(1,3,length.out = n)
    oracle <- 50*((X/3)^2-(X/3)^3)
    Y <- oracle+rnorm(n,sd = sigma)
  } else if (data.type=="continuousWithJump"){
    X <- seq(0,3,length.out = n)
    X.1 <- X[1:round(length(X)/2)]
    X.2 <- X[(round(length(X)/2)+1):length(X)] 
    oracle <- c(50*((X.1/3)^2-(X.1/3)^3),50*((X.2/3)^2-(X.2/3)^3)+5)
    Y <- oracle+rnorm(n,sd = sigma)
  } else if (data.type=="gradualJump"){
    X <- c(seq(0,1,3/n),seq(1+3/n,2,3/n),seq(2+3/n,3,3/(n-1)))
    oracle <- c(rep(1,length(seq(0,1,3/n))),seq(1,5,4/(length(seq(1+3/n,2,3/n))-1)),rep(5,length(seq(2+3/n,3,3/n))))
    Y <- oracle+rnorm(length(oracle),sd = sigma)
  }
  data.frame(X=X,Y=Y,oracle=oracle)
}

## A Monte Carlo replicate written as a function
simulation.ALC <- function(seed,n,data.type="uniformJump",sigma,bw.fixed.value=NULL,
                            repeats=1,ckertype,bwmethod,regtype){
  set.seed(42+seed)## Generate data
  data <- data.generate(n,data.type = data.type,sigma = sigma)
  
  ## Fit data using local constant kernel estimator (LC) 
  bw.lc <- npregbw(Y~X,data=data,regtype=regtype,bwmethod=bwmethod,ckertype=ckertype)
  model.lc <- npreg(bw.lc)
  model.alc <- model.lc #set this temporarily

  ## Use LC \widehat{g}(x) as pilot g(x) in ALC for the first time,
  ## and then ALC for each repeat
  for (i in 1:repeats){
    data <- data.frame(oracle=data$oracle,X=data$X,Y=data$Y,gX=data$oracle) #replace data$oracle with model.alc$mean if you want ALC smoothing 
    ## Refit data using the local constant estimator as an input
    bw.llc <- npregbw(Y~X+gX,data=data,regtype=regtype,bwmethod=bwmethod,ckertype=ckertype)
    if (!is.null(bw.fixed.value)){
      bw.llc$bw <- c(bw.llc$bw[1],bw.fixed.value)
    }
    model.alc <- npreg(bw.llc)
  }

  ## Check fit with oracle function and retern mean squared error
  model.lc.ESE <- 0 #we already have this so we don't bother
  model.alc.ESE <- mean((model.alc$mean-data$oracle)^2)
  c(model.lc.ESE,model.alc.ESE)
}  

## Code to run MC simulations, change number of points, 
## variability of data, and kernel type, and split by number of cores

cores <- detectCores()-1

kernel.type <- "gaussian" #could try other options
bandwidth.selection.method <- "cv.ls" #could try other options
regression.type <- "lc" #lc for local constant, ll for local linear
data.types <- c("uniformJump","continuous","continuousWithJump","gradualJump")
repeats <- 1 #how many times to repeat the ALC smoother, default is 1
bw.fixed.value <- NULL #set this for specifying the ALC bw value

M <- 125 #number of Monte Carlo replicates, 63*2-1 =125 cores
n.seq <- c(100,200,400,800,1600)#,3200,6400,12800)
sigma.seq <- c(0.1,0.5,1,2)

lcr.MSE <- numeric(length(n.seq))
ALC.MSE <- numeric(length(n.seq))
lcr.ESE <- matrix(numeric(M),M,length(n.seq))
ALC.ESE <- matrix(numeric(M),M,length(n.seq))

date.started <- as.character(Sys.Date())
time.started <- proc.time()

for (data.type in data.types){
  for (sigma in sigma.seq){
    for (N in 1:length(n.seq)){
      dum = mclapply(1:M,function(m) simulation.ALC(m,n.seq[N],data.type,
                                                     sigma,bw.fixed.value,
                                                     repeats, kernel.type,
                                                     bandwidth.selection.method,
                                                     regression.type),
                     mc.cores=cores)
      lcr.ESE[,N] <- simplify2array(dum)[1,]
      ALC.ESE[,N] <- simplify2array(dum)[2,]
    }
    
    return.SEs <- data.frame(lcr=lcr.ESE,ALC=ALC.ESE)
    colnames(return.SEs) <- c(paste("lcr",n.seq,sep="."),paste("ALC",n.seq,sep="."))
    ## Write the current data into a file so we can use it in the future
    write.csv(return.SEs,paste("ALCsimulation",date.started,data.type,
              sigma,repeats,"csv",sep = "_"),row.names = FALSE)
  }
}

proc.time()-time.started
date.started
as.character(Sys.Date()) #date finished
