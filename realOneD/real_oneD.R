## Assume the working directory is where it was run from in robotnik-net

## Load packages
require(np)
# options(np.messages=TRUE)
source("load_dataset.R")

### Scanline 42049

scanline <- load.dataset('scanline_42049/scanline_42049.json')
names(scanline)[2] <- "exchangeRate"
# dates_usd_isk <- read.csv(file = "ert_bil_eur_m_1_Data.csv") #see this file for the dates
# usd_isk$Dates <- seq(as.Date("1999/1/1"), as.Date("2019/7/1"), "months")

scanline <- scanline[150:300,]

## Regression settings
ckertype <- "uniform"
bwmethod <- "cv.aic"
regtype <- "lc"

set.seed(11)
## Fit data using local constant kernel estimator (LCKE)
bw.lc <- npregbw(exchangeRate~t,data=scanline,regtype=regtype,bwmethod=bwmethod,ckertype=ckertype)
model.lc <- npreg(bw.lc)
## Set this temporarily...
model.llc <- model.lc
## Use LCKE \widehat{g}(x) as pilot g(x) in LCAS for the first time,
## and then LCAS for each repeat
scanline_local <- data.frame(scanline,gX=model.llc$mean)
bw.llc <- npregbw(exchangeRate~t+gX,data=scanline_local,regtype=regtype,bwmethod=bwmethod,ckertype=ckertype)
model.llc <- npreg(bw.llc)

# bw.llc.a <- npregbw(exchangeRate~t+gX,
#                   data=scanline_local,
#                   regtype=regtype,
#                   bws=c(10,bw.llc$bw[2]),
#                   bandwidth.compute = FALSE,
#                   ckertype=ckertype)
# model.llc.a <- npreg(bw.llc.a)

plot(scanline_local$t,scanline_local$exchangeRate, ylab="Grayscale",xlab="Pixel number",pch=16, cex=1,col="black",las=1)
lines(scanline_local$t,model.llc$mean,col="red",lty = 1,lwd=2)
lines(scanline_local$t,model.lc$mean,col="blue",lty = 2,lwd=2)
# lines(scanline_local$t,model.llc.a$mean,col="blue",lty = 2,lwd=2)
legend("top",
       legend = c("Data","Isotropic","Anisotropic"), 
       col = c("black","blue","red"), 
       pch = c(16,NA,NA),
       lty = c(NA,2,1),lwd=c(NA,2,2))

model.lc$MSE
model.llc$MSE
# model.llc.a$MSE

## Two D

# library(jpeg)
# picture <- readJPEG("scanline_42049/42049.jpg")



### Scanline 126007

scanline <- load.dataset('scanline_126007/scanline_126007.json')
names(scanline)[2] <- "exchangeRate"
# dates_usd_isk <- read.csv(file = "ert_bil_eur_m_1_Data.csv") #see this file for the dates
# usd_isk$Dates <- seq(as.Date("1999/1/1"), as.Date("2019/7/1"), "months")

## Regression settings
ckertype <- "uniform"
bwmethod <- "cv.aic"
regtype <- "lc"

set.seed(11)
## Fit data using local constant kernel estimator (LCKE)
bw.lc <- npregbw(exchangeRate~t,data=scanline,regtype=regtype,bwmethod=bwmethod,ckertype=ckertype)
bw.lc <- npregbw(exchangeRate~t,data=scanline,regtype=regtype,
                 bandwidth.compute=F,
                 bws = 5*bw.lc$bw,
                 ckertype=ckertype)
model.lc <- npreg(bw.lc)
## Set this temporarily...
model.llc <- model.lc

## Use LCKE \widehat{g}(x) as pilot g(x) in LCAS for the first time,
## and then LCAS for each repeat
scanline_local <- data.frame(scanline,gX=model.llc$mean)
bw.llc <- npregbw(exchangeRate~t+gX,data=scanline_local,regtype=regtype,bwmethod=bwmethod,ckertype=ckertype)
model.llc <- npreg(bw.llc)

# bw.llc <- npregbw(exchangeRate~t+gX,
#                   data=scanline_local,
#                   regtype=regtype,
#                   bws=c(25, 25),
#                   bandwidth.compute = FALSE,
#                   ckertype=ckertype)
# model.llc <- npreg(bw.llc)

# bw.lc <- npregbw(exchangeRate~t,data=scanline,regtype=regtype,bwmethod=bwmethod,ckertype=ckertype)
# model.lc <- npreg(bw.lc)


plot(scanline_local$t,scanline_local$exchangeRate, ylab="Grayscale",xlab="Pixel number",pch=16, cex=1,col="black",las=1)
lines(scanline_local$t,model.llc$mean,col="red",lty = 1,lwd=2)
lines(scanline_local$t,model.lc$mean,col="blue",lty = 2,lwd=2)
legend("topleft",
       legend = c("Data","Isotropic","Anisotropic"), 
       col = c("black","blue","red"), 
       pch = c(16,NA,NA),
       lty = c(NA,2,1),lwd=c(NA,2,2))

model.lc$MSE
model.llc$MSE
