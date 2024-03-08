## Load packages
require(np)
source("load_dataset.R") #to work with .json files

### Scanline 42049

scanline <- load.dataset('scanline_42049/scanline_42049.json')
names(scanline)[2] <- "grayscale"
scanline <- scanline[150:300,]

## Regression settings
ckertype <- "uniform"
bwmethod <- "cv.aic"
regtype <- "lc"

set.seed(11)
## Fit data using local constant kernel estimator (LC)
bw.lc <- npregbw(grayscale~t,data=scanline,regtype=regtype,bwmethod=bwmethod,ckertype=ckertype)
model.lc <- npreg(bw.lc)

## Use LC \widehat{g}(x) as pilot g(x) in ALC 
scanline_local <- data.frame(scanline,gX=model.lc$mean)
bw.llc <- npregbw(grayscale~t+gX,data=scanline_local,regtype=regtype,bwmethod=bwmethod,ckertype=ckertype)
model.llc <- npreg(bw.llc)

plot(scanline_local$t,scanline_local$grayscale, ylab="Grayscale",xlab="Pixel number",pch=16, cex=1,col="black",las=1)
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

### Scanline 126007

scanline <- load.dataset('scanline_126007/scanline_126007.json')
names(scanline)[2] <- "grayscale"

## Regression settings
ckertype <- "uniform"
bwmethod <- "cv.aic"
regtype <- "lc"

set.seed(11)
## Fit data using local constant kernel estimator (LC)
bw.lc <- npregbw(grayscale~t,data=scanline,regtype=regtype,bwmethod=bwmethod,ckertype=ckertype)
bw.lc <- npregbw(grayscale~t,data=scanline,regtype=regtype,
                 bandwidth.compute=F,
                 bws = 5*bw.lc$bw,
                 ckertype=ckertype)
model.lc <- npreg(bw.lc)

## Use LC \widehat{g}(x) as pilot g(x) in ALC
scanline_local <- data.frame(scanline,gX=model.lc$mean)
bw.llc <- npregbw(grayscale~t+gX,data=scanline_local,regtype=regtype,bwmethod=bwmethod,ckertype=ckertype)
model.llc <- npreg(bw.llc)

## Plot results
plot(scanline_local$t,scanline_local$grayscale, ylab="Grayscale",xlab="Pixel number",pch=16, cex=1,col="black",las=1)
lines(scanline_local$t,model.llc$mean,col="red",lty = 1,lwd=2)
lines(scanline_local$t,model.lc$mean,col="blue",lty = 2,lwd=2)
legend("topleft",
       legend = c("Data","Isotropic","Anisotropic"), 
       col = c("black","blue","red"), 
       pch = c(16,NA,NA),
       lty = c(NA,2,1),lwd=c(NA,2,2))

model.lc$MSE
model.llc$MSE
