##########################################
### Experiment 2b: Yoked vehicle group ###
### PLOTTING                           ###
##########################################

# #LOADING RAW DATA
# load(file=paste(getwd(), "/Data for R/csacqidx.ridx", sep=""))
# load(file=paste(getwd(), "/Data for R/csacqdata.rdat", sep=""))
# load(file=paste(getwd(), "/Data for R/dataperbin.rdat", sep=""))
# load(file=paste(getwd(), "/Data for R/parseddata.rdat", sep=""))
# load(file=paste(getwd(), "/Data for R/dataPerHalf.rdat", sep=""))
# load(file=paste(getwd(), "/Data for R/YokingData.rdat", sep="))

# Load libraries, install packages
library(matrixStats)
library(ggplot2)
library(plyr)
library(dplyr)
library(RColorBrewer)
install.packages("doBy")
library(doBy)


#############################################
### FIGURE 1.D. KINDS OF CUES PER SESSION ###
#############################################

#Divide data into the 2 groups
BilatAP5cases <- filter(YokingData, group==" BILATAP5")
YokedVEHcases <- filter(YokingData, group==" YOKEDVEH")

# Write a function to draw error bars
MakeErrBars<- function(x, y, errTerm, color="black", lwd=2){
        arrows(x0=x, x1=x, y0=y-errTerm, y1=y+errTerm, code=3, angle=90, length=0.1, col=color, lwd=lwd)
}

#Make plot
plot.new()
maxx=14
maxtop=60
maxy=50
plot.window(xlim=c(0, maxx=14), ylim=c(0, maxtop))

#BILATERAL AP5
dat <- BilatAP5cases

for(i in 1:length(dat$session)){
        selDat <- filter(dat, session==i)
        rect(xleft=i, xright=i+1, ybottom=0, ytop=selDat$meanCSRew, col="gray15", lwd=3)
        rect(xleft=i, xright=i+1, ybottom=selDat$meanCSRew, ytop=selDat$meanCSRew+selDat$meanCSMiss, col="white", lwd=3)
        MakeErrBars(x=i+0.5, y=selDat$meanCSRew, errTerm=selDat$semCSRew)
        MakeErrBars(x=i+0.5, y=selDat$meanCSRew+selDat$meanCSMiss, errTerm=selDat$semCSMiss)
}


#YOKED VEH
dat <- YokedVEHcases

for(i in  1:length(dat$session)){ #(maxx-1):(maxx-length(dat$session))){
        selDat <- filter(dat, session==i)
        panelStart <- maxx/2
        rect(xleft=panelStart+i, xright=panelStart+i+1, ybottom=0, ytop=selDat$meanCSRew, col="gray15", lwd=3)
        rect(xleft=panelStart+i, xright=panelStart+i+1, ybottom=selDat$meanCSRew+selDat$meanCSMiss, ytop=BilatAP5cases$meanCSRew[i]+BilatAP5cases$meanCSMiss[i], border="red", lwd=3)
        rect(xleft=panelStart+i, xright=panelStart+i+1, ybottom=selDat$meanCSRew, ytop=selDat$meanCSRew+selDat$meanCSMiss, col="white", lwd=3)
        MakeErrBars(x=panelStart+i+0.5, y=selDat$meanCSRew+selDat$meanCSMiss, errTerm=selDat$semCSMiss, lwd=2)
        MakeErrBars(x=panelStart+i+0.5, y=BilatAP5cases$meanCSRew[i]+BilatAP5cases$meanCSMiss[i], errTerm=selDat$semextraCues, lwd=2)
}


axis(side=1, at=seq(1, 6)+0.5, labels=seq(1,6), cex.axis=1.5)
axis(side=1, at=maxx/2+seq(1, 6)+0.5, labels=seq(1,6), cex.axis=1.5)

axis(side=2, at=seq(0, maxy, by=10), cex.axis=1.5, pos=0.5, las=2)

text(x=4, y=maxy+5, labels="AP5", cex=2, font=2)
text(x=11, y=maxy+5, labels="Saline", cex=2, font=2)

mtext(text="Day", side=1, line=3, at=4, cex=1.5)
mtext(text="Day", side=1, line=3, at=11, cex=1.5)
mtext(text="Number of trials", side=2, at=25, cex=1.5)



##############################################################################
### EXPERIMENT 2B: TIME SPENT IN RECEPTACLE PRE AND POST CUE DURING TEST   ###
##############################################################################

### I'm looking at the total time (sum) spend in the receptacle in the 10s window after cue end vs. 10s window pre cue end. I looked at the first 40 trials (I want to analyze all trials, but some rats got less and others more trials, but all of them got 40 y pico, so I subset the first 40 trials)

data1 <- csacqidx %>% group_by(group) %>% select(durDiffBin1, durDiffPerTrialBin1, PreCueDurSumBin1, PostCueDurSumBin1) %>% summarize(
        meandurDiffBin1=mean(durDiffBin1),
        meandurDiffPerTrialBin1=mean(durDiffPerTrialBin1),
        PreCueDurSumBin1=mean(PreCueDurSumBin1),
        PostCueDurSumBin1=mean(PostCueDurSumBin1),
        semDurDiff=sd(durDiffBin1)/sqrt(length(durDiffBin1)),
        semDurDiffPerTrialBin1=sd(durDiffPerTrialBin1)/sqrt(length(durDiffPerTrialBin1)),
        semPreCueDurSum=sd(PreCueDurSumBin1)/sqrt(length(PreCueDurSumBin1)),
        semPostCueDurSum=sd(PostCueDurSumBin1)/sqrt(length(PostCueDurSumBin1))
)

colindx = c("red", "blue")

#Plotting a bar plot
png (paste(graphdir, "Exp 2b Mean Sum Time spent in Post cue minus Pre Cue all trials.png", sep=""), width=10, height=10, units="in", res=100)  

plot.new()
plot.window(xlim=c(1, 3), ylim=c(-5, 20))

abline(h=seq(-5, 20, 5), col="gray75")
for(i in 1:length(unique(data1$group))){
        rect(xleft=i, xright=i+1, ybottom=0, ytop=data1[i,]$meandurDiffBin1, col=colindx[i])
        arrows(x0=i+0.5, x1=i+0.5, y0=data1[i,]$meandurDiffBin1, y1=data1[i,]$meandurDiffBin1+data1[i,]$semDurDiff, angle=90, length=0.25)
        arrows(x0=i+0.5, x1=i+0.5, y0=data1[i,]$meandurDiffBin1, y1=data1[i,]$meandurDiffBin1-data1[i,]$semDurDiff, angle=90, length=0.25)
}
axis(side=2, at=seq(-5, 21, 5), cex.axis=2.2, font=2)
axis(side=1, at=seq(1.5, 2.5, 1), labels=c("AP5", "Yoked Saline"), tick=F, cex.axis=2.2, font=2, line=-5)


dev.off()

#### ITI ENTRIES PER SEC
load(paste(getwd(), "/Data for R/csacqidx.ridx", sep=""))
load(paste(getwd(), "/Data for R/dataperbin.rdat", sep=""))

data2 <- group_by(dataperbin, group, bin, add=TRUE) 
data3 <- summarise(data2, 
                   meanITIpersec = mean(ITIperSec, na.rm=T),
                   semITIpersec = sd(ITIperSec, na.rm=T)/sqrt(length(ITIperSec)),
                   meanlatency = mean(latency, na.rm=T),
                   semlatency = sd(latency, na.rm=T)/sqrt(length(latency)),
                   meanTA = mean(latencyCUEvsITI, na.rm=T),
                   semTA = sd(latencyCUEvsITI, na.rm=T)/sqrt(length(latencyCUEvsITI)) )

png (paste(graphdir, "Exp 2b training ITI per sec.png", sep=""), width=11.37, height=11, units="in", res=100)  


plot.new()
plot.window(xlim=c(1, 6), ylim=c(0,0.2))
abline(h=seq(0,0.2,0.05), col="gray75")
colindx = c("#FFBF3D", "#DCF0F0")

axis(side=1, at=seq(1, 6, 1), labels=c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30"), font=2, cex.axis=2.4)
axis(side=2, at=seq(0, .2, .05), font=2, cex.axis=2.4, las=2)


for(i in 1:length(unique(data3$group))){
        a <- unique(data3$group)[i]
        sel <- filter(data3, group==a)
        
        for(j in 1:nrow(sel)){
                arrows(x0=j+((i-1)*0.03), x1=j+((i-1)*0.03), y0=sel$meanITIpersec[j], y1=sel$meanITIpersec[j]+sel$semITIpersec[j], length=0.1, angle=90)
                arrows(x0=j+((i-1)*0.03), x1=j+((i-1)*0.03), y0=sel$meanITIpersec[j], y1=sel$meanITIpersec[j]-sel$semITIpersec[j], length=0.1, angle=90)
        }
        
        lines(sel$meanITIpersec, lwd=3.5, col="black")
        points(sel$meanITIpersec, pch=21, cex=4, col="black", bg=colindx[i])
}


legend (x=2.5, y=0.19, legend=c("Yoked saline", "AP5"), fill=c(colindx[2], colindx[1]), cex=2.1, bg="white", horiz=T)

dev.off()




############# TASK ACCURACY

png (paste(graphdir, "Exp 2b training Task accuracy.png", sep=""), width=11.37, height=11, units="in", res=100)  


plot.new()
plot.window(xlim=c(1, 6), ylim=c(0,7))
abline(h=seq(0,7,1), col="gray75")
colindx = c("#FFBF3D", "#DCF0F0")

axis(side=1, at=seq(1, 6, 1), labels=c("0-5", "5-10", "10-15", "15-20", "20-25", "25-30"), font=2, cex.axis=2.4)
axis(side=2, at=seq(0, 7, 1), font=2, cex.axis=2.4, las=2)


for(i in 1:length(unique(data3$group))){
        a <- unique(data3$group)[i]
        sel <- filter(data3, group==a)
        
        for(j in 1:nrow(sel)){
                arrows(x0=j+((i-1)*0.03), x1=j+((i-1)*0.03), y0=sel$meanTA[j], y1=sel$meanTA[j]+sel$semTA[j], length=0.1, angle=90)
                arrows(x0=j+((i-1)*0.03), x1=j+((i-1)*0.03), y0=sel$meanTA[j], y1=sel$meanTA[j]-sel$semTA[j], length=0.1, angle=90)
        }
        
        lines(sel$meanTA, lwd=3.5, col="black")
        points(sel$meanTA, pch=21, cex=4, col="black", bg=colindx[i])
}


legend (x=1.05, y=7.25, legend=c("Yoked saline", "AP5"), fill=c(colindx[2], colindx[1]), cex=2.1, bg="white", horiz=T)

dev.off()
