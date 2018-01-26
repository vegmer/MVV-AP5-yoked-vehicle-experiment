##########################################
### Experiment 2b: Yoked vehicle group ###
### DATA EXTRACT AND PARSE             ###
##########################################

#Call necessary packages
library("matrixStats")
library("plyr")
library("dplyr")
library("ggplot2")

# Set working directory
# setwd("C:/Users/Mercedes/Dropbox/DISSERTATION/Experiment 2b/Training/")

####################################
### 1. BEHAVIOR DURING TRAINING  ###
####################################

# Load the MedPC files
datafolder = paste(getwd(), "/Raw files/Training/", sep="")
datafiles = list.files(datafolder)

##Use the mpcextract function to extract data from MedPC file###

mpcextract.Rfunc=function(filename){
        rawin = scan(filename, what = "raw", sep = "\n")
        singlevars = grep("^[A-Z]{1}:[[:space:]]{0,}[0-9]", rawin)
        arrayvars = grep("^[A-Z]{1}:$", rawin)
        listout = list(info = as.vector(rawin[c(1:10)]))
        if( length(singlevars) > 0) {
                svarlist = as.numeric(substr(rawin[singlevars], start= 3, stop = 100))
                svarlist = as.list(svarlist)
                names(svarlist) = tolower(substr(rawin[singlevars], start= 1, stop = 1))
        }
        
        if( length(singlevars) == 0){svarlist = list(singlevars = "No single variables found")}
        if(length(arrayvars) > 0){
                for(j in 1:length(arrayvars)){
                        if(j < length(arrayvars)) {iAB = (arrayvars[j]+1):(arrayvars[j+1] - 1)}
                        if(j == length(arrayvars)) {iAB = (arrayvars[j]+1):(length(rawin))}
                        arrayblock = rawin[iAB]
                        avars = substr(arrayblock, start = 8, stop = 100)
                        avars = unlist(strsplit(avars, split = "[[:space:]]{1,}"))
                        avars = as.numeric(avars)
                        avars = avars[!is.na(avars)]
                        if(j == 1 ){avarlist = list(avars)}
                        else {avarlist[[j]] = avars}
                } #End for(j . . .
                names(avarlist) = tolower(substr(rawin[arrayvars], start= 1, stop = 1))
        }# End if(length(arrayvars . . .
        if(length(arrayvars) == 0){avarlist = list(arrayvars = "No array variables found")}
        listout = c(listout, svarlist, avarlist)
        return(listout)
}



#Create empty objects for the idx and the data
csacqidx = data.frame(subject = NA, session = NA, group = NA, filename = datafiles, date = NA)
csacqdata = list()

#Extract data from Med PC
for (i in 1:length(datafiles)){
        toextract = paste(datafolder, datafiles[i], sep = "")
        
        ##EXTRACTION
        rawin =  mpcextract.Rfunc(toextract)
        subject = substring(rawin$info[grep("Subject:",rawin$info)],10, last= 1000000L)
        session = substring(rawin$info[grep("Experiment:",rawin$info)],13, last= 1000000L)
        group = substring(rawin$info[grep("Group:",rawin$info)],7, last= 1000000L)
        filename  = datafiles[i]
        date = substring(rawin$info[grep("Start Date:",rawin$info)],13, last= 1000000L)
        
        medtimes = list(cuestart = rawin[[which(names(rawin) == "y")]],
                        cueend = rawin[[which(names(rawin) == "z")]],
                        reward = rawin[[which(names(rawin) == "t")]],
                        entries = rawin[[which(names(rawin) == "w")]],
                        exits = rawin[[which(names(rawin) == "x")]]
        )
        
        csacqidx$subject[i]= subject
        csacqidx$session[i]= session
        csacqidx$group[i]= group
        csacqidx$date[i]=date
        csacqdata[[i]]=medtimes
}




#####Parse data: create an object with more detailed behavioral data

parseddata = list()
lengthITI = c()

#Establish session bins
sessionbinsize=300
sessionticks=seq(0, 1800, by = sessionbinsize)

#Establish the duration of the pseudolatency window that we'll use to compare the probability/latency of entry during cue. The random ITI is a point during the ITI from the beginning of the ITI to up to 10 sec before the cue starts
wdw=10

#Let's make the dataframe with all the parsed data

for(i in 1:nrow(csacqidx)){
        toanalyze = csacqdata[[i]]
        if(length(toanalyze$cueend)<length(toanalyze$cuestart))(toanalyze$cueend= c(toanalyze$cueend, 1800))
        if(length(toanalyze$exits)<length(toanalyze$ent))(toanalyze$exits= c(toanalyze$exits, 1800))
        
        lengthITI= for (j in 1:length(toanalyze$cuestart)){
                lengthITI=c(lengthITI, toanalyze$cuestart[j+1]-toanalyze$cueend[j])
                lengthITIall=c(toanalyze$cuestart[1], lengthITI[!is.na(lengthITI)])
        }
        
        half <- length(toanalyze$cuestart)/2
        if(half-round(half,0)!=0) {
                divis <- c(rep(1, (floor(half))), rep(2, ceiling(half)))
        }
        if(half-round(half,0)==0) {
                divis <- c(rep(1, half), rep(2, half))
        }
        
        trials = data.frame(bin=findInterval(toanalyze$cuestart, sessionticks), half=divis, trial=c(1:length(toanalyze$cuestart)),start = toanalyze$cuestart,  reward = FALSE, latency = toanalyze$cueend-toanalyze$cuestart, lengthITI = lengthITIall, randomITI= NA, randomITIend= NA, numITIentries=NA)
        
        
        
        trials$reward[which(trials$latency < 10)] = TRUE
        trials$randomITI= toanalyze$cuestart-runif(length(toanalyze$cuestart), min=10, max=lengthITIall)
        trials$randomITIend=trials$randomITI+wdw
        
        
        # Entries data frame
        entries = data.frame(ent = toanalyze$entries, exit = NA, trialno = NA, trialstart = NA)
        entries$exit = toanalyze$exits
        entries$trialno = findInterval(entries$ent,trials$start-15)
        entries$trialno[which(entries$trialno ==0)]=1
        uniquetrial=c(TRUE, entries$trialno!=c(entries$trialno[-1], NA))
        entries$uniquetrial=uniquetrial[1: length(entries$trialno)]
        entries$trialstart = trials$start[entries$trialno]
        entries$diffEnt = entries$ent - entries$trialstart
        entries$diffExit = entries$exit - entries$trialstart
        entries$diffEnt[which(abs(entries$diffEnt)>15)]=NA
        entries$diffExit[which(abs(entries$diffExit)>15)]=NA
        entries$randomITI = trials$randomITI[entries$trialno]
        entries$randomITIend = trials$randomITIend[entries$trialno]
        entries$RandomHE<- entries$ent >= entries$randomITI & entries$ent <= entries$randomITIend & entries$uniquetrial==TRUE
        entries$latITI = entries$ent-entries$randomITI
        entries$latITI [which (entries$RandomHE==FALSE)] = 10
        
        #For previoius data frame
        numITIentries <- c()
        for(j in 1:length(trials$trial)){
                if(trials$reward[j]==F){
                        ITIcount <- length(entries$trialno[entries$trialno==j])
                }
                if(trials$reward[j]==T){
                        ITIcount <- length(entries$trialno[entries$trialno==j])-1  
                }
                numITIentries <- c(numITIentries, ITIcount)
        }
        
        trials$numITIentries=numITIentries
        
        
        # ITI comparison data frame
        ITIcomparison= data.frame (
                bin=findInterval(toanalyze$cuestart, sessionticks),
                half=divis,
                trial= trials$trial, 
                start=trials$start, 
                reward=trials$reward,
                latency=trials$latency,
                RandomHE= FALSE,
                latencyITI= 10,
                TaskAcc=NA)
        trlidx=entries$trialno[which(entries$RandomHE==TRUE)]
        ITIcomparison$RandomHE[trlidx]=TRUE
        ITIcomparison$latencyITI[trlidx]=entries$latITI[entries$RandomHE==TRUE]
        ITIcomparison$TaskAcc <- ITIcomparison$latencyITI-ITIcomparison$latency
        ITIcomparison$ITIentries <- trials$numITIentries
        ITIcomparison$ITIlength <- trials$lengthITI
        
        #ITI entries per second
        TotalITI = data.frame(ITIlengthTotal=NA, ITIentTotal=NA, ITIperSec=NA)
        TotalITI$ITIlengthTotal=sum(trials$lengthITI)
        TotalITI$ITIentTotal=sum(trials$numITIentries)
        TotalITI$ITIperSec=TotalITI$ITIentTotal/TotalITI$ITIlengthTotal
        
        parseddata[[i]]= list(trials = trials, entries = entries, ITIcomparison = ITIcomparison, TotalITI=TotalITI)
}

# ### PROB ENTRY POST VS PRE CUE
# 
# dw <- seq(0, 10, by=0.1)
# 
# IntervalsPreCue <- list()
# IntervalsPostCue <- list()
# IntervalsDF <- list()
# EntryPreCue <- c()
# EntryPostCue <- c()
# probPrePost <- list()
# 
# for(i in 1:nrow(csacqidx)){
#   toanalyze = csacqdata[[i]]
#   if(length(toanalyze$cueend)<length(toanalyze$cuestart))(toanalyze$cueend= c(toanalyze$cueend, 3600))
#   if(length(toanalyze$exits)<length(toanalyze$ent))(toanalyze$exits= c(toanalyze$exits, 3600))
#   
#   #For each of the i sessions in csacqidx, I'm going to make 101 data frames (using a detection window from 0-10s in 0.1s increases gives us 101 intervals to use as detection windows. Each one of these dataframes is going to include: trial (trial number), cue start, Interval Pre cue (the cue start minus the detection window) and interval post cue (cue start plus detection window)). The object with all of these dataframes for each session is going to be called "IntervalsDF" and, in order to subset it, I need a double bracket.
#   
#   probPre <- c()
#   probPost <- c()
#   probPreVSPost <- c()
#   
#   for(j in 1:length(dw)){
#     pre <- toanalyze$cuestart-dw[j]
#     post <- toanalyze$cuestart+dw[j]
#     IntervalsPreCue[j] <- list(pre)
#     IntervalsPostCue[j] <- list(post)
#     IntervalsDF[j] = list(data.frame(
#       trial=1:length(toanalyze$cuestart),
#       cuestart=toanalyze$cuestart,
#       IntervalPreCue=toanalyze$cuestart-dw[j],
#       IntervalPostCue=toanalyze$cuestart+dw[j]
#     ))
#     
#     #Once I have that object (IntervalsDF), I am going to add two columns to it: one saying if there is any entry during the "pre" interval and another one saying if there's any entry during the "post" interval. 
#     
#     testPre2 <- c()
#     testPost2 <- c()
#     
#     for(l in 1:length(IntervalsDF[[j]]$cuestart)){
#       if(length(subset(toanalyze$entries, toanalyze$entries  > IntervalsDF[[j]]$IntervalPreCue[l] & toanalyze$entries <= IntervalsDF[[j]]$cuestart[l]))>=1){testPre <- TRUE} else {testPre <- FALSE}
#       if(length(subset(toanalyze$entries, toanalyze$entries  < IntervalsDF[[j]]$IntervalPostCue[l] & toanalyze$entries >= IntervalsDF[[j]]$cuestart[l]))>=1){testPost <- TRUE} else {testPost <- FALSE}
#       
#       testPre2<- c(testPre2, testPre)
#       testPost2<- c(testPost2, testPost)
#     }
#     
#     IntervalsDF[[j]]$TestPre <- testPre2
#     IntervalsDF[[j]]$TestPost <- testPost2
#     
#     #Now, for each IntervalsDF data frame (remember, I have 101 per session, one per detection window), I'm going to calculate the proportion of trials in which the animal made a head entry given that detection window before (probPre) and after (probPost) the cue.
#     probPre <- c(probPre, length(IntervalsDF[[j]]$TestPre[which(IntervalsDF[[j]]$TestPre==TRUE)])/length(IntervalsDF[[j]]$TestPre))
#     probPost<- c(probPost, length(IntervalsDF[[j]]$TestPost[which(IntervalsDF[[j]]$TestPost==TRUE)])/length(IntervalsDF[[j]]$TestPost))
#   }
#   probPreVSPost <- c(probPreVSPost, probPost-probPre)
#   
#   #Now let's make a data frame with the detection window and the probability of "precue entry" and "postCue entry" for each one of those windows during that session. We'll make a list of data frames (one per session). Note: This might take a while, be patient!  
#   summaryProbPrePost <- data.frame(detectionWindow=dw, probabilityPre=probPre, probabilityPost=probPost, probPreVSPost = probPreVSPost)
#   
#   probPrePost[[i]] <- list(probPrePost = summaryProbPrePost)
# }
# 
# 
# #SaVING THIS OBJECT FOR FUTURE USE (RECOMMENDED SINCE IT TAKES TIME TO COMPUTE)
# save(probPrePost, file = paste(getwd(), "/Data for R/probPrePost.rdat", sep=""))
# 
# 
# #### Now I want to come up with some kind of "Area under the curve" metric. I'm going to add ALL of the values for probPOSTcue-probPREcue and create a vector with those values
# 
# sumProbPreVSPost <- c()
# 
# for(i in 1:length(probPrePost)){
#   sumProbPreVSPost <- c(sumProbPreVSPost, sum(probPrePost[[i]]$probPrePost$probPreVSPost))
# }
# 
# csacqidx$sumProbPreVSPost <- sumProbPreVSPost
# 
# 
# 
# ############################################
# #ANALYSIS 5. 15 LAST MINUTES             ###
# ############################################
# 
# sessionbinsize=300
# sessionticks=seq(0, 1800, by = sessionbinsize)
# findInterval(csacqdata[[1]]$cuestart, sessionticks)
# 
# sessionbybin<- c()
# for(i in 1:length(parseddata)){
#   
#   c(sessionbybin, (data.frame(
#     sessionbin=findInterval(parseddata[[i]]$trials$start, sessionticks), 
#     cuestart=parseddata[[i]]$trials$start,
#     reward=parseddata[[i]]$trials$reward,
#     latency=parseddata[[i]]$trials$latency,
#     lengthITI=parseddata[[i]]$trials$lengthITI,
#     randomITI=parseddata[[i]]$trials$randomITI,
#     randomITIend=parseddata[[i]]$trials$randomITIend
#   )))
# }
# 


##############################################################
#ANALYSIS 1: RATE OF HEAD ENTRIES ACROSS SESSION PER ANIMAL###
##############################################################


#This loop creates the object "rateHE", which is the rate of HE of each item in csacqidx (each session of each animal)
rateHE <- c()
for (i in 1: length(parseddata)){
        rateHE [[i]]<- length(parseddata[[i]]$trials$reward[parseddata[[i]]$trials$reward==TRUE])/length(parseddata[[i]]$trials$reward)
}

#Now I'm adding a column to csacqidx with the rate of HE of each animal each day
csacqidx["RateHE"]<- rateHE



##############################################################################
#ANALYSIS 1.b.: RATE OF HEAD ENTRIES ACROSS SESSION PER ANIMAL LAST 15 MIN!###
##############################################################################

#This allows me to select trials from the half to the end
round(length(parseddata[[1]]$trials$reward)/2):length(parseddata[[1]]$trials$reward)

#Calculate the rate  of HE in those trials
as.numeric(table(parseddata[[1]]$trials$reward[round(length(parseddata[[1]]$trials$reward)/2):length(parseddata[[1]]$trials$reward)])["TRUE"])/(length(parseddata[[1]]$trial$reward)-round(length(parseddata[[1]]$trials$reward)/2))



#now put in a loop
rateHElastHALF <- c()
for (i in 1: length(parseddata)){
        rateHElastHALF [[i]]<- as.numeric(table(parseddata[[i]]$trials$reward[round(length(parseddata[[i]]$trials$reward)/2):length(parseddata[[1]]$trials$reward)])["TRUE"])/(length(parseddata[[i]]$trial$reward)-round(length(parseddata[[1]]$trials$reward)/2))
}

#There's some NA which are animals that didn't make any correct trials in the last half (No TRUE values in parseddata$trials$reward, that's why it comes up as NA)
rateHElastHALF[is.na(rateHElastHALF)] <- 0


#Now I'm adding a column to csacqidx with the rate of HE of each animal each day
csacqidx["RateHElastHALF"]<- rateHElastHALF


######################################################################
#ANALYSIS 2: LATENCY OF CUED HEAD ENTRIES ACROSS SESSION PER ANIMAL###
######################################################################

#Let's make a for loop with that
latencyHE <- c()
latencyHESd <- c()
latencyHESEM <- c()
for (i in 1:length(parseddata)){
        latencyHE [i]<- mean(subset(parseddata[[i]]$trials$latency, parseddata[[i]]$trials$reward==TRUE))
        latencyHESd [i]<- sd(subset(parseddata[[i]]$trials$latency, parseddata[[i]]$trials$reward==TRUE))
        latencyHESEM [i]<- sd(subset(parseddata[[i]]$trials$latency, parseddata[[i]]$trials$reward==TRUE))/sqrt(length(subset(parseddata[[i]]$trials$latency, parseddata[[i]]$trials$reward==TRUE)))
        
}

#Let's add latency and SEM to csacqidx
csacqidx["latencyHE"]<- latencyHE
csacqidx["latencyHESEM"]<- latencyHESEM


#################################################################################
#ANALYSIS 2.b.: LATENCY OF HEAD ENTRIES ACROSS SESSION PER ANIMAL LAST 15 MIN!###
#################################################################################

#This allows me to select trials from the half to the end
round(length(parseddata[[1]]$trials$latency)/2):length(parseddata[[1]]$trials$latency)

#This takes the latency of the last half of the trials, subsets the trials in which the latency was different from 10 (Correct trials) and calculates the mean.
latencyHElastHALF <- mean(parseddata[[1]]$trials$latency[round(length(parseddata[[1]]$trials$latency)/2):length(parseddata[[1]]$trials$latency)][parseddata[[1]]$trials$latency[round(length(parseddata[[1]]$trials$latency)/2):length(parseddata[[1]]$trials$latency)]!=10])


#Let's make a loop
latencyHElastHALF <- c()
latencyHESdlastHALF <- c()
latencyHESEMlastHALF <- c()
for (i in 1:length(parseddata)){
        latencyHElastHALF [i]<- mean(parseddata[[i]]$trials$latency[round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)][parseddata[[i]]$trials$latency[round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)]!=10])
        latencyHESdlastHALF [i]<- sd(parseddata[[i]]$trials$latency[round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)][parseddata[[i]]$trials$latency[round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)]!=10])
        latencyHESEMlastHALF [i]<- sd((parseddata[[i]]$trials$latency[round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)][parseddata[[i]]$trials$latency[round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)]!=10]))/sqrt(length((parseddata[[i]]$trials$latency[round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)][parseddata[[i]]$trials$latency[round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)]!=10])))
}

#Let's add that to csacqidx
csacqidx["latencyHElastHALF"]<- latencyHElastHALF
csacqidx["latencyHESEMlastHALF"]<- latencyHESEMlastHALF



#################################################################
#ANALYSIS 3: LATENCY OF CUED VS. UNCUED ENTRIES               ###
#################################################################



#Let's make a for loop with that
latencyCUEvsITIHE <- c()
latencyCUEvsITIHESd <- c()
latencyCUEvsITIHESEM <- c()
for (i in 1:length(parseddata)){
        latencyCUEvsITIHE [i]<- mean(parseddata[[i]]$ITIcomparison$latencyITI-parseddata[[i]]$ITIcomparison$latency)
        latencyCUEvsITIHESd [i]<- sd(parseddata[[i]]$ITIcomparison$latencyITI-parseddata[[i]]$ITIcomparison$latency)
        latencyCUEvsITIHESEM [i]<- sd(parseddata[[i]]$ITIcomparison$latencyITI-parseddata[[i]]$ITIcomparison$latency)/sqrt(length(parseddata[[i]]$ITIcomparison$latencyITI-parseddata[[i]]$ITIcomparison$latency))
        
}


#Let's add latency and SEM to csacqidx
csacqidx["latencyCUEvsITIHE"]<- latencyCUEvsITIHE
csacqidx["latencyCUEvsITIHESEM"]<- latencyCUEvsITIHESEM


#################################################################
#ANALYSIS 3.b.: LATENCY OF CUED VS. UNCUED ENTRIES LAST HALF  ###
#################################################################

#This allows me to select trials from the half to the end
round(length(parseddata[[1]]$trials$latency)/2):length(parseddata[[1]]$trials$latency)

#This takes the latency of the last half of the trials, subsets the trials in which the latency was different from 10 (Correct trials) and calculates the mean.
latencyCUEvsITIHElastHALF <- mean(parseddata[[1]]$ITIcomparison$latencyITI [round(length(parseddata[[1]]$trials$latency)/2):length(parseddata[[1]]$trials$latency)]-parseddata[[1]]$ITIcomparison$latency[round(length(parseddata[[1]]$trials$latency)/2):length(parseddata[[1]]$trials$latency)])


#Let's make a for loop with that
latencyCUEvsITIHElastHALF <- c()
latencyCUEvsITIHESdlastHALF <- c()
latencyCUEvsITIHESEMlastHALF <- c()
for (i in 1:length(parseddata)){
        latencyCUEvsITIHElastHALF [i]<- mean(parseddata[[i]]$ITIcomparison$latencyITI [round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)]-parseddata[[i]]$ITIcomparison$latency[round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)])
        latencyCUEvsITIHESdlastHALF [i]<- sd(parseddata[[i]]$ITIcomparison$latencyITI [round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)]-parseddata[[i]]$ITIcomparison$latency[round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)])
        latencyCUEvsITIHESEMlastHALF [i]<- sd(parseddata[[i]]$ITIcomparison$latencyITI [round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)]-parseddata[[i]]$ITIcomparison$latency[round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)])/sqrt(length(parseddata[[i]]$ITIcomparison$latencyITI [round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)]-parseddata[[i]]$ITIcomparison$latency[round(length(parseddata[[i]]$trials$latency)/2):length(parseddata[[i]]$trials$latency)]))
        
}


#Let's add latency and SEM to csacqidx
csacqidx["latencyCUEvsITIHElastHALF"]<- latencyCUEvsITIHElastHALF
csacqidx["latencyCUEvsITIHESEMlastHALF"]<- latencyCUEvsITIHESEMlastHALF

######################################################################
#ANALYSIS 4: ITI HEAD ENTRIES                                      ###
######################################################################

#These formulas calculate the AMOUNT of head entries during the ITI (TOTAL AMOUNT OF ENTRIES - REWARDED ENTRIES)

#Let's make a for loop with that
ITIHEnt <- c()

for (i in 1:length(parseddata)){
        ITIHEnt [i]<- length(csacqdata[[i]]$entries)-length(csacqdata[[i]]$reward)
}

#Let's add #of ITI head entries to csacqidx
csacqidx["ITIHEnt"]<- ITIHEnt


## Same thing but per second
ITIHEntPerSec <- c()

for (i in 1:length(parseddata)){
        ITIHEntPerSec[i] <- parseddata[[i]]$TotalITI$ITIperSec
}

#Let's add #of ITI head entries to csacqidx
csacqidx["ITIHEntPerSec"]<- ITIHEntPerSec


#####################################
#ANALYSIS 5: DATA PER BIN         ###
#####################################

rateHEbin <-c()
latencyHEbin <- c()
latencyITIbin <- c()
subjectbin <- c()
sessionbin <- c()
groupbin <- c()

check<- c()

for (j in 1:length(parseddata)){
        for( i in 1:(length(sessionticks)-1)){
                rateHEbin <- c(rateHEbin, as.numeric(table(parseddata[[j]]$trials$reward[parseddata[[j]]$trials$bin==i])["TRUE"])/length(parseddata[[j]]$trial$reward[parseddata[[j]]$trials$bin==i]))
                latencyHEbin <- c(latencyHEbin, mean(subset(parseddata[[j]]$trials$latency[parseddata[[j]]$trials$bin==i], parseddata[[j]]$trials$latency[parseddata[[j]]$trials$bin==i] != (10))))
                latencyITIbin <- c(latencyITIbin, mean(parseddata[[j]]$ITIcomparison$latencyITI[parseddata[[j]]$ITIcomparison$bin==i]))
        }
        subjectbin=c(subjectbin, rep(csacqidx[j,1],6))
        sessionbin=c(sessionbin, rep(csacqidx[j,2],6))
        groupbin=c(groupbin, rep(csacqidx[j,3],6))
}

#There's some NA which are animals that didn't make any correct trials in that bin (No TRUE values in parseddata$trials$reward, that's why it comes up as NA)
rateHEbin[is.na(rateHEbin)] <- 0
latencyHEbin[is.nan(latencyHEbin)] <- 10

#Let's make a data frame with data per bin

dataperbin<- data.frame(
        rat=subjectbin, 
        session=sessionbin, 
        group=groupbin,
        bin=rep(1:6, length(parseddata)), 
        rateHE=rateHEbin, 
        latency=latencyHEbin, 
        latencyITI=latencyITIbin,
        latencyCUEvsITI=latencyITIbin-latencyHEbin)


###############################################
### DATA FOR YOKING ###########################
###############################################

#I want to know, on average, how many CS-Rew pairings and how many CS alone presentations did the Bilateral AP5 group received each day.

idx <- csacqidx$group==" BILATAP5"
BilatAP5cases <- (1:length(csacqidx$group))[idx]

CSRew <- c()
CSMiss <- c()

for(i in 1:length(csacqidx$subject)){
        CSRew[i] <- length(parseddata[[i]]$trials$reward[parseddata[[i]]$trials$reward==T])
        CSMiss[i] <- length(parseddata[[i]]$trials$reward[parseddata[[i]]$trials$reward==F])
}

csacqidx$CSRew <- CSRew
csacqidx$CSMiss <- CSMiss


YokingData <- csacqidx %>% group_by(group, session) %>% 
        summarise(meanCSRew=mean(CSRew), 
                  meanCSMiss=mean(CSMiss), 
                  total=meanCSRew+meanCSMiss,
                  meanRateHE=meanCSRew/total,
                  sem= sd(RateHE)/sqrt(length(RateHE)),
                  n=length(CSRew))

AP5only <- YokingData %>% filter(group==" BILATAP5") 

YokingData$cuealone <- round(AP5only$meanCSMiss, 0)
YokingData$cuerew <- round(AP5only$meanCSRew, 0)
YokingData$extracues <- round(YokingData$cuealone-YokingData$meanCSMiss, 0)


######################################################################
# TOTAL NUMBER OF TRIALS                                           ###
######################################################################

csacqidx$TotalTrials <- csacqidx$CSRew + csacqidx$CSMiss


#############################################################
### LOOKING AT 1st session in halves   (rateHE)           ###
#############################################################

# Number of trials Bilateral AP5 animals on day 1
vec1 <- c()
idx <- 1:length(csacqidx$subject)
sess1BilAP5idx <- idx[csacqidx$session==1 & csacqidx$group==" BILATAP5"]
for(i in sess1BilAP5idx){
        vec1 <- c(vec1, length(parseddata[[i]]$trials$trial))
}

# Number of trials Yoked Veh animals on day 1
vec2 <- c()
idx <- 1:length(csacqidx$subject)
sess1Vehidx <- idx[csacqidx$session==1 & csacqidx$group==" YOKEDVEH"]
for(i in sess1Vehidx){
        vec2 <- c(vec2, length(parseddata[[i]]$trials$trial))
}

# Data in halves

dataPerHalf <- c()
for(i in 1:length(csacqdata)){
        toanalyze <- parseddata[[i]]
        sessHalves <- ddply(toanalyze$ITIcomparison, c("half"), summarise, 
                            N=length(half),
                            meanRateHE=length(reward[which(reward==TRUE)])/length(reward),
                            meanLatency=mean(latency),
                            meanTaskAcc=mean(TaskAcc),
                            meanITIent=mean(ITIentries/ITIlength))
        dataPerHalf[i]=list(sessHalves)
}

HalvesSess1AP5 <- dataPerHalf[sess1BilAP5idx]
HalvesSess1VEH <- dataPerHalf[sess1Vehidx]



#SaVING RAW DATA
save(csacqidx, file=paste(getwd(), "/Data for R/csacqidx.ridx", sep=""))
save(csacqdata, file=paste(getwd(), "/Data for R/csacqdata.rdat", sep=""))
save(dataperbin, file=paste(getwd(), "/Data for R/dataperbin.rdat", sep=""))
save(parseddata, file=paste(getwd(), "/Data for R/parseddata.rdat", sep=""))
save(dataPerHalf,file=paste(getwd(), "/Data for R/dataPerHalf.rdat", sep=""))


# #LOADING RAW DATA
# save(file=paste(getwd(), "/Data for R/csacqidx.ridx", sep=""))
# save(file=paste(getwd(), "/Data for R/csacqdata.rdat", sep=""))
# save(file=paste(getwd(), "/Data for R/dataperbin.rdat", sep=""))
# save(file=paste(getwd(), "/Data for R/parseddata.rdat", sep=""))
# save(file=paste(getwd(), "/Data for R/dataPerHalf.rdat", sep=""))


########## ANOVA #####################

# This is a mixed ANOVA with Group as a between-subjects factor and Session as a within subject factor

aov_group_session_RateHE <- aov (RateHE ~ group*session + Error(subject/session), data=csacqidx)
summary(aov_group_session_RateHE)


