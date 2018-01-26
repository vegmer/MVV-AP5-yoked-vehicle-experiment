#CS PLUS MV ANALYSIS

#Call necessary packages
library("matrixStats")
library("plyr")
library("dplyr")
library("ggplot2")

setwd("C:/Users/mvega/Dropbox/DISSERTATION/")
datafolder = paste(getwd(), "/Experiment 2b/Test/Raw data/", sep="")
datafiles = list.files(datafolder)


##Extract function to extract data from MedPC file###
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
        
        #grep("Subject:",rawin$info) indicates the row in raw$info in which "Subject:" is (i.e. 4)
        #rawin$info[grep("Subject:", rawin$info)] indicates the text in that row
        #rawin$info[grep("Subject:", rawin$info)],10, last=10000000000) indicates to select what is after 10 caracters but before the last 1000000
        
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
        
        #This indicates to create a list in which variable with the name "y" in the object "rawin" is going to be called "cuestart"...
        
        # rm(rawin)
        csacqidx$subject[i]= subject
        csacqidx$session[i]= session
        csacqidx$group[i]= group
        csacqidx$date[i]=date
        csacqdata[[i]]=medtimes
}




#####Parse data
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
        entries$cueend = trials$start[entries$trialno]+trials$latency[entries$trialno]
        entries$diffEnt = entries$ent - entries$trialstart
        entries$diffExit = entries$exit - entries$trialstart
        entries$diffEnt[which(abs(entries$diffEnt)>15)]=NA
        entries$diffExit[which(abs(entries$diffExit)>15)]=NA
        entries$randomITI = trials$randomITI[entries$trialno]
        entries$randomITIend = trials$randomITIend[entries$trialno]
        entries$RandomHE<- entries$ent >= entries$randomITI & entries$ent <= entries$randomITIend & entries$uniquetrial==TRUE
        entries$latITI = entries$ent-entries$randomITI
        entries$latITI [which (entries$RandomHE==FALSE)] = 10
        
        #For previous data frame: # of ITI entries
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
        
        # PreCueDur <- c()
        # PostCueDur <- c()
        # for(j in 1:length(entries$ent)){
        #   if(entries$ent[j]<entries$cueend[j]){PostDur <- 0}
        #   if(entries$ent[j] >= entries$cueend[j] & entries$ent[j]<wdw){
        #     PostDur <- entries$exit[j]-entries$ent[j]
        #     excess <- entries$exit[j]-(entries$cueend[j]+wdw)
        #     if(excess>0){
        #       PostDur <- PostDur-excess
        #     }
        #   }
        #   PostCueDur <- c(PostCueDur, PostDur)
        # }
        # entries$PostCueDur=PostCueDur
        
        
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

### New data frame to analyze the duration of entries during random 10s ITI period vs. after cue end
window <- 10


DurEntries <- list()
for(i in 1:length(csacqdata)){
        dat <- csacqdata[[i]]
        if(length(dat$cueend)<length(dat$cuestart))(dat$cueend= c(dat$cueend, 1800))
        if(length(dat$exits)<length(dat$ent))(dat$exits= c(dat$exits, 1800))
        
        
        ### MAKING A DATA FRAME WITH INFO ABOUT TIME SPENT IN RECEPTACLE BEFORE AND AFTER CUE END
        TimeIn <- data.frame(bin=findInterval(dat$cuestart, sessionticks), trial=1:length(dat$cuestart), cuestart=dat$cuestart, cueend=dat$cueend, postCueInt=dat$cueend+window, preCueInt=dat$cuestart-window, PostCueDur=NA, PostCueDurPerc=NA, PreCueDur=NA, PreCueDurPerc=NA, DurDiff=NA, DurDiffPerc=NA)
        
        #BEFORE CUE END    
        corresp <- findInterval(dat$entries, dat$cueend) #for each element in dat$entries, which "cueend"(trial) they correspond to
        
        DUR <- c()
        for(j in 1:length(dat$entries)){
                if(corresp[j]==0){duration <- 0}
                if(corresp[j]!=0){ 
                        if(dat$entries[j]>=dat$cueend[corresp[j]] & dat$entries[j]<(dat$cueend[corresp[j]]+window)){
                                duration <- dat$exits[j]-dat$entries[j]
                                excess <- dat$exits[j]-(dat$cueend[corresp[j]]+window)
                                if(excess>0){duration <- duration-excess}
                        }
                        else {duration <- 0}}
                DUR <- c(DUR, duration) #Vector with duration of each head entry during the post cue 10s window
        }
        
        PostCueDur <- c()
        for(k in 1:length(dat$cueend)){
                if(sum(corresp==k)>0){
                        D <- sum(DUR[which(corresp==k)])
                } else {D <- 0}
                PostCueDur <- c(PostCueDur, D)
        }
        TimeIn$PostCueDur=PostCueDur
        TimeIn$PostCueDurPerc=PostCueDur/window
        
        
        ### BEFORE CUE START
        
        corresp <- findInterval(dat$entries, dat$cuestart-window) #for each element in dat$entries, which "cueend"(trial) they correspond to
        DUR <- c()
        for(j in 1:length(dat$entries)){
                if(corresp[j]==0){duration <- 0}
                if(corresp[j]!=0){ 
                        if(dat$entries[j]>=(dat$cuestart[corresp[j]]-window) & dat$entries[j]<dat$cuestart[corresp[j]]){
                                duration <- dat$exits[j]-dat$entries[j]
                                excess <- dat$exits[j]-(dat$cuestart[corresp[j]])
                                if(excess>0){duration <- duration-excess}
                        }
                        else {duration <- 0}}
                DUR <- c(DUR, duration) #Vector with duration of each head entry during the pre cue 10s window
        }
        
        PreCueDur <- c()
        for(k in 1:length(dat$cuestart)){
                if(sum(corresp==k)>0){ #If there were entries made during that period
                        D <- sum(DUR[which(corresp==k)])
                } else {D <- 0}
                PreCueDur <- c(PreCueDur, D)
        }
        TimeIn$PreCueDur=PreCueDur
        TimeIn$PreCueDurPerc=PreCueDur/window
        
        TimeIn$DurDiff=PostCueDur-PreCueDur
        TimeIn$DurDiffPerc=TimeIn$PostCueDurPerc-TimeIn$PreCueDurPerc
        
        DurEntries[[i]] <- list(TimeIn=TimeIn)
}

#### CALCULATE SUM OF THE DIFFERENCE BETWEEN TIME SPENT IN COMPT 10S AFTER CUE VS. BEFORE CUE IN THE FIRST BIN OF TRIALS FOR EACH ANIMAL
topBin <- 2 #Bin to analyze (e.g. if topBin=2, only the first bin will)
durDiffEach <- c()
for(i in 1:length(DurEntries)){
        sel <- DurEntries[[i]]$TimeIn %>% filter(bin<topBin)
        calc <- sum(sel$DurDiff)
        durDiffEach <- c(durDiffEach, calc)
}

PreCueEach <- c()
for(i in 1:length(DurEntries)){
        sel <- DurEntries[[i]]$TimeIn %>% filter(bin<topBin)
        calc <- sum(sel$PreCueDur)
        PreCueEach <- c(PreCueEach, calc)
}

PostCueEach <- c()
for(i in 1:length(DurEntries)){
        sel <- DurEntries[[i]]$TimeIn %>% filter(bin<topBin)
        calc <- sum(sel$PostCueDur)
        PostCueEach <- c(PostCueEach, calc)
}

csacqidx$durDiffBin1 <- durDiffEach
csacqidx$PreCueDurSumBin1 <- PreCueEach
csacqidx$PostCueDurSumBin1 <- PostCueEach

#### CALCULATE SUM OF THE DIFFERENCE BETWEEN TIME SPENT IN COMPT 10S AFTER CUE VS. BEFORE CUE IN FIRST 44 TRIALS
minNoTrials <- 44 #This tells us the minimum number of trials all animals got (because some animals got more, and we don't want to count those extra trials since this is a SUM across trials)

durDiffEach <- c()
for(i in 1:length(DurEntries)){
        sel <- DurEntries[[i]]$TimeIn %>% filter(trial<minNoTrials)
        calc <- sum(sel$DurDiff)
        durDiffEach <- c(durDiffEach, calc)
}

PreCueEach <- c()
for(i in 1:length(DurEntries)){
        sel <- DurEntries[[i]]$TimeIn %>% filter(trial<minNoTrials)
        calc <- sum(sel$PreCueDur)
        PreCueEach <- c(PreCueEach, calc)
}

PostCueEach <- c()
for(i in 1:length(DurEntries)){
        sel <- DurEntries[[i]]$TimeIn %>% filter(trial<minNoTrials)
        calc <- sum(sel$PostCueDur)
        PostCueEach <- c(PostCueEach, calc)
}

csacqidx$durDiffSUMALLTRIALS <- durDiffEach
csacqidx$PreCueDurSUMALLTRIALS<- PreCueEach
csacqidx$PostCueDurSUMALLTRIALS <- PostCueEach


### PROB ENTRY POST VS PRE CUE

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
#   if(length(toanalyze$cueend)<length(toanalyze$cuestart))(toanalyze$cueend= c(toanalyze$cueend, 1800))
#   if(length(toanalyze$exits)<length(toanalyze$ent))(toanalyze$exits= c(toanalyze$exits, 1800))
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
# save(probPrePost, file = paste(getwd(), "/Experiment 2b/Test/Data for R/probPrePost.rdat", sep=""))
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
ITIperSecbin <- c()

check<- c()

for (j in 1:length(parseddata)){
        for( i in 1:(length(sessionticks)-1)){
                rateHEbin <- c(rateHEbin, as.numeric(table(parseddata[[j]]$trials$reward[parseddata[[j]]$trials$bin==i])["TRUE"])/length(parseddata[[j]]$trial$reward[parseddata[[j]]$trials$bin==i]))
                latencyHEbin <- c(latencyHEbin, mean(subset(parseddata[[j]]$trials$latency[parseddata[[j]]$trials$bin==i], parseddata[[j]]$trials$latency[parseddata[[j]]$trials$bin==i] != (10))))
                latencyITIbin <- c(latencyITIbin, mean(parseddata[[j]]$ITIcomparison$latencyITI[parseddata[[j]]$ITIcomparison$bin==i], na.rm=T))
                ITIperSecbin <- c(ITIperSecbin, mean(parseddata[[j]]$ITIcomparison$ITIentries[parseddata[[j]]$ITIcomparison$bin==i]/parseddata[[j]]$ITIcomparison$ITIlength[parseddata[[j]]$ITIcomparison$ITIlength]))
        }
        subjectbin=c(subjectbin, rep(csacqidx[j,1],6))
        sessionbin=c(sessionbin, rep(csacqidx[j,2],6))
        groupbin=c(groupbin, rep(csacqidx[j,3],6))
}

#There's some NA which are animals that didn't make any correct trials in that bin (No TRUE values in parseddata$trials$reward, that's why it comes up as NA)
rateHEbin[is.na(rateHEbin)] <- 0

#For latency only analysis, I want NaN in missed trials. But to calculate TAsk accuracy, I need to turn those NaN into 10s
latencyHEbin2 <- latencyHEbin
latencyHEbin2[is.nan(latencyHEbin2)] <- 10

#Let's make a data frame with data per bin

dataperbin<- data.frame(
        rat=subjectbin, 
        session=sessionbin, 
        group=groupbin,
        bin=rep(1:6, length(parseddata)), 
        rateHE=rateHEbin, 
        latency=latencyHEbin, 
        latencyITI=latencyITIbin,
        latencyCUEvsITI=latencyITIbin-latencyHEbin2,
        ITIperSec = ITIperSecbin)

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

YokingData <- ddply(csacqidx, c("group", "session"), summarise, 
                    N=length(CSRew),
                    meanCSRew=round(mean(CSRew), 0),
                    meanCSMiss=round(mean(CSMiss),0))


######################################################################
#ANALYSIS 5b: TOTAL NUMBER OF TRIALS                               ###
######################################################################

csacqidx$TotalTrials <- csacqidx$CSRew + csacqidx$CSMiss


#############################################################
### LOOKING AT 1st session in halves   (rateHE)           ###
#############################################################

# Number of trials Bilateral AP5 animals on day 1
sessionName <- "EXTTTEST"
vec1 <- c()
idx <- 1:length(csacqidx$subject)
sess1BilAP5idx <- idx[csacqidx$session==sessionName & csacqidx$group==" BILATAP5"]
for(i in sess1BilAP5idx){
        vec1 <- c(vec1, length(parseddata[[i]]$trials$trial))
}

# Number of trials Yoked Veh animals on day 1
vec2 <- c()
idx <- 1:length(csacqidx$subject)
sess1Vehidx <- idx[csacqidx$session==sessionName & csacqidx$group==" YOKEDVEH"]
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


#### TEST SESSION BY TRIAL ####

# First, since not all animals got the same number of trials, let's find what was the minimum number of trials and let's only look at that number for all animals

minNoTrials <- min(c(vec1, vec2)) 

DataPerTrial <- c()
for(i in 1:length(parseddata)){
        for(j in 1:minNoTrials){
                DataPerTrial <- rbind(DataPerTrial,
                                      c(csacqidx$subject[i],
                                        csacqidx$group[i],
                                        j,
                                        parseddata[[i]]$trials$reward[parseddata[[i]]$trials$trial==j],
                                        parseddata[[i]]$trials$latency[parseddata[[i]]$trials$trial==j],
                                        parseddata[[i]]$trials$numITIentries[parseddata[[i]]$trials$trial==j],
                                        parseddata[[i]]$trials$lengthITI[parseddata[[i]]$trials$trial==j],
                                        parseddata[[i]]$trials$numITIentries[parseddata[[i]]$trials$trial==j]/parseddata[[i]]$trials$lengthITI[parseddata[[i]]$trials$trial==j]
                                      ))
        }
}

subject<-c()
group <- c()
trial <- c()
reward <- c()
latency <- c()
TaskAcc <- c()
numITIentries <- c()
lengthITI <- c()
ITIperSec <- c()

for(i in 1:length(parseddata)){
        for(j in 1:minNoTrials){
                subject <- c(subject, csacqidx$subject[i])
                group <- c(group,csacqidx$group[i])
                trial <- c(trial, j)
                reward <- c(reward, parseddata[[i]]$trials$reward[parseddata[[i]]$trials$trial==j])
                latency <- c(latency, parseddata[[i]]$trials$latency[parseddata[[i]]$trials$trial==j])
                TaskAcc <- c(TaskAcc, parseddata[[i]]$ITIcomparison$TaskAcc[parseddata[[i]]$trials$trial==j])
                numITIentries <- c(numITIentries, parseddata[[i]]$trials$numITIentries[parseddata[[i]]$trials$trial==j])
                lengthITI <- c(lengthITI, parseddata[[i]]$trials$lengthITI[parseddata[[i]]$trials$trial==j])
                ITIperSec <- c(ITIperSec, parseddata[[i]]$trials$numITIentries[parseddata[[i]]$trials$trial==j]/parseddata[[i]]$trials$lengthITI[parseddata[[i]]$trials$trial==j])
        }
}

DataPerTrial <- data.frame(subject, group, trial, reward, latency, TaskAcc, numITIentries, lengthITI, ITIperSec)

DataPerTrial$latency[which(DataPerTrial$latency==10)] <- NaN


DataPerTrialSumm <- ddply(DataPerTrial, c("group", "trial"), summarise,
                          N=length(subject),
                          meanITIperSec=mean(ITIperSec),
                          sdITI=sd(ITIperSec),
                          seITI=sdITI/sqrt(N),
                          meanRateHE=length(reward[reward==T])/N,
                          sdRateHE=sd(mean(length(reward[reward==T]))),
                          seRateHE=sdRateHE/sqrt(N),
                          meanlatency=mean(latency[which(latency!=10)], na.rm=T),
                          sdLat=sd(latency, na.rm=T),
                          seLat=sdLat/sqrt(N), 
                          meanTaskAcc=mean(TaskAcc),
                          sdTA=sd(TaskAcc),
                          seTA=sdTA/sqrt(N)
)

############################################################################
### DATA FOR PERIEVENT HISTOGRAM (ENTRIES AROUND TIME OF CUE ONSET) ########
############################################################################

# I want a dataframe in which:
## ## - Subject
## ## - Session
## ## - Trial
## ## - Cue onset
## ## - 40 tickmarks of 0.5ms around cue onset (20 before, 20 after)
## ## - Entry Y/N: was an entry made by that animal on that bin before/after cue onset?
## The number of rows is going to be 40*Number of trials*Sessions


binlength <- 1 #in sec
window <- 15 #in sec
tickmarks<- seq(-window, +window, by=binlength)
#This gives us 41 values. We'll drop the last one because 10s is not a bin, is the end of the last bin.
tickmarks<- tickmarks[-length(tickmarks)]


#Let's make our data.frame
subject<- c()
group <- c()
session<- c()
trial<- c()
cueStart<- c()
bin<- c()
entry<- c()

for(i in 1:length(csacqdata)){
        for(j in 1:length(csacqdata[[i]]$cuestart)){
                for(k in 1:length(tickmarks)){
                        
                        subject<- c(subject,csacqidx$subject[i])
                        group <- c(group, csacqidx$group[i])
                        session<- c(session, csacqidx$session[i])
                        trial <- c(trial, j)
                        cueStart <- c(cueStart, parseddata[[i]]$trials$start[j])
                        bin <- c(bin, tickmarks[k])
                        if(length(csacqdata[[i]]$entries[csacqdata[[i]]$entries>(csacqdata[[i]]$cuestart[j]+tickmarks[k]) & csacqdata[[i]]$entries<(csacqdata[[i]]$cuestart[j]+tickmarks[k]+binlength)])==0){check<- FALSE} else {check<-TRUE}
                        entry <- c(entry, check) 
                }
        }
}

PETH <- data.frame(subject=subject, group=group, session=session, trial=trial, cueStart=cueStart, bin=bin, entry=entry)



PETHsumm <- PETH %>% filter(trial>=1 & trial<=5) %>% group_by(group, bin) %>% summarise(ratio=sum(entry)/length(trial))

ggplot(PETHsumm, aes(x=bin, y=ratio, col=group))+geom_line()

#SaVING RAW DATA
save(csacqidx, file="C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/csacqidx.ridx")
save(csacqdata, file="C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/csacqdata.rdat")
save(dataperbin, file="C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/dataperbin.rdat")
save(parseddata, file="C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/parseddata.rdat")
save(dataPerHalf,file="C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/dataPerHalf.rdat" )
save(DataPerTrial,file="C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/DataPerTrial.rdat" )
save(DataPerTrialSumm,file="C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/DataPerTrialSumm.rdat" )
save(PETH,file="C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/PETH.rdat" )


#LOADING RAW DATA
load("C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/csacqidx.ridx")
load("C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/csacqdata.rdat")
load("C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/parseddata.rdat")
load("C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/dataperbin.rdat")
load("C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/dataPerHalf.rdat")
load("C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/DataPerTrial.rdat")
load("C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/DataPerTrialSumm.rdat")
load("C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 2b/Test/Data for R/PETH.rdat")



########## ANOVA #####################

# This is a mixed ANOVA with GROUP as a between-subjects factor and BIN as a within subject factor

# Effects of group and bin on Rate of head entries per bin
aov_group_bin_RateHE <- aov (ITIperSec ~ group*bin + Error(rat/bin), data=dataperbin)
summary(aov_group_bin_RateHE)

# Effects of group and bin on ITI ENTRIES/SEC per bin
aov_group_bin_ITIperSec <- aov (ITIperSec ~ group*bin + Error(rat/bin), data=dataperbin)
summary(aov_group_bin_ITIperSec)

# Effects of group and bin on Latency Post cue vs. Pre Cue per bin
aov_group_bin_TaskAcc <- aov (latencyCUEvsITI ~ group*bin + Error(rat/bin), data=dataperbin)
summary(aov_group_bin_TaskAcc)

# Effects of group on total time spent in receptacle 10s after cue vs. 10s before the cue throughout the whole session

#Test homoscedastity
install.packages("lawstat")
library(lawstat)

with(csacqidx, levene.test(durDiffSUMALLTRIALS, as.factor(group))) #Equal variances TRUE

ttest_group_SumDurEntryDiffALL <- t.test(durDiffSUMALLTRIALS ~ group, paired=FALSE, data=csacqidx, alternative="less", var.equal=TRUE)
ttest_group_SumDurEntryDiffALL


# Effects of group on latency to enter receptacle 10s after cue vs. 10s before the cue (random ITI window) throughout the whole session
with(csacqidx, levene.test(latencyCUEvsITIHE, as.factor(group))) #Equal variances TRUE

ttest_group_TASKACCL <- t.test(latencyCUEvsITIHE ~ group, paired=FALSE, data=csacqidx, var.equal=TRUE)
ttest_group_TASKACCL

########### ADDITIONAL GRAPHS FOR POSTER ###########
colindx1=c("#ffb521", "#dcf0f0")
colindx2=c("#eb7f00", "#539e9e")

#### Sum of diff of time spent inside compartment post cue vs. pre cue (10s window)
plot.new()
plot.window(xlim=c(0.5, 2.5), ylim=c(-10, 40))

abline(h=seq(-10, 40, 10), col="gray75")
axis(side=1, at=c(1, 2), labels=c("AP5/VEH", "VEH/VEH"), font=2, cex.axis=2.2, lwd=0)
axis(side=2, at=seq(-10, 40, 10), font=2, cex.axis=2.2, las=2)

for(i in 1:length(unique(csacqidx$group))){
        sel <- filter(csacqidx, group==unique(csacqidx$group)[[i]])
        rect(xleft = i-0.5, xright = i+0.5, ybottom=0, ytop=mean(sel$durDiffSUMALLTRIALS), col=colindx1[i])
        points(x=rep(i+0.1, nrow(sel)), y=sel$durDiffSUMALLTRIALS, pch=21, col="black", bg=colindx2[i], cex=1.2)
        arrows(x0=i, x1=i, y0=mean(sel$durDiffSUMALLTRIALS), y1=mean(sel$durDiffSUMALLTRIALS)+sd(sel$durDiffSUMALLTRIALS)/sqrt(nrow(sel)), length=0.25, angle=90)
        arrows(x0=i, x1=i, y0=mean(sel$durDiffSUMALLTRIALS), y1=mean(sel$durDiffSUMALLTRIALS)-sd(sel$durDiffSUMALLTRIALS)/sqrt(nrow(sel)), length=0.25, angle=90)
}


for(i in 1:length(unique(csacqidx$group))){
        sel=filter(csacqidx, group==unique(csacqidx$group)[i])
        points(x=rep(i, nrow(sel)), y=sel$durDiffSUMALLTRIALS, pch=21, bg=colindx2[i], col="black", cex=1.3)
}


#### BOXPLOT OF TASK ACCURACY DURING LAST SESSION
pdf (paste(graphdir, "Exp 2b extinction test ACC LAT barplot.pdf", sep=""), width=8, height=7) 

plot.new()
plot.window(xlim=c(0.5, 2.5), ylim=c(-0.4, 2.8))

abline(h=seq(-0.4, 2.8, 0.4), col="gray75")
axis(side=1, at=c(1, 2), labels=c("AP5/VEH", "VEH/VEH"), font=2, cex.axis=2.2, lwd=0)
axis(side=2, at=seq(-0.4, 2.8, 0.4), font=2, cex.axis=2.2, las=2)

boxplot(latencyCUEvsITIHE ~ group, data=csacqidx, col=colindx1, axes=F, add=T)
 
colindx2=c("#eb7f00", "#539e9e")
for(i in 1:length(unique(csacqidx$group))){
         sel=filter(csacqidx, group==unique(csacqidx$group)[i])
         points(x=rep(i, nrow(sel)), y=sel$latencyCUEvsITIHE, pch=21, bg=colindx2[i], col="black", cex=1.3)
 }



for(i in 1:length(unique(csacqidx$group))){
        sel <- filter(csacqidx, group==unique(csacqidx$group)[[i]])
        rect(xleft = i-0.5, xright = i+0.5, ybottom=0, ytop=mean(sel$latencyCUEvsITIHE), col=colindx1[i])
        points(x=rep(i+0.1, nrow(sel)), y=sel$latencyCUEvsITIHE, pch=21, col="black", bg=colindx2[i], cex=1.2)
        arrows(x0=i, x1=i, y0=mean(sel$latencyCUEvsITIHE), y1=mean(sel$latencyCUEvsITIHE)+sd(sel$latencyCUEvsITIHE)/sqrt(nrow(sel)), length=0.25, angle=90)
        arrows(x0=i, x1=i, y0=mean(sel$latencyCUEvsITIHE), y1=mean(sel$latencyCUEvsITIHE)-sd(sel$latencyCUEvsITIHE)/sqrt(nrow(sel)), length=0.25, angle=90)
}
dev.off()
