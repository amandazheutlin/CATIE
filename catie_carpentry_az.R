# Preprocessing Script for Machine Learning CATIE Trial Data
# Authors: D. Siljak and A.M. Chekroud, May-Jul 2015


#### Housekeeping
# Data manipulation and plotting
libs <- c("DMwR","stringr","plyr","dplyr","ggplot2","RColorBrewer",
          "reshape2","rtf","caret","psych","mi","doMC")
lapply(libs, require, character.only = TRUE)

# Parallel Computing (default:use n/2 cores during model build)
registerDoMC(12)

### Directories and seeds
set.seed(1) # Reading from the same sheet
setwd("/data/swe_gwas/ABZ/CATIE/")
workDir <- getwd()
dataDir <- paste0(workDir, "/data/") 

### Functions to help processing the data frames

getplot <- function(column){
  # plots histogram of given column in data frame
  # column: the column we want to plot, needs to be of type numeric
  #error handling
  if(!(is.numeric(column))){
    stop("Cannot plot, column not in correct class")
  }
  r   <- range(column, na.rm=TRUE)[2]-range(column, na.rm=TRUE)[1]
  w   <- r/30
  bw  <- 0
  if(w<1){
    bw <- 1
  }
  else{
    bw <- w
  }
  qplot(column, geom="histogram", binwidth=bw)
}

typeconvert <- function(frame, converts, transf = as.numeric){
  # converts the type of specific columns in a frame
  # Args: 
  #frame: data frame in which columns are to be converted
  #converts: vector of names of columns to be converted
  #transf: data type to transform the column into. Default is set to numeric
  # Returns:
  # df
  frame[converts] = lapply(frame[converts], transf)
  frame
}

updateDict <- function(frame, dict){
  # Extract variables and their explanations into the dictionary
  # Args:
  # frame: data frame with the variables
  # dict: data frame that acts as a dictionary
  # Returns:
  # dict: the updated dictionary
  for(nm in names(frame)){
    vec <- data.frame(variable = nm, meaning = frame[1, nm], stringsAsFactors = FALSE)
    dict = rbind(dict, vec)
  }
  dict
}

defineMe <- function(var, dict){
  # Get a variable's explanation from the dictionary
  # Args:
  # var: string representing the variable 
  # dict: data frame used as the dictionary
  # Returns:
  # the entry in dict that represents the variable's explanation
  if(var %in% dict$variable){
    index=which(dict$variable==var)
    return(dict[index, "meaning"])
  }
  # error handling
  else stop("Error: variable not found in dictionary")
}

mergecols <- function(frame, colnames){
  #helper function for combinecols
  #frame: data frame where the columns are to be combined
  #colnames: vector of names of specific columns
  print(colnames)
  if(length(colnames)==1){
    return(frame)
  }
  newcolname <- paste(str_sub(colnames[1], start=1, end=-1), "-", str_sub(rev(colnames)[1], start=-1, end=-1), sep="")
  frame[newcolname] <- rep(NA, nrow(frame))
  for(i in 1:dim(frame)[1]){
    
    for(nm in rev(colnames)){
      print(nm)
      if(is.na(frame[i,nm])){
        next
      }
      else{
        frame[i, newcolname] <- frame[i, nm]
        break
      }
    }
  }
  target <- which(names(frame) == colnames[1])
  #print(target)
  #frame  <- frame[, c(names(frame)[c(1:target)], newcolname, names(frame)[c((target+1):(dim(frame)[2]-1))])]
  frame[colnames] <- list(NULL)
  frame
}

combinecols <- function(frame, name, start=0, end, step){
  #combines columns in a df that was recoded from long to wide
  #frame: df in which we are combining columns
  #name: name of the column that was recoded
  #start: the minimum value of the variable that the df was recoded according to
  #end: maximum value
  #step: how many columns to combine at a time
  range <- end-start
  mod <- (range+1)%%step
  if(mod!=0){
    range <- range-mod
  }
  sequence  <- seq(start, end, step)
  for(i in sequence){
    tocombine <- NULL
    for(j in 0:(step-1)) {
      tocombine <- c(tocombine, paste(name, "_", toString(i+j), sep=""))
      
    }
    #print(tocombine)
    frame <- mergecols(frame, tocombine)
  }
  
  tocombine <- NULL
  if(mod!=0){
    print("hi")
    for(i in (range+1):(range+mod)){
      #print(i)
      tocombine <- c(tocombine, paste(name, "_", toString(i), sep=""))
    }
    frame <- mergecols(frame, tocombine)
  }
  frame
}

extractdiags <- function(subs){
  # CATIE has many awkwardly coded spreadsheets
  # Extracts diagnoses from the rows of a data frame and their frequencies
  # Args:
  # subs: a data frame that is a subset containing the relevant diagnosis columns
  # Returns:
  # diagstable: a data frame containing the diagnoses in one column and their frequencies in the other
  n <- dim(subs)[2]
  m <- dim(subs)[1]
  diagstable <- data.frame(diagnosis=character(), freq=numeric(), stringsAsFactors=FALSE)
  for (i in 1:n){
    for(j in 1:m){
      if (subs[j, i] %in% diagstable$diagnosis){
        index <- which(diagstable$diagnosis == subs[j,i])
        diagstable[index,] <- c(subs[j,i], as.numeric(diagstable$freq[index]) + 1)
        next
      }
      else{
        diagstable <- rbind(diagstable, data.frame(diagnosis=subs[j,i], freq=1, stringsAsFactors=FALSE))
      }
    }
    
  }
  return(diagstable)
}

getsummary<-function(sourcef, destinationf){
  # Automatically calculate summary stats for all columns in a frame
  # Function outputs summary into a new destination frame
  # Args:
  # sourcef: to-be-summarized data frame 
  # destinationf: output df to put all the summary statistics in
  
  for(i in 1:dim(sourcef)[2]){
    if(!is.numeric(sourcef[,i])){
      next
    }
    su <- (summary(sourcef[,i]))
    destinationf[i, "min"]         <- su[1]
    destinationf[i, "Q1"]          <- su[2]
    destinationf[i, "median"]      <- su[3]
    destinationf[i, "mean"]        <- su[4]
    destinationf[i, "Q3"]          <- su[5]
    destinationf[i, "max"]         <- su[6]
    destinationf[i, "NAs"]         <- su[7]
    destinationf[i, "Prop_NAs"]    <- round((su[7]/14.6), digits=1)
    variance   <-  var(sourcef[,i], na.rm=TRUE)
    destinationf[i, "variance"] <- variance
    ran        <- range(sourcef[,i], na.rm=TRUE)
    destinationf[i, "range"]    <- ran[2]-ran[1]
  }
  return(destinationf)
}

# Intialize empty df for the data dictionary
dict <- data.frame(x=character(), y=character(), stringsAsFactors=FALSE)








### Reading in demographic data about each patient from demo01.csv
### Notes
# Removed columns that contained indicators of whether patient is black vs. not black, etc. 
# since we have a race column with the relevant details. 
# Also took out das1ms because we have a more simplified marital status column
# Removed blanks/junk: demo4a, demo5a, soc, ran002, dema6, dema7, dema9, living1, strata, trt, 
# fseqno, daysrz, version_form, collection_title 
# Converted data to numeric and added variable definitions into dictionary
# Many rows were duplicated, only unique rows kept
# Recoded 'race' column into numeric type. 1=American Indian/Alaska Native, 2=Asian, 3=Black or African American,
# 4=Hawaiian or Pacific Islander 5=More than one race 6= Unknown or not reported 7=White
# Recoded 'curreduccompleted' into numeric. 1=Advanced degree completed [e.g. Ph.D.], 2=Advanced degree courses, not graduated 
# 3=College graduate 4=College graduate and some Master's level 5=Community college or technical school 
# 6=Did not complete high school 7=GED/High school diploma 8=Master's degree completed 9=Some college, did not graduate 
demo <- read.csv("data/demo01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
# Exclude some columns
demo <- subset(demo, select=-c(collection_id, visit, gender, dataset_id, subjectkey, interview_date, white, 
                               black, native, asian, pacific, race_c, hispanic, any_mhx,
                               b1_panss, scid17a, scid18a, scid19, scid20, das1ms, demo4a, demo5a, soc, 
                               prf_term, me2, me3, me3a, ran002, dema6, dema7, dema9, living1, strata, trt, 
                               fseqno, daysrz, version_form, collection_title))
dict <- updateDict(demo, dict)
demo <- demo[2:dim(demo)[1],]
demo$psysoc_89 <- ifelse(demo$psysoc_89 == "Y", 1, 2)

demo$race <- as.factor(demo$race)
levels(demo$race) <- c(6,1,2,3,4,5,6,7) 
# 1:7 - indian, asian, black, hawaiian / PI, multiple races, other/missing, white
demo$race <- as.character(demo$race) %>% as.numeric()

demo$curreduccompleted <- as.factor(demo$curreduccompleted)
levels(demo$curreduccompleted) <- c(9,5,3,7,6,1,4,8,2)
# 1:9 - did not complete HS, some college, college grad, GED/HS diploma
# advanced degree courses, comm college / tech, college grad + some MA, 
# MA completed, other (including advanced degree completed)
demo$curreduccompleted <- as.character(demo$curreduccompleted) %>% as.numeric()

demo$sitetype <- as.factor(demo$sitetype)
levels(demo$sitetype) <- c(2,3,14,4,5,6,7,8,9,10,11,12,13,1,15)
# 1:15: UC;VA, ALL, PN, PP:MC;SH, PP;RO, RO, SH, SH;RO, UC,
# UC;MC, UC;PN; UC;PP, UC;SH, PP, VA
demo$sitetype <- as.character(demo$sitetype) %>% as.numeric()

demo <- typeconvert(demo, names(demo))
demo <- unique(demo)
demo$me0a<-ifelse(is.na(demo$me0a),0,demo$me0a)

### NB: cataract (too few rows) and adverse events data (poorly coded) are unusable.



### Reading in the drug attitude inventory data from dai01.csv
### Notes
# Deleted administrative fields + blanks; converted to numeric; updated dictionary

dai <- read.csv("data/dai01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
dai <- subset(dai, select = -c(collection_id, collection_title, visitid, truncvis, visit, 
                               visday, dataset_id, subjectkey, interview_age, interview_date, fseqno, version_form, gender))
dict <- updateDict(dai, dict)
dai  <- dai[dai$phase_ct=="Pre-Rand",]
dai["phase_ct"] <- list(NULL)
dai <- typeconvert(dai, c("src_subject_id",names(dai)[2:12]))


### Reading in ecg patient data from ecg01.csv
### Notes
# Deleted administrative fields + blanks; converted to numeric; updated dictionary
# Df is messy - requires some work
# First we select the 10 columns containing diagnoses
# We only kept ECG diagnoses that were sufficiently common (5% of patients)
# Create new binary regressors for presence/absence of a specific diagnosis (anywhere in the 10 diagnosis fields)

ecg  <- read.csv("data/ecg01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
ecg  <- subset(ecg, select = -c(collection_title, dataset_id, subjectkey, interview_age, interview_date, gender))
dict <- updateDict(ecg, dict)
ecg  <- ecg[ecg$phase_ct == "Pre-Rand",]
diags <- subset(ecg, select=c(diag1, diag2, diag3, diag4, diag5, diag6, diag7, diag8, diag9, diag10))
alldiags <- extractdiags(diags)
alldiags$freq <- as.numeric(alldiags$freq)
for(i in 1:dim(alldiags)[1]){
  alldiags$freq[i] <- alldiags$freq[i]/15
}
alldiags <- alldiags[alldiags$freq>5,]
alldiags[which(alldiags[,"diagnosis"]=="", arr.ind=TRUE), 1] <- "none"
ecgnew <- subset(ecg, select=-c(diag1, diag2, diag3, diag4, diag5, diag6, diag7, diag8, diag9, diag10))
for(nm in alldiags$diagnosis){
  ecgnew[,nm] <- 0
}
for(nm in names(diags)){
  for(i in 1:dim(ecgnew)[1]){
    if(ecg[i, nm] %in% alldiags$diagnosis){
      ecgnew[i, ecg[i,nm]]=1
    }
  }
}

# ecgnew is the recoded version of ecg with the new diagnosis columns
ecgnew <- typeconvert(ecgnew, c("src_subject_id"))
todel  <- which(ecgnew$visitid==201)
ecgnew <- ecgnew[-todel,]
ecgnew["truncvis"] <- list(NULL)
ecgnew <- ecgnew[order(ecgnew$src_subject_id),]
dupids <- ecgnew[duplicated(ecgnew$src_subject_id), 2]
# Delete rows that were multiple measurements/entries for same patient. Keep the first one
ecgnew <- ecgnew[-c(63, 376, 443, 636, 640, 716, 723, 734, 759, 870, 1215, 1259),]
todel  <- names(ecgnew)[c(3:7, 24:101)]
ecgnew[todel] <- list(NULL)
ecgnew <- typeconvert(ecgnew, names(ecgnew))



### Reading in aims global severity score data from aims01.csv.
### Contains patients' evaluations on the severity of abnormal movements scale
### Notes:
# A lot of junk in this df. Removed *many* blanks variables, and some admin fields

aims <- read.csv("data/aims01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
aims <- subset(aims, select = -c(collection_title, neuro_tonguemove, aims_extrem_score_date1, aims_trunk_score_date1, aims_observ_move_ts, aims_dental13_date1, aims_14_date1, neuro_tonguemovedesc, dataset_id, subjectkey, interview_date, interview_age, gender, food_mouth, rev_teethdentalab, rev_teethdentalabdesc,
                                 neuro_involunmove,neuro_involunmovedesc, neuro_bodymove, neuro_bodymovedesc, neuro_handsmove, neuro_handsmovedesc, neuro_tonguerest,
                                 neuro_tonguerestdesc, neuro_facelegmove, neuro_facelegmovedesc, neuro_bodyprof, neuro_bodyprofdesc, neuro_trunk, neuro_trunkdesc, neuro_gait, 
                                 neuro_gaitdescribe, c_score, c_index, c_obj, c_gca, c_eps, site, trt_grp, etype, respond_detail_oth_spec, relationship, aimsmed_ttl, aimston1,
                                 aimston2, aimston3, aimston4, aimsgait1, aimsgait2, aimsgait3, aimstrem1, aimscw1, aimscw2, aimscw3, aimscw4, aimsmf))


aims_del <- c(names(aims)[c(35:41, 69:75)],"aims_facial_score_date1")
aims[aims_del] <- list(NULL)
dict <- updateDict(aims, dict)
aims <- aims[aims$phase_ct=="Pre-Rand",]
aims <- subset(aims, select=-c(visitid, visday, visit, base1, truncvis, collection_id, phase_ct))
aims <- typeconvert(aims, names(aims))




### Reading in the clinician global impressions file from cgis01.csv.
### The file contains the CGI evaluations for each patient
### This is one of the DVs for the CATIE trial (and thus is converted from long->wide)
### Notes:
# Deleted administrative fields + blanks; converted to numeric; updated dictionary
# Filtered for data only collected during "Pre-Rand" and "Phase1/1A" stages of the study.

cgis   <- read.csv("data/cgis01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
cgis   <- subset(cgis, select=-c(collection_id, visitid, visit, visday, dataset_id, gender, 
                                 subjectkey, collection_title, interview_age, interview_date, 
                                 last1, site, dsmvers, date_updated, session_id, completed))
dict   <- updateDict(cgis, dict)
cgisiv <- cgis[cgis$phase_ct == "Pre-Rand",]
cgisdv <- cgis[cgis$phase_ct == "Phase 1/1A", ]
cgisiv["phase_ct"] <- list(NULL)
cgisdv["phase_ct"] <- list(NULL)
cgis_conv <- names(cgisiv)[1:35]
cgisiv    <- typeconvert(cgisiv, cgis_conv)
cgisdv    <- typeconvert(cgisdv, cgis_conv)
cgisdv["base1"]    <- list(NULL)
cgisiv["base1"]    <- list(NULL)
cgisiv["truncvis"] <- list(NULL)

#write.csv(cgisiv, file = paste(workDir, "cgis_new.csv", sep = ""))

# Convert long -> wide
cgisdv_wide <- dcast(cgisdv, src_subject_id ~ truncvis, mean, value.var="cs01")
for(i in 2:23){
  names(cgisdv_wide)[i] <- paste0("cs01_", toString(i-2))
  dict<-rbind(dict, c(names(cgisdv_wide)[i], defineMe("cs01", dict)))
}

for(nm in names(cgisdv)[4:34]){
  cgisdv_wide_new <- dcast(cgisdv, src_subject_id ~ truncvis, mean, value.var=nm)
  for(i in 2:23){
    names(cgisdv_wide_new)[i] <- paste0(nm,"_",toString(i-2))
    dict <- rbind(dict, c(names(cgisdv_wide_new)[i], defineMe(nm, dict)))
  }
  cgisdv_wide <- full_join(cgisdv_wide, cgisdv_wide_new, by="src_subject_id")
}
rm(cgisdv_wide_new)
cgisdv_wide["truncvis"] <- list(NULL)
# Reorder cgisdv by subject id just for checking that previous code worked
cgisdv <- cgisdv[order(cgisdv$src_subject_id),]




### Reading in the calgary depression score file from clgry01.csv.
### The file contains evaluations of all patients on the calgary depression scale
# Deleted administrative fields + blanks; converted to numeric; kept pre-rand only; updated dictionary

clgry <- read.csv("data/clgry01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
clgry <- subset(clgry, select = -c(collection_title, gender, visitid, visit, 
                                   truncvis, visday, collection_id, dataset_id, subjectkey, 
                                   interview_age, interview_date, last1, c1_calg, base1))
dict  <- updateDict(clgry, dict)
clgry <- clgry[clgry$phase_ct=="Pre-Rand",]
clgry["phase_ct"] <- list(NULL)
clgry_conv <- names(clgry)
clgry <- typeconvert(clgry, clgry_conv)


### Reading in the file dosage from dosage01.csv.
### Contains information on the medications and doses received by each patient
# Deleted administrative fields + blanks; converted to numeric; updated dictionary
# Converted long -> wide since dosage varied throughout stage 1

dosecomp <- read.csv("data/dosecomp01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
dosecomp <- subset(dosecomp, select = -c(collection_title, visitid, visday, visit, 
                                         collection_id, dataset_id, subjectkey, interview_age, interview_date, 
                                         gender, medad11,medad13, max2))
dict     <- updateDict(dosecomp, dict)
dosecomp <- dosecomp[dosecomp$phase_ct=="Phase 1/1A",]
dosecomp["phase_ct"] <- list(NULL)
dosecomp["medad01"]  <- list(NULL)
dosecomp <- typeconvert(dosecomp, c("src_subject_id", "truncvis", names(dosecomp)[3:28]))
dosecomp_wide <- dcast(dosecomp, src_subject_id~truncvis, mean, value.var="medad03")
for(i in 2:23){
  names(dosecomp_wide)[i] <- paste0("medad03_", toString(i-2))
  dict<-rbind(dict, c(names(dosecomp_wide)[i], defineMe("medad03", dict)))
}
for(nm in names(dosecomp)[4:28]){
  dosecomp_wide_new <- dcast(dosecomp, src_subject_id ~ truncvis, mean, value.var=nm)
  for(i in 2:23){
    names(dosecomp_wide_new)[i] <- paste0(nm, "_", toString(i-2))
    dict<-rbind(dict, c(names(dosecomp_wide_new)[i], defineMe(nm, dict)))
  }
  dosecomp_wide <- full_join(dosecomp_wide, dosecomp_wide_new, by="src_subject_id")
}
rm(dosecomp_wide_new)


### Reading in the endphase file from endphase01.csv.
### The file contains reasons why each patient ended a phase, any Rx side-effects, and Rx guesses by both clinicians and patients
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary

endphase <- read.csv("data/endphase01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)

endphase <- subset(endphase, select = -c(collection_title, visitid, visday, visit, protocol, 
                                         collection_id, dataset_id, subjectkey, interview_age, 
                                         interview_date,gender, truncvis))
dict <- updateDict(endphase, dict)
endphase <- endphase[endphase$phase_ct=="Phase 1/1A",]
endphase["phase_ct"] <- list(NULL)
endphase <- typeconvert(endphase, names(endphase))




### Reading in "fint", from fint01.csv.
### This file contains the family/caretaker interviews conducted for each patient
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary

fint <- read.csv("data/fint01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
fint <- subset(fint, select = -c(collection_title, visitid, truncvis, visday, visit, gender, 
                                 collection_id, dataset_id, subjectkey, interview_age, interview_date, protocol))
dict <- updateDict(fint, dict)
fint <- fint[fint$phase_ct=="Pre-Rand",]
fint["phase_ct"] <- list(NULL)
fint <- typeconvert(fint, names(fint))


### Reading in hair drug test file from hair01.csv.
### File contains results of hair drig test
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary
# Where there were QNS for cocaine and THC tests, we imputed the mean for these two.
# Beta values for these variables may not be accurate, but non-parametric models should be fine

hair   <- read.csv("data/hair01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
hair   <- subset(hair, select = -c(collection_title, visit, visitid, visday, truncvis,gender, protocol, 
                                   collection_id, dataset_id, subjectkey, interview_date, interview_age))
dict   <- updateDict(hair, dict)
hair   <- hair[hair$phase_ct=="Pre-Rand",]
hair["phase_ct"] <- list(NULL)
hair   <- typeconvert(hair, names(hair))
i      <- which(hair$src_subject_id==1929)[2] #remove duplicate entry
hair   <- hair[-c(i),]
hair$thc     <- ifelse(!is.na(hair$meth) & is.na(hair$thc), mean(hair$thc, na.rm=TRUE), hair$thc) #impute the mean
hair$cocaine <- ifelse(!is.na(hair$meth) & is.na(hair$cocaine), mean(hair$cocaine, na.rm=TRUE), hair$cocaine) #impute the mean

### Reading in the Insight & Treatment Attitudes Questionnaire file from itaq01.csv.
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary

itaq <- read.csv("data/itaq01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
itaq <- subset(itaq, select = -c(collection_title, visitid, visit, truncvis, visday, gender, collection_id, 
                                 dataset_id, subjectkey, interview_date, interview_age, phase, days_baseline, 
                                 site, substanc, attitude, stress, blue, illness_emo, treatmt, prefer, week))
dict <- updateDict(itaq, dict)
itaq <- itaq[itaq$phase_ct == "Pre-Rand",]
itaq["phase_ct"] <- list(NULL)
itaq <- typeconvert(itaq, names(itaq))

### Read in information on different phases and when the patients ended
# Deleted administrative fields+blank+columns whose contents were combined with another column into a more complete one
# +cols not related to phase1/1a. converted to numeric; updated dictionary
# separated dependent and independent variables
keyvars  <- read.csv("data/keyvars01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
keyvars[names(keyvars)[c(1:3, 5:12, 16:19, 21, 24:25, 28:56, 58:59, 61:62, 64:65, 67:80, 85)]]<-list(NULL)
dict     <- updateDict(keyvars, dict)
keyvars  <- keyvars[2:dim(keyvars)[1], ]
keyvars  <- typeconvert(keyvars, names(keyvars))
dkeyvars <- subset(keyvars, select = c(src_subject_id, e1_day, comp1_1a, dcr1_1, comp_s, es_day, dcr_s))
keyvars  <- subset(keyvars, select = - c(e1_day, comp1_1a, dcr1_1, comp_s, es_day, dcr_s))


### Read in lab results from lab01.csv.
### NB: Required some tweaking, which is explained inline
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary

lab <- read.csv("data/lab01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
lab <- subset(lab, select = -c(collection_title, rstatus, visitid, chgval1, baseval1, result, cvhinorm, cvunits, cvlonorm, cvresult, 
                               colldy, truncvis, visday, visit, collection_id, dataset_id, subjectkey, interview_date, interview_age,
                               gender, testgrp, testcode, collectm, base1bf, last1bf, basval1b, chgval1b, base2f,
                               last2f, baseval2, chgval2, base3f, last3f, baseval3, chgval3, lmealtm, ldosetm, ldosedy, lmealdy, last1f))
dict <- updateDict(lab, dict)
lab["repeat."] <- list(NULL)
lab <- lab[lab$phase_ct == "Pre-Rand",]
lab <- typeconvert(lab, c("src_subject_id", "hinorm", "lonorm", "resn", "base1f"))
# extract tests
tests <- unique(lab$testname)
# remove empty string from test name
tests <- tests[tests!=""]
# extract all the patient ids
ids <- unique(lab$src_subject_id)
# initialize new data frame to contain just the patient IDs as rows and test names as columns
# Recoded the results because many IDs were repeated for each different test result and flag
# very inconsistently so we could not conventionally turn the data frame from long into wide
newlab <- NULL
newlab <- data.frame(src_subject_id = ids)

# fill out the new data frame
for(t in tests){
  newlab[,t] <- NA
}

for(i in 1:dim(lab)[1]){
  ind <- which(newlab$src_subject_id == lab[i,1])
  if(lab[i, "testname"]!="" & !is.na(lab[i, "resn"])){
    newlab[ind, lab[i, "testname"]] <- lab[i,"resn"]
  }
}

#delete columns that are all NA from the new data frame
newlab["Is Pregnancy Test Required?"] <- list(NULL)
emptycols <- character()
for(nm in names(newlab)){
  if(all(is.na(newlab[nm]))){
    emptycols<-c(emptycols, nm)
  }
}
newlab[emptycols] <- list(NULL)


### Read in macarthur assesment file from macvlnce01.csv
### Contains patients' responses to an assesment on violence committed both towards the patient and by the patient
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary

macvlnce <- read.csv("data/macvlnce01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
macvlnce <- subset(macvlnce, select = -c(collection_title, visitid, visday, visit, truncvis, gender, 
                                         collection_id, protocol, dataset_id, subjectkey,
                                         interview_date, interview_age))
dict <- updateDict(macvlnce, dict)
macvlnce <- macvlnce[macvlnce$phase_ct=="Pre-Rand",]
macvlnce["phase_ct"] <- list(NULL)
macvlnce <- typeconvert(macvlnce, names(macvlnce))

### Read in medication info file from med01.csv for each patient.
### Contains information on all the medications patients took in their lifetime
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary
# Required a lot of tweaking
# A lot of drugs were named inconsistently (c.f. sertraline and sertraline hydrochloride). These were combined manually

med <- read.csv("data/med01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
med <- subset(med, select = -c(collection_title, collection_id, dataset_id, subjectkey, interview_date, interview_age, gender, drgterm2, shrtcode, 
                               atc2lev3, atc1lev3, startdy, endday, othmed4, othmed4a, psy_med, daysrz, othmed2, othmed2a, othmed2b, othmed2c, othmed2e))
med[names(med)[37:98]] <- list(NULL)
med[names(med)[37:75]] <- list(NULL)
dict   <- updateDict(med, dict)
# Remove first row which contains the variable defitions
med    <- med[2:dim(med)[1], ]
med    <- typeconvert(med, c(names(med)[1:2], names(med)[4:36]))
# make new frame containing all the columns from the med df except for drugterm and medcode
newmed <- subset(med, select = -c(drugterm, medcode))
# remove duplicates
newmed <- unique(newmed)
# extract vector containing all the patient IDs and sort it
newmedids <- unique(newmed$src_subject_id)
newmedids <- sort(newmedids)
# initialize a df with same column names as newmed
medframe <- as.data.frame(matrix(nrow = length(newmedids), ncol = dim(newmed)[2], dimnames = list(NULL, names(newmed))))
# sort newmed by ID
newmed <- newmed[order(newmed[,1]),]
# fill in the medframe data frame. I created this new data frame because there were duplicate entries in med coded
# in such a way that unique() could not get rid of them. Only option was recode the information in a new df
medframe["src_subject_id"] <- newmedids
for(j in 2:dim(newmed)[2]){
  c <- 1
  for(id in newmedids){
    index <- which(newmed$src_subject_id == id)
    if(length(index)>1){
      val <- NA
      for(el in index){
        if(!is.na(newmed[el,j])){
          val <- newmed[el,j]
        }
        else{
          next
        }
      }
      medframe[c,j] <- val
    }
    else{
      medframe[c,j] <- newmed[index,j]
    }
    c <- c+1
  }
}
# We manually selected medications of interest (only psychotropic meds), listed here as a vector of medication codes 
medcodes  <- c(1011401, 1143201, 1270901, 1319401, 1487001, 228501, 273201, 27401, 285201, 582601, 700501, 724401)
# Extract the name of each drug code  
medsnames <- vector()
# Add a column for each code, named appropriately
for(code in medcodes){
  ind <- which(med$medcode == code)
  medsnames <- c(medsnames, med[ind, "drugterm"])
}
medsnames <- unique(medsnames)
for(nm in medsnames){
  medframe[,nm]<-0
}
# Fill in medframe with binary regressors for each medname
for(i in 1:dim(medframe)[1]){
  index   <- which(med$src_subject_id == medframe[i,1])
  allmeds <- med$drugterm[index]
  for(nm in medsnames){
    if(nm %in% allmeds){
      medframe[i,nm] <- 1
    }
  }
}
medframe <- plyr::rename(medframe, c("SERTRALINE HYDROCHLORIDE"="SERTRALINE_HYDROCHLORIDE",  "ZIPRASIDONE HYDROCHLORIDE"= "ZIPRASIDONE_HYDROCHLORIDE",  
                                     "VALPROATE SEMISODIUM" =  "VALPROATE_SEMISODIUM",  "VALPROIC ACID" = "VALPROIC_ACID", "CITALOPRAM HYDROBROMIDE"="CITALOPRAM_HYDROBROMIDE",
                                     "FLUOXETINE HYDROCHLORIDE"="FLUOXETINE_HYDROCHLORIDE", "AMFEBUTAMONE HYDROCHLORIDE"="AMFEBUTAMONE_HYDROCHLORIDE"))
medframe$fluoxetine   <- ifelse(medframe$FLUOXETINE_HYDROCHLORIDE == 1 | medframe$FLUOXETINE == 1, 1, 0)
medframe$sertraline   <- ifelse(medframe$SERTRALINE_HYDROCHLORIDE == 1 | medframe$SERTRALINE == 1, 1, 0)
medframe$valproate    <- ifelse(medframe$VALPROATE_SEMISODIUM == 1 | medframe$VALPROIC_ACID == 1, 1, 0)
medframe$citalopram   <- ifelse(medframe$CITALOPRAM_HYDROBROMIDE == 1 | medframe$CITALOPRAM == 1, 1, 0)
medframe$amfebutamone <- ifelse(medframe$AMFEBUTAMONE_HYDROCHLORIDE == 1 | medframe$AMFEBUTAMONE == 1, 1, 0)
medframe$ziprasidone  <- ifelse(medframe$ZIPRASIDONE_HYDROCHLORIDE == 1 | medframe$ZIPRASIDONE == 1, 1, 0)

medframe[c("SERTRALINE_HYDROCHLORIDE", "SERTRALINE", "ZIPRASIDONE_HYDROCHLORIDE", "ZIPRASIDONE", "VALPROATE_SEMISODIUM", "VALPROIC_ACID",
           "CITALOPRAM_HYDROBROMIDE", "CITALOPRAM", "AMFEBUTAMONE_HYDROCHLORIDE", "AMFEBUTAMONE", "FLUOXETINE_HYDROCHLORIDE", "FLUOXETINE")] <- list(NULL)

# Replace NA's in flag columns with 0's
flags<-names(medframe)[c(3:19, 24:26)]
for(flag in flags){
  medframe[,flag] <- ifelse(is.na(medframe[,flag]), 0, medframe[,flag])
}






### Read in medication dispensation information for the study from meddispn01.csv
### This df tells you the drug that the patient was randomized to
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary

meddispn <- read.csv("data/meddispn01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
meddispn <- subset(meddispn, select = -c(collection_title, collection_id, gender, dataset_id, protocol,        
                                         subjectkey, interview_date, interview_age, copyid))
dict <- updateDict(meddispn, dict)
meddispn <- meddispn[meddispn$phase_ct=="Pre-Rand",]
meddispn <- typeconvert(meddispn, c(names(meddispn)[6:13],"src_subject_id","truncvis", "visitid"))
meddispn <- unique(meddispn)

# ###Read in subj df from ndar_subject-1.csv
# subj <- read.csv("data/ndar_subject01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
# #delete empty/irrelevant fields
# subj <- subset(subj, select = -c(collection_title, collection_id, gender, dataset_id, subjectkey, interview_date, interview_age,
#                                twins_study, sibling_study, family_study, family_user_def_id, subjectkey_mother, src_mother_id,
#                                subjectkey_father, src_father_id, subjectkey_sibling1, src_sibling1_id, sibling_type1, subjectkey_sibling2,
#                                src_sibling2_id, sibling_type2, subjectkey_sibling3, src_sibling3_id, sibling_type3, subjectkey_sibling4,
#                                src_sibling4_id, sibling_type4, zygosity, patient_id_biorepository, cell_id_original, cell_id_biorepository,
#                                agre_subject_id))
# subj_del <- names(subj)[11:20]
# subj[subj_del] <- list(NULL)
# dict <- updateDict(subj, dict)
# #remove first row
# subj <- subj[2:dim(subj)[1],]
# #convert columns to numeric
# subj <- typeconvert(subj, c("src_subject_id", "ethnic_group", "sample_id_biorepository"))



### Read in neurobatt file from neurobatt01.csv.
### Contains results of neurocognitive tests (details in documentation)
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary

neurobatt <- read.csv("data/neurobatt01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
neurobatt <- subset(neurobatt, select = -c(collection_title, visitid, visit, visday, truncvis, protocol, collection_id, gender, dataset_id, subjectkey, interview_date, interview_age))
dict      <- updateDict(neurobatt, dict)
neurobatt <- neurobatt[neurobatt$phase_ct=="Pre-Rand",]
neurobatt["phase_ct"] <- list(NULL)
neurobatt <- typeconvert(neurobatt, names(neurobatt))
#delete one duplicate entry for a patient
i         <- which(neurobatt$src_subject_id==1466)[2]
neurobatt <- neurobatt[-c(i),]

### Read in PANSS file from panss01.csv
### This will be the basis of the classification target
### Therefore it is converted from long -> wide
# Deleted administrative fields + blanks; converted to numeric; updated dictionary

panss     <- read.csv("data/panss01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
panss     <- subset(panss, select=-c(collection_id, gender, dataset_id, subjectkey, interview_date, interview_age))
panss_del <- names(panss)[2:227]
panss[panss_del] <- list(NULL)
panss_del <- names(panss)[59:78]
panss[panss_del] <- list(NULL)
panss     <- subset(panss, select=-c(visit, visitid, visday, last1))
dict      <- updateDict(panss, dict)
#Filter for phases we wanted
panssdv   <- panss[panss$phase_ct=="Phase 1/1A",]
panss     <- panss[panss$phase_ct=="Pre-Rand",] 
#panss is the df that contains information from the "Pre-Rand" phase (independent variables) and panssdv is the dependent variable one 
panss[c("phase_ct", "base4")] <- list(NULL)
panssdv[c("phase_ct", "base4")] <- list(NULL)
panss_conv  <- names(panss)
panss       <- typeconvert(panss, panss_conv)
panssdv     <- typeconvert(panssdv, panss_conv)

write.csv(panss, file = paste(dataDir, "panss_new.csv", sep = ""))

# move around columns so that truncvis is the second column (more convenient for when we join all the dfs)
col_idx <- grep("truncvis", names(panssdv))
panssdv <- panssdv[, c(1,col_idx,2:(col_idx-1),(col_idx+1):dim(panssdv)[2])]

panssdv_wide <- NULL
panssdv_wide <- dcast(panssdv, src_subject_id~truncvis, mean, value.var="panss_general")
for(i in 2:22){
  names(panssdv_wide)[i] <- paste0("panss_general_", toString(i-2))
  dict<-rbind(dict, c(names(panssdv_wide)[i], defineMe("panss_general", dict)))
}

for(nm in names(panssdv)[4:52]){
  panssdv_wide_new <- dcast(panssdv, src_subject_id ~ truncvis, mean, value.var=nm)
  for(i in 2:22){
    names(panssdv_wide_new)[i] <- paste0(nm, "_", toString(i-2))
    dict<-rbind(dict, c(names(panssdv_wide_new)[i], defineMe(nm, dict))) 
  }
  panssdv_wide <- full_join(panssdv_wide, panssdv_wide_new, by="src_subject_id")
}
rm(panssdv_wide_new)

# Create a new data frame containing the subject ID's of each patient as well as the last day (i.e largest truncvis number)
# That patient was evaluated on the panss scale during this time period, in order to use that as our time variable
idvec <- unique(panssdv$src_subject_id)
lastTime  <- data.frame(src_subject_id=idvec, truncvis = c(rep(0,length(idvec))))
for(id in idvec){
  index1 <- which(panssdv$src_subject_id==id)
  max    <- 0
  for(ind in index1){
    if(panssdv$truncvis[ind]>max){
      max <- panssdv$truncvis[ind]
    }
  }
  index2  <- which(lastTime$src_subject_id==id)
  lastTime[index2,2] <- max
}

panss["truncvis"]   <- list(NULL)
panssdv["truncvis"] <- list(NULL)




### Read in the quality of life questionnaire file from qol01.csv. 
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary
qol   <- read.csv("data/qol01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
qol   <- subset(qol, select = -c(collection_title, visitid, c1_qol, c1_qrel, c1_qins, c1_qpsy, visit, visday, truncvis, dataset_id, subjectkey, interview_date, interview_age, gender, last1))
dict  <- updateDict(qol, dict)
qol   <- qol[qol$phase_ct=="Pre-Rand",]
qol["phase_ct"] <- list(NULL)
qol   <- typeconvert(qol, names(qol))



### Read in the scid/psychiatric history df from scid_ph01.csv.
### Contains participants' responses to the Structured Clinical Interview
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary
# Manually recoded NAs as 0's where they depended on a previous question (answer not NA, answer == 0)

scid <- read.csv("data/scid_ph01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
scid <- subset(scid, select=-c(collection_title, collection_id, visitid, visit, visday, truncvis, dataset_id, 
                               subjectkey, interview_date, interview_age, gender, protocol, aescode, bafinfo,
                               comments_mics, ld_othrs, md_othrs, mda_othrs, q048_sza_life, fseqno, scda15, 
                               scda17, scda13, scda14, version_form,
                               daysrz, scda1, lda_othrs, scda2, scda16, yrs_ause))
dict <- updateDict(scid, dict)
scid <- scid[scid$phase_ct=="Pre-Rand",]
scid["phase_ct"] <- list(NULL)
scid <- typeconvert(scid, names(scid))
scid$scid02  <- ifelse(is.na(scid$scid02) & scid$scid01 == 0, 0, scid$scid02)
scid$scid04  <- ifelse(is.na(scid$scid04) & scid$scid03 == 0, 0, scid$scid04)
scid$scid06  <- ifelse(is.na(scid$scid06) & scid$scid05 == 0, 0, scid$scid06)
scid$scid08  <- ifelse(is.na(scid$scid08) & scid$scid07 == 0, 0, scid$scid08)
scid$scid10  <- ifelse(is.na(scid$scid10) & scid$scid09 == 0, 0, scid$scid10)
scid$scid12  <- ifelse(is.na(scid$scid12) & scid$scid11 == 0, 0, scid$scid12)
scid$scid14  <- ifelse(is.na(scid$scid14) & scid$scid13 == 0, 0, scid$scid14)
scid$scid16  <- ifelse(is.na(scid$scid16) & scid$scid15 == 0, 0, scid$scid16)
scid$scid21a <- ifelse(is.na(scid$scid21a) & scid$scid21 == 0, 0, scid$scid21a)
scid$scid17a <- ifelse(is.na(scid$scid17a) & scid$scid17 == 0, 0, scid$scid17a)

### Read in screen data frame from screen01.csv.
### Contains results of screening for comorbid disorders and TD (Tardive dyskinesia)
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary

screen <- read.csv("data/screen01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
screen <- subset(screen, select = -c(collection_title, visitid, visit, truncvis, visday, protocol, collection_id, dataset_id, 
                                     subjectkey, interview_date, interview_age, gender,
                                     missed, attended, discon, level, week, crcid, clinid, clinicid, update, hrsd_r, qccur_r, qscur_r, 
                                     menop, medication1_name, medication1_dosage, medication2_name, medication2_dosage, medication3_name, 
                                     medication3_dosage, medication4_name, medication4_dosage, atmtsuic, risksuic, fnaid, stsd1,
                                     stpdate1, stsd2, stpdate2, stsd3, stpdate3, stsd4, stpdate4, mpdat))
dict    <- updateDict(screen, dict)
screen  <- screen[screen$phase_ct == "Pre-Rand",]
screen["phase_ct"] <- list(NULL)
screen  <- typeconvert(screen, names(screen))



###Read in the short form health survey data frame from sf1201.csv
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary

sf   <- read.csv("data/sf1201.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
#delete c1_sfphy, c1_sfmtl, last1 since they were indicators for end of phase1
sf   <- subset(sf, select=-c(collection_title, collection_id, visitid, truncvis, visday, visit, 
                             dataset_id, subjectkey, interview_date, interview_age, gender, 
                             c1_sfphy, c1_sfmtl, base1, last1))
dict <- updateDict(sf, dict)
sf   <- sf[sf$phase_ct=="Pre-Rand",]
sf["phase_ct"] <- list(NULL)
sf   <- typeconvert(sf, names(sf))








### Read in service utilization data frame from surf01.csv.
### Contains information of how/if participants use services such as hospitals, doctors, nursing homes etc.
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary

surf <- read.csv("data/surf01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
surf <- subset(surf, select = -c(collection_title, collection_id, dataset_id, protocol, subjectkey, interview_date, fseqno,
                                 interview_age, gender, typecode, avgmin, numvisit))
surf_del <- names(surf)[c(2,3,4,6, 14:128)]
surf[surf_del] <- list(NULL)
dict <- updateDict(surf, dict)
surf <- surf[surf$phase_ct=="Pre-Rand",]
surf$phase_ct <- NULL
surf <- unique(surf)
surf <- typeconvert(surf, c(names(surf)))




### Read in surfq data frame from surfq01.csv.
### Contains details of participants' financial situation (wages, benefits, whether supported by family), 
### and any history of legal troubles, details of their medical plan, coverage
# Deleted administrative fields + blanks; converted to numeric; kept only Pre-rand; updated dictionary

surfq <- read.csv("data/surfq01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
surfq <- subset(surfq, select = -c(protocol, visitid, visit, truncvis, visday, collection_title, collection_id, dataset_id, subjectkey, interview_date, interview_age, gender, fseqno, version_form))
dict  <- updateDict(surfq, dict)
surfq <- surfq[surfq$phase_ct=="Pre-Rand",]
surfq["phase_ct"] <- list(NULL)
surfq["base1"]    <- list(NULL)
surfq <- typeconvert(surfq, names(surfq))




### Read in timeto df from timeto01.csv.
### Contains details on after how long a phase was finished (if at all), why it was finished etc.
### Very interesting data!
# Deleted administrative fields + blanks; converted to numeric; updated dictionary

timeto <- read.csv("data/timeto01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
#delete empty/unnecessary fields as well as any fields relating to phases other than phase1/1A
timeto <- subset(timeto, select = -c(collection_title, collection_id, dataset_id, subjectkey, interview_date, interview_age, gender,
                                     tmdisc1b, tmdisc2, tmdisc3, censor1b, censor2, trncdr1b, trncdur2, dcr_eff2, dcr_tae2, dcr_pat2,
                                     dcr_pgo2, dcr_pst2, b1b_day, b2_day, b3_day, enddy_1b, enddy_2, enddy_3))
dict   <- updateDict(timeto, dict)
timeto <- timeto[2:dim(timeto)[1],]
timeto <- typeconvert(timeto, names(timeto)[1:14])
timeto <- timeto[order(timeto[,1]),] #reorder





### Read in patients' vitals df from vitals01.csv
### Contains details patients general medical checkup information (height/weight/bp/bmi etc.)
# Deleted administrative fields + blanks; converted to numeric; updated dictionary

vitals      <- read.csv("data/vitals01.csv", header = TRUE, sep = ",", quote = "\"'", stringsAsFactors = FALSE)
vitals      <- subset(vitals, select = -c(collection_title, collection_id, dataset_id, subjectkey, interview_date, interview_age, gender))
vitals_del  <- names(vitals)[2:21]
vitals[vitals_del] <- list(NULL)
vitals_del  <- names(vitals)[10:16]
vitals[vitals_del] <- list(NULL)
vitals_del  <- names(vitals)[22:28]
vitals[vitals_del] <- list(NULL)
vitals_del  <- names(vitals)[25:138]
vitals[vitals_del] <- list(NULL)
vitals_del  <- c("c_wt", "p_gain")
vitals[vitals_del] <- list(NULL)
vitals_del  <- "vital_pulse"
vitals[vitals_del] <- list(NULL)
dict        <- updateDict(vitals, dict)
vitals      <- vitals[vitals$phase_ct=="Pre-Rand",]
vitals_del  <- c("visitid", "visit", "truncvis", "visday", "phase_ct", "base1")
vitals[vitals_del] <- list(NULL)
vitals      <- typeconvert(vitals, names(vitals))






############## End of data carpentry
############## Next we prepare to export the data dictionary





# Remove duplicate entries from dictionary
dict       <- unique(dict)
dictdups   <- dict[duplicated(dict$variable), 1]
for(dup in dictdups){
  i<-which(dict$variable==dup)[2]
  dict<-dict[-c(i),]
}


#generate rtf of the finalized dictionary
rtfdoc     <- RTF("DZ_data_dictionary.rtf", width=10, height=11, font.size=10, omi=c(.5,.5,.5,.5))
addTable(rtfdoc, dict, col.widths=c(3,6), font.size=10, row.names=F, NA.string="-")
done(rtfdoc)

### Initialize the df which we will form by joining all the data frames we will be using for analysis
wholeset <- NULL

### Join the dfs
wholeset <- full_join(panssdv_wide, cgisdv_wide, by="src_subject_id")
wholeset <- full_join(wholeset, lastTime, by="src_subject_id")

#  move columns around so truncvis is the second column
col_idx  <- grep("truncvis", names(wholeset))
wholeset <- wholeset[, c(col_idx, (1:ncol(wholeset))[-col_idx])]
col_idx  <- grep("src_subject_id", names(wholeset))
wholeset <- wholeset[, c(col_idx, (1:ncol(wholeset))[-col_idx])]
wholeset <- full_join(wholeset, dkeyvars, by="src_subject_id")
wholeset <- full_join(wholeset, timeto, by="src_subject_id")
wholeset <- full_join(wholeset, endphase, by="src_subject_id")

# when we join data frames, sometimes there are variables named the same
# R will then rename the two into Column.X and Column.Y
# This is a little error handling tool we use after each join
# If same-name columns are ever joined, this loop will flag it to us.

for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}

wholeset <- full_join(wholeset, keyvars, by="src_subject_id")

wholeset <- full_join(wholeset, panss, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, medframe, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x", nm) | grepl("\\.y", nm))
    print(nm)
}
wholeset <- full_join(wholeset, newlab, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, aims, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, cgisiv, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, dai, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, clgry, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, demo, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, dosecomp_wide, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, ecgnew, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, hair, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, itaq, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, macvlnce, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, neurobatt, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, scid, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, screen, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, sf, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, surfq, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}
wholeset <- full_join(wholeset, vitals, by="src_subject_id")
for(nm in names(wholeset)){
  if(grepl("\\.x\\b", nm) | grepl("\\.y\\b", nm))
    print(nm)
}


# Remove a few columns that were all NA but slipped through initial screening and a few unnecessary columns
wholeset[c("base1", "rtd13", "rtd14", "b2_score", "b2_index", "b2_obj", "b2_gca", "b2_eps", "collection_id", "base2_ecg", "base3_ecg",
           "b2_qt", "b3_qt", "colldy", "bat01", "bat02", "bat03", "bat04", "bat05", "bat06", "bat07", "bat08", "bat09", "bat10", 
           "bat11", "bat12", "white", "black", "hispanic", "native", "rescreen", "b2_wt")]<-list(NULL)

for(nm in names(wholeset)){
  if(all(is.na(wholeset[nm]))){
    wholeset[nm] <- list(NULL)
  }
}


#remove some columns
for(i in c(0:8)){
  nm <- paste("panss_total_", toString(i), sep="")
  wholeset[nm] <- list(NULL)
}

for(i in c(0, 2:20)){
  nm <- paste("b1_panss_", toString(i), sep="")
  wholeset[nm] <- list(NULL)
}

#combine panss columns
wholeset <- combinecols(wholeset, "panss_total", 9, 20, 12)
wholeset$Olanzapine   <- ifelse(wholeset$treat11a==1, 1, 0)
wholeset$Quetiapine   <- ifelse(wholeset$treat11a==2, 1, 0)
wholeset$Risperidone  <- ifelse(wholeset$treat11a==3, 1, 0)
wholeset$Ziprasidone  <- ifelse(wholeset$treat11a==4, 1, 0)
wholeset$Perphenazine <- ifelse(wholeset$treat11a==5, 1, 0)

#### Write the data frame as a CSV 
write.csv(wholeset, file = paste(dataDir, "CATIE_proc.csv", sep = ""))

# initialize and generate data summary sheet
summarystats <- data.frame(nm=names(wholeset), meaning=rep(NA, ncol(wholeset)), min=rep(0, ncol(wholeset)),
                           Q1=rep(0, ncol(wholeset)), median=rep(0, ncol(wholeset)),
                           mean=rep(0, ncol(wholeset)), Q3=rep(0, ncol(wholeset)),
                           max=rep(0, ncol(wholeset)),NAs=rep(0, ncol(wholeset)),
                           variance=rep(0, ncol(wholeset)), range=rep(0, ncol(wholeset)), stringsAsFactors=FALSE)
summarystats <- getsummary(wholeset, summarystats)
summarystats[1, c("NAs","Prop_NAs")] <- 0


sharedNames <- character()
summarystats$NAs      <- ifelse(is.na(summarystats$NAs),0,summarystats$NAs)
summarystats$Prop_NAs <- ifelse(is.na(summarystats$Prop_NAs),0,summarystats$Prop_NAs)

for(nm in summarystats$nm){
  if(!(inherits(try(defineMe(nm, dict)), 'try-error'))){
    i <- which(summarystats$nm == nm)
    summarystats[i, "meaning"] <- defineMe(nm, dict)
    if(summarystats[i, "Prop_NAs"] < 20){
      sharedNames <- c(sharedNames, nm)
    }
  }
}


write.csv(summarystats, file = paste(workDir, "/CATIE_summarystats.csv", sep = ""))

#filter only 15% and less of NAs
#remove dependent variables before we impute
subs <- summarystats[summarystats$Prop_NAs<15,]$nm
good_vars <-  wholeset[,subs]

#generate summary stats of the variables with <15% NAs with describe
describe_summary_stats <- describe(good_vars)

#write to a CSV
write.csv(describe_summary_stats, file = paste(workDir, "/CATIEdescribe_summary_stats.csv", sep = ""))

#remove near-zero variance
nzv <- nearZeroVar(good_vars)
good_vars <- good_vars[-nzv]
dim(good_vars)

independent_names <- names(good_vars)[c(93, 100:512)]
dependent_names <- setdiff(names(good_vars), independent_names)

ind_vars <- good_vars[independent_names]
dep_vars <- good_vars[dependent_names]



#numNAs <- apply(ind_vars, 1, function(z) sum(is.na(z))/dim(ind_vars)[1])
#ind_vars <- ind_vars[numNAs<.5,]
contnms <- character()

for(x in names(ind_vars)) {
  if(dim(table(ind_vars[x]))>2){
    contnms <- c(contnms, x)
  }
}

contnms <- contnms[contnms!="cut_pans"]
continous <- ind_vars[,contnms]
categorical <- ind_vars
categorical[contnms] <- list(NULL)
categorical$src_subject_id <- wholeset$src_subject_id
continous$src_subject_id <- wholeset$src_subject_id

for(x in names(categorical)){
  categorical[,x] <- ifelse(is.na(categorical[,x]), -1, categorical[,x])
}

continous <- typeconvert(continous)

#mean imputation
for(nm in names(continous)){
  continous[,nm] <- ifelse(is.na(continous[,nm]), mean(continous[,nm], na.rm=TRUE), continous[,nm])
}

imputed_data <- merge(continous, categorical, by.x="src_subject_id", by.y="src_subject_id")



### write csv of imputed data with under 15% NAs and the dependent vars (not imputed)
write.csv(imputed_data, file = paste(workDir, "/CATIE_imputed_data.csv", sep = ""))
write.csv(dep_vars, file = paste(workDir, "/CATIE_dep_vars.csv", sep = ""))
### No point sharing full data dictionary when many variables have 90% missing data
### Here we make a smaller one for sharing with collaborators


shareableSummary <- summarystats[summarystats$Prop_NAs < 15,]

write.csv(shareableSummary, file = paste(workDir, "/CATIE_indSummary.csv", sep = ""))

ind_vars_summary <- shareableSummary[shareableSummary$nm %in% names(ind_vars),]

write.csv(ind_vars_summary, file = paste(workDir, "/CATIE_indSummary.csv", sep = ""))

shareableDictionary <- NULL
shareableDictionary <- dict[dict$variable %in% sharedNames,]

# Write the summary stats spreadsheet
write.csv(shareableSummary, file = paste(workDir, "/CATIE_shareableSummary.csv", sep = ""))
# Write the data dictionary for variables with less than 20% missing data
shareDictRTF     <- RTF("shareable_dictionary.rtf", width=10, height=11, font.size=10, omi=c(.5,.5,.5,.5))
addTable(shareDictRTF, shareableDictionary, col.widths=c(3,6), font.size=10, row.names=F, NA.string="-")
done(shareDictRTF)

