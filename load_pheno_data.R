#Louise CHORFI, 25/05/23
#load pheotyping datas 
#run with: source("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell numeration panels/load_pheno_data.R")

setwd("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell numeration panels")
input_dir="//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell numeration panels/SEP - tableau par tube"

options(stringsAsFactors = FALSE)
options(error = traceback)
options(show.error.locations = TRUE)
library("plyr")       
library("dplyr")        
library(tidyverse)

# --------- --Import datas from excel ------------
#### numerations (TBNK) ####
num=data.frame()
num <- read.table("SEP-numerations.csv", sep=";",dec=",",stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                  check.names = FALSE, fileEncoding = "UTF-8")
names(num)[names(num)=="ID"]="Sample_name"

num <- num[,c(1,2,3,5,7,9,11)]

#### Tube 1 : panel for phenotyping main immune cells populations ####
T1=data.frame()
T1 <- read.table(paste0(input_dir, "/SEP_C3_T1utf8.csv"), sep=";",dec=",",stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                 check.names = FALSE, fileEncoding = "UTF-8")
#View(T1)

names(T1)[names(T1)=="Sample name"]="Sample_name"
T1$Sample_name=gsub("J0$","J00",T1$Sample_name)
T1$Sample_name=gsub("M1","M01",T1$Sample_name)
T1$Visit=gsub("J00", "J0", T1$Visit)
T1$Visit=gsub("M01", "M1", T1$Visit)

#### Tube 2 : panel for deep phenotyping B cells ####
T2=data.frame()
T2 <- read.table(paste0(input_dir, "/SEP_C3_T2utf8.csv"), sep=";",dec=",",stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                 check.names = FALSE, fileEncoding = "UTF-8")

names(T2)[names(T2)=="Sample name"]="Sample_name"
T2$Sample_name=gsub("J0$","J00",T2$Sample_name)
T2$Sample_name=gsub("M1","M01",T2$Sample_name)
T2$Visit=gsub("J00", "J0", T2$Visit)
T2$Visit=gsub("M01", "M1", T2$Visit)

#### Tube 4 : panel for deep phenotyping T cells ####
T4=data.frame()
T4 <- read.table(paste0(input_dir, "/SEP_C3_T4utf8.csv"), sep=";",dec=",",stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                 check.names = FALSE, fileEncoding = "UTF-8")

names(T4)[1]="Sample_name"
T4$Sample_name=gsub("J0$","J00",T4$Sample_name)
T4$Sample_name=gsub("M1","M01",T4$Sample_name)
T4$Visit=gsub("J00", "J0", T4$Visit)
T4$Visit=gsub("M01", "M1", T4$Visit)

#### Tube 9 : panel for phenotyping T cells activation markers ####
T9=data.frame()
T9 <- read.table(paste0(input_dir, "/SEP_C3_T9utf8.csv"), sep=";",dec=",", stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                 check.names = FALSE, fileEncoding = "UTF-8")

names(T9)[names(T9)=="Sample name"]="Sample_name"
T9$Sample_name=gsub("J0$","J00",T9$Sample_name)
T9$Sample_name=gsub("M1","M01",T9$Sample_name)
T9$Visit=gsub("J00", "J0", T9$Visit)
T9$Visit=gsub("M01", "M1", T9$Visit)


#### Tregs ####
Treg=data.frame()
Treg <- read.table(paste0(input_dir, "/SEP_C3_Tregutf8.csv"), sep=";",dec=",", stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                   check.names = FALSE, fileEncoding = "UTF-8")

names(Treg)[names(Treg)=="id"]="Sample_name"
names(Treg)[names(Treg)=="VISIT"]="Visit"
Treg$Sample_name=gsub("J0$","J00",Treg$Sample_name)
Treg$Sample_name=gsub("M1","M01",Treg$Sample_name)
Treg$Visit=gsub("J00", "J0", Treg$Visit)
Treg$Visit=gsub("M01", "M1", Treg$Visit)


#### Treatments : list of treatment each patient receive####
treatlist <- read.table("Patient treatments.csv", sep=";",dec=",",stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                        check.names = FALSE, fileEncoding = "UTF-8")
names(treatlist)<- c("Sample_name","treatment")

#### total : merging all dataframes ####
total=merge(num, T9, by= c("Sample_name", "Visit"), all=TRUE)
total=merge(total, T4, by= c("Sample_name", "Visit"), all=TRUE)
total=merge(total, T2,  by= c("Sample_name", "Visit"), all=TRUE)
total=merge(total, T1,  by= c("Sample_name", "Visit"), all=TRUE)
total=merge(total, Treg,  by= c("Sample_name", "Visit"), all=TRUE)

#Removing duplicated columns, fluorescences and absolute numbers datas. I will only use percentages.
total=total[-grep(c("F"),names(total))]
total=total[-grep(c("N"),names(total))]
total=total[-grep(c(".x"), names(total))]
total=total[-grep(c(".y"), names(total))]

for (n in 2:length(names(total))){
  TEMP_char <- names(total)[n]
  TEMP_char <- gsub("/","_",TEMP_char, fixed = TRUE)
  names(total)[n]<-TEMP_char
}

#Mysterious column called Var.number with NA values to remove.
total <-total[,!grepl("^Var",colnames(total))]
#total %>% dplyr::select(!starts_with("Var."))

total$Visit[total$Visit == "j28"] <- "J28"
names(total) <- gsub("Visit","time_point", colnames(total))
names(total)
View(total)

#-------------------- Renaming time points according to treatments : -----------------
# aCD20 patients have       J0 -     - J56  -  J84 
# OtherDMTs patients have   J0 - J28 - J56  -
# reminder : sampling are done before vaccination shot
total <- merge(treatlist, total, by = "Sample_name")
aCD20 <- total$treatment == "Ocrelizumab" | total$treatment == "Rituximab" 

total$time_point[aCD20 == TRUE] <- gsub(pattern = "M1", replacement = "J84", 
                                       x = total$time_point[aCD20], 
                                       fixed = TRUE)
total$time_point[aCD20 == F] <-  gsub("M1","J56",total$time_point[aCD20 == F], fixed = TRUE)

#------------- Dataframes for patients under aCD20 (Bdep = "B depleted") and Oter DMTs treatments -----------
Bdep <- total %>% 
  subset((total$treatment == "Ocrelizumab" | total$treatment == "Rituximab"))

Other <- total %>%
  subset(total$treatment != "Ocrelizumab" & total$treatment != "Rituximab") 


# -------------- Useful vectors for analysis : ---------------
#measured parameters (for loops), patients names (for loops), time_points (for loops), names of all the variables (for plotting)

param=names(total)[4:length(names(total))]
param=param[param !="Cal.Factor"]

patients=as.character(unique(total$Sample_name))

time_points=as.character(unique(total$time_point))

treatments = unique(treatlist$treatment)

nomscols=names(total)[-c(1:3)]
nomscols=gsub("P ","% ",nomscols)
nomscols=gsub("P_","% ",nomscols)
nomscols=gsub("N ","Abs. Number of ",nomscols)

#--------------- cleaning environment -----------------
rm(T1,T2,T4, T9, Treg, num, treatlist, aCD20, TEMP_char, n, input_dir)
