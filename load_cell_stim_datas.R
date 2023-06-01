##Louise CHORFI, 05/05/2023
##Calling for my datas from excel and removing annoying characters.

## use following line in other files to run script :
## source("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell stim/Covivac stim/R/load_datas.R")

setwd("Y:/LCH/Documents/COVIVAC - cell stim")
library(tidyverse)

# --------- --Import datas from excel ------------
#data_N : table of data in absolute numbers
#data_P : table of data in percentages

data_P=data.frame()
data_N=data.frame()
TEMP=data.frame()

for (j in list.files("./excels/cell stim data")){
  TEMP=read.table(paste0("./excels/cell stim data/", j), sep=";",dec=",",stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                  check.names = FALSE)
  names(TEMP) <-  gsub(" ","_",names(TEMP), fixed = TRUE)
  #print(TEMP)
  TEMP_N <- TEMP[,1:23]
  TEMP_P <- TEMP[,24:46]
  TEMP_P <- TEMP_P %>% relocate(c(SAMPLE_ID, time_point, patient_ID, condition, treatment))
  TEMP_N <- TEMP_N %>% relocate(c(SAMPLE_ID, time_point, patient_ID, condition, treatment))
  data_P=rbind(data_P, TEMP_P)
  data_N=rbind(data_N, TEMP_N)
}
# View(data_P)
# View(data_N)


# --------- remove datas that are not exploitable (invalid pos or neg controls, no cells) --------
exploitable=read.table("./excels/exploitable.csv", sep=";",dec=",",stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                       fileEncoding="UTF-8",check.names = FALSE)
data_P <- merge(exploitable, data_P, by = "patient_ID", all.x = F)
data_P <- data_P %>% subset(subset = data_P$include == "oui")

data_N <- merge(exploitable, data_N, by = "patient_ID", all.x = F)
data_N <- data_N %>% subset(subset = data_N$include == "oui")

#------------- renaming properly ---------------
#### column names ####
for (n in 6:23){
  TEMP_char <- names(data_P)[n]
  TEMP_char <- gsub("[","",TEMP_char, fixed = TRUE)
  TEMP_char <- gsub("]","",TEMP_char, fixed = TRUE)
  TEMP_char <- gsub("(","",TEMP_char, fixed = TRUE)
  TEMP_char <- gsub(")","",TEMP_char, fixed = TRUE)
  TEMP_char <- gsub("+",".1",TEMP_char, fixed = TRUE)
  TEMP_char <- gsub("-",".0",TEMP_char, fixed = TRUE)
  TEMP_char <- gsub(" ","_",TEMP_char, fixed = TRUE)
  TEMP_char <- gsub("/",".",TEMP_char, fixed = TRUE)
  TEMP_char <- paste0("P_", TEMP_char)
  names(data_P)[n]<-TEMP_char
  TEMP_char <- gsub("P", "N", TEMP_char)
  names(data_N)[n]<-TEMP_char
}

#### conditions ####
data_P$condition <- gsub("_newAg", "", data_P$condition, ignore.case = TRUE)
data_N$condition <- gsub("_newAg", "", data_N$condition, ignore.case = TRUE)
data_P$condition <- gsub(" new Ag", "", data_P$condition, ignore.case = TRUE)
data_N$condition <- gsub(" new Ag", "", data_N$condition, ignore.case = TRUE)
data_P$condition <- gsub("pos", "+", data_P$condition, ignore.case = TRUE)
data_N$condition <- gsub("pos", "+", data_N$condition, ignore.case = TRUE)
data_P$condition <- gsub("neg", "-", data_P$condition, ignore.case = TRUE)
data_N$condition <- gsub("neg", "-", data_N$condition, ignore.case = TRUE)
data_N$condition <- gsub("/cmv", "", data_N$condition, ignore.case = TRUE)

data_P$condition <- gsub("/cmv", "", data_P$condition, ignore.case = TRUE)
data_N$condition <- gsub("-cmv", "", data_N$condition, ignore.case = TRUE)
data_P$condition <- gsub("-cmv", "", data_P$condition, ignore.case = TRUE)

data_P$condition <- gsub("cmv", "cmh", data_P$condition, ignore.case = TRUE)
data_P$condition <- gsub("CMH", "cmh", data_P$condition, ignore.case = F)
data_N$condition <- gsub("CMH", "cmh", data_N$condition, ignore.case = F)
data_N$condition <- gsub("cmv", "cmh", data_N$condition, ignore.case = TRUE)

data_N$condition <- gsub("WT", "wt", data_N$condition, ignore.case = TRUE)
data_P$condition <- gsub("WT", "wt", data_P$condition, ignore.case = TRUE)

#### Renaming time points according to treatments : ####
# aCD20 patients have       J0 -     - J56  -  J84 
# OtherDMTs patients have   J0 - J28 - J56  -
# reminder : sampling are done before vaccination shot

aCD20 <- data_P$treatment == "Ocrelizumab" | data_P$treatment == "Rituximab" 

data_P$time_point[aCD20 == TRUE] <- gsub(pattern = "M1", replacement = "J84", 
                        x = data_P$time_point[aCD20], fixed = TRUE)
data_P$time_point[aCD20 == F] <-  gsub("M1","J56",data_N$time_point[aCD20 == F], fixed = TRUE)


data_N$time_point[aCD20 == TRUE] <- gsub(pattern = "M1", replacement = "J84", 
                                       x = data_P$time_point[aCD20], fixed = TRUE)
data_N$time_point[aCD20 == F] <-  gsub("M1","J56",data_N$time_point[aCD20 == F], fixed = TRUE)
 