##Louise CHORFI, 05/05/2023
##Calling for the humoral datas from excel and cleaning annoying characters.

## use following code in other files to new script :
## source("//contrebasse.transimmunom.org/Users/lch/Documents/COVIVAC - cell stim/Covivac stim/R/load_humoral_datas.R")

setwd("Y:/LCH/Documents/COVIVAC - cell stim")
library(tidyverse)


# --------- --Import datas from excel ------------
humoral=data.frame()
TEMP=data.frame()

for (j in list.files("./excels/humoral data")){
  TEMP=read.table(paste0("./excels/humoral data/", j), sep=";",dec=",",stringsAsFactors=FALSE, 
                  strip.white=TRUE, header=TRUE,
                  check.names = FALSE)
  names(TEMP) <-  gsub(" ","_",names(TEMP), fixed = TRUE)
  #print(TEMP)
  humoral=rbind(humoral, TEMP)
  
}
View(humoral)

treatments=read.table("./excels/Patient treatments.csv", sep=";",dec=",",stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                       fileEncoding="UTF-8",check.names = FALSE)

humoral <- merge(treatments, humoral, by = "patient_ID")

# --------- remove datas that are not exploitable (invalid pos or neg controls, no cells) --------
exploitable=read.table("./excels/exploitable.csv", sep=";",dec=",",stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE,
                       fileEncoding="UTF-8",check.names = FALSE)
humoral <- merge(exploitable, humoral, by = "patient_ID", all.x = F)
humoral <- humoral %>% subset(subset = humoral$include == "oui")

#------------- renaming properly ---------------
for (n in names(humoral)){
  TEMP_char <- n
  TEMP_char <- gsub("(","",TEMP_char, fixed = TRUE)
  TEMP_char <- gsub(")","",TEMP_char, fixed = TRUE)
  TEMP_char <- gsub("/",".",TEMP_char, fixed = TRUE)
  names(humoral)[names(humoral) == n] <-TEMP_char
  #T1_P[n] <- gsub("N/A", NA, T1_P[n])
}
