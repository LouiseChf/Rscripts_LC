rm(list = ls(all = TRUE))
options(stringsAsFactors = FALSE)
options(error = traceback)
options(show.error.locations = TRUE)

library("plyr")       
library("dplyr")        
library("readr")
library("stringr")
library(Hmisc)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(tidyverse)
library(rstatix)

setwd("//DRUM/facs/Clinical_Trials/FC_Analysis/Studies/Covivac/Cytométrie")

input_dir="//DRUM/facs/Clinical_Trials/FC_Analysis/Studies/Covivac/Cytométrie/Export"
output_dir="//DRUM/facs/Clinical_Trials/FC_Analysis/Studies/Covivac/Cytométrie/Graphs"

if (dir.exists(output_dir)==FALSE){
  dir.create(output_dir)
}


###T1
T1=data.frame()
for (j in list.files("Export/T1")){
  TEMP=read.delim(paste0("Export/T1/", j), sep=",",dec=".",stringsAsFactors=FALSE, strip.white=TRUE)
  T1=rbind(T1, TEMP)
}

names(T1)[names(T1)=="ï..Sample.name"]="Sample_name"

err2=read.delim(paste0("Export/Erreurs_noms_T2.txt"), sep=";",dec=".",stringsAsFactors=FALSE, strip.white=TRUE, header=F)


for (i in T1$Sample_name){
  if (i %in% err2$V1){
    T1$Sample_name[T1$Sample_name==i]=err2$V2[err2$V1==i]
    # print(i)
  }
}


T1$Sample_name=gsub("J0$","J00",T1$Sample_name)
T1$Sample_name=gsub("M1","M01",T1$Sample_name)

T1$Visit=substr(T1$Sample_name, 13, 15)
T1$Visit=gsub("J00", "J0", T1$Visit)
T1$Visit=gsub("M01", "M1", T1$Visit)

T1$Sample_name=substr(T1$Sample_name, 1, 8)
# T1[10,1]="007-0004"
# T1=T1[-50,]

###T2
T2=data.frame()
for (j in list.files("Export/T2")){
  TEMP=read.delim(paste0("Export/T2/", j), sep=",",dec=".",stringsAsFactors=FALSE, strip.white=TRUE)
  T2=rbind(T2, TEMP)
}

names(T2)[names(T2)=="ï..Sample.name"]="Sample_name"

err2=read.delim(paste0("Export/Erreurs_noms_T2.txt"), sep=";",dec=".",stringsAsFactors=FALSE, strip.white=TRUE, header=F)


for (i in T2$Sample_name){
  if (i %in% err2$V1){
    T2$Sample_name[T2$Sample_name==i]=err2$V2[err2$V1==i]
    # print(i)
  }
}


T2$Sample_name=gsub("J0$","J00",T2$Sample_name)
T2$Sample_name=gsub("M1","M01",T2$Sample_name)

T2$Visit=substr(T2$Sample_name, 13, 15)
T2$Visit=gsub("J00", "J0", T2$Visit)
T2$Visit=gsub("M01", "M1", T2$Visit)

T2$Sample_name=substr(T2$Sample_name, 1, 8)
T2[10,1]="007-0004"
# T2=T2[-50,]

###T4
T4=data.frame()
for (j in list.files("Export/T4")){
  TEMP=read.delim(paste0("Export/T4/", j), sep=",",dec=".",stringsAsFactors=FALSE, strip.white=TRUE)
  names(TEMP)[names(TEMP)%in%c("ï..Name.of.the.data.file.containing.the.data.set","Name.of.the.data.file.containing.the.data.set")]="Sample_name"
  T4=rbind(T4, TEMP)
  # print(dim(TEMP))
}


err4=read.delim(paste0("Export/Erreurs_noms_T4.txt"), sep=";",dec=".",stringsAsFactors=FALSE, strip.white=TRUE, header=F)


for (i in T4$Sample_name){
  if (i %in% err4$V1){
    T4$Sample_name[T4$Sample_name==i]=err4$V2[err4$V1==i]
    # print(i)
  }
}

T4$Sample_name=str_split(T4$Sample_name, "-COVIVAC-", simplify = T)[,1]
err4$V1=str_split(err4$V1, "-COVIVAC-", simplify = T)[,1]

# T4$Sample_name=gsub("017-0001","007-0001",T4$Sample_name)
# T4$Sample_name=gsub("007-0027-ML$","007-0027-ML_M1",T4$Sample_name)
# T4$Sample_name=gsub("007-0006-AR$","007-0006-AR_J0",T4$Sample_name)
T4$Sample_name=gsub("J0$","J00",T4$Sample_name)
T4$Sample_name=gsub("M1","M01",T4$Sample_name)

T4$Visit=substr(T4$Sample_name, 13, 15)
T4$Visit=gsub("J00", "J0", T4$Visit)
T4$Visit=gsub("M01", "M1", T4$Visit)

T4$Sample_name=substr(T4$Sample_name, 1, 8)
# T4[10,1]="007-0004"
# T4=T4[-53,]

###T9
T9=data.frame()
for (j in list.files("Export/T9")){
  TEMP=read.delim(paste0("Export/T9/", j), sep=",",dec=".",stringsAsFactors=FALSE, strip.white=TRUE)
  T9=rbind(T9, TEMP)
}

names(T9)[names(T9)=="ï..Sample.name"]="Sample_name"

err9=read.delim(paste0("Export/Erreurs_noms_T9.txt"), sep=";",dec=".",stringsAsFactors=FALSE, strip.white=TRUE, header=F)

for (i in T9$Sample_name){
  if (i %in% err9$V1){
    T9$Sample_name[T9$Sample_name==i]=err9$V2[err9$V1==i]
    # print(i)
  }
}

# T9$Sample_name=gsub("017-0001","007-0001",T9$Sample_name)
# T9$Sample_name=gsub("007-0014-GG","007-0014-DJ_J0",T9$Sample_name)
T9$Sample_name=gsub("J0$","J00",T9$Sample_name)
T9$Sample_name=gsub("M1","M01",T9$Sample_name)

T9$Visit=substr(T9$Sample_name, 13, 15)
T9$Visit=gsub("J00", "J0", T9$Visit)
T9$Visit=gsub("M01", "M1", T9$Visit)

T9$Sample_name=substr(T9$Sample_name, 1, 8)

###total

total=merge(T2, T4, by= c("Sample_name", "Visit"), all=TRUE)
total=merge(total, T9,  by= c("Sample_name", "Visit"), all=TRUE)
total=merge(total, T1,  by= c("Sample_name", "Visit"), all=TRUE)
total=total[-grep(c("F\\."),names(total))]
total=total[-grep(c("N\\."),names(total))]
total=total[total$Visit != "J56",]

for (i in unique(total$Sample_name)){
  print(i)
  print(total$Visit[total$Sample_name==i])
  for (j in unique(total$Visit))
  if (!j %in% total$Visit[total$Sample_name==i]){
    # total[total$Sample_name==i & total$Visit=="M1", c(3:length(names(total)))]=NA
    total=rbind(total, c(i, j, rep(NA, length(names(total))-2)))
  }
}

# 
# for(i in names(total)[-c(1,2)]){
#   names(total)[names(total)==i]=nomscols[i]
# }
write.table(total, file="//DRUM/facs/Clinical_Trials/FC_Analysis/Studies/Covivac/Cytométrie/tableau_007.csv",
            row.names = F, sep=",")

# 
# total=total[total$Sample_name != "007-0016",]
# total=total[total$Sample_name != "007-0023",]
# total=total[total$Sample_name != "007-0017",]


## graphs

param=names(total)[3:length(names(total))]
param=param[param !="Cal.Factor"]
patients=as.character(unique(total$Sample_name))
visites=as.character(unique(total$Visit))

toplot=melt(total, id.vars=c("Sample_name","Visit"))
toplot$variable=as.character(toplot$variable)
toplot$Visit=as.factor(toplot$Visit)
toplot$value=as.numeric(gsub(",","",toplot$value))

my_comparisons <- list(c("J0","J28"), c("J28","M1"),c("J0","M1"))


nomscolsT2=names(read.delim(paste0("Export/T2/",list.files("Export/T2")[1]), sep=",",dec=".",stringsAsFactors=FALSE,
                            strip.white=TRUE, header = TRUE, check.names=F))[-c(1:2)]
nomscolsT4=names(read.delim(paste0("Export/T4/",list.files("Export/T4")[1]), sep=",",dec=".",stringsAsFactors=FALSE,
                            strip.white=TRUE, header = TRUE, check.names=F))[-1]
nomscolsT9=names(read.delim(paste0("Export/T9/",list.files("Export/T9")[1]), sep=",",dec=".",stringsAsFactors=FALSE,
                            strip.white=TRUE, header = TRUE, check.names=F))[-1]
nomscolsT1=names(read.delim(paste0("Export/T1/",list.files("Export/T1")[1]), sep=",",dec=".",stringsAsFactors=FALSE,
                            strip.white=TRUE, header = TRUE, check.names=F))[-c(1:2)]
nomscols=c(nomscolsT2,nomscolsT4,nomscolsT9,nomscolsT1)
nomscols=nomscols[grep("P ", nomscols)]
names(nomscols)=param
nomscols=gsub("P ","% ",nomscols)
nomscols=gsub("N ","Abs. Number of ",nomscols)


sapply(param, function(something){
 
  pj=position_jitter(width = 0.2, height=0, seed = 9)
  
  toplot2=toplot
  for (i in unique(toplot2$Sample_name)){
    if (nrow(na.omit(toplot2[toplot2$variable==something & toplot2$Sample_name==i,])) != 3){
      toplot2=toplot2[!toplot2$Sample_name==i,]
    }
  }
  
  ggplot(na.omit(toplot2[toplot2["variable"]==something,]), aes(x=Visit, y=value,  width=0.5))+
    geom_boxplot(notch=F, position=position_dodge(),aes(fill=Visit), alpha=0.5, outlier.colour = "white")+
    geom_line(position=pj,
              aes_string(group="Sample_name"), color="peachpuff3")+
    geom_jitter(position = pj,
                aes(fill=Visit), pch=21, size=3, alpha=1, show.legend = FALSE) +
    # facet_wrap(as.formula(paste("~","Visit")), ncol=4)+
    stat_compare_means(comparisons = my_comparisons, paired = TRUE, label= "p.signif", hide.ns = T)+
    # labs(caption = get_pwc_label(pwc))+
    theme_pubr()+
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(axis.text.x = element_text(face="bold", size=12, angle=90, vjust = 0.5), axis.text.y = element_text(face="bold", size=10))+
    ylab(something) +
    scale_fill_manual(values=c("darkorange2","tomato2","sienna","red3"))+
    scale_color_manual(values=c("darkorange2","tomato2","sienna","red3"))+
    theme(strip.text.x = element_text(size=12, color="black", face="bold"),
          strip.text.y = element_text(size=12, color="red",face="bold.italic"))+
    theme(strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid"))+
    theme(axis.title.x = element_blank())+
    theme(axis.title.y = element_blank())+
    # theme(axis.text.x = element_blank())+
    # scale_x_discrete(labels=c("ILT-101 N", "ILT-101 Y", "Placebo N", "Placebo Y"))+
    # scale_y_continuous(trans = "log2",breaks=c(seq(0,1,0.2),seq(1,10,1)))+
    theme(legend.position="none")+
    theme(legend.title=element_blank())+
    ggtitle(nomscols[something])
  
  ggsave(file = paste0(output_dir,"/",something,".png"), units = "cm", dpi = 150)
  
   
})


