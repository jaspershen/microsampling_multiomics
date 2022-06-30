masstools::setwd_project()


library(tidyverse)
library(data.table)
library(cowplot)
library(plotly)
library(ComplexHeatmap)
library(ggforce)
library(lubridate)
library(scales)
library(MetaCycle)
library(pheatmap)
library(ggpubr)
library(dtwclust)
library(jsonify)
library(jsonlite)
library(readxl)

rm(list = ls())

# WEARABLES
steps = read.csv("~/Box/Microsampling/Analysis/RK/05_2019/iPOP1_steps_05_2019", sep="\t", header=TRUE)
hr = read.csv("~/Box/Microsampling/Analysis/RK/05_2019/iPOP1_HR_05_2019", sep="\t", header=TRUE)
glu = read.csv("~/Box/Microsampling/Analysis/RK/CGM_data.txt", sep="\t", header=TRUE) 
samples = read.csv("~/Box/Microsampling/Analysis/RK/sample_registration.csv", header=T, colClasses = c("integer","factor","character","character","POSIXct","character","character","character"))

# sleep processing
a=stream_in(file("~/Box/Microsampling/Analysis/RK/05_2019/Sleep.json"))
sleep = rbindlist((a[[1]][[1]])[[7]][[1]])
for (i in 2:length(a[[1]])) {
  sleep = rbind(sleep, rbindlist( (a[[1]][[i]])[[7]][[1]], fill=TRUE))
}

s2 = a[[1]][[1]][1:4]
for (i in 2:length(a[[1]])) {
  s2 = rbind(s2, a[[1]][[i]][1:4]) 
}

sleep$type="Sleep"
steps$type="Steps"
hr$type="HR"
glu$type="CGM"

sleep$MolClass="Sleep"
steps$MolClass="Steps"
hr$MolClass="HR"
glu$MolClass="CGM"

samples$type = ""
samples[!is.na(samples$comment) & samples$comment != "","type"] = "Comment"
samples[!is.na(samples$food) & samples$food != "","type"] = "Food"
samples[!is.na(samples$activity) & samples$activity != "","type"] = "Activity"
samples[!is.na(samples$sample.ID) & samples$sample.ID != "","type"] = "Sample"

names(steps)[names(steps) == "steps"] = "value"
names(hr)[names(hr) == "heart_rate"] = "value"
names(glu)[names(glu) == "GlucoseValue"] = "value"
names(glu)[names(glu) == "GlucoseDisplayTime"] = "DT"
sleep = dplyr::rename(sleep, DT = dateTime, value = level)
steps = dplyr::rename(steps, DT = date_time)
hr = dplyr::rename(hr, DT = date_time)
samples = dplyr::rename(samples, DT = date_time)


steps$DT = as.POSIXct(steps$DT, tz = "America/Los_Angeles")
hr$DT = as.POSIXct(hr$DT, tz = "America/Los_Angeles")
glu$DT = as.POSIXct(glu$DT, tz = "America/Los_Angeles")
samples$DT = as.POSIXct(samples$DT, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
sleep$DT = as_datetime(sleep$DT, tz="America/Los_Angeles")

d1 = as.POSIXct("2019-04-29 00:00:00", tz = "America/Los_Angeles")
d2 = as.POSIXct("2019-05-07 23:59:59", tz = "America/Los_Angeles")

steps = steps[steps$DT >= d1 & steps$DT <= d2, ]
hr = hr[hr$DT >= d1 & hr$DT <= d2, ]
glu = glu[glu$DT >= d1 & glu$DT <= d2, ]
sleep = sleep[sleep$DT >= d1 & sleep$DT <= d2, ]

steps$time_window = cut(steps$DT, breaks="15 min", labels=FALSE)
hr$time_window = cut(hr$DT, breaks="15 min", labels=FALSE)
glu$time_window = cut(glu$DT, breaks="15 min", labels=FALSE)
sleep$time_window = cut(sleep$DT, breaks="15 min", labels=FALSE)


sleep = sleep %>% mutate(DT = DT, Day = date(DT), Time = hms::as.hms(DT)) %>% dplyr::rename(Intensity = value, MolName = type)
steps = steps %>% mutate(DT = DT, Day = date(DT), Time = hms::as.hms(DT)) %>% dplyr::rename(Intensity = value, MolName = type)
hr = hr %>% mutate(DT = DT, Day = date(DT), Time = hms::as.hms(DT)) %>% dplyr::rename(Intensity = value, MolName = type)
glu = glu %>% mutate(DT = DT, Day = date(DT), Time = hms::as.hms(DT)) %>% dplyr::rename(Intensity = value, MolName = type)

wearables = bind_rows(steps,hr,glu)
wearables = select(wearables,MolName,MolClass,DT,Intensity)
wearables = wearables %>% mutate(Day = date(DT), Time = hms::as.hms(DT), Hr = hour(DT))




# FOOD PROFILE
TimeAnno = read.csv("~/Box/Microsampling/Analysis/RK/TimeAnno.csv")
FoodList = read.csv("~/Box/Microsampling/Analysis/RK/Food_list.csv")
FoodTime = TimeAnno %>% separate_rows(food_entry, food_multiplier) 
FoodTime$food_entry = as.integer(FoodTime$food_entry)
#FoodTime = FoodTime %>% select(date_time,activity,food,food_entry,food_multiplier)
FoodTime = left_join(FoodTime, FoodList, by = c("food_entry" = "EntryNum")) %>%
  subset(!is.na(food_entry))

FoodTimeL = FoodTime %>% 
  gather("Parameter","Value","Calories":"Alcohol_g") %>% 
  select(-sample.ID,-finger,-comment,-activity, -time,-date, -timef, -datef)

FoodTimeL$DT = as.POSIXct(FoodTimeL$date_time, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
FoodTimeL = FoodTimeL %>% mutate(Day = date(DT), Time = hms::as.hms(DT), Hr = hour(DT))
FoodTimeL$food_multiplier = as.integer(FoodTimeL$food_multiplier)

# temp = 
# FoodTimeL %>% 
#   plyr::dlply(.variables = .(DT, Day, Time, food_multiplier, Parameter))

FoodTimeL_sum = FoodTimeL %>% 
  dplyr::group_by(DT, Day, Time, food_multiplier, Parameter) %>% 
  dplyr::summarise(Value = sum(Value))

FoodTimeL_sum = FoodTimeL_sum %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(MolName = Parameter, Intensity = Value) 

FoodTimeL_sum %>% 
  dplyr::filter(is.na(Intensity))


# FoodTimeL_sum %>% 
#   dplyr::group_by(DT, MolName) %>% 
#   dplyr::summarise(n = n()) %>% 
#   dplyr::filter(n > 1)


###this is a method to avoid the duplicated variables
FoodTimeL_sum$MolName = 
  paste(FoodTimeL_sum$food_multiplier, FoodTimeL_sum$MolName, sep = "_")

FoodTimeL_sum$MolClass = "Food"

class(FoodTimeL_sum)

library(plyr)


# 
# all_omes = all_omes %>% mutate(Day = date(DT), Time = hms::as.hms(DT), Hr = hour(DT))
# 
# 
# ggplot(subset(FoodTimeL_sum,Parameter %in% c("Carbs_g","Fat_g","Protein_g","Alcohol_g")), aes(x = Time, y = food_multiplier*Value, color=Parameter, fill = Parameter, group=Parameter)) + 
# #  geom_area() +
#   geom_point() + geom_line() + 
#   geom_line(data = subset(all_omes,MolName=="INSULIN"), mapping = aes(x=Time, y=Intensity)) #+ 
#     #facet_wrap(~Day) + theme(axis.text.x = element_text(angle=90))
# 
# 
# 
# ggplot() + 
#   geom_line(data = subset(FoodTimeL_sum,Parameter %in% c("Carbs_g","Protein_g")), mapping = aes(x = DT, y = zscore(food_multiplier*Value), color=Parameter, fill = Parameter, group=Parameter)) +
#  geom_line(data = subset(all_omes,MolName=="CHEX4"), mapping = aes(x=DT, y=zscore(Intensity), color=MolName)) #+ 
#     #facet_wrap(~Day) + theme(axis.text.x = element_text(angle=90))
# 
# # # CHEX4 correlations 
# # ss = subset(all_omes, MolClass=="MetabolicPanel")
# # CHEX4 = ss[ss$MolName=="CHEX4_MP",c("DT","Intensity")] %>% dplyr::rename(CHEX4 = Intensity)
# # sdf = left_join(ss,CHEX4)
# # ggplot(sdf, aes(CHEX4, Intensity)) + facet_wrap(~MolName,scales="free") + geom_point() + geom_smooth()
# # ss = subset(all_omes, MolClass=="Cytokine_41Plex")
# # CHEX4 = ss[ss$MolName=="CHEX4",c("DT","Intensity")] %>% dplyr::rename(CHEX4 = Intensity)
# # sdf = left_join(ss,CHEX4)
# # ggplot(sdf, aes(CHEX4, Intensity)) + facet_wrap(~MolName,scales="free") + geom_point() + geom_smooth()
#  
# ggplot(subset(all_omes,grepl("CHEX4", all_omes$MolName)), aes(DT, Intensity, group=MolName, color=MolName)) + geom_line() + geom_smooth()
# 
# ggplot(subset(FoodTimeL,Parameter %in% c("Carbs_g","Fat_g","Protein_g","Alcohol_g")), aes(x = DT,y = Value, fill = Parameter, group=Parameter)) + geom_bar(stat="identity",position="stack",width = 1000) #geom_point() + geom_line()
# 
# 
# ggplot()
# 
# subset(FoodTimeL,Parameter=="Carbs_g" && Day=="2019-05-01")



# INTEGRATE


load("~/Box/Microsampling/Analysis/Dan_2019-08-17/MOmics_01.RData")

#TimeAnno$datef = format(TimeAnno$date_time,"%a %b %d")
#TimeAnno$datef = factor(TimeAnno$datef,levels = unique(TimeAnno$datef))
#TimeAnno$timef = as.POSIXct(format(TimeAnno$date_time, "2000-01-01 %H:%M:%S"))
library(tidyverse)
TimeAnno = read.csv("~/Box/Microsampling/Analysis/RK/TimeAnno.csv")
TimeAnno = TimeAnno %>% dplyr::select(-time, -date) %>% dplyr::rename(DT = date_time)
events = gather(TimeAnno,"Type","Details",5:6)
events = events[!is.na(events$Details),]

start = as_datetime("2019-04-29 03:00:00", tz="America/Los_Angeles")
events$DT = as.POSIXct(events$DT, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
events = events %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))

# lipids
lipids = DF247_lipids2[,c(1,2,4,5,10,13)]
lipids$SampleID = as.numeric(gsub("24_7_prepsample","",lipids$SampleID))
lipids = lipids[!is.na(lipids$SampleID),]
lipids$MolClass = "Lipid"
lipids = 
  dplyr::rename(lipids, MolName = LipidSpecies, MolSubclass = LipidClass, Intensity = log_sample_nmol_per_g_concentration, DT = CollectionTime_format)
lipids = dplyr::select(lipids, -collector)
lipids$SampleID = as.character(lipids$SampleID)

lipids %>% 
  dplyr::filter(MolName == "CE.12.0.") %>% 
  dplyr::pull(Intensity) %>% 
  plot()

# proteins - plasma

proteins <- data.table(fread("~/Box/Microsampling/proteins/merged_library/allproteins_overtime_annotated_Mike247plasma.csv"), sep='\t')
uniprot = read.csv("~/Box/Microsampling/Analysis/RK/uniprot_match.csv")
uniprot$Entry_name = gsub("_HUMAN","",uniprot$Entry_name)
proteinsL = gather(proteins[,c(1:313,321,325)], Sample, Intensity, 1:313)
proteinsL = left_join(proteinsL,select(uniprot,Entry,Entry_name), by=c("Sample"="Entry"))
proteinsL = subset(proteinsL,!is.na(Entry_name))
proteinsL = unite(proteinsL,"MolName",Sample,Entry_name)
proteinsL = proteinsL %>% dplyr::rename(DT = date_time, SampleID = SampleIndex)
proteinsL$MolClass = "Protein"
proteinsL$MolSubclass = ""
proteinsL$SampleID = as.character(proteinsL$SampleID)
proteinsL$DT = as.POSIXct(proteinsL$DT, tz = "America/Los_Angeles", format = "%Y-%m-%d %H:%M:%S")
#proteinsL = proteinsL %>% group_by(SampleID) %>% subset(!duplicated(MolName)) %>% ungroup() #mutate(MolName = make.unique(as.character(MolName)))

head(proteinsL)
length(unique(proteinsL$SampleID))
length(unique(proteinsL$MolName))

proteinsL %>% 
  dplyr::filter(MolName == "B9A064_IGLL5") %>% 
  dplyr::pull(Intensity) %>% 
  plot()

# metabolites
#metab_ = read.csv("~/Box/Microsampling/metabolites/MS_247/Data/metabolomics/Microsampling_247_val_mets.csv")
metab_all = read.csv("~/Box/Microsampling/metabolites/MS_247/Data/metabolomics/Microsampling_247_combined_all.csv")
length(unique(paste(metab_all$Mode, metab_all$Compound)))

# ### metID
# pR = read.csv("~/Box/Microsampling/metabolites/MS_247/metID_Microsampling_247/RPLC/pos/identification.table.new_pRPLC.csv")
# nR = read.csv("~/Box/Microsampling/metabolites/MS_247/metID_Microsampling_247/RPLC/neg/identification.table.new_nRPLC.csv")
# pH = read.csv("~/Box/Microsampling/metabolites/MS_247/metID_Microsampling_247/HILIC/pos/identification.table.new_pHILIC.csv")
# nH = read.csv("~/Box/Microsampling/metabolites/MS_247/metID_Microsampling_247/HILIC/neg/identification.table.new_nHILIC.csv")
# pR$mode="pR"
# nR$mode="nR"
# pH$mode="pH"
# nH$mode="nH"
# all_metID = rbind(pR,nR,pH,nH)
# 
# #all_metID = all_metID[!duplicated(all_metID$name),]
# mets_join = left_join(metab_all,all_metID,by=c("Compounds_ID"="name"))
# write.table(mets_join,file="24-7_joined_110820b.csv",sep=",",row.names=F)
# 
# 
# mets_join %>% ungroup %>% group_by(mode) %>% summarise(count = n())
# 
metab_all = subset(metab_all, !Metabolite_val %in% c("0","","NA")) %>% select(Metabolite_val:CAS_val,Super.pathway:Sub.pathway,X1:X106)
# #metab_all = select(metab_all, Compounds_ID, Metabolite_val:CAS_val,Super.pathway:Sub.pathway,X1:X106)
# 
# #left_join(metab_all, select(all_metID,name,Compound.name:Database), by=c("name"="Compounds_ID"))
# 
# j=inner_join(metab_all, select(all_metID,name,Compound.name:Database), by=c("Compounds_ID"="name"))
# 
# j=select(j,Compound.name:Database,everything())
# j=j[!duplicated(j$Compounds_ID),]
# 
# write.table(j,file="24-7_joined.csv",sep=",",row.names=F)

inj_order = read_excel("~/Box/Microsampling/metabolites/MS_247/Data/metabolomics/M-Sampling-Brittany.xlsx","247 ome")
inj_order = dplyr::rename(inj_order, SampleID = "PrepIndex (SampleName)")
inj_order$AcqOrder = paste0("X", inj_order$AcqOrder)
#inj_order$AcqOrder = gsub("\\..","",inj_order$AcqOrder)

metab = gather(metab_all,AcqOrder,Intensity,-(Metabolite_val:Sub.pathway))
metab = dplyr::rename(metab,MolName = Metabolite_val, HMDB_ID = HMDB_val)
metab = left_join(metab,select(inj_order,AcqOrder,SampleID)) %>% select(-AcqOrder)

# median norm
metab = metab %>% group_by(SampleID) %>% mutate(Intensity = log10(Intensity/median(Intensity)))
#metab$SampleID = as.numeric(gsub("X24_7_prepsample","",metab$SampleID))
metab = metab[!is.na(metab$SampleID),]
metab$SampleID = as.character(metab$SampleID)

anno = read.csv("~/Box/Microsampling/Analysis/annotations/24_7_Stability_MasterIDReferenz_DH4_KE2_rk.csv")
metab = left_join(metab,anno[,c("PrepID","CollectionTime", "collector")], by=c("SampleID" = "PrepID"))
metab$DT = as.POSIXct(metab$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
metab = select(metab, -collector, -CollectionTime)
metab$MolClass = "Metabolite"
metab$MolSubclass = NA
metab = subset(metab, !is.na(DT))
#metab$MolName = make.unique(as.character(metab$MolName))
metab = metab %>% group_by(SampleID) %>% mutate(MolName = make.unique(as.character(MolName)))

# weird sample #50
metab = metab[metab$SampleID != 50,]
metab=subset(metab, !(DT=="2019-05-05 09:11:00" & metab$MolName=="Caffeine"))


length(unique(metab$SampleID))
metab$MolName


# metabolic protein panel
MetaData=read.csv("~/Box/Microsampling/Analysis/annotations/Copy of 24_7_Stability_MasterIDReferenz_DH4_KE2.csv") 

MetaData = MetaData[MetaData$Study == "F0" & MetaData$Original_SampleID %in% 1:97,]
MetaData$Name = paste("DH-",MetaData$Original_SampleID,sep="")

MP = read.csv("~/Box/Microsampling/cytokines/24-7 omics-MetabolicHormone-Ryan Metabolic Plate-1.csv")
MP = MP[2:97,]

MP = left_join(MP,MetaData[,c("Name","CollectionTime","PrepIndex")])
MP$DT = as.POSIXct(MP$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
MPL = gather(MP,"Cytokine","MFI",5:21)
MPL$Intensity = log10(as.numeric(MPL$MFI))
MPL = dplyr::rename(MPL, MolName = Cytokine, SampleID=PrepIndex)
MPL$MolClass = "MetabolicPanel"
MPL$MolSubclass = NA
MPL = select(MPL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)
MPL[MPL$MolName=="TNFA",]$MolName = "TNFA_MP"
MPL[MPL$MolName=="MCP1",]$MolName = "MCP1_MP"
MPL[MPL$MolName=="CHEX1",]$MolName = "CHEX1_MP"
MPL[MPL$MolName=="CHEX2",]$MolName = "CHEX2_MP"
MPL[MPL$MolName=="CHEX3",]$MolName = "CHEX3_MP"
MPL[MPL$MolName=="CHEX4",]$MolName = "CHEX4_MP"
MPL[MPL$MolName=="IL6",]$MolName = "IL6_MP"

length(unique(MPL$SampleID))
length(unique(MPL$MolName))

MPL %>% 
  dplyr::filter(MolName == "C.Peptide") %>% 
  dplyr::pull(Intensity) %>% 
  plot()


# cortisol 
MetaData=read.csv("~/Box/Microsampling/Analysis/annotations/Copy of 24_7_Stability_MasterIDReferenz_DH4_KE2.csv") 

MetaData = MetaData[MetaData$Study == "F0" & MetaData$Original_SampleID %in% 1:97,]
MetaData$Name = paste("DH-",MetaData$Original_SampleID,sep="")

cort = read.csv("~/Box/Microsampling/cytokines/RYAN CORTISOL 8-28-19.csv")
cort = cort[10:102,]

cort = left_join(cort,MetaData[,c("Name","CollectionTime","PrepIndex")])
cort$DT = as.POSIXct(cort$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
cort = dplyr::rename(cort, SampleID=PrepIndex)

cortL = gather(cort,"Cytokine","MFI",5)
cortL$Intensity = log10(as.numeric(cortL$MFI))
cortL$MolClass = "CortisolEnzymatic"
cortL$MolSubclass = NA
cortL$MolName = "Cortisol"
cortL = select(cortL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)

length(unique(cortL$SampleID))
length(unique(cortL$MolName))

cortL %>% 
  dplyr::filter(MolName == "Cortisol") %>% 
  dplyr::pull(Intensity) %>% 
  plot()


#41-plex
cyto = read.csv("~/Box/Microsampling/cytokines/24-7 omics-HumanLuminexMAG42plex-Mitra H41 Plex-1.csv")
cyto = cyto[10:105,]

cyto = left_join(cyto,MetaData[,c("Name","CollectionTime","PrepIndex")])
cyto$DT = as.POSIXct(cyto$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")

cytoL = gather(cyto,"Cytokine","MFI",5:49)
cytoL$Intensity = log10(as.numeric(cytoL$MFI))
cytoL$MolClass = "Cytokine_41Plex"
cytoL$MolSubclass = NA
cytoL = dplyr::rename(cytoL, SampleID = PrepIndex, MolName = Cytokine)
cytoL = select(cytoL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)


length(unique(cytoL$SampleID))
length(unique(cytoL$MolName))

cytoL %>% 
  dplyr::filter(MolName == "EGF") %>% 
  dplyr::pull(Intensity) %>% 
  plot()


all_omes = bind_rows(metab, lipids, proteinsL, MPL, cortL, cytoL)

# add MolNum
temp = tibble(MolName = unique(all_omes$MolName)) #, MolIndex = 1:n(MolName))
temp$MolNum = 1:nrow(temp)
all_omes = left_join(all_omes, temp)

all_omes_wear = bind_rows(all_omes, wearables)
all_omes_wear_food = bind_rows(all_omes, wearables, FoodTimeL_sum)

all_omes = all_omes %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))
all_omes = all_omes[all_omes$DT < as_date("2019-05-07"),]

all_omes_wear = all_omes_wear %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))
all_omes_wear_food = all_omes_wear_food %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))


grid = seq(as.POSIXct(floor_date(min(all_omes$DT),"hour")),as.POSIXct(max(all_omes$DT)),by="1 hour")

save(all_omes_wear_food, file = "data/24_7_study/all_omes_wear_food")

# all_omes_hr = 
#   all_omes %>% 
#   dplyr::select(-KEGG_val, -HMDB_ID, -CAS_val, -Super.pathway, -Sub.pathway) %>% 
#   dplyr::group_by(MolName,MolClass,MolSubclass,MolNum) %>% 
#   dplyr::mutate(MI = mean(Intensity), SDI = sd(Intensity), 
#                 Intensity_z = (Intensity - mean(Intensity,na.rm=T))/sd(Intensity,na.rm=T)) %>% 
#   tidyr::nest() %>%
#   dplyr::mutate(hr = purrr::map(data, function(x) tibble::as_tibble(approx(x$DT, x$Intensity_z, xout=grid, rule=1)))) %>%
#   unnest(hr) %>% 
#   dplyr::rename(DT = x, Intensity = y)
# 
# all_omes_hr = all_omes_hr[!is.na(all_omes_hr$Intensity),]
# all_omes_hr = all_omes_hr %>% 
#   group_by(MolName) %>% 
#   mutate(Hr = as.factor(as.numeric(difftime(DT, first(DT)), units="hours")), Hr_day = hour(DT), AM = am(DT))
# all_omes_hr = all_omes_hr[all_omes_hr$DT < as_date("2019-05-07"),]
# all_omes_hr = all_omes_hr %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))
load("data/24_7_study/all_omes_hr")

# all_omes_wear = all_omes_wear %>% 
#   select(-KEGG_val, -HMDB_ID, -CAS_val, -Super.pathway, -Sub.pathway) %>% 
#   group_by(MolName,MolClass,MolSubclass,MolNum) %>% mutate(MI = mean(Intensity), SDI = sd(Intensity), Intensity_z = (Intensity - mean(Intensity,na.rm=T))/sd(Intensity,na.rm=T)) %>% nest() 
# all_omes_wear = all_omes_wear %>%
#   mutate(hr = map(data, function(x) as_tibble(approx(x$DT, x$Intensity_z, xout=grid, rule=1)))) %>%
#   unnest(hr) %>% dplyr::rename(DT = x, Intensity = y)
# all_omes_wear = all_omes_wear[!is.na(all_omes_wear$Intensity),]
# all_omes_wear = all_omes_wear %>% group_by(MolName) %>% mutate(Hr = as.factor(as.numeric(difftime(DT, first(DT)), units="hours")), Hr_day = hour(DT), AM = am(DT))
# all_omes_wear = all_omes_wear[all_omes_wear$DT > "2019-04-29 03:00:00" & all_omes_wear$DT < as_date("2019-05-07"),]
# all_omes_wear = all_omes_wear %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))
# 
# # 
# # all_omes_wear_food = all_omes_wear_food %>% group_by(MolName,KEGG_val, HMDB_ID, CAS_val, Super.pathway, Sub.pathway, MolClass,MolSubclass,MolNum) %>% mutate(MI = mean(Intensity), SDI = sd(Intensity), Intensity_z = (Intensity - mean(Intensity,na.rm=T))/sd(Intensity,na.rm=T)) %>% nest() %>%
# #   mutate(hr = map(data, function(x) as_tibble(approx(x$DT, x$Intensity_z, xout=grid, rule=1)))) %>%
# #   unnest(hr) %>% dplyr::rename(DT = x, Intensity = y)
# # all_omes_wear_food = all_omes_wear_food[!is.na(all_omes_wear_food$Intensity),]
# # all_omes_wear_food = all_omes_wear_food %>% group_by(MolName) %>% mutate(Hr = as.factor(as.numeric(difftime(DT, first(DT)), units="hours")), Hr_day = hour(DT), AM = am(DT))
# # all_omes_wear_food = all_omes_wear_food[all_omes_wear_food$DT > "2019-04-29 03:00:00" & all_omes_wear_food$DT < as_date("2019-05-07"),]
# # all_omes_wear_food = all_omes_wear_food %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))
# # 
# # 
# # grid2 = seq(as.POSIXct(floor_date(min(all_omes$DT),"hour")),as.POSIXct(max(all_omes$DT)),by="2 hour")
# # #all_omes %>% group_by(MolName,MolClass,MolSubclass) %>% mutate(Intensity_z = (Intensity - mean(Intensity)/sd(Intensity)))
# # 
# # all_omes_2hr = all_omes %>% group_by(MolName,MolClass,MolSubclass,MolNum) %>% mutate(Intensity_z = (Intensity - mean(Intensity))/sd(Intensity)) %>% nest() %>%
# #   mutate(hr = map(data, function(x) as_tibble(approx(x$DT, x$Intensity_z, xout=grid2, rule=1)))) %>%
# #   unnest(hr) %>% dplyr::rename(DT = x, Intensity = y)
# # all_omes_2hr = all_omes_2hr[!is.na(all_omes_2hr$Intensity),]
# # all_omes_2hr = all_omes_2hr %>% group_by(MolName) %>% mutate(Hr = as.factor(as.numeric(difftime(DT, first(DT)), units="hours")), Hr_day = hour(DT), AM = am(DT))
# # all_omes_2hr = all_omes_2hr[all_omes_2hr$DT < as_date("2019-05-07"),]
# # all_omes = all_omes[all_omes$DT < as_date("2019-05-07"),]
# # 
# # all_omes_2hr = all_omes_2hr %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))
# 
# # metab = metab[-c(which(metab$DT=="2019-05-05 09:11:00" & metab$MolName=="Caffeine")),]
# # m2=subset(metab, !(DT=="2019-05-05 09:11:00" & metab$MolName=="Caffeine"))
# # p=ggplot(subset(metab,MolName=="Caffeine"),aes(DT,Intensity)) + geom_line()
# # ggplotly(p)
# 
# #ggplot(subset(all_omes_wear, MolName %in% c("HR")), aes(DT,Intensity,group=MolClass, color=MolName)) + geom_line()

load("data/24_7_study/all_omes_wear")


# BIG ALL-DATA FIG
all_omes_hr_d = all_omes_hr %>% 
  mutate(Day = day(DT))

anno_t = as.data.frame(select(subset(ungroup(all_omes_hr_d), MolNum == 1), DT, Day, Hr, Hr_day, AM))
# anno_shading = anno_t %>% group_by(Day) %>% subset(AM) %>% summarize(shading_start = first(DT), shading_end = last(DT))
anno_shading = anno_t %>%
  group_by(Day) %>%
  summarise(shading_start = DT[(Hr_day ==
                                  6)], shading_end = DT[(Hr_day == 18)])
events = subset(events, DT < "2019-05-07")


# all_omes = all_omes %>% mutate(Day = date(DT), Time = hms::as.hms(DT), Hr = hour(DT))
# sleep$date_time = as.POSIXct(sleep$date_time,tz="America/Los_Angeles")
# ggplot() + geom_rect(data=sleep, aes(xmin=Time, xmax=Time+seconds, ymin=0, ymax=1, fill=value)) + 
#   facet_wrap(~Day) + 
#   geom_line(data=subset(all_omes,MolName=="Cortisol"), aes(Time,Intensity))

mytheme = theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.line.x = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                panel.background = element_rect("grey95"),
                strip.background = element_blank(),
                strip.text.x = element_blank(),
                strip.text.y = element_blank()
) 


plot_labels = data.frame(MolName = "Sleep", Label = "Sleep")
p1 = ggplot() +
  geom_rect(
    data = anno_shading,
    mapping = aes(
      xmin = shading_start,
      xmax = shading_end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "lightyellow"
  ) +
  #geom_vline(data = subset(events, grepl("milk",events$Details)), mapping=aes(xintercept=DT), color="gray",size=1) +
  geom_rect(data = sleep,
            aes(
              xmin = DT,
              xmax = DT + seconds,
              ymin = 0,
              ymax = 1,
              fill = Intensity
            )) + labs(fill = "") +
  scale_x_datetime(
    breaks = date_breaks("4 hour"),
    date_labels = "%a %b-%d %H:%M",
    limits = c(min(all_omes$DT), max(all_omes$DT)),
    timezone = "America/Los_Angeles"
  ) +
  scale_fill_brewer(palette = "Set1") +
  mytheme + theme(
    legend.position = "top",
    legend.box = "horizontal",
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  ) + guides(fill = guide_legend(nrow = 1)) +
  facet_wrap( ~ MolClass, ncol = 1, strip.position = "right") +
  geom_label(
    data = plot_labels,
    mapping = aes(
      x = as.POSIXct("2019-04-29 04:00:00", tz = "America/Los_Angeles"),
      y = Inf,
      label = Label
    ),
    color = "black",
    size = 5,
    hjust = 0,
    vjust = 1,
    nudge_y = -10,
    nudge_x = -60
  ) +
  mytheme

p1
# ggplot(data = steps, aes(x=DT,y=value,color=type,group=1)) + geom_line() +
#   geom_vline(data = subset(events, grepl("milk",events$Details)), mapping=aes(xintercept=DT), color="gray",size=1) +
# scale_x_datetime(breaks=date_breaks("2 hour"), date_labels = "%H:%M") + theme(axis.text.x = element_blank()) +
# geom_line(data = hr, aes(x=DT,y=value, color=type, group=1)) + geom_line(data = glu, aes(x=DT, y=value, color=type, group=1)) + theme_void() + labs(color="",size=5) +
# scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %b-%d %H:%M", limits=c(min(all_omes$DT),max(all_omes$DT)), timezone = "America/Los_Angeles") 

plot_labels = data.frame(MolName = c("Steps","HR","CGM"), Label = c("Steps","HR","CGM"))
p2 = ggplot() + 
  geom_rect(data = anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="lightyellow") +
  geom_line(data = wearables, aes(DT, Intensity)) + 
  facet_wrap(~MolName, ncol=1, scales="free_y", strip.position="right") + 
  scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %b-%d %H:%M", limits=c(min(all_omes$DT),max(all_omes$DT)), timezone = "America/Los_Angeles") + 
  #theme(axis.text.x=element_text(angle=90)) +
  geom_label(data = plot_labels, mapping = aes(x=as.POSIXct("2019-04-29 04:00:00", tz="America/Los_Angeles"), y=Inf, label = Label), color="black", size=5, hjust=0, vjust=1, nudge_y=-10, nudge_x=-60) +
  mytheme 
p2
# p3 = 
#   ggplot() + geom_rect(data = anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="lightyellow") +
#   #geom_line(data=subset(all_omes_hr,MolName=="Cortisol"), aes(DT,Intensity)) + 
#   geom_line(data=subset(all_omes,MolName %in% c("Cortisol")), aes(DT,Intensity,color=MolName)) + labs(color="") +
#   geom_point(data=subset(all_omes,MolName %in% c("Cortisol")), aes(DT,Intensity,color=MolName)) + theme(axis.line.y = element_blank(), axis.text.y=element_blank(), axis.title.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x = element_text(angle=90)) + xlab("") +
#   scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %b-%d %H:%M", limits=c(min(all_omes$DT),max(all_omes$DT))) 

# 
# p5 = 
#   ggplot() + geom_rect(data = anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="lightyellow") +
#   #geom_line(data=subset(all_omes_hr,MolName=="Cortisol"), aes(DT,Intensity)) + 
#   geom_line(data=subset(all_omes,MolName %in% c("GLP1")), aes(DT,Intensity,color=MolName)) + theme(axis.text.x = element_blank()) + theme_void() + labs(color="") +
#   geom_point(data=subset(all_omes,MolName %in% c("GLP1")), aes(DT,Intensity,color=MolName)) + theme(axis.text.x = element_blank()) + labs(color="") 

#ggplot(data=subset(all_omes_2hr,MolName %in% c("Cortisol","INSULIN","LEPTIN","GHRELIN","GLP1"))) + geom_tile(aes(x=DT,y=MolName,fill=Intensity))

FoodTimeL_sum$DT = as.POSIXct(FoodTimeL_sum$DT,tz="America/Los_Angeles")
FoodTimeL_sum$MolName = factor(FoodTimeL_sum$MolName, levels = c("Calcium_PDV","Calories","Carbs_g","Fat_g","Alcohol_g","Fiber_g","Iron_PDV","Phosphorus_mg","Potassium_mg","Protein_g","Saturated_fat_g","Sodium_mg","Sugars_g", "Vitamin_A_PDV", "Vitamin_C_PDV"))
plot_labels = data.frame(MolClass = "Food", Label = "Food")
p4 = 
  ggplot() +
  geom_rect(
    data = anno_shading,
    mapping = aes(
      xmin = shading_start,
      xmax = shading_end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "lightyellow"
  ) +
  geom_point(
    data = subset(
      FoodTimeL_sum,
      Intensity > 0 &
        MolName %in% c("Carbs_g", "Fat_g", "Protein_g", "Alcohol_g")
    ),
    aes(DT, food_multiplier * Intensity, color = MolName),
    shape = "square",
    size = 3,
    alpha = .75
  ) +
  #theme_void() +
  theme(panel.background = element_rect("grey95")) + mytheme +
  labs(color = "") + scale_x_datetime(
    breaks = date_breaks("4 hour"),
    date_labels = "%a %b-%d %H:%M",
    limits = c(min(all_omes$DT), max(all_omes$DT)),
    timezone = "America/Los_Angeles"
  ) +
  facet_wrap( ~ MolClass, ncol = 1, strip.position = "right") + #+ theme(axis.text = element_blank()) #, axis.line.x = element_blank())
  geom_label(
    data = plot_labels,
    mapping = aes(
      x = as.POSIXct("2019-04-29 04:00:00", tz = "America/Los_Angeles"),
      y = Inf,
      label = Label
    ),
    color = "black",
    size = 5,
    hjust = 0,
    vjust = 1,
    nudge_y = -10,
    nudge_x = -60
  )
p4
mols = c(4,53,267,271,1415,1407,1405)
plot_labels = data.frame(unique(all_omes[all_omes$MolNum %in% mols, "MolName"]), unique(all_omes[all_omes$MolNum %in% mols, "MolName"])) %>% dplyr::rename(Label = MolName.1)
all_omes$MolName = factor(all_omes$MolName, levels = rev(unique(all_omes$MolName))) # Hack to get cortisol on top
p5 = ggplot() + 
  geom_rect(data = anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="lightyellow") +
  geom_vline(data = subset(events, grepl("milk",events$Details)), mapping=aes(xintercept=DT), color="gray",size=1) +
  geom_line(data = subset(all_omes, MolNum %in% mols), aes(DT,Intensity,color=MolName, group=MolNum)) + 
  geom_point(data = subset(all_omes, MolNum %in% mols), aes(DT,Intensity,color=MolName, group=MolNum)) + 
  facet_wrap(~MolName,ncol=1,scales="free_y", strip.position = "right") + theme(panel.background = element_rect("grey95")) +
  scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %H:%M", timezone = "America/Los_Angeles", limits=c(min(all_omes$DT),max(all_omes$DT))) +
  theme(axis.text.x = element_text(angle=90), legend.position="none") + labs(color="") + xlab("") + ylab("") + 
  geom_label(data = plot_labels, mapping = aes(x=as.POSIXct("2019-04-29 04:00:00", tz="America/Los_Angeles"), y=Inf, color=MolName, label = Label), size=5, hjust=0, vjust=1, nudge_y=-10, nudge_x=-60) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank())

p5

# ggplot() + 
# geom_vline(data = subset(events, grepl("milk",events$Details)), mapping=aes(xintercept=DT), color="gray",size=1) +
# geom_line(data = subset(all_omes, MolClass %in% "MetabolicPanel" & !grepl("CHEX", MolName)), aes(DT,Intensity,color=MolName, group=MolNum)) + 
# geom_point(data = subset(all_omes, MolClass %in% "MetabolicPanel" & !grepl("CHEX", MolName)), aes(DT,Intensity,color=MolName, group=MolNum)) + 
# facet_wrap(~MolName,ncol=2,scales="free_y") +
# scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%H:%M", timezone = "America/Los_Angeles", limits=c(min(all_omes$DT),max(all_omes$DT))) +
# theme(axis.text.x = element_text(angle=90), legend.position="none") + labs(color="")


#  
# plot(p1)
# plot(p2)
# plot(p3)
# ggplotly(p4)

#plot_grid(p2,p1, p4,p3, ncol=1, axis="tblr", align="v", rel_heights = c(2,1,5,5))
#plot_grid(p1,p2,p4,p3, ncol=1 , axis="tblr", align="v", rel_heights = c(1,5,2,13))
#plot_grid(p1,p2,p4,p3, ncol=1 , axis="tblr", align="v", rel_heights = c(1,5,16,13))

plot_grid(p1,p2,p4, p5, ncol=1, axis="tblr", align="v", rel_heights = c(2,4,2,12))

# ggplot(all_omes_hr, aes(DT,MolName)) + geom_raster(aes(fill=zscore(Intensity))) + facet_grid(MolClass~., scales="free",space="free") + theme(axis.text.y=element_blank()) 
# #line 323
# all_omes_hr_clust = left_join(all_omes_hr,pdat[,c("MolNum","ClusterName")])
# all_omes_hr_clust = all_omes_hr_clust[!duplicated(all_omes_hr_clust),]
# all_omes_hr_clust = all_omes_hr_clust %>% mutate(Day = date(DT), Time = hms::as.hms(DT), Hr = hour(DT))
# 
# all_omes_hr_clust_avg = subset(all_omes_hr_clust,!is.na(ClusterName)) %>% group_by(DT, ClusterName) %>% summarize(clust_avg = mean(Intensity))
# 
# ggplot() + geom_rect(data = anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="lightyellow") + geom_line(data=all_omes_hr_clust_avg, mapping=aes(DT, clust_avg, group=ClusterName)) +
#   facet_grid(ClusterName~.,space="free", scales="free") + theme(legend.position="top") + xlab("Hour of day") + ylab("Intensity (z-score)") + labs(color = "Molecule Type:") 
# 
# p4 = ggplot(all_omes_hr_clust_avg, aes(DT, Intensity, color=MolClass, group=MolNum)) + geom_line(alpha=.5) +
#   geom_smooth(se=T, span=.1, color="black", aes(group=MolNum)) +  #geom_smooth(span=0.01,color="black",se=F, aes(group=Cluster)) +
#   facet_grid(ClusterName~.,space="free", scales="free") + theme(legend.position="top") + xlab("Hour of day") + ylab("Intensity (z-score)") + labs(color = "Molecule Type:")
# 
# 
# p4







# Salicylic acid

ggplot(subset(metab, MolName=="Salicylic acid"), aes(DT,Intensity)) + geom_line() + geom_point() 




# PDF PLOTS

library(ggdendro)

all_omes_z = all_omes %>% select(-KEGG_val, -HMDB_ID, -CAS_val, -Super.pathway, -Sub.pathway) %>% group_by(MolName,MolClass,MolSubclass,MolNum) %>% mutate(Intensity_z = (Intensity - mean(Intensity))/sd(Intensity)) 
all_omes_z = all_omes_z[!is.na(all_omes_z$Intensity_z),]
all_omes_z = all_omes_z[all_omes_z$DT < as_date("2019-05-07"),]

all_omes_zw = all_omes_z %>% ungroup() %>% subset(!duplicated(MolName)) %>% select(MolName,SampleID,Intensity_z) %>% spread(SampleID,Intensity_z)

all_omes_hr_w = all_omes_hr %>% ungroup() %>% select(MolName,DT,Intensity) %>% spread(DT,Intensity)

all = column_to_rownames(all_omes_hr_w,"MolName")
dd <- dist(all, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
#plot(hc)
#plot(as.phylo(hc), type="fan")

dend = as.dendrogram(hc)
dend_data = dendro_data(dend,type="rectangle")

ord = label(dend_data) %>% dplyr::rename(MolName = label)

allo = all_omes %>% arrange(factor(MolName, levels = ord$MolName))
allo$MolName = factor(allo$MolName, levels=ord$MolName)
allo$MolNumName = paste0(allo$MolNum,"_",allo$MolName)

i = 1
for (i in 1:300) {
  ggplot() + 
    geom_vline(data = subset(events, grepl("milk",events$Details)), mapping=aes(xintercept=DT), color="gray",size=1) +
    geom_line(allo, mapping=aes(x = DT, y = Intensity, group=MolNumName, color=MolNumName)) +
    geom_point(allo, mapping=aes(x = DT, y = Intensity, group=MolNumName, color=MolNumName)) +
    #geom_line(allo, mapping=aes(x = DT, y = Intensity), color="red", size=3) +   
    #geom_smooth(color="red", aes(group=Hr_day), se=T) + geom_line(alpha=0.5, aes(group=day)) + 
    theme(legend.position="none") +
    scale_x_datetime(breaks=date_breaks("8 hour"), date_labels = "%H:%M") +
    #    theme_bw() + 
    #facet_wrap(~MolName, scales="free", ncol=3, nrow=9)
    facet_wrap_paginate(~MolNumName, scales="free", ncol=2, nrow=9,page=i)
  
  
  ggsave(paste("Plot_Similarity_3/",i,".png",sep=""), width=28, height=16, units="in")
  
  print(i)
  
  
}





# MUSCLE MILK ANALYSIS


metab = metab %>% mutate(day = date(DT), hour_of_day=hour(DT), tod = hms::as.hms(DT, tz="America/Los_Angeles"))

subset(metab,grepl(""))

ggplot() + 
  geom_vline(data = subset(events, grepl("milk",events$Details)), mapping=aes(xintercept=tod), color="red",size=2) +
  geom_line(data = subset(all_omes, MolNum %in% c(4,183,187)), aes(tod,Intensity,color=MolName, group=MolNum)) + 
  geom_point(data = subset(all_omes, MolNum %in% c(4,183,187)), aes(tod,Intensity,color=MolName, group=MolNum)) + 
  #geom_line(data = subset(all_omes, MolName=="Caffeine"), aes(tod,Intensity,color=MolName)) +
  facet_wrap(~day, ncol=1) 


ggplot() + 
  geom_vline(data = subset(events, grepl("milk",events$Details)), mapping=aes(xintercept=DT), color="gray",size=1) +
  geom_line(data = subset(all_omes_2hr, MolNum %in% c(4,53,267,271,1437)), aes(DT,Intensity,color=MolName, group=MolNum)) + 
  geom_point(data = subset(all_omes_2hr, MolNum %in% c(4,53,267,271,1437)), aes(DT,Intensity,color=MolName, group=MolNum)) + facet_wrap(~MolName,ncol=1,scales="free_y")


ggplot(data = steps, aes(x=Time,y=value,color=type,group=1)) + geom_line() +
  #scale_x_datetime(breaks=date_breaks("2 hour"), date_labels = "%H:%M") + theme(axis.text.x = element_blank()) +
  geom_line(data = hr, aes(x=Time,y=value, color=type, group=1)) + geom_line(data = glu, aes(x=Time, y=value, color=type, group=1)) + theme_void() + labs(color="") +
  #scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %b-%d %H:%M", limits=c(min(all_omes$DT),max(all_omes$DT)), timezone = "America/Los_Angeles") 
  
  facet_wrap(~Day, ncol=1) 

subset(all_omes, grepl("Salicylic",MolName))$MolNum #4
unique(subset(all_omes, grepl("Hydroxyphenyllactic",MolName))$MolNum) #37
subset(all_omes, grepl("benzenetriol",MolName))$MolNum #183
subset(all_omes, grepl("Caffeine",MolName))$MolNum #187
unique(subset(all_omes, grepl("Cortisol",MolName))$MolNum) #1415
unique(subset(all_omes, grepl("^PP",MolName))$MolNum) #1407
unique(subset(all_omes, grepl("LEPTIN",MolName))$MolNum) #1407
unique(subset(all_omes, grepl("TNFA",MolName))$MolNum) #1454


d=subset(metab, grepl("L-Acetylcarnitine",MolName)) #37





# cytokine correction

# metabolic protein panel
MetaData=read.csv("~/Box/Microsampling/Analysis/annotations/Copy of 24_7_Stability_MasterIDReferenz_DH4_KE2.csv") 

MetaData = MetaData[MetaData$Study == "F0" & MetaData$Original_SampleID %in% 1:97,]
MetaData$Name = paste("DH-",MetaData$Original_SampleID,sep="")

# total protein
tp = read.csv("~/Box/Microsampling/cytokines/27-4CytokineTotalProtein.csv")
tp = left_join(tp, MetaData[,c("Name","CollectionTime","PrepIndex")])
tp$DT = as.POSIXct(tp$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
tp = dplyr::rename(tp, SampleID=PrepIndex, Intensity=Total_protein)
tpL = tp # gather(tp,"Cytokine","MFI",5)
tpL$Intensity = zscore(tpL$Intensity)
#tpL$Intensity = log10(as.numeric(tpL$MFI))
tpL$MolClass = "BCA"
tpL$MolSubclass = NA
tpL$MolName = "TotalProtein"
tpL = select(tpL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)


# metabolic protein panel
MP = read.csv("~/Box/Microsampling/cytokines/24-7 omics-MetabolicHormone-Ryan Metabolic Plate-1.csv")
MP = MP[2:97,]
MP = left_join(MP,MetaData[,c("Name","CollectionTime","PrepIndex")])
MP$DT = as.POSIXct(MP$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
MPL = gather(MP,"Cytokine","MFI",5:21)
MPL$Intensity = log10(as.numeric(MPL$MFI))
MPL = dplyr::rename(MPL, MolName = Cytokine, SampleID=PrepIndex)
MPL$MolClass = "MetabolicPanel"
MPL$MolSubclass = NA
MPL = select(MPL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)
MPL[MPL$MolName=="TNFA",]$MolName = "TNFA_MP"
MPL[MPL$MolName=="MCP1",]$MolName = "MCP1_MP"
MPL[MPL$MolName=="CHEX1",]$MolName = "CHEX1_MP"
MPL[MPL$MolName=="CHEX2",]$MolName = "CHEX2_MP"
MPL[MPL$MolName=="CHEX3",]$MolName = "CHEX3_MP"
MPL[MPL$MolName=="CHEX4",]$MolName = "CHEX4_MP"
MPL[MPL$MolName=="IL6",]$MolName = "IL6_MP"


# cortisol 
cort = read.csv("~/Box/Microsampling/cytokines/RYAN CORTISOL 8-28-19.csv")
cort = cort[10:102,]
cort = left_join(cort,MetaData[,c("Name","CollectionTime","PrepIndex")])
cort$DT = as.POSIXct(cort$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
cort = dplyr::rename(cort, SampleID=PrepIndex)
cortL = gather(cort,"Cytokine","MFI",5)
cortL$Intensity = log10(as.numeric(cortL$MFI))
cortL$MolClass = "CortisolEnzymatic"
cortL$MolSubclass = NA
cortL$MolName = "Cortisol"
cortL = select(cortL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)

#41-plex
cyto = read.csv("~/Box/Microsampling/cytokines/24-7 omics-HumanLuminexMAG42plex-Mitra H41 Plex-1.csv")
cyto = cyto[10:105,]
cyto = left_join(cyto,MetaData[,c("Name","CollectionTime","PrepIndex")])
cyto$DT = as.POSIXct(cyto$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
cytoL = gather(cyto,"Cytokine","MFI",5:49)
cytoL$Intensity = log10(as.numeric(cytoL$MFI))
cytoL$MolClass = "Cytokine_41Plex"
cytoL$MolSubclass = NA
cytoL = dplyr::rename(cytoL, SampleID = PrepIndex, MolName = Cytokine)
cytoL = select(cytoL, SampleID, DT, MolName, Intensity, MolClass, MolSubclass)

all = bind_rows(MPL, cortL, cytoL, tpL)

tpL = dplyr::rename(tpL, TP = Intensity)
allTP = left_join(all, tpL[,c("SampleID","TP")])


ggplot(allTP,aes(Intensity,TP)) + geom_point() + facet_wrap(~MolName, scales="free") + geom_smooth(method="lm")

all = subset(all,!(SampleID %in% c(106,25)))
ggplot(all, aes(DT,Intensity)) + geom_line() + geom_point() + facet_wrap(~MolName, ncol=1, scales="free_y")

ggplot(all, aes(DT,Intensity)) + geom_boxplot() + facet_wrap(~MolName, scales="free_y")

#all %>% group_by(MolName) %>% mutate 








# CORRELATION ANALYSIS

library(d3heatmap)

all_omes_wear_w = all_omes_wear %>% select(MolName,MolClass,DT,Intensity) %>% spread(DT, Intensity)
anno_rc = data.frame(Class = all_omes_wear_w$MolClass)
rownames(anno_rc) = all_omes_wear_w$MolName
all_omes_wear_w = select(all_omes_wear_w,-MolClass)
all = column_to_rownames(all_omes_wear_w, "MolName")
all_cor = cor(t(all))
d3heatmap(all_cor)
#pheatmap(all_cor, annotation_row = anno_rc, annotation_col = anno_rc)

get_cor = as.data.frame(t(all_cor[c("CGM","Steps","HR","Cortisol","LEPTIN","INSULIN","TNFA"),])) %>% rownames_to_column("MolName")
target = c("CGM","Steps","HR","Cortisol","LEPTIN","INSULIN","TNFA")
for (i in 1:length(target)) {
  # Positive correlation 
  get_cor = subset(get_cor, !grepl("CHEX", MolName) & !grepl(target[i], MolName))
  top = get_cor[order(get_cor[,target[i]], decreasing=T),"MolName"]
  a=subset(all_omes_wear, MolName %in% top[1:10])
  dat = left_join(a,subset(all_omes_wear, MolName==target[i]) %>% ungroup %>% select("DT","Intensity"), by="DT") 
  dat$MolName = factor(dat$MolName, levels = top)
  ggplot() + 
    geom_line(dat, mapping=aes(DT,Intensity.y,group=MolName),color="green") + 
    geom_line(dat, mapping=aes(DT,Intensity.x,group=MolName),color="red") + 
    facet_wrap(~MolName,ncol=1,scales="free_y") +
    theme(strip.background = element_rect(fill="white")) + 
    theme(axis.title.y = element_text(color="green")) +
    theme(strip.text = element_text(colour = 'red')) + xlab("Date") + ylab(target[i])
  ggsave(paste0(target[i],"_pos_corr.png"), width=9, height=12, units="in")
  
  # Negative correlation 
  get_cor = subset(get_cor, !grepl("CHEX", MolName) & !grepl(target[i], MolName))
  top = get_cor[order(get_cor[,target[i]], decreasing=F),"MolName"]
  a=subset(all_omes_wear, MolName %in% top[1:10])
  dat = left_join(a,subset(all_omes_wear, MolName==target[i]) %>% ungroup %>% select("DT","Intensity"), by="DT") 
  dat$MolName = factor(dat$MolName, levels = top)
  ggplot() + 
    geom_line(dat, mapping=aes(DT,Intensity.y,group=MolName),color="green") + 
    geom_line(dat, mapping=aes(DT,Intensity.x,group=MolName),color="red") + 
    facet_wrap(~MolName,ncol=1,scales="free_y") +
    theme(strip.background = element_rect(fill="white")) + 
    theme(axis.title.y = element_text(color="green")) +
    theme(strip.text = element_text(colour = 'red')) + xlab("Date") + ylab(target[i])
  ggsave(paste0(target[i],"_neg_corr.png"), width=9, height=12, units="in")
  print(i)
}

# Positive correlation with Cortisol
get_cor = subset(get_cor, !grepl("CHEX", MolName) & !grepl(target[i], MolName))
top = get_cor[order(get_cor[,target[i]], decreasing=T),"MolName"]
a=subset(all_omes_wear, MolName %in% top[1:10])
dat = left_join(a,subset(all_omes_wear, MolName==target[i]) %>% ungroup %>% select("DT","Intensity"), by="DT") 
dat$MolName = factor(dat$MolName, levels = top)
ggplot() + 
  geom_line(dat, mapping=aes(DT,Intensity.y,group=MolName),color="green") + 
  geom_line(dat, mapping=aes(DT,Intensity.x,group=MolName),color="red") + 
  facet_wrap(~MolName,ncol=1,scales="free_y") +
  theme(strip.background = element_rect(fill="white")) + 
  theme(axis.title.y = element_text(color="green")) +
  theme(strip.text = element_text(colour = 'red')) + xlab("Date") + ylab(target[i])
ggsave(paste0(target[i],"_pos_corr.png"), width=9, height=12, units="in")

# Negative correlation with Cortisol
get_cor = subset(get_cor, !grepl("CHEX", MolName) & !grepl(target[i], MolName))
top = get_cor[order(get_cor[,target[i]], decreasing=F),"MolName"]
a=subset(all_omes_wear, MolName %in% top[1:10])
dat = left_join(a,subset(all_omes_wear, MolName==target[i]) %>% ungroup %>% select("DT","Intensity"), by="DT") 
dat$MolName = factor(dat$MolName, levels = top)
ggplot() + 
  geom_line(dat, mapping=aes(DT,Intensity.y,group=MolName),color="green") + 
  geom_line(dat, mapping=aes(DT,Intensity.x,group=MolName),color="red") + 
  facet_wrap(~MolName,ncol=1,scales="free_y") +
  theme(strip.background = element_rect(fill="white")) + 
  theme(axis.title.y = element_text(color="green")) +
  theme(strip.text = element_text(colour = 'red')) + xlab("Date") + ylab(target[i])
ggsave(paste0(target[i],"_neg_corr.png"), width=9, height=12, units="in")

### WITH FOOD
# all_omes_wear_food_w = all_omes_wear_food %>% select(MolName,MolClass,DT,Intensity) %>% spread(DT, Intensity)
# anno_rc = data.frame(Class = all_omes_wear_food_w$MolClass)
# rownames(anno_rc) = all_omes_wear_food_w$MolName
# all_omes_wear_food_w = select(all_omes_wear_food_w,-MolClass)
# all = column_to_rownames(all_omes_wear_food_w, "MolName")
# all_cor = cor(t(all))
# pheatmap(all_cor, annotation_row = anno_rc, annotation_col = anno_rc)
# 
# get_cor = as.data.frame(t(all_cor[c("CGM","Steps","HR","Cortisol","LEPTIN","INSULIN","TNFA","Carbs_g","Fat_g","Protein_g","Calories"),]))


#ccc = correlate(t(all))
# library(tidyquant)
# library(corrr)
# network_plot(all_cor)



# PCA

library(ggfortify)
library(cluster)
library(viridis)


all_omes_pca = all_omes %>% ungroup() %>% subset(!MolClass %in% "Protein") %>% select(-MolClass, -MolNum, -MolSubclass, -KEGG_val, -HMDB_ID, -CAS_val, -Sub.pathway, -Super.pathway) 
all_omes_pca = all_omes_pca %>% spread(MolName, Intensity)
all_omes_pca = subset(all_omes_pca, !SampleID %in% c(50,52))

all_omes_pca = all_omes_pca %>% arrange(factor(DT, levels = sort(DT)))
all_omes_pca[is.na(all_omes_pca)] = 0


autoplot(prcomp(all_omes_pca[,-c(1,2)], center=T, scale=T), data=all_omes_pca, colour="DT", palette = "viridis") + labs(colour = "Date")

p = prcomp(all_omes_pca[,-c(1,2)], center=T, scale=T)

autoplot(clara(as.data.frame(p$x),3))






# Compare duplicated cytokines


ggplot(subset(all_omes_hr, MolName %in% c("TNFA","TNFA_MP")), aes(x=Hr, y=Intensity, group=MolName)) + geom_line() +
  geom_vline(events, mapping = aes(xintercept = Hr, color = Type, linetype=Type))

ggplot(subset(all_omes_hr, MolName %in% c("MCP1","MCP1_MP")), aes(x=Hr, y=Intensity, group=MolName, color=MolName)) + geom_line() +
  geom_vline(events, mapping = aes(xintercept = Hr, color = Type, linetype=Type))

ggplot(subset(all_omes_hr, MolName %in% c("IL6","IL6_MP")), aes(x=Hr, y=Intensity, group=MolName, color=MolName)) + geom_line()

s1 = subset(all_omes_hr, MolName %in% c("TNFA","MCP1", "IL6"))
s2 = subset(all_omes_hr, MolName %in% c("TNFA_MP","MCP1_MP", "IL6_MP"))

s3 = s1 %>% select(MolName, Hr, Intensity) %>% dplyr::rename(MolName1 = MolName)
s3$Intensity2 = s2$Intensity

ggplot(s3,aes(Intensity,Intensity2, color=MolName1)) + geom_point() + facet_wrap(~MolName1, ncol=1) + stat_cor(aes(label=paste0(..r.label..)), method = "pearson", label.x = -1.5, label.y = 3) + theme(legend.position = "none") + xlab("Metabolic Panel") + ylab("41-Plex Panel")

ggscatter(s3, x="Intensity", y = "Intensity2", add="reg.line", conf.int=T)




# IDENTIFY CIRCADIAN FLUCTUATIONS

# metab_hr_w = all_omes_hr %>% subset(MolClass=="Metabolite") %>% 
#   select("MolName","Hr","Intensity") %>% spread(Hr, Intensity)
# metab_hr_w = metab_hr_w[,order(as.numeric(colnames(metab_hr_w)))]
# metab_hr_w = metab_hr_w %>% select(MolName, everything())
# write.csv(metab_hr_w,"metab_hr_w.csv", row.names=F)

meta2d(infile="metab_hr_w.csv", filestyle="csv", timepoints="Line1", outdir = "metab")
res = read.csv("metab/meta2d_metab_hr_w.csv")

all_omes_hr_w = all_omes_hr %>% ungroup() %>%
  select("MolName","Hr","Intensity") %>% spread(Hr, Intensity)
all_omes_hr_w = all_omes_hr_w[,order(as.numeric(colnames(all_omes_hr_w)))]
all_omes_hr_w = all_omes_hr_w %>% select(MolName, everything())
write.csv(all_omes_hr_w,"all_omes_hr_w_9feb20.csv", row.names=F)

meta2d(infile="all_omes_hr_w_9feb20.csv", filestyle="csv", timepoints="Line1", outdir = "all_omes_cyc")




res = read.csv("all_omes_cyc/meta2d_all_omes_hr_w.csv")
sel = res[res$JTK_BH.Q < 1e-5,"CycID"]

#all_omes_hr = all_omes_hr %>% mutate(Intensity_smooth = smooth(Intensity))
#all_omes_hr = all_omes_hr %>% group_by(MolName,Hr_day) %>% mutate(Day_avg = mean(Intensity))
all_omes_hr = all_omes_hr %>% group_by(MolName,Hr_day) %>% mutate(Day_avg = mean(Intensity))

ggplot(subset(all_omes_hr, MolName %in% sel), aes(x=Hr, y=Intensity, group=MolName)) + geom_line() + facet_wrap(~MolName)

ggplot(subset(all_omes_hr, MolName %in% sel), aes(x=Hr, y=MolName, fill = Intensity)) + geom_tile() + scale_fill_gradient(low="white", high="blue")

ggplot(subset(all_omes_hr, MolName %in% sel), aes(x=Hr_day, y=MolName, fill = Day_avg)) + geom_tile() + scale_fill_gradient(low="white", high="blue")

ggplot(subset(all_omes_hr, MolName %in% sel), aes(x=Hr_day, y=Intensity, group=MolName)) + geom_smooth() + facet_wrap(~MolName)

sel = c(
  "Hydroxyphenyllactic.acid",
  "LysoPE.20.2..1",
  "LysoPE.20.1..1",
  "C20.4.DC.FA.1",
  "X1.2.3.benzenetriol.sulfate", #noncyc:
  "C12H14O5_pos",
  "C10.1.DC.FA..Decenedicarboxylic.acid.",
  "C10H19NO3",
  "C15H16N2O2")

subset(res, res$CycID %in% sel)


# cluster and plot 24hr averages CIRCADIAN


sel = res[res$JTK_BH.Q < 1e-3,"CycID"]

all_omes_hr_avg = all_omes_hr %>% group_by(MolName,MolNum,MolClass,MolSubclass,Hr_day) %>% summarize(Hr_avg = as.numeric(mean(Intensity)))

#s = select(all_omes_hr_avg,-c(MolNum,MolClass,MolSubclass)) %>% subset(MolName=="TNFA") %>% group_by(MolName,Hr_day) %>% mutate(Hr_avg_smooth = smooth(Hr_avg))

all_omes_hr_avg_w = subset(all_omes_hr_avg, (MolName %in% sel) & !(MolSubclass %in% c("TAG"))) %>% spread(Hr_day,Hr_avg) %>% ungroup

d = select(all_omes_hr_avg_w, -c("MolName","MolClass","MolSubclass","MolNum"))
m = select(all_omes_hr_avg_w, c("MolName","MolClass","MolSubclass","MolNum"))


pc <- tsclust(d, type = "partitional", k = 30L, # preproc = zscore, 
              distance = "dtw_basic", centroid = "pam", 
              seed = 3247L, trace = TRUE,
              args = tsclust_args(dist = list(window.size = 20L)))

pc <- tsclust(d, type = "h", k = 20L, preproc = zscore, 
              distance = "sbd", centroid = shape_extraction, 
              seed = 3247L, control = hierarchical_control(method = "average"))

pc_dtw <- tsclust(
  d,
  k = 20L,
  distance = "dtw_basic",
  centroid = "dba",
  trace = TRUE,
  seed = 8,
  norm = "L2",
  window.size = 20L,
  args = tsclust_args(cent = list(trace = TRUE))
)

fc <- tsclust(d, type = "f", k = 20L,
              preproc = zscore, distance = "L2",
              seed = 42)

plot(pc)
plot(pc, type = "sc")
plot(fc)

all_omes_hr_avg_w$Cluster = pc@cluster

pdat = gather(all_omes_hr_avg_w,"Hr_day", "Hr_avg", -c("MolName","MolClass","MolSubclass","MolNum","Cluster"))
pdat$Hr_day = as.numeric(pdat$Hr_day)
#pdat_med = pdat %>% group_by(Cluster, Hr) %>% summarize(Int_med = median(Intensity))

pdat[pdat$Cluster %in% c(3,5,12,17,23,22),"ClusterName"] = "Meal Responsive"
pdat[pdat$Cluster %in% c(8,24,29),"ClusterName"] = "Sharp Mid-day Dip"
pdat[pdat$Cluster %in% c(13),"ClusterName"] = "Late Peak"
pdat[pdat$Cluster %in% c(9,14,16,19,30),"ClusterName"] = "Mid-day Dip"
pdat[pdat$Cluster %in% c(1,2),"ClusterName"] = "Mid-day Peak"
pdat[pdat$Cluster %in% c(27),"ClusterName"] = "Early Peak"


ggplot(subset(pdat,!is.na(ClusterName)), aes(Hr_day, Hr_avg, color=MolClass, group=MolNum)) + geom_line(alpha=.5) + 
  geom_smooth(se=T, color="black", aes(group=ClusterName)) +  #geom_smooth(span=0.01,color="black",se=F, aes(group=Cluster)) + 
  facet_wrap(~ClusterName) + theme(legend.position="top") + xlab("Hour of day") + ylab("Intensity (z-score)") + labs(color = "Molecule Type:")

# #sets = {c(1,2,4,5,6,7,8,9,11,12,13,14,)}
# ggplot(subset(pdat,T), aes(Hr_day, Hr_avg, color=MolClass, group=MolNum)) + geom_line() + 
#   #geom_smooth(se=F) +  #geom_smooth(span=0.01,color="black",se=F, aes(group=Cluster)) + 
#   facet_wrap(~Cluster) + theme(legend.position="bottom") 
# ggplotly(p)
# 
# ggplot(subset(pdat,MolClass=="Lipid" & MolSubclass=="TAG"), aes(Hr_day, Hr_avg, color=MolClass, group=MolNum)) + geom_line() + #geom_smooth(span=0.01,color="black",se=F, aes(group=Cluster)) + 
#   facet_wrap(~Cluster) + theme(legend.position="bottom") 
# 
# ggplot(subset(pdat,MolName=="Cortisol" | MolName=="GLP1"), aes(as.numeric(Hr), Intensity, color=MolName, group=MolNum)) + geom_line() + #geom_smooth(span=0.01,color="black",se=F, aes(group=Cluster)) + 
#   facet_wrap(~Cluster) + theme(legend.position="bottom") 
# 
# x <- rbinom(100, size=2, p=0.1)
# y <- rbinom(100, size=2, p=0.1)
# 
# ## try these in an R markdown document for best results
# kTable(x)
# my_table <- kTable(x, y, top.left.cell="foo", left.label="bar", top.label="baz")
# pxt( my_table )

pdc = select(pdat,MolName,ClusterName)
all_omes_clust = left_join(all_omes, pdc[!duplicated(pdc),], by="MolName")



#pdat[pdat$Cluster %in% c(3,5,12,17,23,22),"ClusterName"] = "Meal Responsive"
pdat[pdat$Cluster %in% c(8,24,29),"ClusterName"] = "Sharp Mid-day Dip"
pdat[pdat$Cluster %in% c(13),"ClusterName"] = "Late Peak"
pdat[pdat$Cluster %in% c(5,6,7),"ClusterName"] = "Mid-day Dip"
pdat[pdat$Cluster %in% c(26,27,11),"ClusterName"] = "Mid-day Peak"
pdat[pdat$Cluster %in% c(25,29,30,3,14, 9),"ClusterName"] = "Early Peak"




# 24-HOUR PLOTS


library(viridis)


i = 1
for (i in 101:200) {
  ggplot() + geom_line(all_omes_hr, mapping=aes(x = Hr_day, y = Intensity, group=day, color=day)) +
    geom_line(all_omes_hr_avg, mapping=aes(x=Hr_day, y=Hr_avg), color="red", size=3) +   
    #geom_smooth(color="red", aes(group=Hr_day), se=T) + geom_line(alpha=0.5, aes(group=day)) + 
    
    theme_bw() + 
    facet_wrap_paginate(~MolName, scales="free", ncol=2, nrow=7, page=i)
  
  ggsave(paste("all_omes_day_plots/all_omes_day_",i,".png",sep=""), width=8, height=12, units="in")
  
  print(i)
  
  
}


i = 1
for (i in 1:100) {
  ggplot() + 
    geom_line(all_omes_hr, mapping=aes(x = Hr_day, y = Intensity, group=day, color=day, alpha=0.5)) +
    geom_smooth(all_omes_hr, mapping=aes(x = Hr_day, y = Intensity, group=MolName), color="red", span=0.5) +
    theme_bw() + 
    facet_wrap_paginate(~MolName, scales="free", ncol=6, nrow=6, page=i)
  
  ggsave(paste("all_omes_day_plots/all_omes_day_smooth_",i,".png",sep=""), width=20, height=12, units="in")
  
  print(i)
  
  
}






# PLOT CYTOKINES

MetaData=read.csv("~/Box/Microsampling/Analysis/annotations/Copy of 24_7_Stability_MasterIDReferenz_DH4_KE2.csv") 

MetaData = MetaData[MetaData$Study == "F0" & MetaData$Original_SampleID %in% 1:97,]
MetaData$Name = paste("DH-",MetaData$Original_SampleID,sep="")

MP = read.csv("~/Box/Microsampling/cytokines/24-7 omics-MetabolicHormone-Ryan Metabolic Plate-1.csv")
MP = MP[2:97,]

MP = left_join(MP,MetaData[,c("Name","CollectionTime","PrepIndex")])
MP$DT = as.POSIXct(MP$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")
MPL = gather(MP,"Cytokine","MFI",5:21)
MPL$Intensity = log10(as.numeric(MPL$MFI))
MPL = dplyr::rename(MPL, MolName = Cytokine, SampleID=PrepIndex)
MPL$MolClass = "MetabolicPanel"
MPL$MolSubclass = NA
MPL = MPL[MPL$DT < as_date("2019-05-07"),]

i = 1
for (i in 1:100) {
  
  
  ggplot(MPL, aes(x = DT, y = as.numeric(MFI))) + geom_point() + geom_line() + facet_wrap_paginate(~MolName, scales="free", ncol=3, nrow=8, page=i) 
  
  ggsave(paste("MetabolicPanel_",i,".png",sep=""), width=22, height=11, units="in")
  
  print(i)
  
}


#41-plex
cyto = read.csv("~/Box/Microsampling/cytokines/24-7 omics-HumanLuminexMAG42plex-Mitra H41 Plex-1.csv")
cyto = cyto[10:105,]

cyto = left_join(cyto,MetaData[,c("Name","CollectionTime","PrepIndex")])
cyto$DT = as.POSIXct(cyto$CollectionTime, tz = "America/Los_Angeles", format = "%m/%d/%y %H:%M")

cytoL = gather(cyto,"Cytokine","MFI",5:49)
cytoL$Intensity = log10(as.numeric(cytoL$MFI))
cytoL$MolClass = "Cytokine_41Plex"
cytoL$MolSubclass = NA
cytoL = dplyr::rename(cytoL, SampleID = PrepIndex, MolName = Cytokine)
cytoL = cytoL[cytoL$DT < as_date("2019-05-07"),]
cytoL$MFI = as.numeric(cytoL$MFI)
#cytoL = cytoL %>% group_by(MolName) %>% mutate(m = mean(MFI), csd_l = mean(MFI) - 3*sd(MFI), csd_u = mean(MFI) + 3*sd(MFI), MFI = ifelse(MFI < csd_l | MFI > csd_u, NA, MFI))

mod = lm(TNFA ~ CHEX4, cyto)
tidy(mod)
r = as.data.frame(data = mod$residuals)
#r = dplyr::rename
#ggplot(r, aes(rowna)

#MFI = 
ifelse(cytoL$MFI < quantile(cytoL$MFI,0.000001,na.rm=TRUE) || cytoL$MFI < quantile(cytoL$MFI, 0.999999,na.rm=TRUE), NA, cytoL$MFI)

ggplot(subset(cytoL, MolName=="RANTES"), aes(x = DT, y = as.numeric(MFI))) + geom_point() + geom_line() + geom_hline(aes(yintercept = m))

i = 1
for (i in 1:100) {
  
  
  ggplot(cytoL, aes(x = DT, y = as.numeric(MFI))) + geom_point() + geom_line() + facet_wrap_paginate(~MolName, scales="free", ncol=3, nrow=8, page=i) 
  
  ggsave(paste("Cytokine41-plex_",i,".png",sep=""), width=22, height=11, units="in")
  
  print(i)
  
}

# CD40L, EOTAXIN, IL1A, IL1B
# MCP1 -end


# PLOT CIRCADIAN 


library(ggforce)

res = read.csv("all_omes_cyc/meta2d_all_omes_hr_w.csv")
sel = res[res$JTK_BH.Q < 1e-10,"CycID"]

all_omes_hr_d = all_omes_hr %>% mutate(Day = day(DT))
all_omes_hr_d$Hr = as.numeric(all_omes_hr_d$Hr)
anno_t = as.data.frame(select(subset(ungroup(all_omes_hr_d),MolNum==1), DT, Day, Hr, Hr_day, AM))
anno_shading = anno_t %>% group_by(Day) %>% subset(AM) %>% summarize(shading_start = first(Hr), shading_end = last(Hr))

pdata = subset(all_omes_hr, MolName %in% sel)
pdata$Hr = as.numeric(pdata$Hr)
i = 1
for (i in 1:100) {
  
  ggplot() + geom_rect(anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), color=NA, fill="lightyellow") +
    geom_line(pdata, mapping = aes(x=Hr, y = Intensity, group=MolNum)) + geom_point(pdata, mapping = aes(x=Hr, y = Intensity, group=MolNum)) + scale_x_continuous(breaks=seq(0,180,4)) + 
    facet_wrap_paginate(~MolNum, scales="free", ncol=3, nrow=8, page=i) 
  
  
  ggsave(paste("Circ_",i,".png",sep=""), width=22, height=11, units="in")
  
  print(i)
  
}

ggplot() + geom_rect(anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), color=NA, fill="lightyellow") + 
  geom_line(pdata, mapping = aes(x=Hr, y = Intensity, group=MolNum), alpha=.1) + scale_x_continuous(breaks=seq(0,180,4)) + geom_smooth(pdata, mapping = aes(x=Hr, y = Intensity, group=MolNum), stat="summary", fun.y="median", size=2, color="red") + 
  geom_vline(events, mapping = aes(xintercept = Hr, color = Type))



s = pdata %>% ungroup() %>% select(-DT, -Hr_day, -AM, -day, -hour_of_day, -t) %>% spread(Hr,Intensity)

#all_omes_w = s %>% ungroup() %>% select(-DT) %>% spread(Hr,Intensity)

d = select(s, -c("MolName","MolClass","MolSubclass","MolNum"))
m = select(s, c("MolName","MolClass","MolSubclass","MolNum"))
d = d[,order(as.numeric(colnames(d)))]

dd=as.matrix(d)
df = as.data.frame(d)
rownames(df) = m$MolName 

pheatmap(df, cluster_rows = T, cluster_cols = F, annotations_col = anno_t) #, annotation_row = m) #labels_col = colnames(d), annotation_col = t(tt["AM",]))


ha1 = HeatmapAnnotation(Hr = anno_t$Hr_day)
Heatmap(df, cluster_columns = F, top_annotation=ha1, show_column_names = F)











#temp
mols = c(4,271,1415)
plot_labels = data.frame(unique(all_omes[all_omes$MolNum %in% mols,"MolName"]), unique(all_omes[all_omes$MolNum %in% mols,"MolName"])) %>% dplyr::rename(Label = MolName.1)
all_omes$MolName = factor(all_omes$MolName, levels = rev(unique(all_omes$MolName))) # Hack to get cortisol on top
ggplot() + 
  geom_rect(data = anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="lightyellow") +
  #geom_vline(data = subset(events, grepl("milk",events$Details)), mapping=aes(xintercept=DT), color="gray",size=1) +
  geom_line(data = subset(all_omes, MolNum %in% mols), aes(DT,Intensity,color=MolName, group=MolNum)) + 
  geom_point(data = subset(all_omes, MolNum %in% mols), aes(DT,Intensity,color=MolName, group=MolNum)) + 
  facet_wrap(~MolName,ncol=1,scales="free_y", strip.position = "right") + theme(panel.background = element_rect("grey95")) +
  scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %H:%M", timezone = "America/Los_Angeles", limits=c(min(all_omes$DT),max(all_omes$DT))) +
  theme(axis.text.x = element_text(angle=90), legend.position="none") + labs(color="") + xlab("") + ylab("") + 
  geom_label(data = plot_labels, mapping = aes(x=as.POSIXct("2019-04-29 04:00:00", tz="America/Los_Angeles"), y=Inf, color=MolName, label = Label), size=5, hjust=0, vjust=1, nudge_y=-10, nudge_x=-60) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank())



# INFLAMMATORY EVENT BIG ALL-DATA FIG

all_omes_hr_d = all_omes_hr %>% mutate(Day = day(DT))
anno_t = as.data.frame(select(subset(ungroup(all_omes_hr_d),MolNum==1), DT, Day, Hr, Hr_day, AM))
# anno_shading = anno_t %>% group_by(Day) %>% subset(AM) %>% summarize(shading_start = first(DT), shading_end = last(DT))
anno_shading = anno_t %>% group_by(Day) %>% summarise(shading_start = DT[(Hr_day==6)], shading_end = DT[(Hr_day == 18)])
events = subset(events, DT < "2019-05-07")


mytheme = theme(axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.line.x = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                panel.background = element_rect("grey95"),
                strip.background = element_blank(),
                strip.text.x = element_blank(),
                strip.text.y = element_blank()
) 


plot_labels = data.frame(MolName = "Sleep", Label = "Sleep")
p1 = ggplot() + 
  geom_rect(data = anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="lightyellow") +
  #geom_vline(data = subset(events, grepl("milk",events$Details)), mapping=aes(xintercept=DT), color="gray",size=1) +
  geom_rect(data=sleep, aes(xmin=DT, xmax=DT+seconds, ymin=0, ymax=1, fill=Intensity)) + labs(fill="") +
  scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %b-%d %H:%M", limits=c(min(all_omes$DT),max(all_omes$DT)), timezone = "America/Los_Angeles") + 
  scale_fill_brewer(palette = "Set1") + 
  mytheme + theme(legend.position = "top", legend.box = "horizontal", axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank()) + guides(fill = guide_legend(nrow=1)) +
  facet_wrap(~MolClass, ncol=1, strip.position = "right") +
  geom_label(data = plot_labels, mapping = aes(x=as.POSIXct("2019-04-29 04:00:00", tz="America/Los_Angeles"), y=Inf, label = Label), color="black", size=5, hjust=0, vjust=1, nudge_y=-10, nudge_x=-60) + 
  mytheme

plot_labels = data.frame(MolName = c("Steps","HR","CGM"), Label = c("Steps","HR","CGM"))
p2 = ggplot() + 
  geom_rect(data = anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="lightyellow") +
  geom_line(data = wearables, aes(DT, Intensity)) + 
  facet_wrap(~MolName, ncol=1, scales="free_y", strip.position="right") + 
  scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %b-%d %H:%M", limits=c(min(all_omes$DT),max(all_omes$DT)), timezone = "America/Los_Angeles") + 
  #theme(axis.text.x=element_text(angle=90)) +
  geom_label(data = plot_labels, mapping = aes(x=as.POSIXct("2019-04-29 04:00:00", tz="America/Los_Angeles"), y=Inf, label = Label), color="black", size=5, hjust=0, vjust=1, nudge_y=-10, nudge_x=-60) +
  mytheme 

FoodTimeL_sum$DT = as.POSIXct(FoodTimeL_sum$DT,tz="America/Los_Angeles")
FoodTimeL_sum$MolName = factor(FoodTimeL_sum$MolName, levels = c("Calcium_PDV","Calories","Carbs_g","Fat_g","Alcohol_g","Fiber_g","Iron_PDV","Phosphorus_mg","Potassium_mg","Protein_g","Saturated_fat_g","Sodium_mg","Sugars_g", "Vitamin_A_PDV", "Vitamin_C_PDV"))
plot_labels = data.frame(MolClass = "Food", Label = "Food")
p4 = 
  ggplot() + 
  geom_rect(data = anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="lightyellow") +
  geom_point(data = subset(FoodTimeL_sum, Intensity > 0 & MolName %in% c("Carbs_g","Fat_g","Protein_g","Alcohol_g")), aes(DT,food_multiplier*Intensity, color=MolName), shape="square",size=3,alpha=.75) + 
  #theme_void() + 
  theme(panel.background = element_rect("grey95")) + mytheme +
  labs(color="") + scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %b-%d %H:%M", limits=c(min(all_omes$DT),max(all_omes$DT)), timezone = "America/Los_Angeles") +
  facet_wrap(~MolClass, ncol=1, strip.position = "right") + #+ theme(axis.text = element_blank()) #, axis.line.x = element_blank())
  geom_label(data = plot_labels, mapping = aes(x=as.POSIXct("2019-04-29 04:00:00", tz="America/Los_Angeles"), y=Inf, label = Label), color="black", size=5, hjust=0, vjust=1, nudge_y=-10, nudge_x=-60)


mols = c(1454, 1432, 971, 99, 1418, 1164, 254, 164)
mol_names = subset(all_omes, SampleID==70 & MolNum %in% mols)$MolName
mols = rev(mols)
dat1 = all_omes[all_omes$MolNum %in% mols,]
dat1$MolName = factor(dat1$MolName, levels = rev(as.matrix(unique(dat1[dat1$MolNum %in% mols,"MolName"]))))
plot_labels = data.frame(unique(all_omes[all_omes$MolNum %in% mols,"MolName"]), unique(all_omes[all_omes$MolNum %in% mols,"MolName"]),  unique(all_omes[all_omes$MolNum %in% mols,"MolNum"])) %>% dplyr::rename(Label = MolName.1) 
p5 = ggplot() + 
  geom_rect(data = anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="lightyellow") +
  geom_line(data = dat1, aes(DT,Intensity,color=MolName, group=MolNum)) + 
  geom_point(data = dat1, aes(DT,Intensity,color=MolName, group=MolNum)) + 
  #geom_smooth(data = dat1, aes(DT,Intensity,color=MolName, group=MolNum), span = .2, se=F) + 
  facet_wrap(~MolName,ncol=1,scales="free_y", strip.position = "right") + theme(panel.background = element_rect("grey95")) +
  scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %H:%M", timezone = "America/Los_Angeles", limits=c(min(all_omes$DT),max(all_omes$DT))) +
  theme(axis.text.x = element_text(angle=90), legend.position="none") + labs(color="") + xlab("") + ylab("") + 
  geom_label(data = plot_labels, mapping = aes(x=as.POSIXct("2019-04-29 04:00:00", tz="America/Los_Angeles"), y=Inf, color=MolName, label = Label), size=5, hjust=0, vjust=1, nudge_y=-10, nudge_x=-60) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank()) 

plot_grid(p1,p2,p4, p5, ncol=1, axis="tblr", align="v", rel_heights = c(2,4,2,12))





# EXERCISE ZOOM IN

# all_omes_hr_d = all_omes_hr %>% mutate(Day = day(DT))
# anno_t = as.data.frame(select(subset(ungroup(all_omes_hr_d),MolNum==1), DT, Day, Hr, Hr_day, AM))
# # anno_shading = anno_t %>% group_by(Day) %>% subset(AM) %>% summarize(shading_start = first(DT), shading_end = last(DT))
# anno_shading = anno_t %>% group_by(Day) %>% summarise(shading_start = DT[(Hr_day==6)], shading_end = DT[(Hr_day == 18)])
# events = subset(events, DT < "2019-05-07")


# mytheme = theme(axis.text.x = element_blank(),
#                 axis.ticks.x = element_blank(),
#                 axis.line.x = element_blank(),
#                 axis.title.x = element_blank(),
#                 axis.title.y = element_blank(),
#                 panel.background = element_rect("grey95"),
#                 strip.background = element_blank(),
#                 strip.text.x = element_blank(),
#                 strip.text.y = element_blank()
#                 ) 

d = ymd(20190431)

plot_labels = data.frame(MolName = "Sleep", Label = "Sleep")
p1 = ggplot() + 
  geom_rect(data=subset(sleep,Time > hms::as.hms(5*60*60) & Time < hms::as.hms(11*60*60) & Day == d), aes(xmin=DT, xmax=DT+seconds, ymin=0, ymax=1, fill=Intensity)) + labs(fill="") +
  scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %b-%d %H:%M", limits=c(min(all_omes$DT),max(all_omes$DT)), timezone = "America/Los_Angeles") + 
  scale_fill_brewer(palette = "Set1") + 
  mytheme + theme(legend.position = "top", legend.box = "horizontal", axis.text.y = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank()) + guides(fill = guide_legend(nrow=1)) +
  facet_wrap(~MolClass, ncol=1, strip.position = "right") +
  geom_label(data = plot_labels, mapping = aes(x=as.POSIXct("2019-04-29 04:00:00", tz="America/Los_Angeles"), y=Inf, label = Label), color="black", size=5, hjust=0, vjust=1, nudge_y=-10, nudge_x=-60) + 
  mytheme

plot_labels = data.frame(MolName = c("Steps","HR","CGM"), Label = c("Steps","HR","CGM"))
p2 = ggplot() + 
  geom_line(data=subset(wearables,Time > hms::as.hms(5*60*60) & Time < hms::as.hms(11*60*60) & Day == d), aes(DT, Intensity)) + 
  facet_wrap(~MolName, ncol=1, scales="free_y", strip.position="right")  
#scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %b-%d %H:%M", limits=c(min(all_omes$DT),max(all_omes$DT)), timezone = "America/Los_Angeles") + 
#theme(axis.text.x=element_text(angle=90)) +
#geom_label(data = plot_labels, mapping = aes(x=as.POSIXct("2019-04-29 04:00:00", tz="America/Los_Angeles"), y=Inf, label = Label), color="black", size=5, hjust=0, vjust=1, nudge_y=-10, nudge_x=-60) +
#mytheme 


mols = c(1454, 1432, 971, 99, 1418, 1164, 254, 164)
mol_names = subset(all_omes, SampleID==70 & MolNum %in% mols)$MolName
mols = rev(mols)
dat1 = all_omes[all_omes$MolNum %in% mols,]
dat1 = subset(dat1, tod > hms::as.hms(5*60*60) & tod < hms::as.hms(11*60*60) & day == d)
dat1$MolName = factor(dat1$MolName, levels = rev(as.matrix(unique(dat1[dat1$MolNum %in% mols,"MolName"]))))
plot_labels = data.frame(unique(all_omes[all_omes$MolNum %in% mols,"MolName"]), unique(all_omes[all_omes$MolNum %in% mols,"MolName"]),  unique(all_omes[all_omes$MolNum %in% mols,"MolNum"])) %>% dplyr::rename(Label = MolName.1) 
p5 = ggplot() + 
  geom_line(data = dat1, aes(DT,Intensity,color=MolName, group=MolNum)) + 
  geom_point(data = dat1, aes(DT,Intensity,color=MolName, group=MolNum)) + 
  #geom_smooth(data = dat1, aes(DT,Intensity,color=MolName, group=MolNum), span = .2, se=F) + 
  facet_wrap(~MolName,ncol=1,scales="free_y", strip.position = "right") + theme(panel.background = element_rect("grey95")) 
#scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %H:%M", timezone = "America/Los_Angeles", limits=c(min(all_omes$DT),max(all_omes$DT))) +
#theme(axis.text.x = element_text(angle=90), legend.position="none") + labs(color="") + xlab("") + ylab("") + 
#geom_label(data = plot_labels, mapping = aes(x=as.POSIXct("2019-04-29 04:00:00", tz="America/Los_Angeles"), y=Inf, color=MolName, label = Label), size=5, hjust=0, vjust=1, nudge_y=-10, nudge_x=-60) +
#theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank()) 

plot_grid(p1,p2, p5, ncol=1, axis="tblr", align="v", rel_heights = c(2,4,12))

plot(p1)



# MULTIDAY PATTERNS

# search = c("TNFA","CD40L","PDGFAA","LPE.18.0.","LPE.16.0","PC.18.0.18.3.","TAG54.7.FA16.1","TAG54.7.FA22.6","PC.16.0.16.1.","C8H12N2O3.1","C7H10O6.5","N-(1-Deoxy-1-fructosyl)valine.1","APOL1")
# for (i in 1:length(search)) {
#   print(unique(subset(all_omes, grepl(search[i],MolName),c("MolNum","MolName","MolClass"))))
#   
# }
# 
# search = c("C15H16N2O2","Q14151","C18:2.DC FA","C6H11NO4","Acetylglycine","(CMPF)","C12H14O5 pos","APOA1","AHSP","Sebacic acid","APOA2","AFAM","C14:2.2OH FA","C13H22O4.1","EOTAXIN")
# for (i in 1:length(search)) {
#   print(unique(subset(all_omes, grepl(search[i],MolName),c("MolNum","MolName","MolClass"))))
#   
# }

mols = c(1454, 1432, 971, 99, 1418, 1164, 254, 164)

p1 = c(1454,1436,1432,1113,413,461,971,978,441,83,99)
mol_names = subset(all_omes, SampleID==70 & MolNum %in% p1)$MolName
mols = rev(p1)
dat1 = all_omes[all_omes$MolNum %in% mols,]
dat1$MolName = factor(dat1$MolName, levels = rev(as.matrix(unique(dat1[dat1$MolNum %in% mols,"MolName"]))))
plot_labels = data.frame(unique(all_omes[all_omes$MolNum %in% mols,"MolName"]), unique(all_omes[all_omes$MolNum %in% mols,"MolName"]),  unique(all_omes[all_omes$MolNum %in% mols,"MolNum"])) %>% dplyr::rename(Label = MolName.1) 
#all_omes$MolName = factor(all_omes$MolName, levels = unique(all_omes[all_omes$MolNum %in% mols,"MolName"])) # Hack to get cortisol on top
ggplot() + 
  geom_rect(data = anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="lightyellow") +
  geom_line(data = dat1, aes(DT,Intensity,color=MolName, group=MolNum)) + 
  geom_point(data = dat1, aes(DT,Intensity,color=MolName, group=MolNum)) + 
  #geom_smooth(data = dat1, aes(DT,Intensity,color=MolName, group=MolNum), span = .2, se=F) + 
  facet_wrap(~MolName,ncol=1,scales="free_y", strip.position = "right") + theme(panel.background = element_rect("grey95")) +
  scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %H:%M", timezone = "America/Los_Angeles", limits=c(min(all_omes$DT),max(all_omes$DT))) +
  theme(axis.text.x = element_text(angle=90), legend.position="none") + labs(color="") + xlab("") + ylab("") + 
  geom_label(data = plot_labels, mapping = aes(x=as.POSIXct("2019-04-29 04:00:00", tz="America/Los_Angeles"), y=Inf, color=MolName, label = Label), size=5, hjust=0, vjust=1, nudge_y=-10, nudge_x=-60) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank()) 

p2 = c(188, 1360, 164, 110, 1418, 254, 266, 1164,1393,243,1166,1303,195,201,1418)
mol_names = subset(all_omes, SampleID==70 & MolNum %in% p1)$MolName
mols = rev(p2)
dat1 = all_omes[all_omes$MolNum %in% mols,]
dat1$MolName = factor(dat1$MolName, levels = rev(as.matrix(unique(dat1[dat1$MolNum %in% mols,"MolName"]))))
plot_labels = data.frame(unique(all_omes[all_omes$MolNum %in% mols,"MolName"]), unique(all_omes[all_omes$MolNum %in% mols,"MolName"]),  unique(all_omes[all_omes$MolNum %in% mols,"MolNum"])) %>% dplyr::rename(Label = MolName.1) 
ggplot() + 
  geom_rect(data = anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="lightyellow") +
  geom_line(data = dat1, aes(DT,Intensity,color=MolName, group=MolNum)) + 
  geom_point(data = dat1, aes(DT,Intensity,color=MolName, group=MolNum)) + 
  #geom_smooth(data = dat1, aes(DT,Intensity,color=MolName, group=MolNum), span = .2, se=F) + 
  facet_wrap(~MolName,ncol=1,scales="free_y", strip.position = "right") + theme(panel.background = element_rect("grey95")) +
  scale_x_datetime(breaks=date_breaks("4 hour"), date_labels = "%a %H:%M", timezone = "America/Los_Angeles", limits=c(min(all_omes$DT),max(all_omes$DT))) +
  theme(axis.text.x = element_text(angle=90), legend.position="none") + labs(color="") + xlab("") + ylab("") + 
  geom_label(data = plot_labels, mapping = aes(x=as.POSIXct("2019-04-29 04:00:00", tz="America/Los_Angeles"), y=Inf, color=MolName, label = Label), size=5, hjust=0, vjust=1, nudge_y=-10, nudge_x=-60) +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), strip.text.y = element_blank()) 












TimeAnno$datef = format(TimeAnno$date_time,"%a %b %d")
TimeAnno$datef = factor(TimeAnno$datef,levels = unique(TimeAnno$datef))
TimeAnno$timef = as.POSIXct(format(TimeAnno$date_time, "2000-01-01 %H:%M:%S"))

write.csv(TimeAnno,"TimeAnno.csv", row.names=F)


library(scales)

#TimeAnno = TimeAnno[!is.na(TimeAnno$c("activity","food","datef","timef")]

events = gather(TimeAnno,"Type","Details",7:8)
events = events[!is.na(events$Details),]
days = unique(dataL$datef)


sal = 229
trend = c(133, 211, 212, 214,239,243,73,9)  
metab$Time = hms::as.hms(metab$DT, tz="America/Los_Angeles")
metab$Date = as_date(metab$DT, tz="America/Los_Angeles")

a = subset(metab, metab$met_num %in% sal)
ggplot(a, aes(x=Time, y = Intensity, group=1)) + geom_line()  + 
  #geom_text_repel(data = subset(events, Type=="activity"), aes(x=date_time,y=approx(a$DT,a$Intensity,date_time)$y,label=Details)) + 
  facet_wrap(~Date, ncol=2)


geom_rect(aes(xmin=as.POSIXct("2019-04-29 00:00:00", tz="America/Los_Angeles"), xmax=as.POSIXct("2019-04-29 12:00:00", tz="America/Los_Angeles"), ymin=-Inf, ymax=Inf), color=NA, fill="lightyellow") + 
  geom_rect(aes(xmin=as.POSIXct("2019-4-30 00:00:00", tz="America/Los_Angeles"), xmax=as.POSIXct("2019-04-30 12:00:00", tz="America/Los_Angeles"), ymin=-Inf, ymax=Inf), color=NA, fill="lightyellow") + 
  geom_rect(aes(xmin=as.POSIXct("2019-05-01 00:00:00", tz="America/Los_Angeles"), xmax=as.POSIXct("2019-05-01 12:00:00", tz="America/Los_Angeles"), ymin=-Inf, ymax=Inf), color=NA, fill="lightyellow") + 
  geom_rect(aes(xmin=as.POSIXct("2019-05-02 00:00:00", tz="America/Los_Angeles"), xmax=as.POSIXct("2019-05-02 12:00:00", tz="America/Los_Angeles"), ymin=-Inf, ymax=Inf), color=NA, fill="lightyellow") + 
  geom_rect(aes(xmin=as.POSIXct("2019-05-03 00:00:00", tz="America/Los_Angeles"), xmax=as.POSIXct("2019-05-03 12:00:00", tz="America/Los_Angeles"), ymin=-Inf, ymax=Inf), color=NA, fill="lightyellow") + 
  geom_rect(aes(xmin=as.POSIXct("2019-05-04 00:00:00", tz="America/Los_Angeles"), xmax=as.POSIXct("2019-05-04 12:00:00", tz="America/Los_Angeles"), ymin=-Inf, ymax=Inf), color=NA, fill="lightyellow") + 
  geom_rect(aes(xmin=as.POSIXct("2019-05-05 00:00:00", tz="America/Los_Angeles"), xmax=as.POSIXct("2019-05-05 12:00:00", tz="America/Los_Angeles"), ymin=-Inf, ymax=Inf), color=NA, fill="lightyellow") + 
  geom_rect(aes(xmin=as.POSIXct("2019-05-06 05:00:00", tz="America/Los_Angeles"), xmax=as.POSIXct("2019-05-06 12:00:00", tz="America/Los_Angeles"), ymin=-Inf, ymax=Inf), color=NA, fill="lightyellow") + 
  geom_rect(aes(xmin=as.POSIXct("2019-05-07 00:00:00", tz="America/Los_Angeles"), xmax=as.POSIXct("2019-05-07 12:00:00", tz="America/Los_Angeles"), ymin=-Inf, ymax=Inf), color=NA, fill="lightyellow") +
  geom_point() + geom_line() + facet_wrap(~met_ID, ncol=2, scales="free")

ggsave(paste("test.png",sep=""), width=30, height=15, units="in")





load("~/Box/Microsampling/Analysis/Dan_2019-08-17/MOmics_01.RData")


i = 1
for (i in 1:100) {
  
  ggplot(DF247_lipids2, aes(x=CollectionTime_format, y = log_sample_nmol_per_g_concentration, group=LipidSpecies)) + 
    geom_point() + geom_line() +
    facet_wrap_paginate(~LipidSpecies, scales="free", ncol=3, nrow=8, page=i)
  
  ggsave(paste("Lipids_",i,".png",sep=""), width=22, height=11, units="in")
  
  print(i)
  
}




library(lubridate)
grid = seq(as.POSIXct(floor_date(min(all_omes$DT),"hour")),as.POSIXct(max(all_omes$DT)),by="1 hour")
#all_omes %>% group_by(MolName,MolClass,MolSubclass) %>% mutate(Intensity_z = (Intensity - mean(Intensity)/sd(Intensity)))

unique(all_omes_hr$MolClass)

all_omes_hr = all_omes %>% 
  group_by(MolName,MolClass,MolSubclass,MolNum) %>% 
  mutate(Intensity_z = (Intensity - mean(Intensity))/sd(Intensity)) %>% 
  nest() %>% 
  mutate(hr = map(data, function(x) as_tibble(approx(x$DT, x$Intensity_z, xout=grid, rule=1)))) %>%
  unnest(hr) %>% dplyr::rename(DT = x, Intensity = y)


all_omes_hr = all_omes_hr[!is.na(all_omes_hr$Intensity),]
all_omes_hr = all_omes_hr %>% group_by(MolName) %>% mutate(Hr = as.factor(as.numeric(difftime(DT, first(DT)), units="hours")), Hr_day = hour(DT), AM = am(DT))
all_omes_hr = all_omes_hr[all_omes_hr$DT < as_date("2019-05-07"),]

all_omes_w = all_omes_hr %>% ungroup() %>% select(-DT) %>% spread(Hr,Intensity)

d = select(all_omes_w, -c("MolName","MolClass","MolSubclass","MolNum","Hr_day","AM"))
m = select(all_omes_w, c("MolName","MolClass","MolSubclass","MolNum","Hr_day","AM"))

library(dtwclust)
pc <- tsclust(d, type = "partitional", k = 20L, # preproc = zscore, 
              distance = "dtw_basic", centroid = "pam", 
              seed = 3247L, trace = TRUE,
              args = tsclust_args(dist = list(window.size = 20L)))
plot(pc)

all_omes_w$Cluster = pc@cluster

pdat = gather(all_omes_w,"Hr", "Intensity", -c("MolName","MolClass","MolSubclass","MolNum","Cluster"))

#pdat_med = pdat %>% group_by(Cluster, Hr) %>% summarize(Int_med = median(Intensity))
ggplot(subset(pdat,MolClass=="Lipid" & MolSubclass=="TAG"), aes(as.numeric(Hr), Intensity, color=MolClass, group=MolNum)) + geom_line() + #geom_smooth(span=0.01,color="black",se=F, aes(group=Cluster)) + 
  facet_wrap(~Cluster) + theme(legend.position="bottom") 

ggplot(subset(pdat,MolName=="Cortisol" | MolName=="GLP1"), aes(as.numeric(Hr), Intensity, color=MolName, group=MolNum)) + geom_line() + #geom_smooth(span=0.01,color="black",se=F, aes(group=Cluster)) + 
  facet_wrap(~Cluster) + theme(legend.position="bottom") 







anno_t = as.data.frame(select(subset(ungroup(all_omes_hr),MolNum==1), Hr, Hr_day, AM))
rownames(anno_t) = anno_t$Hr
anno_t$Hr = anno_t$Hr>12
#anno_t = anno_t[,-1]

tt = t(anno_t)


s = subset(all_omes_hr,MolClass=="MetabolicPanel" & !grepl("CHEX",all_omes_hr$MolName)) %>% ungroup() %>% select(-DT, -Hr_day, -AM) %>% spread(Hr,Intensity)

#all_omes_w = s %>% ungroup() %>% select(-DT) %>% spread(Hr,Intensity)

d = select(s, -c("MolName","MolClass","MolSubclass","MolNum"))
m = select(s, c("MolName","MolClass","MolSubclass","MolNum"))
d = d[,order(as.numeric(colnames(d)))]

dd=as.matrix(d)
df = as.data.frame(d)
rownames(df) = m$MolName 

pheatmap(df, cluster_rows = T, cluster_cols = F, annotations_col = anno_t) #, annotation_row = m) #labels_col = colnames(d), annotation_col = t(tt["AM",]))

ha1 = HeatmapAnnotation(Hr = anno_t$Hr_day)
Heatmap(df, cluster_columns = F, top_annotation=ha1, show_column_names = F)

ggplot(subset(all_omes, MolClass=="MetabolicPanel"), aes(DT,Intensity)) + geom_point() + geom_line() + geom_vline()
facet_wrap(~MolName, ncol=1, scales="free")






anno_t = as.data.frame(select(subset(ungroup(all_omes_hr),MolNum==1), Hr, Hr_day, AM))
rownames(anno_t) = anno_t$Hr
anno_t$Hr = anno_t$Hr>12
#anno_t = anno_t[,-1]

tt = t(anno_t)

s = subset(all_omes_hr,MolClass=="MetabolicPanel" & !grepl("CHEX",all_omes_hr$MolName)) %>% ungroup() %>% select(-DT, -Hr_day, -AM) %>% spread(Hr,Intensity)

s = subset(all_omes_hr,MolClass=="Lipid") %>% ungroup() %>% select(-DT,-Hr_day, -AM) %>% spread(Hr,Intensity)
d = select(s, -c("MolName","MolClass","MolSubclass","MolNum"))
m = select(s, c("MolName","MolClass","MolSubclass","MolNum"))
d = d[,order(as.numeric(colnames(d)))] 
df = as.data.frame(d)

ha1 = HeatmapAnnotation(Hr = anno_t$Hr_day)
p4=Heatmap(df, cluster_columns = F, top_annotation=ha1, show_column_names = F)




all_omes_hr_d = all_omes_hr %>% mutate(Day = day(DT))

anno_t = as.data.frame(select(subset(ungroup(all_omes_hr_d),MolNum==1), DT, Day, Hr, Hr_day, AM))
anno_shading = anno_t %>% group_by(Day) %>% subset(AM) %>% summarize(shading_start = first(DT), shading_end = last(DT))



circ = c(62,81,116,149,151,157,170,174)
s = subset(all_omes_hr, MolSubclass=="TAG") %>% ungroup() %>% select(-DT, -Hr_day, -AM) %>% spread(Hr,Intensity)
d = select(s, -c("MolName","MolClass","MolSubclass","MolNum"))
m = select(s, c("MolName","MolClass","MolSubclass","MolNum"))
d = d[,order(as.numeric(colnames(d)))] 
df = as.data.frame(d)

ha1 = HeatmapAnnotation(Hr = anno_t$AM)
Heatmap(df, cluster_columns = F, top_annotation=ha1, show_row_names=T)

superheat(df)



# plot with AM/PM shading and food/activity annotation 


sel = subset(all_omes_hr,MolClass=="Lipid") #%>% ungroup() %>% select(-DT, -Hr_day, -AM) %>% spread(Hr,Intensity)

i = 1
for (i in 1:100) {
  
  ggplot() + geom_rect(anno_shading, mapping = aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), color=NA, fill="lightyellow") +
    geom_line(sel, mapping = aes(x=DT, y = Intensity, group=MolNum)) + geom_point(sel, mapping = aes(x=DT, y = Intensity, group=MolNum)) +
    geom_vline(events, mapping=aes(xintercept = date_time, color=Type)) + 
    facet_wrap_paginate(~MolNum, scales="free", ncol=3, nrow=8, page=i) 
  
  ggsave(paste("Lipids_MolNum_",i,".png",sep=""), width=22, height=11, units="in")
  
  print(i)
  
}

ggplot(anno_shading, aes(DT, )) + geom_point() geom_vline(aes(xintercept=shading_start))
#geom_rect(aes(xmin=shading_start, xmax=shading_end, ymin=-Inf, ymax=Inf), fill="black")

dfr = data.frame(r1 = c(1,5), r2 = c(2,6))
dft = data.frame(xx=1:10, yy=21:30)

ggplot() + geom_line(dft, mapping=aes(xx,yy)) + geom_rect(dfr, mapping = aes(xmin=r1, xmax=r2, ymin=-Inf, ymax=Inf))
ggplot(dft,aes(x=x,y=y, group=1)) + geom_line() + geom_vline(dfr, mapping = aes(xintercept=r1))



