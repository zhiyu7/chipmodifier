gc()
rm(list=ls())

library(dplyr)
library(tidyverse)
library(data.table)
library(ukbtools) 
library(survival)
library(survminer)
library(ggrepel)
library(gridExtra)
library(plinkQC)
library(tableone)
library(ggforestplot)

###################################################################
############################read in data############################
#############################################################################

############################read in score by Yunfeng############################
rs=fread("/medpop/esp2/zyu/chip_protemoics/data/eQTL_score/20220910_UKB500k_eachracebest_FHSMESA_Yunfengversion_PRSCSandPT_ukbafreq_inflammatorygenes.txt.gz")
rs=rs[!duplicated(rs)]

colnames(rs) = gsub("PRS", "score", colnames(rs))

##########################read in pheno by Maryam##############################
pheno = fread("/medpop/esp2/mzekavat/UKBB/ukbb_PhenoFile.ALL_500k.UpdatedIncdPhenos_202020.txt")
pheno = data.frame(pheno)
pheno$IDs_toRemove_SampleQCv3 = ifelse((!(pheno$Submitted_Gender == pheno$Inferred_Gender) |  pheno$Non_Consented== 1),1,0) ## N:
pheno = pheno[which(pheno$IDs_toRemove_SampleQCv3==0),]

######################read in chip curated by Mesbah###################################
newchip_overall=get(load("/medpop/esp2/mesbah/projects/Meta_GWAS/MetaGWAS_650k/chip_call/ukb450k_mgbb53k.CHIP_vars.for_phewas.rda"))
newchip_overall[which(newchip_overall$CHIP==1),][is.na(newchip_overall[which(newchip_overall$CHIP==1),])] <- 0

newchip_vaf=fread("/medpop/esp2/zyu/chip_protemoics/data/200K_CHIP/UKB200k_CHIP_Variants.11Aug2021.csv")
newchip_250k_vaf=fread("/medpop/esp2/mesbah/datasets/CHIP/UKBB/450k/chipVars_CV_AB.ukb250k_all.csv.gz")
newchip_250k_vaf=newchip_250k_vaf%>%rename(eid_7089=Broad_ID)

newchip_vaf_CHIP = newchip_vaf%>% group_by(eid_7089) %>% summarise(CHIP_VAF=sum(AF)) %>% mutate(CHIP_VAF=10*CHIP_VAF)
newchip_250k_vaf_CHIP = newchip_250k_vaf%>% group_by(eid_7089) %>% summarise(CHIP_VAF=sum(AF)) %>% mutate(CHIP_VAF=10*CHIP_VAF)
newchip_vaf_CHIP=rbind(newchip_vaf_CHIP, newchip_250k_vaf_CHIP)
newchip_overall=merge(newchip_overall, newchip_vaf_CHIP, by="eid_7089", all.x=T)
newchip_overall=newchip_overall%>%mutate(CHIP_VAF=ifelse(CHIP==0, 0, ifelse(CHIP==1, CHIP_VAF, NA)))
newchip_vaf_JAK2=newchip_vaf%>%filter(Gene=="JAK2")
newchip_vaf_JAK2=newchip_vaf_JAK2%>%select(eid_7089, AF) %>% group_by(eid_7089) %>% summarise(JAK2_VAF=sum(AF))%>%mutate(JAK2_VAF=10*JAK2_VAF)
newchip_250k_vaf_JAK2=newchip_250k_vaf%>%filter(Gene.refGene=="JAK2")
newchip_250k_vaf_JAK2=newchip_250k_vaf_JAK2%>%select(eid_7089, AF) %>% group_by(eid_7089) %>% summarise(JAK2_VAF=sum(AF))%>%mutate(JAK2_VAF=10*JAK2_VAF)
newchip_vaf_JAK2=rbind(newchip_vaf_JAK2, newchip_250k_vaf_JAK2)
newchip_overall=merge(newchip_overall, newchip_vaf_JAK2, by="eid_7089", all.x=T)
newchip_overall=newchip_overall%>%mutate(JAK2_VAF=ifelse(JAK2==0, 0, ifelse(JAK2==1, JAK2_VAF, NA)))
newchip_vaf_DNMT3A=newchip_vaf%>%filter(Gene=="DNMT3A")
newchip_vaf_DNMT3A=newchip_vaf_DNMT3A%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(DNMT3A_VAF=sum(AF))%>%mutate(DNMT3A_VAF=10*DNMT3A_VAF)
newchip_250k_vaf_DNMT3A=newchip_250k_vaf%>%filter(Gene.refGene=="DNMT3A")
newchip_250k_vaf_DNMT3A=newchip_250k_vaf_DNMT3A%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(DNMT3A_VAF=sum(AF))%>%mutate(DNMT3A_VAF=10*DNMT3A_VAF)
newchip_vaf_DNMT3A=rbind(newchip_vaf_DNMT3A, newchip_250k_vaf_DNMT3A)
newchip_overall=merge(newchip_overall, newchip_vaf_DNMT3A, by="eid_7089", all.x=T)
newchip_overall=newchip_overall%>%mutate(DNMT3A_VAF=ifelse(DNMT3A==0, 0, ifelse(DNMT3A==1, DNMT3A_VAF, NA)))
newchip_vaf_TET2=newchip_vaf%>%filter(Gene=="TET2")
newchip_vaf_TET2=newchip_vaf_TET2%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(TET2_VAF=sum(AF))%>% mutate(TET2_VAF=10*TET2_VAF)
newchip_250k_vaf_TET2=newchip_250k_vaf%>%filter(Gene.refGene=="TET2")
newchip_250k_vaf_TET2=newchip_250k_vaf_TET2%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(TET2_VAF=sum(AF))%>% mutate(TET2_VAF=10*TET2_VAF)
newchip_vaf_TET2=rbind(newchip_vaf_TET2, newchip_250k_vaf_TET2)
newchip_overall=merge(newchip_overall, newchip_vaf_TET2, by="eid_7089", all.x=T)
newchip_overall=newchip_overall%>%mutate(TET2_VAF=ifelse(TET2==0, 0, ifelse(TET2==1, TET2_VAF, NA)))
newchip_vaf_ASXL1=newchip_vaf%>%filter(Gene=="ASXL1")
newchip_vaf_ASXL1=newchip_vaf_ASXL1%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(ASXL1_VAF=sum(AF))%>%mutate(ASXL1_VAF=10*ASXL1_VAF)
newchip_250k_vaf_ASXL1=newchip_250k_vaf%>%filter(Gene.refGene=="ASXL1")
newchip_250k_vaf_ASXL1=newchip_250k_vaf_ASXL1%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(ASXL1_VAF=sum(AF))%>%mutate(ASXL1_VAF=10*ASXL1_VAF)
newchip_vaf_ASXL1=rbind(newchip_vaf_ASXL1, newchip_250k_vaf_ASXL1)
newchip_overall=merge(newchip_overall, newchip_vaf_ASXL1, by="eid_7089", all.x=T)
newchip_overall=newchip_overall%>%mutate(ASXL1_VAF=ifelse(ASXL1==0, 0, ifelse(ASXL1==1, ASXL1_VAF, NA)))
newchip_vaf_TP53=newchip_vaf%>%filter(Gene=="TP53")
newchip_vaf_TP53=newchip_vaf_TP53%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(TP53_VAF=sum(AF))%>%mutate(TP53_VAF=10*TP53_VAF)
newchip_250k_vaf_TP53=newchip_250k_vaf%>%filter(Gene.refGene=="TP53")
newchip_250k_vaf_TP53=newchip_250k_vaf_TP53%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(TP53_VAF=sum(AF))%>%mutate(TP53_VAF=10*TP53_VAF)
newchip_vaf_TP53=rbind(newchip_vaf_TP53, newchip_250k_vaf_TP53)
newchip_overall=merge(newchip_overall, newchip_vaf_TP53, by="eid_7089", all.x=T)
newchip_overall=newchip_overall%>%mutate(TP53_VAF=ifelse(TP53==0, 0, ifelse(TP53==1, TP53_VAF, NA)))
newchip_vaf_PPM1D=newchip_vaf%>%filter(Gene=="PPM1D")
newchip_vaf_PPM1D=newchip_vaf_PPM1D%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(PPM1D_VAF=sum(AF))%>%mutate(PPM1D_VAF=10*PPM1D_VAF)
newchip_250k_vaf_PPM1D=newchip_250k_vaf%>%filter(Gene.refGene=="PPM1D")
newchip_250k_vaf_PPM1D=newchip_250k_vaf_PPM1D%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(PPM1D_VAF=sum(AF))%>%mutate(PPM1D_VAF=10*PPM1D_VAF)
newchip_vaf_PPM1D=rbind(newchip_vaf_PPM1D, newchip_250k_vaf_PPM1D)
newchip_overall=merge(newchip_overall, newchip_vaf_PPM1D, by="eid_7089", all.x=T)
newchip_overall=newchip_overall%>%mutate(PPM1D_VAF=ifelse(PPM1D==0, 0, ifelse(PPM1D==1, PPM1D_VAF, NA)))
newchip_vaf_SF3B1=newchip_vaf%>%filter(Gene=="SF3B1")
newchip_vaf_SF3B1=newchip_vaf_SF3B1%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(SF3B1_VAF=sum(AF))%>%mutate(SF3B1_VAF=10*SF3B1_VAF)
newchip_250k_vaf_SF3B1=newchip_250k_vaf%>%filter(Gene.refGene=="SF3B1")
newchip_250k_vaf_SF3B1=newchip_250k_vaf_SF3B1%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(SF3B1_VAF=sum(AF))%>%mutate(SF3B1_VAF=10*SF3B1_VAF)
newchip_vaf_SF3B1=rbind(newchip_vaf_SF3B1, newchip_250k_vaf_SF3B1)
newchip_overall=merge(newchip_overall, newchip_vaf_SF3B1, by="eid_7089", all.x=T)
newchip_overall=newchip_overall%>%mutate(SF3B1_VAF=ifelse(SF3B1==0, 0, ifelse(SF3B1==1, SF3B1_VAF, NA)))
newchip_vaf_SRSF2=newchip_vaf%>%filter(Gene=="SRSF2")
newchip_vaf_SRSF2=newchip_vaf_SRSF2%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(SRSF2_VAF=sum(AF))%>%mutate(SRSF2_VAF=10*SRSF2_VAF)
newchip_250k_vaf_SRSF2=newchip_250k_vaf%>%filter(Gene.refGene=="SRSF2")
newchip_250k_vaf_SRSF2=newchip_250k_vaf_SRSF2%>%select(eid_7089, AF)%>% group_by(eid_7089) %>% summarise(SRSF2_VAF=sum(AF))%>%mutate(SRSF2_VAF=10*SRSF2_VAF)
newchip_vaf_SRSF2=rbind(newchip_vaf_SRSF2, newchip_250k_vaf_SRSF2)
newchip_overall=merge(newchip_overall, newchip_vaf_SRSF2, by="eid_7089", all.x=T)
newchip_overall=newchip_overall%>%mutate(SRSF2_VAF=ifelse(SRSF2==0, 0, ifelse(SRSF2==1, SRSF2_VAF, NA)))

newchip=merge(newchip_overall,rs,all.x=T, all.y=T, by.x="eid_7089", by.y="id")
newchip=newchip%>%rename(eid=eid_7089)%>%select(-SEX)

pheno = pheno%>%select(id, in_white_British_ancestry_subset, NewCHIP,IL6R_rs2228145_C, TwoHundredk_hasCHIP, age, age2, Sex, contains("smok"), contains("Coronary_Artery_Disease"), PC1:PC10, 
genotyping_array, contains("composite_mi_cad_stroke"), contains("heme"), BMI, contains("Diabetes_Type_2"),  contains("Hypercholesterolemia"), contains("Hypertension"),
 White.blood.cell..leukocyte..count, Red.blood.cell..erythrocyte..count, Haemoglobin.concentration, Haematocrit.percentage, Mean.corpuscular.volume, 
Mean.corpuscular.haemoglobin, Mean.corpuscular.haemoglobin.concentration, Red.blood.cell..erythrocyte..distribution.width, Platelet.count, Platelet.crit, 
Mean.platelet..thrombocyte..volume, Platelet.distribution.width, Lymphocyte.count, Monocyte.count, Neutrophill.count, Eosinophill.count, Basophill.count, 
Nucleated.red.blood.cell.count, Lymphocyte.percentage, Monocyte.percentage, Neutrophill.percentage, Eosinophill.percentage, Basophill.percentage, 
Nucleated.red.blood.cell.percentage, Reticulocyte.percentage, Reticulocyte.count, Mean.reticulocyte.volume, Mean.sphered.cell.volume, Immature.reticulocyte.fraction, 
High.light.scatter.reticulocyte.percentage, High.light.scatter.reticulocyte.count, White.blood.cell..leukocyte..count_NoOutlierBySex, Red.blood.cell..erythrocyte..count_NoOutlierBySex, 
Haemoglobin.concentration_NoOutlierBySex, Haematocrit.percentage_NoOutlierBySex, Mean.corpuscular.volume_NoOutlierBySex, Mean.corpuscular.haemoglobin_NoOutlierBySex, 
Mean.corpuscular.haemoglobin.concentration_NoOutlierBySex, Red.blood.cell..erythrocyte..distribution.width_NoOutlierBySex, Platelet.count_NoOutlierBySex, Platelet.crit_NoOutlierBySex, 
Mean.platelet..thrombocyte..volume_NoOutlierBySex, Platelet.distribution.width_NoOutlierBySex, Lymphocyte.count_NoOutlierBySex, Monocyte.count_NoOutlierBySex, Neutrophill.count_NoOutlierBySex, 
Eosinophill.count_NoOutlierBySex, Basophill.count_NoOutlierBySex, Nucleated.red.blood.cell.count_NoOutlierBySex, Lymphocyte.percentage_NoOutlierBySex, Monocyte.percentage_NoOutlierBySex, 
Neutrophill.percentage_NoOutlierBySex, Eosinophill.percentage_NoOutlierBySex, Basophill.percentage_NoOutlierBySex, Nucleated.red.blood.cell.percentage_NoOutlierBySex, 
Reticulocyte.percentage_NoOutlierBySex, Reticulocyte.count_NoOutlierBySex, Mean.reticulocyte.volume_NoOutlierBySex, Mean.sphered.cell.volume_NoOutlierBySex, Immature.reticulocyte.fraction_NoOutlierBySex, 
High.light.scatter.reticulocyte.percentage_NoOutlierBySex, High.light.scatter.reticulocyte.count_NoOutlierBySex, contains("C.reactive.protein"), contains("SBP"), contains("DBP"), contains("C.reactive.protein"),
Albumin.x:Urate, fvc, fev1, fev1_fvc_ratio, Total.Cholesterol, HDL.Cholesterol,Total.Triglycerides,LDL.Cholesterol)

pheno_newchip=merge(pheno,newchip,by.x="id",by.y="eid")
pheno_newchip$FiftyK_indicator=ifelse(is.na(pheno_newchip$NewCHIP), 0, 1)
pheno_newchip=pheno_newchip[!is.na(pheno_newchip$AIM2_score),]
pheno_newchip=pheno_newchip[!is.na(pheno_newchip$CHIP),]
#pheno_newchip=pheno_newchip[which(pheno_newchip$CHIP_Batch=="UKB200k"),]

#bd_related <- fread("/medpop/esp2/pradeep/UKBiobank/v3data/ukb7089_rel_s488363.dat")
#bd_related_try=bd_related[which(bd_related$ID1 %in% pheno_newchip$id & bd_related$ID2 %in% pheno_newchip$id),]
#bd_related_try <- relatednessFilter(bd_related_try, otherCriterion = NULL, relatednessTh = 0.0884, otherCriterionTh = NULL, relatednessIID1 = "ID1", relatednessIID2 = "ID2",
#                                relatednessFID1 = NULL,
#                                relatednessFID2 = NULL,
#                                relatednessRelatedness = "Kinship",
#                                #otherCriterionIID = "IID",
#                                otherCriterionMeasure = NULL,
#                                verbose = FALSE)
#bd_related_excluded <- bd_related_try$failIDs
#write.table(bd_related_excluded,file="/medpop/esp2/zyu/chip_protemoics/listofremovefor450K.txt", sep="\t", quote=FALSE, row.names=FALSE)
bd_related_excluded=fread("/medpop/esp2/zyu/chip_protemoics/listofremovefor450K.txt")
pheno_newchip_unrelated=pheno_newchip[!pheno_newchip$id %in% bd_related_excluded$IID,]
#pheno_newchip_unrelated=merge(pheno_newchip_unrelated,rs,by.x="id",by.y="FID")
pheno_newchip=pheno_newchip_unrelated

######code IL6R SNP
pheno_newchip$rs2228145_C_bi=ifelse(pheno_newchip$IL6R_rs2228145_C==0, 0, 1)

#########################################################
########################Table 1#############################
###############################################################
factorVars <- c("DNMT3A", "TET2", "ASXL1", "JAK2", "TP53", "PPM1D", "SF3B1", "SRSF2", "DTA", "DDR", "Splice", "in_white_British_ancestry_subset", "Sex", "ever_smoked", 
"Prev_Diabetes_Type_2", "Prev_Coronary_Artery_Disease", "Prev_Hypercholesterolemia","Prev_Hypertension","expandedCHIP", "expandedDNMT3A", "expandedTET2", 
"expandedASXL1", "expandedJAK2", "expandedPPM1D", "expandedSF3B1", "expandedTP53", "expandedSRSF2", "expandedDTA", "expandedDDR","expandedSplice")

vars <- c("DNMT3A", "TET2", "ASXL1", "JAK2", "TP53","PPM1D", "SF3B1", "SRSF2", "DTA", "DDR", "Splice", "AGE_assessment", "in_white_British_ancestry_subset", 
"age", "age2", "Sex", "ever_smoked","BMI","Prev_Diabetes_Type_2", "Prev_Coronary_Artery_Disease", "Prev_Hypercholesterolemia","Prev_Hypertension","expandedCHIP", 
"expandedDNMT3A", "expandedTET2", "expandedASXL1", "expandedJAK2", "expandedPPM1D", "expandedSF3B1", "expandedTP53", "expandedSRSF2", "expandedDTA", "expandedDDR",
"expandedSplice")

tableOne <- CreateTableOne(vars = vars, data = pheno_newchip[which(pheno_newchip$Prev_Coronary_Artery_Disease==0 & pheno_newchip$Prev_AnyHemeCa==0),], strata="CHIP", factorVars = factorVars)
write.table(print(tableOne),"/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_ukb450K_tableone_unrelated_bychip.txt", sep="\t", quote=FALSE, row.names=TRUE)

tableOne <- CreateTableOne(vars = vars, data = pheno_newchip[which(pheno_newchip$Prev_Coronary_Artery_Disease==0 & pheno_newchip$Prev_AnyHemeCa==0),], factorVars = factorVars)
write.table(print(tableOne),"/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_ukb450K_tableone_unrelated.txt", sep="\t", quote=FALSE, row.names=TRUE)

table(newchip_vaf[which(newchip_vaf$Gene=="JAK2"&newchip_vaf$eid_7089 %in% pheno_newchip[which(pheno_newchip$Prev_Coronary_Artery_Disease==0 & pheno_newchip$Prev_AnyHemeCa==0),]$id),]$Protein_Change)
table(newchip_250k_vaf[which(newchip_250k_vaf$Gene.refGene=="JAK2"& newchip_250k_vaf$eid_7089 %in% pheno_newchip[which(pheno_newchip$Prev_Coronary_Artery_Disease==0 & pheno_newchip$Prev_AnyHemeCa==0),]$id),]$NonsynOI)

#############################################################################
#########################Modeling###################################
###############################################################################

#########################Numbers for text#####################################
table(pheno_newchip[which(pheno_newchip$Prev_composite_mi_cad_stroke_death==0 & pheno_newchip$Prev_AnyHemeCa==0),]$Incd_composite_mi_cad_stroke_death)
#     0      1
#372608  44962
table(pheno_newchip[which(pheno_newchip$Prev_composite_mi_cad_stroke_death==0 ),]$Incd_composite_mi_cad_stroke_death)
#     0      1
#373781  45450

summary(pheno_newchip[which(pheno_newchip$Prev_composite_mi_cad_stroke_death==0 & pheno_newchip$Prev_AnyHemeCa==0),]$FollowUp_composite_mi_cad_stroke_death)
#median 11.030555

#######################Correlations################################
correlation_test <- function(df, target_var, test_vars, method=c("pearson", "spearman")) {
  results <- data.frame()
  
  for (var in test_vars) {
    for (m in method) {
      test_result <- cor.test(df[[target_var]], df[[var]], method=m)
      
      results <- rbind(results, data.frame(
        Target = target_var,
        Variable = var,
        Method = m,
        Correlation = test_result$estimate,
        P.Value = test_result$p.value
      ))
    }
  }
  
  return(results)
}

# Subsetting the dataframe
df <- pheno_newchip[which(pheno_newchip$ASXL1==1),]

# Perform correlation tests
results <- correlation_test(df, 'ASXL1_VAF', TRS, method=c("pearson", "spearman"))

# Print results
print(results)

###############no interaction##################################

chip_list=c("CHIP", "DNMT3A", "TET2", "ASXL1", "JAK2", "PPM1D", "TP53", "SF3B1", "SRSF2")
chip_vaf_list=c("CHIP_VAF", "DNMT3A_VAF", "TET2_VAF", "ASXL1_VAF", "JAK2_VAF", "TP53_VAF", "PPM1D_VAF", "SF3B1_VAF", "SRSF2_VAF")
chip_large_list=c("expandedCHIP", "expandedDNMT3A", "expandedTET2", "expandedASXL1", "expandedJAK2","expandedPPM1D","expandedTP53", "expandedSF3B1", "expandedSRSF2", "expandedDTA", "expandedDDR", "expandedSplice")

outcome_list=c("composite_mi_cad_stroke_death", "Coronary_Artery_Disease_SOFT", "Coronary_Artery_Disease_INTERMEDIATE", "Coronary_Artery_Disease_HARD")

chip_fulllist=c(chip_list, chip_vaf_list, chip_large_list)

h<-data.frame(chip=character(), outcome=character(), n=double(), hr=double(), se=double(), z=double(), pval=double(), lci=double(), uci=double(), stringsAsFactors=FALSE)      

for(m in 1:length(outcome_list)) {
	for(i in 1:length(chip_fulllist)) {
      		h[i+(m-1)*length(chip_fulllist),"chip"]<-chip_fulllist[i]
		h[i+(m-1)*length(chip_fulllist),"outcome"]<-outcome_list[m]
		Incd_outcome<-paste("Incd_", outcome_list[m], sep="")
      		fmla<-as.formula(paste("Surv(as.numeric(FollowUp_", outcome_list[m], "),Incd_", outcome_list[m], ") ~ ", chip_fulllist[i],"+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+Prev_Diabetes_Type_2+ever_smoke",sep=""))
      		if(!is.na(Incd_outcome)){
			f<-coxph(fmla,data=pheno_newchip[which(pheno_newchip$Prev_AnyHemeCa==0& pheno_newchip$Prev_composite_mi_cad_stroke_death==0),])
      			h[i+(m-1)*length(chip_fulllist),"n"]<-summary(f)$n
      			h[i+(m-1)*length(chip_fulllist),"hr"]<-summary(f)$coefficient[1,2]
      			h[i+(m-1)*length(chip_fulllist),"se"]<-summary(f)$coefficient[1,3]
			h[i+(m-1)*length(chip_fulllist),"z"]<-summary(f)$coefficient[1,4]
      			h[i+(m-1)*length(chip_fulllist),"pval"]<-summary(f)$coefficient[1,5]
      			h[i+(m-1)*length(chip_fulllist),"lci"]<-exp(confint(f)[1,1])
      			h[i+(m-1)*length(chip_fulllist),"uci"]<-exp(confint(f)[1,2])
		} 	
	}
}
write.table(h, "/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_Updated1030_CHIP_simplelargevaf_CADoutcomes_nointeraction_adjustforsmoke.txt", sep="\t", quote=FALSE, row.names=FALSE)	


#######################chip association with score############################
h<-data.frame(chip=character(), outcome=character(), n=double(), beta=double(), se=double(), z=double(), pval=double(), lci=double(), uci=double(), stringsAsFactors=FALSE)      

TRS=c("AIM2_score", "CASP1_score", "CASP5_score", "IFNGR2_score","IL10_score", "IL18BP_score", "IL18RAP_score", "IL1B_score",
"IL1R1_score", "IL1R2_score", "IL6_score", "IL6ST_score","JAK2_score", "NEK7_score", "NLRC4_score", "NLRP3_score",
"TNF_score", "TYK2_score", "CARD8_score", "IFNGR1_score","IL18_score", "IL18R1_score", "IL1RAP_score", 
 "JAK3_score", "STAT4_score", "STAT6_score")

for(m in 1:length(TRS)) {
	for(i in 1:length(chip_list)) {
      		h[i+(m-1)*length(chip_list),"chip"]<-chip_list[i]
		h[i+(m-1)*length(chip_list),"outcome"]<-TRS[m]
      		fmla<-as.formula(paste(TRS[m], " ~ ", chip_list[i],"+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+Prev_Diabetes_Type_2+ever_smoke",sep=""))
      		if(!is.na(Incd_outcome)){
			f<-glm(fmla,family = gaussian, data=pheno_newchip[which(pheno_newchip$Prev_AnyHemeCa==0 & pheno_newchip$Prev_composite_mi_cad_stroke_death==0 ),])
      			h[i+(m-1)*length(chip_list),"n"]<-summary(f)$df.null+1
      			h[i+(m-1)*length(chip_list),"beta"]<-summary(f)$coefficient[2,1]
      			h[i+(m-1)*length(chip_list),"se"]<-summary(f)$coefficient[2,2]
			h[i+(m-1)*length(chip_list),"z"]<-summary(f)$coefficient[2,3]
      			h[i+(m-1)*length(chip_list),"pval"]<-summary(f)$coefficient[2,4]
      			h[i+(m-1)*length(chip_list),"lci"]<-summary(f)$coefficient[2,1]-1.96*summary(f)$coefficient[2,2]
      			h[i+(m-1)*length(chip_list),"uci"]<-summary(f)$coefficient[2,1]+1.96*summary(f)$coefficient[2,2]
		} 	
	}
}
write.table(h, "/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_Updated1030_CHIP_scoreonly.txt", sep="\t", quote=FALSE, row.names=FALSE)	

#####################stratified analysis and interaction testing##################################
pheno_newchip$Incd_composite_mi_cad_stroke_death=ifelse(pheno_newchip$Prev_composite_mi_cad_stroke_death==1, NA, pheno_newchip$Incd_composite_mi_cad_stroke_death)

chip_var=c("CHIP","DNMT3A", "TET2","ASXL1","JAK2")

TRS=c("AIM2_score", "CASP1_score", "CASP5_score", "IFNGR2_score","IL10_score", "IL18BP_score", "IL18RAP_score", "IL1B_score",
"IL1R1_score", "IL1R2_score", "IL6_score", "IL6ST_score","JAK2_score", "NEK7_score", "NLRC4_score", "NLRP3_score",
"TNF_score", "TYK2_score", "CARD8_score", "IFNGR1_score","IL18_score", "IL18R1_score", "IL1RAP_score", 
 "JAK3_score", "STAT4_score", "STAT6_score")

for(j in 1:length(chip_var)) {

	h_chip<-data.frame(exposure =character(), chip=character(), n_nochip=double(), hr_nochip=double(), se_nochip=double(), z_nochip=double(), pval_nochip=double(), lci_nochip=double(), uci_nochip=double(), 
	n_chip=double(), hr_chip=double(), se_chip=double(), z_chip=double(), pval_chip=double(),lci_chip=double(), uci_chip=double(), beta_inter=double(), se_inter=double(), z_inter=double(), pval_inter=double(), stringsAsFactors=FALSE)      

	for(i in 1:length(TRS)) {
		print(paste(TRS[i], chip_var[j], sep="_"))
      		h_chip[i,"exposure"]<-TRS[i]
      		h_chip[i,"chip"]<-chip_var[j]
            	temp=pheno_newchip[which(pheno_newchip$Prev_composite_mi_cad_stroke_death==0 & pheno_newchip$Prev_AnyHemeCa==0),]
      		fmla<-as.formula(paste("Surv(FollowUp_composite_mi_cad_stroke_death, Incd_composite_mi_cad_stroke_death) ~ ", TRS[i],"+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+Prev_Diabetes_Type_2+ever_smoked",sep=""))#
      		f1<-coxph(fmla,data=temp[which(temp[which(colnames(temp)==chip_var[j])]==0),])
      		f2<-coxph(fmla,data=temp[which(temp[which(colnames(temp)==chip_var[j])]==1),])
		fmla_inter=as.formula(paste("Surv(FollowUp_composite_mi_cad_stroke_death, Incd_composite_mi_cad_stroke_death) ~ ", TRS[i], "*", chip_var[j], "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+Prev_Diabetes_Type_2+ever_smoked",sep=""))
		f3<-coxph(fmla_inter,data=temp)
      		h_chip[i,"n_nochip"]<-summary(f1)$n
      		h_chip[i,"hr_nochip"]<-summary(f1)$coefficient[1,2]
      		h_chip[i,"se_nochip"]<-summary(f1)$coefficient[1,3]
		h_chip[i,"z_nochip"]<-summary(f1)$coefficient[1,4]
      		h_chip[i,"pval_nochip"]<-summary(f1)$coefficient[1,5]
      		h_chip[i,"lci_nochip"]<-exp(confint(f1)[1,1])
      		h_chip[i,"uci_nochip"]<-exp(confint(f1)[1,2])
      		h_chip[i,"n_chip"]<-summary(f2)$n
      		h_chip[i,"hr_chip"]<-summary(f2)$coefficient[1,2]
      		h_chip[i,"se_chip"]<-summary(f2)$coefficient[1,3]
		h_chip[i,"z_chip"]<-summary(f2)$coefficient[1,4]
      		h_chip[i,"pval_chip"]<-summary(f2)$coefficient[1,5]
      		h_chip[i,"lci_chip"]<-exp(confint(f2)[1,1])
      		h_chip[i,"uci_chip"]<-exp(confint(f2)[1,2])
		h_chip[i,"beta_inter"]<-summary(f3)$coefficient[19,2]
     		h_chip[i,"se_inter"]<-summary(f3)$coefficient[19,3]
		h_chip[i,"z_inter"]<-summary(f3)$coefficient[19,4]
      		h_chip[i,"pval_inter"]<-summary(f3)$coefficient[19,5]
	}
	
	#write.table(h_chip, file=paste0("/medpop/esp2/zyu/chip_protemoics/output/aim2/FFinalfinal_20221030_stratify_PTPRSCS_ukbafreq_", chip_var[j], ".txt",sep=""), sep="\t", quote = FALSE,row.names = FALSE)
	
	if (j==1){
		h_chip_all=h_chip	
	}else{
		h_chip_all=rbind(h_chip_all, h_chip)
	}
}

write.table(h_chip_all,file="/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_20221030_stratify_PTPRSCS_ukbafreq_overall.txt",sep="\t",quote = FALSE,row.names = FALSE)

####################################plot##################################
h=fread("/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_20221030_stratify_PTPRSCS_ukbafreq_overall.txt")
h_nochip=h[,1:9]
h_chip=h[,c(1:2,10:16)]
colnames(h_nochip)=c("exposure","chip","n","hr","se","z","p","lci","uci")
h_nochip$chip=paste("No ", h_nochip$chip, sep="")
h_nochip$beta=log(h_nochip$hr)
colnames(h_chip)=c("exposure","chip","n","hr","se","z","p","lci","uci")
h_chip$beta=log(h_chip$hr)
df<-rbind(h_chip, h_nochip)
df$exposure=gsub(".{6}$","",df$exposure)

pdf(file = '/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_450K_stratify_hr_5genes.pdf',width=15, height=9) 
a=ggforestplot::forestplot(df = df[which(df$chip=="CHIP"|df$chip=="No CHIP"),],name = exposure,estimate = beta, 
se = se,pvalue = p,psignif = 0.05,colour = chip, title = "CHIP", logodds = TRUE)+ggplot2::scale_color_manual(values=c("#999999", "#404040"))+ggplot2::theme(legend.title = element_blank(),legend.position="bottom",legend.text = element_text(size=10),axis.title.y = element_text(face = "italic"))+ggplot2::xlim(0.75, 1.50)

b=ggforestplot::forestplot(df = df[which(df$chip=="DNMT3A"|df$chip=="No DNMT3A"),],name = exposure,estimate = beta, 
se = se,pvalue = p,psignif = 0.05,colour = chip,title = "DNMT3A", logodds = TRUE)+ggplot2::scale_color_manual(values=c("#E41A1C", "#404040"))+ggplot2::theme(legend.title = element_blank(),legend.position="bottom",legend.text = element_text(size=10),axis.title.y = element_text(face = "italic"))+ggplot2::xlim(0.75, 1.50)

c=ggforestplot::forestplot(df = df[which(df$chip=="TET2"|df$chip=="No TET2"),]%>% mutate(chip = fct_relevel(chip, "TET2", "No TET2")), name = exposure,estimate = beta, 
se = se,pvalue = p,psignif = 0.05,colour = chip, title = "TET2", logodds = TRUE)+ggplot2::scale_color_manual(values=c("#377EB8", "#404040"))+ggplot2::theme(legend.title = element_blank(),legend.position="bottom",legend.text = element_text(size=10),axis.title.y = element_text(face = "italic"))+ggplot2::xlim(0.75, 1.50)

d=ggforestplot::forestplot(df = df[which(df$chip=="ASXL1"|df$chip=="No ASXL1"),],name = exposure,estimate = beta, 
se = se,pvalue = p,psignif = 0.05,colour = chip, title = "ASXL1", logodds = TRUE)+ggplot2::scale_color_manual(values=c("#4DAF4A", "#404040"))+ggplot2::theme(legend.title = element_blank(),legend.position="bottom",legend.text = element_text(size=10),axis.title.y = element_text(face = "italic"))+ggplot2::xlim(0.75, 1.50)

e=ggforestplot::forestplot(df = df[which(df$chip=="JAK2"|df$chip=="No JAK2"),],name = exposure,estimate = beta, 
se = se,pvalue = p,psignif = 0.05,colour = chip, title = "JAK2", logodds = TRUE)+ggplot2::scale_color_manual(values=c("#FF7F00", "#404040"))+ggplot2::theme(legend.title = element_blank(),legend.position="bottom",legend.text = element_text(size=10),axis.title.y = element_text(face = "italic"))

figure <- ggarrange(a,b,c,d,e, ncol = 5, nrow = 1)
figure
dev.off()
#title: "HR for incident CVD events (95% CI)\nper 1 SD increment in TRS of each gene"

p0.05chip=unique(h_chip[p<0.05,]$exposure)
h=h%>%filter(exposure %in% p0.05chip) %>%mutate(z_inter =ifelse(pval_chip<0.05, z_inter, 0)) %>%select(exposure, chip, z_inter) %>% spread(chip, z_inter)%>%select(exposure, CHIP, DNMT3A, TET2, ASXL1, JAK2)
write.table(h,file=paste("/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_zforp0.05interaction.txt",sep=""),, sep="\t", quote=FALSE, row.names=F)

######################################blood panel modifier scan########################################
blood_trait_list=c( "White.blood.cell..leukocyte..count", "Red.blood.cell..erythrocyte..count", "Haemoglobin.concentration", "Haematocrit.percentage", "Mean.corpuscular.volume", "Mean.corpuscular.haemoglobin", "Mean.corpuscular.haemoglobin.concentration", 
"Red.blood.cell..erythrocyte..distribution.width", "Platelet.count", "Platelet.crit", "Mean.platelet..thrombocyte..volume", "Platelet.distribution.width", "Lymphocyte.count", "Monocyte.count", "Neutrophill.count", "Eosinophill.count", "Basophill.count", "Nucleated.red.blood.cell.count", 
"Lymphocyte.percentage", "Monocyte.percentage", "Neutrophill.percentage", "Eosinophill.percentage", "Basophill.percentage", "Nucleated.red.blood.cell.percentage", "Reticulocyte.percentage", "Reticulocyte.count", "Mean.reticulocyte.volume", "Mean.sphered.cell.volume", 
"Immature.reticulocyte.fraction", "High.light.scatter.reticulocyte.percentage", "High.light.scatter.reticulocyte.count")

score_list=c("AIM2_score","IL1RAP_score")

h_jak2=data.frame(blood_trait=character(), score=character(), beta=double(), se=double(), z=double(), p=double(), n=double())
for (j in 1:length(score_list)) {
	for (i in 1:length(blood_trait_list)) {
		temp=pheno_newchip[which(pheno_newchip$Prev_composite_mi_cad_stroke_death==0 & pheno_newchip$Prev_AnyHemeCa==0 &pheno_newchip$JAK2==1 ),]
		#temp[blood_trait_list[i]]=scale(temp[blood_trait_list[i]])
		temp[blood_trait_list[i]]=scale(log2(temp[blood_trait_list[i]]+1))
		#temp=temp[is.infinite(temp[blood_trait_list[i]])]
		fmla=as.formula(paste(blood_trait_list[i]," ~ ", score_list[j], "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+Prev_Diabetes_Type_2+ever_smoked", sep=""))
		f=glm(fmla, data=temp, family = "gaussian", na.action=na.exclude)
		h_jak2[i+(j-1)*length(blood_trait_list), "blood_trait"]=blood_trait_list[i]
		h_jak2[i+(j-1)*length(blood_trait_list), "score"]= score_list[j]
		h_jak2[i+(j-1)*length(blood_trait_list), "beta"]=summary(f)$coefficient[2, 1]
		h_jak2[i+(j-1)*length(blood_trait_list), "se"]=summary(f)$coefficient[2, 2]
		h_jak2[i+(j-1)*length(blood_trait_list), "z"]=summary(f)$coefficient[2, 1]/summary(f)$coefficient[2, 2]
		h_jak2[i+(j-1)*length(blood_trait_list), "p"]=summary(f)$coefficient[2, 4]
		h_jak2[i+(j-1)*length(blood_trait_list), "n"]=summary(f)$df.null+1
		print(paste(blood_trait_list[i], score_list[j]))
	}
	h_jak2$p_FDR_singlescore<-p.adjust(h_jak2$p, method = "fdr")
	h_jak2$chip="JAK2"
}

score_list=c("AIM2_score","IL18RAP_score")
h_asxl1=data.frame(blood_trait=character(), score=character(), beta=double(), se=double(), z=double(), p=double(), n=double())

for (j in 1:length(score_list)) {
	for (i in 1:length(blood_trait_list)) {
		temp=pheno_newchip[which(pheno_newchip$Prev_composite_mi_cad_stroke_death==0 & pheno_newchip$Prev_AnyHemeCa==0 &pheno_newchip$ASXL1==1 ),]
		#temp[blood_trait_list[i]]=scale(temp[blood_trait_list[i]])
		temp[blood_trait_list[i]]=scale(log2(temp[blood_trait_list[i]]+1))
		#temp=temp[is.infinite(temp[blood_trait_list[i]])]
		fmla=as.formula(paste(blood_trait_list[i]," ~ ", score_list[j], "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+Prev_Diabetes_Type_2+ever_smoked", sep=""))
		f=glm(fmla, data=temp, family = "gaussian", na.action=na.exclude)
		h_asxl1[i+(j-1)*length(blood_trait_list), "blood_trait"]=blood_trait_list[i]
		h_asxl1[i+(j-1)*length(blood_trait_list), "score"]= score_list[j]
		h_asxl1[i+(j-1)*length(blood_trait_list), "beta"]=summary(f)$coefficient[2, 1]
		h_asxl1[i+(j-1)*length(blood_trait_list), "se"]=summary(f)$coefficient[2, 2]
		h_asxl1[i+(j-1)*length(blood_trait_list), "z"]=summary(f)$coefficient[2, 1]/summary(f)$coefficient[2, 2]
		h_asxl1[i+(j-1)*length(blood_trait_list), "p"]=summary(f)$coefficient[2, 4]
		h_asxl1[i+(j-1)*length(blood_trait_list), "n"]=summary(f)$df.null+1
		print(paste(blood_trait_list[i], score_list[j]))
	}
	h_asxl1$p_FDR_singlescore<-p.adjust(h_asxl1$p, method = "fdr")
	h_asxl1$chip="ASXL1"
}

score_list=c("IL1RAP_score")
h_dnmt3a=data.frame(blood_trait=character(), score=character(), beta=double(), se=double(), z=double(), p=double(), n=double())

for (j in 1:length(score_list)) {
	for (i in 1:length(blood_trait_list)) {
		temp=pheno_newchip[which(pheno_newchip$Prev_composite_mi_cad_stroke_death==0 & pheno_newchip$Prev_AnyHemeCa==0 &pheno_newchip$DNMT3A==1 ),]
		#temp[blood_trait_list[i]]=scale(temp[blood_trait_list[i]])
		temp[blood_trait_list[i]]=scale(log2(temp[blood_trait_list[i]]+1))
		#temp=temp[is.infinite(temp[blood_trait_list[i]])]
		fmla=as.formula(paste(blood_trait_list[i]," ~ ", score_list[j], "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+Prev_Diabetes_Type_2+ever_smoked", sep=""))
		f=glm(fmla, data=temp, family = "gaussian", na.action=na.exclude)
		h_dnmt3a[i+(j-1)*length(blood_trait_list), "blood_trait"]=blood_trait_list[i]
		h_dnmt3a[i+(j-1)*length(blood_trait_list), "score"]= score_list[j]
		h_dnmt3a[i+(j-1)*length(blood_trait_list), "beta"]=summary(f)$coefficient[2, 1]
		h_dnmt3a[i+(j-1)*length(blood_trait_list), "se"]=summary(f)$coefficient[2, 2]
		h_dnmt3a[i+(j-1)*length(blood_trait_list), "z"]=summary(f)$coefficient[2, 1]/summary(f)$coefficient[2, 2]
		h_dnmt3a[i+(j-1)*length(blood_trait_list), "p"]=summary(f)$coefficient[2, 4]
		h_dnmt3a[i+(j-1)*length(blood_trait_list), "n"]=summary(f)$df.null+1
		print(paste(blood_trait_list[i], score_list[j]))
	}
	h_dnmt3a$p_FDR_singlescore<-p.adjust(h_dnmt3a$p, method = "fdr")
	h_dnmt3a$chip="DNMT3A"
}

score_list=c("IL1RAP_score")
h_chip=data.frame(blood_trait=character(), score=character(), beta=double(), se=double(), z=double(), p=double(), n=double())

for (j in 1:length(score_list)) {
	for (i in 1:length(blood_trait_list)) {
		temp=pheno_newchip[which(pheno_newchip$Prev_composite_mi_cad_stroke_death==0 & pheno_newchip$Prev_AnyHemeCa==0 &pheno_newchip$CHIP==1 ),]
		temp[blood_trait_list[i]]=scale(log2(temp[blood_trait_list[i]]+1))
		fmla=as.formula(paste(blood_trait_list[i]," ~ ", score_list[j], "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+Prev_Diabetes_Type_2+ever_smoked", sep=""))
		f=glm(fmla, data=temp, family = "gaussian", na.action=na.exclude)
		h_chip[i+(j-1)*length(blood_trait_list), "blood_trait"]=blood_trait_list[i]
		h_chip[i+(j-1)*length(blood_trait_list), "score"]= score_list[j]
		h_chip[i+(j-1)*length(blood_trait_list), "beta"]=summary(f)$coefficient[2, 1]
		h_chip[i+(j-1)*length(blood_trait_list), "se"]=summary(f)$coefficient[2, 2]
		h_chip[i+(j-1)*length(blood_trait_list), "z"]=summary(f)$coefficient[2, 1]/summary(f)$coefficient[2, 2]
		h_chip[i+(j-1)*length(blood_trait_list), "p"]=summary(f)$coefficient[2, 4]
		h_chip[i+(j-1)*length(blood_trait_list), "n"]=summary(f)$df.null+1
		print(paste(blood_trait_list[i], score_list[j]))
	}
	h_chip$p_FDR_singlescore<-p.adjust(h_chip$p, method = "fdr")
	h_chip$chip="CHIP"
}

h=rbind(h_jak2, h_dnmt3a, h_chip, h_asxl1)

h$p_FDR_overall<-p.adjust(h$p, method = "fdr")
write.table(h, file="/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_chiponly_stratified_bloodtraits.txt", sep="\t", quote=FALSE, row.names=FALSE)

h_z=h%>%separate(score, c("score", NA))%>%mutate(inter=paste(score, chip, sep="_"))%>%select(blood_trait, inter, z)
h_z=reshape(h_z, idvar = "blood_trait", timevar = "inter", direction = "wide")
h_z=h_z%>%select(blood_trait, z.IL1RAP_DNMT3A, z.IL1RAP_CHIP,  z.AIM2_ASXL1, z.IL18RAP_ASXL1, z.IL1RAP_JAK2, z.AIM2_JAK2)     

write.table(h_z, file="/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_Zscore_chiponly_stratified_bloodtraits.txt", sep="\t", quote=FALSE, row.names=FALSE)

######################################blood panel association########################################

chip_list=chip_var

h=data.frame(blood_trait=character(), chip=character(), beta=double(), se=double(), z=double(), p=double(), n=double())
for (j in 1:length(chip_list)) {
	for (i in 1:length(blood_trait_list)) {
		temp=pheno_newchip[which(pheno_newchip$Prev_AnyHemeCa==0 & pheno_newchip$Prev_composite_mi_cad_stroke_death==0),]
		temp[blood_trait_list[i]]=scale(log2(temp[blood_trait_list[i]]+1))
		fmla=as.formula(paste(blood_trait_list[i]," ~ ", chip_list[j], "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+ever_smoke+BMI+Prev_Diabetes_Type_2+ever_smoked", sep=""))
		f=glm(fmla, data=temp, family = "gaussian", na.action=na.exclude)
		h[i+(j-1)*length(chip_list), "blood_trait"]=blood_trait_list[i]
		h[i+(j-1)*length(chip_list), "chip"]= chip_list[j]
		h[i+(j-1)*length(chip_list), "beta"]=summary(f)$coefficient[2, 1]
		h[i+(j-1)*length(chip_list), "se"]=summary(f)$coefficient[2, 2]
		h[i+(j-1)*length(chip_list), "z"]=summary(f)$coefficient[2, 1]/summary(f)$coefficient[2, 2]
		h[i+(j-1)*length(chip_list), "p"]=summary(f)$coefficient[2, 4]
		h[i+(j-1)*length(chip_list), "n"]=summary(f)$df.null+1
		print(paste(blood_trait_list[i], chip_list[j]))
	}
	h$p_FDR_singlechip<-p.adjust(h$p, method = "fdr")
}
h$p_FDR_overall<-p.adjust(h$p, method = "fdr")
write.table(h, file="/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_chip_bloodtraits.txt", sep="\t", quote=FALSE, row.names=TRUE)


#Supplemental table 6
tab1=CreateTableOne(vars=blood_trait_list, data=pheno_newchip)
tab1_chip=CreateTableOne(vars=blood_trait_list, strata="CHIP", data=pheno_newchip)

median=blood_trait_list

tab1_export=print(tab1, nonnormal=median, exact="stage", quote=FALSE, noSpaces=TRUE, printTOggle=FALSE)
tab1_chip_export=print(tab1_chip, nonnormal=median, exact="stage", quote=FALSE, noSpaces=TRUE, printTOggle=FALSE)

write.csv(tab1_export,file="/medpop/esp2/zyu/chip_protemoics/output/aim2/blood_12052021.csv")
write.csv(tab1_chip_export,file="/medpop/esp2/zyu/chip_protemoics/output/aim2/blood_chip_stratum_12052021.csv")

######################################biomarker modifier scan########################################
biomarker_trait_list=c("C.reactive.protein", "Cholesterol", "HDL.cholesterol", "LDL.direct","Triglycerides")

score_list=c("AIM2_score","IL1RAP_score")

h_jak2=data.frame(biomarker_trait=character(), score=character(), beta=double(), se=double(), z=double(), p=double(), n=double())
for (j in 1:length(score_list)) {
	for (i in 1:length(biomarker_trait_list)) {
		temp=pheno_newchip[which(pheno_newchip$Prev_composite_mi_cad_stroke_death==0 & pheno_newchip$Prev_AnyHemeCa==0 &pheno_newchip$JAK2==1 ),]
		#temp[biomarker_trait_list[i]]=scale(temp[biomarker_trait_list[i]])
		temp[biomarker_trait_list[i]]=scale(log2(temp[biomarker_trait_list[i]]+1))
		#temp=temp[is.infinite(temp[biomarker_trait_list[i]])]
		fmla=as.formula(paste(biomarker_trait_list[i]," ~ ", score_list[j], "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+Prev_Diabetes_Type_2+ever_smoked", sep=""))
		f=glm(fmla, data=temp, family = "gaussian", na.action=na.exclude)
		h_jak2[i+(j-1)*length(biomarker_trait_list), "biomarker_trait"]=biomarker_trait_list[i]
		h_jak2[i+(j-1)*length(biomarker_trait_list), "score"]= score_list[j]
		h_jak2[i+(j-1)*length(biomarker_trait_list), "beta"]=summary(f)$coefficient[2, 1]
		h_jak2[i+(j-1)*length(biomarker_trait_list), "se"]=summary(f)$coefficient[2, 2]
		h_jak2[i+(j-1)*length(biomarker_trait_list), "z"]=summary(f)$coefficient[2, 1]/summary(f)$coefficient[2, 2]
		h_jak2[i+(j-1)*length(biomarker_trait_list), "p"]=summary(f)$coefficient[2, 4]
		h_jak2[i+(j-1)*length(biomarker_trait_list), "n"]=summary(f)$df.null+1
		print(paste(biomarker_trait_list[i], score_list[j]))
	}
	h_jak2$p_FDR_singlescore<-p.adjust(h_jak2$p, method = "fdr")
	h_jak2$chip="JAK2"
}

score_list=c("AIM2_score","IL18RAP_score")
h_asxl1=data.frame(biomarker_trait=character(), score=character(), beta=double(), se=double(), z=double(), p=double(), n=double())

for (j in 1:length(score_list)) {
	for (i in 1:length(biomarker_trait_list)) {
		temp=pheno_newchip[which(pheno_newchip$Prev_composite_mi_cad_stroke_death==0 & pheno_newchip$Prev_AnyHemeCa==0 &pheno_newchip$ASXL1==1 ),]
		#temp[biomarker_trait_list[i]]=scale(temp[biomarker_trait_list[i]])
		temp[biomarker_trait_list[i]]=scale(log2(temp[biomarker_trait_list[i]]+1))
		#temp=temp[is.infinite(temp[biomarker_trait_list[i]])]
		fmla=as.formula(paste(biomarker_trait_list[i]," ~ ", score_list[j], "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+Prev_Diabetes_Type_2+ever_smoked", sep=""))
		f=glm(fmla, data=temp, family = "gaussian", na.action=na.exclude)
		h_asxl1[i+(j-1)*length(biomarker_trait_list), "biomarker_trait"]=biomarker_trait_list[i]
		h_asxl1[i+(j-1)*length(biomarker_trait_list), "score"]= score_list[j]
		h_asxl1[i+(j-1)*length(biomarker_trait_list), "beta"]=summary(f)$coefficient[2, 1]
		h_asxl1[i+(j-1)*length(biomarker_trait_list), "se"]=summary(f)$coefficient[2, 2]
		h_asxl1[i+(j-1)*length(biomarker_trait_list), "z"]=summary(f)$coefficient[2, 1]/summary(f)$coefficient[2, 2]
		h_asxl1[i+(j-1)*length(biomarker_trait_list), "p"]=summary(f)$coefficient[2, 4]
		h_asxl1[i+(j-1)*length(biomarker_trait_list), "n"]=summary(f)$df.null+1
		print(paste(biomarker_trait_list[i], score_list[j]))
	}
	h_asxl1$p_FDR_singlescore<-p.adjust(h_asxl1$p, method = "fdr")
	h_asxl1$chip="ASXL1"
}

score_list=c("IL1RAP_score")
h_dnmt3a=data.frame(biomarker_trait=character(), score=character(), beta=double(), se=double(), z=double(), p=double(), n=double())

for (j in 1:length(score_list)) {
	for (i in 1:length(biomarker_trait_list)) {
		temp=pheno_newchip[which(pheno_newchip$Prev_composite_mi_cad_stroke_death==0 & pheno_newchip$Prev_AnyHemeCa==0 &pheno_newchip$DNMT3A==1 ),]
		#temp[biomarker_trait_list[i]]=scale(temp[biomarker_trait_list[i]])
		temp[biomarker_trait_list[i]]=scale(log2(temp[biomarker_trait_list[i]]+1))
		#temp=temp[is.infinite(temp[biomarker_trait_list[i]])]
		fmla=as.formula(paste(biomarker_trait_list[i]," ~ ", score_list[j], "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+Prev_Diabetes_Type_2+ever_smoked", sep=""))
		f=glm(fmla, data=temp, family = "gaussian", na.action=na.exclude)
		h_dnmt3a[i+(j-1)*length(biomarker_trait_list), "biomarker_trait"]=biomarker_trait_list[i]
		h_dnmt3a[i+(j-1)*length(biomarker_trait_list), "score"]= score_list[j]
		h_dnmt3a[i+(j-1)*length(biomarker_trait_list), "beta"]=summary(f)$coefficient[2, 1]
		h_dnmt3a[i+(j-1)*length(biomarker_trait_list), "se"]=summary(f)$coefficient[2, 2]
		h_dnmt3a[i+(j-1)*length(biomarker_trait_list), "z"]=summary(f)$coefficient[2, 1]/summary(f)$coefficient[2, 2]
		h_dnmt3a[i+(j-1)*length(biomarker_trait_list), "p"]=summary(f)$coefficient[2, 4]
		h_dnmt3a[i+(j-1)*length(biomarker_trait_list), "n"]=summary(f)$df.null+1
		print(paste(biomarker_trait_list[i], score_list[j]))
	}
	h_dnmt3a$p_FDR_singlescore<-p.adjust(h_dnmt3a$p, method = "fdr")
	h_dnmt3a$chip="DNMT3A"
}

score_list=c("IL1RAP_score")
h_chip=data.frame(biomarker_trait=character(), score=character(), beta=double(), se=double(), z=double(), p=double(), n=double())

for (j in 1:length(score_list)) {
	for (i in 1:length(biomarker_trait_list)) {
		temp=pheno_newchip[which(pheno_newchip$Prev_composite_mi_cad_stroke_death==0 & pheno_newchip$Prev_AnyHemeCa==0 &pheno_newchip$CHIP==1 ),]
		temp[biomarker_trait_list[i]]=scale(log2(temp[biomarker_trait_list[i]]+1))
		fmla=as.formula(paste(biomarker_trait_list[i]," ~ ", score_list[j], "+in_white_British_ancestry_subset+age+Sex+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+BMI+Prev_Diabetes_Type_2+ever_smoked", sep=""))
		f=glm(fmla, data=temp, family = "gaussian", na.action=na.exclude)
		h_chip[i+(j-1)*length(biomarker_trait_list), "biomarker_trait"]=biomarker_trait_list[i]
		h_chip[i+(j-1)*length(biomarker_trait_list), "score"]= score_list[j]
		h_chip[i+(j-1)*length(biomarker_trait_list), "beta"]=summary(f)$coefficient[2, 1]
		h_chip[i+(j-1)*length(biomarker_trait_list), "se"]=summary(f)$coefficient[2, 2]
		h_chip[i+(j-1)*length(biomarker_trait_list), "z"]=summary(f)$coefficient[2, 1]/summary(f)$coefficient[2, 2]
		h_chip[i+(j-1)*length(biomarker_trait_list), "p"]=summary(f)$coefficient[2, 4]
		h_chip[i+(j-1)*length(biomarker_trait_list), "n"]=summary(f)$df.null+1
		print(paste(biomarker_trait_list[i], score_list[j]))
	}
	h_chip$p_FDR_singlescore<-p.adjust(h_chip$p, method = "fdr")
	h_chip$chip="CHIP"
}

h=rbind(h_jak2, h_dnmt3a, h_chip, h_asxl1)

h$p_FDR_overall<-p.adjust(h$p, method = "fdr")
write.table(h, file="/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_chiponly_stratified_biomarkertraits.txt", sep="\t", quote=FALSE, row.names=FALSE)

h_z=h%>%separate(score, c("score", NA))%>%mutate(inter=paste(score, chip, sep="_"))%>%select(biomarker_trait, inter, z)
h_z=reshape(h_z, idvar = "biomarker_trait", timevar = "inter", direction = "wide")
h_z=h_z%>%select(biomarker_trait, z.IL1RAP_DNMT3A, z.IL1RAP_CHIP,  z.AIM2_ASXL1, z.IL18RAP_ASXL1, z.IL1RAP_JAK2, z.AIM2_JAK2)     

write.table(h_z, file="/medpop/esp2/zyu/chip_protemoics/output/aim2/Finalfinal_Zscore_chiponly_stratified_biomarkertraits.txt", sep="\t", quote=FALSE, row.names=FALSE)

