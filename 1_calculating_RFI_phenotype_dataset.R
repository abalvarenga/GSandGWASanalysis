#' ---
#' title: "Estimating metrics (RFI) for ovine"
#' author: "Amanda Botelho Alvarenga"
#' date: "Sep 10, 2018"
#' ---

#---------------------------------------------------------------------------------------------------------------
#                                   Organizing dataset and estimating metrics RFI
#---------------------------------------------------------------------------------------------------------------

# Importing datasets
setwd("C:\\Users\\Amanda\\Desktop\\Artigos_em_desenvolvimento\\GWAS-GS\\New_analysis_genomic_2018\\Metrics\\Data_input")
phen=read.table("1phenotypes.txt",header=T);head(phen) # all phenotypes reported in experiment- crude phenotypes
a=read.table("2AOL_lote4.txt",header=T);head(a) # Muscle depth information for 4 contemporaneos group
GC_info=read.table("3GC_phen.txt",header=T);head(GC_info) # Information of contemporaneous group

#Merging phen and a
phen_1=merge(phen,a,by= intersect("lote_Id_amostra","lote_Id_amostra"),all=TRUE);head(phen_1)
phen_1$AOL=phen_1$AOL.x
phen_1$AOL[!is.na(phen_1$AOL.y)] = phen_1$AOL.y[!is.na(phen_1$AOL.y)];head(phen_1)
head(GC_info)

#Merging phen_1 with GC_info
colnames(GC_info)=c("Id_amostra","Baia","Lote","Trat","CMS_kg","lote_Id_amostra","GC_trat")
phen_final=merge(phen_1,GC_info,by= c("lote_Id_amostra", "lote_Id_amostra"),all=FALSE);head(phen_final)
nrow(phen_final)

#____________________________Calculating RFI from differents methods ___________________________________________________ 
# The models for RFI was applied to sheep in Knott 2008 (doi:10.1016/j.anifeedsci.2007.05.013)

# Average daily gain: ADG
phen_final$ADG=(phen_final$Peso_abate-phen_final$Peso_Inicial)/phen_final$DIAS_CONF

# Feed intake: FI
phen_final$FI=phen_final$CMS_kg

# Feed conversion rate: FCR
phen_final$FCR=phen_final$ADG/phen_final$FI
fcr1=phen_final[,c("ID_genotipagem","FCR","GC_trat")]
fcr=fcr1[order(fcr1$ID_genotipagem),]

# Mid-test of liveweight
phen_final$mt=(phen_final$Peso_abate+phen_final$Peso_Inicial)/2

# MMWT
phen_final$mmwt=((phen_final$Peso_abate+phen_final$Peso_Inicial)/2)^0.73

# MWT
phen_final$mwt=(phen_final$Peso_Inicial)+((phen_final$Peso_abate-phen_final$Peso_Inicial)/2)

# 1) Estimating the regression for FI by Francois
reg=lm(phen_final$FI~phen_final$GC_trat+phen_final$ADG+phen_final$mt+phen_final$AOL + phen_final$EGS)
summary(reg)
coef=as.matrix(reg$coefficients)

# 2) Estimating by Australian beef cattle
phen_final$FI_t=phen_final$FI*phen_final$DIAS_CONF
reg2=lm(phen_final$FI_t~phen_final$GC_trat+phen_final$ADG+phen_final$mmwt)
summary(reg2)
coef2=t(as.matrix(reg2$coefficients));coef2

# 3) Estimating by 1963 models of residual feed intake and residual liveweight gain
reg3=lm(phen_final$FI~phen_final$GC_trat+phen_final$ADG+phen_final$mwt)
summary(reg3)
coef3=as.matrix(reg3$coefficients)

# Estimating RFI
# The methodology "3) RFI by Australian beef cattle" was selected.
#		It method had a simular r² among methods, however, to calculate it is more easier than other
#			due to use commonly traits reported.
phen_final$FE=coef2[,1]+coef2[,2]*phen_final$GC_trat+coef2[,3]*phen_final$ADG+coef2[,3]*phen_final$mmwt
phen_final$RFI=phen_final$CMS_kg-phen_final$FE
head(phen_final)

rfi=cbind(phen_final$ID_genotipagem,phen_final$RFI,phen_final$GC_trat)
colnames(rfi)=c("SNP_name","rfi","gc$");head(rfi)
write.table(rfi,"gensel_rfi.txt", quote=F, col.names=T, row.names=F,sep="\t")

#__________________ Saving dataset _______________________________________________________________________
# FCR
fcr=cbind(phen_final$ID_genotipagem,phen_final$FCR,phen_final$GC_trat)
colnames(fcr)=c("SNP_name","fcr","gc$");head(fcr)
write.table(fcr,"gensel_fcr.txt", quote=F, col.names=T, row.names=F,sep="\t")

# FI total
phen_final$FI_t=phen_final$FI*phen_final$DIAS_CONF
fi=cbind(phen_final$ID_genotipagem,phen_final$FI_t,phen_final$GC_trat)
colnames(fi)=c("SNP_name","FI","gc$");head(fi)
write.table(fi,"gensel_fi.txt", quote=F, col.names=T, row.names=F,sep="\t")

# ADG
adg=cbind(phen_final$ID_genotipagem,phen_final$ADG,phen_final$GC_trat)
colnames(adg)=c("SNP_name","ADG","gc$");head(adg)
write.table(adg,"gensel_adg.txt", quote=F, col.names=T, row.names=F,sep="\t")

# MMWT
mmwt=cbind(phen_final$ID_genotipagem,phen_final$mmwt,phen_final$GC_trat)
colnames(mmwt)=c("SNP_name","mmwt","gc$");head(mmwt)
write.table(mmwt,"gensel_mmwt.txt", quote=F, col.names=T, row.names=F,sep="\t")

# Dataset with all phenotypes
phen_final$FI_total=phen_final$CMS_kg*phen_final$DIAS_CONF
phen=phen_final[,c(5,20,10,21,23,25,28,29)]
head(phen)
write.table(phen,"data_all.txt", quote=F, col.names=T, row.names=F,sep="\t")
