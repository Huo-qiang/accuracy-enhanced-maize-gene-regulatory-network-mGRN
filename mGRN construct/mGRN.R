## Integration of GRNboost2 and PCC Networks
library(dplyr)
GRNboost2_grn<-read.table("maize_mGRN.tsv",header = T,sep = "\t")
pcc<-read.table("PCCs.unique.txt",header = T,sep = "\t")
grnboost2_pcc = inner_join(GRNboost2_grn,pcc,by = c("TF","Target"))
###
atac_gene<-read.csv("anno_atac.csv"，header=T)
grnboost2_pcc_atac<-inner_join(grnboost2_pcc,atac_gene,by = "Target")
###
cisbp_grn<-read.csv("cisbp_grn.csv",header = T)
grnboost2_pcc_atac_cisbp<-inner_join(grnboost2_pcc_atac,cisbp_grn,by = c("TF","Target"))
##
dap_grn<-read.csv("dap_grn",header= T)
grnboost2_pcc_atac_cisbp_dap<-inner_join(grnboost2_pcc_atac_cisbp,dap_grn,by = c("TF","Target"))
write.csv(grnboost2_pcc_atac_cisbp_dap,"all_mGRN.csv")
