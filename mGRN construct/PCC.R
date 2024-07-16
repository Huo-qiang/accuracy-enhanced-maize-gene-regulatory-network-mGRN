library(reshape2)

TPM_expression <- 'RNA_fpkm.txt'
TFs <- 'TFgenesallmaize.txt'
candidate_num <- 1600 # dynamic number

data <- read.delim(TPM_expression,stringsAsFactors = F)
data <- t(log2(data+1))

TFs <- read.delim(TFs,header = F,stringsAsFactors = F)
TFs <- TFs[which(TFs$V1%in%colnames(data)),1]

for(i in 1:length(TFs)){
  tmp_TF <- TFs[i]
  tmp_PCC <- data.frame(cor(data[,tmp_TF],data))
  rownames(tmp_PCC) <- tmp_TF
  tmp_PCC <- melt(tmp_PCC)
  tmp_PCC$TF <- tmp_TF
  tmp_PCC$value <- abs(tmp_PCC$value)
  tmp_PCC <- tmp_PCC[order(tmp_PCC$value,decreasing = T),]
  tmp_PCC2 <- tmp_PCC[1:(candidate_num+1),]
  if(i==1){
    total_PCC <- tmp_PCC2
  }else{
    total_PCC <- rbind(total_PCC,tmp_PCC2)
  } 
}

write.table(total_PCC,'PCCs.txt',row.names = F,sep = '\t',quote = F)
system("cat PCCs.txt | gawk -vOFS='\t' '{if($1!=$3)print $3,$1,$2}' > PCCs.unique.txt")
