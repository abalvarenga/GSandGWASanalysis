#reading genotype:
geno=as.data.frame(read.table("data.new",header=F))
indid=as.data.frame(geno[1,])
geno=as.data.frame(geno[-1,])
n=nrow(geno)
geno=geno[order(geno[,3]),]
#SNPs/chr
SNP_chr=table(geno[,2])
write.table(SNP_chr,"SNP_per_chr.txt", quote=F, col.names=T, row.names=F,sep=" ")

new_id=as.data.frame(seq(1:n)) #creating a sequence to rename the SNP
geno2=data.frame(new_id,geno)
a="new_id"
indid2=data.frame(a,indid)

#SNP_map
map=as.data.frame(geno2[,c(1,3,4,2)])
write.table(map,"map.txt", quote=F, col.names=F, row.names=F,sep=" ")

#SNP_file
file=as.data.frame(geno2[,-c(1,2,3,4)])
indfile=as.data.frame(indid2[,-c(1,2,3,4)])
file=as.data.frame(ifelse(file=="AA",0, ifelse(file=="BB",2,ifelse(file=="AB",1,ifelse(file=="BA",1,5)))))

tfile=as.data.frame(t(file))
tindfile=as.data.frame(t(indfile))
genofinal=data.frame(tindfile,tfile)
write.table(genofinal,"geno_pre.txt", quote=F, col.names=F, row.names=F,sep=" ")


