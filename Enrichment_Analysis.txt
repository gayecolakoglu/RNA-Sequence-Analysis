library(clusterProfiler)
library(GSEABase)
library(org.Hs.eg.db)

#R Codes for finding ENTREZID for DE Top20
hs <- org.Hs.eg.db
my.symbols<-c("FTX","CCN2","POFUT1","EEF1A1P38","PIP4K2A","ZMYM3","PTPRA","FBLN5","MIR1244-1","ANO1","C1S","LMOD1","RPLP0P6","C7orf50","PTP4A2","RBM23","MRPL42","APMAP","RPL10AP2","DHRS7B","HSPA1B","BRD2")
select(hs, keys = my.symbols,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")

ENTREZID:
	DE genes: 100302692,1490,23509,124199,5305,9203,5786,10516,100302285,55107,716,25802,220717,84310,8073,55147,28977,57136,253986,25979

#R Codes for finding ENTREZID Annovar Top20 variants gene
	hs <- org.Hs.eg.db
	my.symbols<-c("ZC3H12A","HSPA6","DUSP1","HSPA1A","HSPA1B","BRD2","CXCL8","HSPH1")
	select(hs, keys = my.symbols,columns = c("ENTREZID", "SYMBOL"),keytype = "SYMBOL")
ENTREZID:
	Annovar genes: 80149, 3310, 1843, 3303, 3304, 6046, 3576, 10808 

#GSEA for 6DE genes
filename <- "c7.all.v7.1.entrez.gmt" 
gmtfile <- system.file(filename)
c6 <- read.gmt(gmtfile)
yourEntrezIdList<- c(100302692,1490,5786,10516,100302285,55107) #ENTREZID of DE genes 
ImmunSigEnrich <- enricher(yourEntrezIdList, TERM2GENE=c6, pvalueCutoff = 0.05)
ImmunSigEnrich <- setReadable(ImmunSigEnrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
#write.csv(ImmunSigEnrich,"MyImmunePathwayRelatedGenes6DE.csv")

goEnrich<-enrichGO(gene= yourEntrezIdList,OrgDb= org.Hs.eg.db, ont= "ALL",pAdjustMethod="BH",pvalueCutoff = 0.05,readable= TRUE)
#write.csv(goEnrich,"MyGORelatedGenes6DE.csv")

keggEnrich<-enrichKEGG(gene= yourEntrezIdList,organism= "hsa",pAdjustMethod="BH", pvalueCutoff = 0.05)
#write.csv(keggEnrich,"MyKEGGRelatedGenes6DE.csv")


#Annovar genes
filename <- "c7.all.v7.1.entrez.gmt"
gmtfile <- system.file(filename)
c6 <- read.gmt(gmtfile)
yourEntrezIdList<- c(80149,3310,1843,3303,3304,6046,3576,10808) #ENTREZID of Annovar genes
ImmunSigEnrich <- enricher(yourEntrezIdList, TERM2GENE=c6, pvalueCutoff = 0.05)
ImmunSigEnrich <- setReadable(ImmunSigEnrich, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
#write.csv(ImmunSigEnrich,"MyImmunePathwayRelatedGenesAnn.csv")

goEnrich<-enrichGO(gene= yourEntrezIdList,OrgDb= org.Hs.eg.db, ont= "ALL",pAdjustMethod="BH",pvalueCutoff = 0.05,readable= TRUE)
#write.csv(goEnrich,"MyGORelatedGenesAnn.csv")

keggEnrich<-enrichKEGG(gene= yourEntrezIdList,organism= "hsa",pAdjustMethod="BH", pvalueCutoff = 0.05)
#write.csv(keggEnrich,"MyKEGGRelatedGenesAnn.csv")

#Create outfile for putput plots
outfile="/media/data04/gaye/workspace/rnaseq/entrez/GSEA_DE_output.pdf"
#outfile="/media/data04/gaye/workspace/rnaseq/entrez/GSEA_Annovar_output.pdf"
pdf(file=outfile)
yourEntrezIdList <- sort(yourEntrezIdList,decreasing=TRUE)

library(enrichplot)TRUE)
barplot(goEnrich, showCategory=20)
dotplot(goEnrich, showCategory=30)

barplot(keggEnrich, showCategory=20)
dotplot(keggEnrich, showCategory=30)
cnetplot(keggEnrich) #for plotting KEGG pathways 

cnetplot(ImmunSigEnrich) #for plotting ImmunSigEnrich pathways 

dev.off()

# Exit the R session
quit(save="no")
