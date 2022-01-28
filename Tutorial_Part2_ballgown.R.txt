library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# Define a path for the output PDF to be written
outfile="/media/data04/gaye/workspace/rnaseq/de/ballgown/ref_only/Tutorial_Part2_ballgown_output.pdf"

# Load phenotype data
pheno_data = read.csv("NORMAL_vs_INFECTED.csv")

# Display the phenotype data
pheno_data

# Load the ballgown object from file
load('bg.rda')

# The load command, loads an R object from a file into memory in our R session. 
# You can use ls() to view the names of variables that have been loaded
ls()

# Print a summary of the ballgown object
bg

# Open a PDF file where we will save some plots. We will save all figures and then view the PDF at the end
pdf(file=outfile)

# Extract FPKM values from the 'bg' object
fpkm = texpr(bg,meas="FPKM")

# View the last several rows of the FPKM table
tail(fpkm)

# Transform the FPKM values by adding 1 and convert to a log2 scale
fpkm = log2(fpkm+1)

# View the last several rows of the transformed FPKM table
tail(fpkm)

# Create boxplots to display summary statistics for the FPKM values for each sample
boxplot(fpkm,col=as.numeric(pheno_data$type)+1,las=2,ylab='log2(FPKM+1)')

#function for finding the top gene row in trancsript
for (val in 1: 227818)
{
    if (geneNames(bg)[val] == "FTX") {
	print(val)
	}
}

# Display the transcript ID for a single row of data
ballgown::transcriptNames(bg)[223535]

# Display the gene name for a single row of data 
ballgown::geneNames(bg)[223535]

# Create a BoxPlot comparing the expression of a single gene for all replicates of both conditions
plot(fpkm[223535,] ~ as.factor(pheno_data$type), border=c(2,3), main=paste(ballgown::geneNames(bg)[223535],' : ', ballgown::transcriptNames(bg)[223535]),pch=19, xlab="Type", ylab='log2(FPKM+1)')

# Add the FPKM values for each sample onto the plot 
points(fpkm[223535,] ~ jitter(as.numeric(pheno_data$type)), col=as.numeric(pheno_data$type)+1, pch=16)

# Create a plot of transcript structures observed in each replicate and color transcripts by expression level
plotTranscripts(ballgown::geneIDs(bg)[223535], bg, main=c('Gene in sample SRR12424255'), sample=c('SRR12424255'))
plotTranscripts(ballgown::geneIDs(bg)[223535], bg, main=c('Gene in sample SRR12424256'), sample=c('SRR12424256'))
plotTranscripts(ballgown::geneIDs(bg)[223535], bg, main=c('Gene in sample SRR12424257'), sample=c('SRR12424257'))
plotTranscripts(ballgown::geneIDs(bg)[223535], bg, main=c('Gene in sample SRR12424243'), sample=c('SRR12424243'))
plotTranscripts(ballgown::geneIDs(bg)[223535], bg, main=c('Gene in sample SRR12424244'), sample=c('SRR12424244'))
plotTranscripts(ballgown::geneIDs(bg)[223535], bg, main=c('Gene in sample SRR12424245'), sample=c('SRR12424245'))

#plotMeans('FTX',bg,groupvar="type",legend=FALSE)

# Close the PDF device where we have been saving our plots
dev.off()

# Exit the R session
quit(save="no")
