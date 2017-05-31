# Example script on how to perform a WGCNA analysis
# Gokul Ramaswami 5-26-2017
# Input dataset: Brainspan RNA-seq samples from 6 brain regions spanning fetal development through old age 
# Network: RNA editing levels measured at ~2.5 million editing sites from RADAR v2 (www.rnaedit.com)
# Objective: Generate a Co-Editing network throughout brain development

# set up R and load libraries
rm(list=ls());
options(stringsAsFactors = FALSE);
library(WGCNA);
library(flashClust);
library(gplots)

# Output directory for figures
outputDir <- "../figures/";
# Output directory for data
outputDirData <- "../processed_data/";

# Load Metadata
metaData <- read.table("../Metadata.txt",sep=" ",header=TRUE)
rownames(metaData) = paste(metaData[,"Braincode"],metaData[,"Region"],sep=".")

# Define 6 gross regions: Neocortex (NCX), Amygdala (AMY), Striatum (STR), Hippocampus (HIP), Cerebellum (CBC), Thalamus (MD)
regions <-       c("NCX","NCX","NCX","AMY","NCX","NCX","NCX","STR","NCX","HIP","NCX","NCX","NCX","NCX","CBC","MD")
names(regions) = c("S1C","V1C","A1C","AMY","OFC","VFC","MFC","STR","ITC","HIP","STC","M1C","DFC","IPC","CBC","MD")

# Define developmental time periods: from Kang et al. 2011
periods <- c(2,2,3,4,5,5,6,6,6,7,7,8,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,14)
names(periods) = c("8PCW","9PCW","12PCW","13PCW","16PCW","17PCW","19PCW","21PCW","22PCW","25PCW","26PCW","4M","6M","10M","12M","2Y","3Y","4Y", "8Y","11Y","13Y","15Y","18Y","19Y","21Y","23Y","30Y","36Y","37Y","40Y") 

# convert age into period in metaData
metaData[,"period"] <- apply(metaData,1,function(x) periods[names(periods)==x["Age"]])

# Load data matrix - editing level measurements
dataMatrix <- read.table(gzfile("../data/dataMatrix.txt.gz"),sep="\t",header=TRUE,row.name=1)

# Get rid of rows with maximum value of zero
dataMatrix <- t(dataMatrix[apply(dataMatrix,1,function(x) max(x[which(!is.na(x))])>0),])

# use WGCNA function to filter data by missingness entries
if (!file.exists(paste0(outputDirData,"GSG.RData"))) {
	gsg = goodSamplesGenes(dataMatrix, minNGenes=1, verbose = 3);
	save(gsg, file=paste0(outputDirData,"GSG.RData"))
} else {
	load(paste0(outputDirData,"GSG.RData"))
}
dataMatrix_removeOutlier = dataMatrix[gsg$goodSamples, gsg$goodGenes]

# Filter metadata to only rows present in data matrix
metaData = metaData[rownames(metaData) %in% rownames(dataMatrix_removeOutlier),]

# convert region into Gross Region in metaData
metaData[,"GrossRegion"] <- apply(metaData,1,function(x) regions[names(regions)==x["Region"]])

# look for outliers using sample-sample connectivity
if (!file.exists(paste0(outputDirData,"zConnectivity.RData"))) {
	sdout <- 2; normadj <- (0.5+0.5*bicor(t(dataMatrix_removeOutlier), use='pairwise.complete.obs'))^2
	netsummary <- fundamentalNetworkConcepts(normadj); 
	K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
	C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
	outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
	print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(t(dataMatrix_removeOutlier))[outliers]); print(table(outliers))
	pdf(paste0(outputDir,"Outliers_ZConnect.pdf"))
	plot(Z.K, col = as.numeric(as.factor(metaData[,"GrossRegion"])), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
	abline(h=-2, lty=2)
	dev.off()
	save(Z.K, outliers, file=paste0(outputDirData,"zConnectivity.RData"))
} else {
	load(paste0(outputDirData,"zConnectivity.RData"))
}
dataMatrix_removeOutlier_keep <- dataMatrix_removeOutlier[!outliers,]
metaData = metaData[!outliers,]

# cluster with traits to check for associations
sampleTree2 = hclust(dist(dataMatrix_removeOutlier_keep), method = "average")
# filter metadata to kept samples
metaData_filt <- metaData[rownames(dataMatrix_removeOutlier_keep),c("Hemisphere","RIN","Dissectionscore","period","GrossRegion")]
# convert metaData to colors
metaData_filt_colors <- metaData_filt
metaData_filt_colors$RIN <- numbers2colors(metaData_filt_colors$RIN)
metaData_filt_colors$Dissectionscore <- numbers2colors(metaData_filt_colors$Dissectionscore)
metaData_filt_colors$period <- numbers2colors(metaData_filt_colors$period)
colorPalette <- c("blue","red","orange","green","pink","brown")
names(colorPalette) = c(1:6)
metaData_filt_colors$Hemisphere <- colorPalette[as.numeric(as.factor(metaData_filt$Hemisphere))]
metaData_filt_colors$GrossRegion <- colorPalette[as.numeric(as.factor(metaData_filt$GrossRegion))]
# Plot the sample dendrogram and the colors underneath.
pdf(file = paste0(outputDir,"sampleClustering_withTraits.pdf"), width = 12, height = 9)
plotDendroAndColors(sampleTree2, metaData_filt_colors,
                    groupLabels = names(metaData_filt_colors), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# soft-thresholding
datExpr <- dataMatrix_removeOutlier_keep
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
if (!file.exists(paste0(outputDirData,"signed.sft.thresh.RData"))) {
	# Call the network topology analysis function
	sft = pickSoftThreshold(datExpr, networkType="signed",corFnc="bicor", powerVector = powers, verbose = 5)
	save(sft,file=paste0(outputDirData,"signed.sft.thresh.RData"))
	# Plot the results:
	pdf(file = paste0(outputDir,"softThreshold.pdf"), width = 9, height = 5)
	par(mfrow = c(1,2));
	cex1 = 0.9;
	# Scale-free topology fit index as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     	xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    	main = paste("Scale independence"));
	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     	labels=powers,cex=cex1,col="red");
	# this line corresponds to using an R^2 cut-off of h
	abline(h=0.80,col="red")
	# Mean connectivity as a function of the soft-thresholding power
	plot(sft$fitIndices[,1], sft$fitIndices[,5],
     	xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     	main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
	dev.off()
} else {
	load(paste0(outputDirData,"signed.sft.thresh.RData"))
}

# TOM matrix
softPower = 9; # based on graph from above
if (!file.exists(paste0(outputDirData,"TOM_SP",softPower,".RData"))) {
	adjacency = adjacency(datExpr, power = softPower, type = "signed",corFnc="bicor");

	TOM = TOMsimilarity(adjacency);
	dissTOM = 1-TOM
	# Plot the resulting clustering tree (dendrogram)
	save(TOM,dissTOM,adjacency,softPower,file=paste0(outputDirData,"TOM_SP",softPower,".RData"))
} else {
	load(paste0(outputDirData,"TOM_SP",softPower,".RData"))
}

# Make the resulting clustering tree (dendrogram)
geneTree = flashClust(as.dist(dissTOM), method = "average");

# Relating dendrogram with traits
datTraits <- metaData[rownames(datExpr),c("Hemisphere","RIN","Dissectionscore","period","GrossRegion")]
traitmat=as.data.frame(cbind(as.factor(datTraits[,1]),as.numeric(datTraits[,2]),as.numeric(datTraits[,3]),as.numeric(datTraits[,4]),datTraits[,5]))# convert categorical variables in factor and numeric as numeric
rownames(traitmat)=rownames(datTraits)
colnames(traitmat)=c("Hemisphere","RIN","Dissectionscore","period","GrossRegion")

if (!file.exists(paste0(outputDirData,"DendroTraits_SP",softPower,".RData"))) {

	geneSigs=matrix(NA,nrow=10,ncol=ncol(datExpr)) # create a vector to hold the data
	for(i in 1:ncol(geneSigs)) {
		# calculate r correlation value for numeric variables
		# calculate adjusted R^2s square-root for categorical variables (factors)
		exprvec=as.numeric(datExpr[,i]) # get the expression vector for ith gene
		hemir=sqrt(max(summary(lm(exprvec~as.factor(traitmat[,1])))$adj.r.squared,0))
		rinr=bicor(traitmat[,2],exprvec,use="pairwise.complete.obs")
		dissectionr=bicor(traitmat[,3],exprvec,use="pairwise.complete.obs")
		periodr=bicor(traitmat[,4],exprvec,use="pairwise.complete.obs")
		
		NCX_r=sqrt(max(summary(lm(exprvec~as.numeric(traitmat[,5]=="NCX")))$adj.r.squared,0))
		AMY_r=sqrt(max(summary(lm(exprvec~as.numeric(traitmat[,5]=="AMY")))$adj.r.squared,0))
		STR_r=sqrt(max(summary(lm(exprvec~as.numeric(traitmat[,5]=="STR")))$adj.r.squared,0))
		HIP_r=sqrt(max(summary(lm(exprvec~as.numeric(traitmat[,5]=="HIP")))$adj.r.squared,0))
		CBC_r=sqrt(max(summary(lm(exprvec~as.numeric(traitmat[,5]=="CBC")))$adj.r.squared,0))
		MD_r=sqrt(max(summary(lm(exprvec~as.numeric(traitmat[,5]=="MD")))$adj.r.squared,0))
		regionr = max(NCX_r,AMY_r,STR_r,HIP_r,CBC_r,MD_r)
	
		geneSigs[,i]=c(hemir, rinr, dissectionr, periodr, NCX_r, AMY_r, STR_r, HIP_r, CBC_r, MD_r)
	}
	colnames(geneSigs)=colnames(datExpr)
	rownames(geneSigs)=c("Hemisphere","RIN","Dissectionscore","period","NCX","AMY","STR","HIP","CBC","MD")
	# convert to colors
	# For categorical variables we do not want values, thus lim=c(0,1), and signed and centered=F
	geneSigsColor=matrix(NA,nrow=nrow(geneSigs),ncol=ncol(datExpr)) # create a vector to hold the data
	geneSigsColor[1,] =numbers2colors(as.numeric(geneSigs[1,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
	geneSigsColor[2,] =numbers2colors(as.numeric(geneSigs[2,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1)) 
	geneSigsColor[3,] =numbers2colors(as.numeric(geneSigs[3,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
	geneSigsColor[4,] =numbers2colors(as.numeric(geneSigs[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
	geneSigsColor[5,] =numbers2colors(as.numeric(geneSigs[5,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1)) 
	geneSigsColor[6,] =numbers2colors(as.numeric(geneSigs[6,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
	geneSigsColor[7,] =numbers2colors(as.numeric(geneSigs[7,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
	geneSigsColor[8,] =numbers2colors(as.numeric(geneSigs[8,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
	geneSigsColor[9,] =numbers2colors(as.numeric(geneSigs[9,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
	geneSigsColor[10,] =numbers2colors(as.numeric(geneSigs[10,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1))
	rownames(geneSigsColor)=c("Hemisphere","RIN","Dissectionscore","period","NCX","AMY","STR","HIP","CBC","MD")
	colnames(geneSigsColor)=colnames(datExpr)
	# Try out tree cutting parameters
	mColorh <- mLabelh <- colorLabels <- NULL  
 	for (minModSize in c(20,50,100)) {
   		for (dthresh in c(0.1,0.2,0.25)) {
     		for (ds in c(2,4)) {
       			print("Trying parameters:")
       			print(c(minModSize,dthresh,ds))
       			tree = cutreeHybrid(dendro = geneTree, pamStage=FALSE, minClusterSize = minModSize, cutHeight = 0.9999, deepSplit = ds, distM = as.matrix(dissTOM))
       	 		merged <- mergeCloseModules(exprData = datExpr,colors = tree$labels,
                                   cutHeight = dthresh)
      	 		mColorh <- cbind(mColorh,labels2colors(merged$colors))
      	 		mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
     		}
   		}
	} 
	mColorh1=cbind(mColorh,geneSigsColor[1,],geneSigsColor[2,],geneSigsColor[3,],geneSigsColor[4,],geneSigsColor[5,],geneSigsColor[6,],geneSigsColor[7,],geneSigsColor[8,],geneSigsColor[9,],geneSigsColor[10,])
	mLabelh1=c(mLabelh,rownames(geneSigsColor)[c(1:10)])
	save(mColorh1,mLabelh1,geneSigs,geneSigsColor,file=paste0(outputDirData,"DendroTraits_SP",softPower,".RData"))
} else {
	load(paste0(outputDirData,"DendroTraits_SP",softPower,".RData"))
}

pdf(paste0(outputDir,"Signed_Dendro_parameters_SP",softPower,".pdf"),height=25,width=20)
plotDendroAndColors(geneTree,mColorh1,groupLabels=mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Different Module Cutting Parameters")
dev.off()

# make final cut based on parameters chosen
if (!file.exists(paste0(outputDirData,"modules.RData"))) {
	mms=50
	ds =4
	dthresh=0.1
	tree = cutreeHybrid(dendro = geneTree, pamStage=F, minClusterSize =mms, cutHeight = 0.9999, deepSplit = ds, distM = as.matrix(dissTOM))

	merged <- mergeCloseModules(exprData = datExpr,colors = tree$labels, cutHeight = dthresh)
	mColorh <- cbind(labels2colors(merged$colors),geneSigsColor[1,],geneSigsColor[2,],geneSigsColor[3,],geneSigsColor[4,],geneSigsColor[5,],geneSigsColor[6,],geneSigsColor[7,],geneSigsColor[8,],geneSigsColor[9,],geneSigsColor[10,])
	mLabelh <- c("Merged Colors",rownames(geneSigsColor)[c(1:10)])

	pdf(paste0(outputDir,"Signed_Dendro_parameters_SP",softPower,"_Final.pdf"),height=10,width=16)
	plotDendroAndColors(geneTree, mColorh, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power = ",softPower,"mms=",mms,"ds=",ds,"dthresh=",dthresh));
	dev.off()

	mergedColors = labels2colors(merged$colors);

	# Eigengenes of the new merged modules:
	MEList=moduleEigengenes(datExpr, colors = mergedColors,softPower= softPower, nPC=1)
	MEs=MEList$eigengenes
	MEs=orderMEs(MEs)
	rownames(MEs) = rownames(datExpr)
	moduleColors = mergedColors
	names(moduleColors) <- colnames(datExpr)
	save(geneTree,moduleColors,MEs,file=paste0(outputDirData,"modules.RData"))
} else {
	load(paste0(outputDirData,"modules.RData"))
}

# For modules - correlation module eigengenes with traits
Pmat <- matrix(NA,nrow=ncol(MEs),ncol=10)
colnames(Pmat) = c("Hemisphere","RIN","DissectionScore","Period","NCX","AMY","STR","HIP","CBC","MD")
rownames(Pmat) = colnames(MEs)
Bmat <- matrix(NA,nrow=ncol(MEs),ncol=10)
for (i in 1:ncol(MEs)) {
	hemi_lm = lm(MEs[,i]~as.factor(metaData[rownames(MEs),"Hemisphere"]))
	rin_lm = lm(MEs[,i]~metaData[rownames(MEs),"RIN"])
	dissection_lm = lm(MEs[,i]~metaData[rownames(MEs),"Dissectionscore"])
	period_lm = lm(MEs[,i]~metaData[rownames(MEs),"period"])
		
	NCX_lm = lm(MEs[,i]~as.numeric(metaData[rownames(MEs),"GrossRegion"] == "NCX"))
	AMY_lm = lm(MEs[,i]~as.numeric(metaData[rownames(MEs),"GrossRegion"] == "AMY"))
	STR_lm = lm(MEs[,i]~as.numeric(metaData[rownames(MEs),"GrossRegion"] == "STR"))
	HIP_lm = lm(MEs[,i]~as.numeric(metaData[rownames(MEs),"GrossRegion"] == "HIP"))
	CBC_lm = lm(MEs[,i]~as.numeric(metaData[rownames(MEs),"GrossRegion"] == "CBC"))
	MD_lm = lm(MEs[,i]~as.numeric(metaData[rownames(MEs),"GrossRegion"] == "MD"))
	
	Pmat[i,] = c(summary(hemi_lm)$coefficients[2,4],summary(rin_lm)$coefficients[2,4],summary(dissection_lm)$coefficients[2,4],summary(period_lm)$coefficients[2,4],summary(NCX_lm)$coefficients[2,4],summary(AMY_lm)$coefficients[2,4],summary(STR_lm)$coefficients[2,4],summary(HIP_lm)$coefficients[2,4],summary(CBC_lm)$coefficients[2,4],summary(MD_lm)$coefficients[2,4])
	Bmat[i,] = c(summary(hemi_lm)$coefficients[2,1],summary(rin_lm)$coefficients[2,1],summary(dissection_lm)$coefficients[2,1],summary(period_lm)$coefficients[2,1],summary(NCX_lm)$coefficients[2,1],summary(AMY_lm)$coefficients[2,1],summary(STR_lm)$coefficients[2,1],summary(HIP_lm)$coefficients[2,1],summary(CBC_lm)$coefficients[2,1],summary(MD_lm)$coefficients[2,1])	
}

pdf(paste0(outputDir,"MEtraits.pdf"),height=10,width=8)
## Plot the module-level summaries
rmMod = colnames(MEs) %in% "MEgrey" # remove grey module

BFmat <- Pmat*(nrow(Pmat)-1) # Bonferroni correction of P-values
BFmat[BFmat>1] <- 1

textBFmat <- signif(BFmat,2)
textBFmat[BFmat > 0.05] <- ""
dispMat <- -log10(Pmat)*sign(Bmat)
dispMat[BFmat > 0.05] <- 0
dispMat <- dispMat[!rmMod,]
textBFmat <- textBFmat[!rmMod,]
dispMat[dispMat > 10] <- 10
dispMat[dispMat < -10] <- -10
rownames(dispMat) <- colnames(MEs)[!rmMod]

maxval <- 10
labeledHeatmap(Matrix=dispMat,
               textMatrix=textBFmat,
               yLabels=colnames(MEs)[!rmMod],
               yColorLabels=TRUE,
               ySymbols=gsub("ME","M.",colnames(MEs)[!rmMod]),
               xLabels=colnames(Pmat),
               colors=blueWhiteRed(100),
               cex.lab.x=1,
               cex.lab.y=0.7,
               cex.text=0.5,
               zlim=c(-maxval,maxval),
               main="ME - trait associations")
dev.off()
