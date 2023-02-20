### example of doing enrichment test on SomaScan data
### User assign these 4 variables
myPathIn1 <- "where is your protein filtering file"
myPathIn2 <- "where is your protein expression data"
myPathIn3 <- "where is your Reactom downloaded file" 
databaseFile <- "downloaded file"
ranksOriginal <- "your ranked metrics, with seq names"
### Then run the following, and View(enrichedPath), which should be your final results for Enrichment test based on Reactom

###____________________________________________________________________________________________________________________________________________________
### main body
###____________________________________________________________________________________________________________________________________________________
source("EnrichmentTest.SomaScan.R")

RFUstart ="seq.10000.28"
### clear up protein meta data table
### For the purpose of extract information for multiple genes encoding complex proteins, tidy up the EntrezGeneID, EntrezGeneSymbol as well as UniProt which have multiple spaces
proFil <-  read.csv(paste0(myPathIn1,'somascan/QC/patch1/Protein_filters_metadata.csv'),row.names = 1)

for(Counter in 1:nrow(proFil)){
  proFil$UniProt[Counter] <- gsub("\\s+","\\|",proFil$UniProt[Counter]) ### remove space only leave "|" as delimiter
  proFil$EntrezGeneID[Counter] <- gsub("\\s+","\\|",proFil$EntrezGeneID[Counter]) ### merge multiple spaces as one single space
  proFil$EntrezGeneSymbol[Counter] <- gsub("\\s+","\\|",proFil$EntrezGeneSymbol[Counter])
}

# ### find proteins with multiple EntrezGenes; some multiple EntrezGeneID have the same EntrezGeneSymbol
dupRowIdA <- grep("\\|",proFil$EntrezGeneID)

# ### a data frame showing proteins with multiple EntrezGenes. proComplexRef is used to replace element in the reference gene sets.
proComplexRef = proFil[dupRowIdA,c("SomaId","EntrezGeneID","EntrezGeneSymbol","UniProt")]

### resolve the problem of duplicated genes
### multiple Somamer match one gene. Mean paired correlation coefficients table
SecondarySet1Raw <- read.csv(paste0(myPathIn2,"somascan/csv/SomaScan_rel2_secondary1.csv"),row.names = 2)
expDat1Raw <- SecondarySet1Raw[,which(colnames(SecondarySet1Raw)==RFUstart):ncol(SecondarySet1Raw)]
multiSomaSameGeneCorTable1 <- GetMultiSomaSameGeneCorTable(ProMeta,expDat1Raw,"EntrezGeneSymbol")
multiSomaSameGeneCorTable2 <- GetMultiSomaSameGeneCorTable(ProMeta,expDat1Raw,"EntrezGeneID")

### select the Representative Somamer for the gene.
# ranksX <- ReplaceRanksX(ProMeta,ranksOriginal,"EntrezGeneID",multiSomaSameGeneCorTable2)
ranksX <- ReplaceRanksX(ProMeta,ranksOriginal,"EntrezGeneSymbol",multiSomaSameGeneCorTable1)

### prepare the Reactom genesets, replace the protein complexies with SomaID
GeneSets <- getPath(myPathIn3,databaseFile,proComplexRef,"EntrezGeneSymbol")

### Perform fgsea enrichment test
enrichedPath <- getEnrichPath(GeneSets,ranksX)
