### This script includes enrichment test functions specialized for SomaScan data, which tackle with two extra issues beyond normal enrichment testing.
### 1. Protein complexes -- one Somamer correpsonds to different genes.
### 2. Multiple somamers correspond to the same gene.

###____________________________________________________________________________________________________________________________________________________
### function section
###____________________________________________________________________________________________________________________________________________________
### function for fgsea pathway enrichment analysis
### input: pathBaseFile-- pathway database local file; myPathIn -- file path directs to pathway databse file; proComplexRef: protein complexy frame
### output: all_gene_sets -- pathway gene sets, after protein complexes replacement.
getPath <- function(myPathIn,pathBaseFile,proComplexRef,whatEntrez){
  ### read in pathway gene sets from database
  all_gene_setsPre <- fgsea::gmtPathways(paste0(myPathIn,pathBaseFile))
  ### change EntrezGeneIDs for complex proteins to unique SomaIds
  all_gene_sets <- replaceGenes(all_gene_setsPre,proComplexRef,whatEntrez) ### only use GeneID to anonate genes)
  return(all_gene_sets)
}

# function to tackle with protein complexies
# input: all_gene_setsPre -- gene sets before replacing protein complexies; proComplexRef -- table showing protein complexies and corresponding protein meta data; whichCriteria -- indicate use EntrezID or Symbol for replacement in the gene sets.
# output: all_gene_sets -- gene sets after replacing protein complexies with consistent somaID.
replaceGenes <- function(all_gene_setsPre,proComplexRef,whichCriteria){
  ### replace genes which encode complex proteins in the pathways gene sets with the unique SomaId
  all_gene_sets = all_gene_setsPre
  setK = 1
  
  ### pre-genesets may include entrez names or entrez ids. define which entrez item type is in the pre genesets.
  criteriaColumn = which(colnames(proComplexRef)==whichCriteria)
  
  ### replace protein complexies with the consensus symbol.
  for (geneSet in all_gene_setsPre){
    for(complexPro in 1:nrow(proComplexRef)){
      toReplace = strsplit(proComplexRef[,criteriaColumn][complexPro],"\\|")[[1]]
      replaceGene <- unlist(sapply(toReplace,function(x){which(geneSet==x)})) ### define the positions of the protein complixies in the current geneSet
      if(length(replaceGene)!=0){
        geneSet[replaceGene] = rep(proComplexRef$SomaId[complexPro],length(replaceGene)) ### if current gene set include protein complexies, replace all of them with the unique label (here our case is somaID)
        all_gene_sets[[setK]] = unique(geneSet) ### only keep the unique entries within the gene set
      }else{next}
    }
    setK = setK +1}
  
  return(all_gene_sets)
}

### function to generate a table, where multiple somamer per EntrezGeneSymbol/EntrezGeneID
# input: ProMeta -- protein meta table; expDat -- protein intensity data; whichEntrez -- indicate use EntrezID or Symbol for replacement in the gene sets.
# output: ; multiSomaSameGeneCorTable -- correlation table indicate Somamer representative per gene.
GetMultiSomaSameGeneCorTable <- function(ProMeta,expDat,whichEntrez){
  dupID <- which(duplicated(ProMeta[,whichEntrez]))
  uniDup <- unique(ProMeta[dupID,whichEntrez]) 
  uniDup <- uniDup[which(uniDup!="")]
  
  multiSomaSameGeneCorTable = vector(mode="list",length=length(uniDup))
  names(multiSomaSameGeneCorTable) = uniDup 
  
  for(SomaGroupC in uniDup){
    whichSomas <- rownames(ProMeta)[which(ProMeta[,whichEntrez]==SomaGroupC)]
    corVector <- apply(cor(expDat[,whichSomas]),1,mean)
    names(corVector) <- whichSomas
    multiSomaSameGeneCorTable[[SomaGroupC]] <- corVector
  }
  return(multiSomaSameGeneCorTable)
}

### function to select which seqID corresponding the same gene for enrichment test
# input: ProMeta -- protein meta table; rankPro -- ranking metric for proteins; whichCriteria -- indicate use EntrezID or Symbol for replacement in the gene sets; multiSomaSameGeneCorTable -- correlation table indicate Somamer representative per gene.
# output: ranksX.uniqueGene -- ranked gene with unique gene name
ReplaceRanksX <- function(ProMeta,rankPro,whichEntrez,multiSomaSameGeneCorTable){
  ranksX <- rankPro
  names(ranksX) <- ProMeta[names(rankPro),whichEntrez]
  
  if(any(duplicated(names(ranksX)))){
    EntrezReplace <- vector("numeric",length=length(rankPro))
    for(EntrezC in 1:length(names(ranksX))){
      Entrez <- names(ranksX)[EntrezC]
      somaSeq <-names(rankPro)[which(names(ranksX)==Entrez)]
      if(Entrez %in% names(multiSomaSameGeneCorTable)){
        EntrezReplace[EntrezC] <- names(sort(multiSomaSameGeneCorTable[[Entrez]][somaSeq],decreasing=TRUE))[1]
      }else{EntrezReplace[EntrezC] <- names(rankPro)[EntrezC]}
    }
    ranksX.uniqueGene <- rankPro[unique(EntrezReplace)]
    names(ranksX.uniqueGene) <- ProMeta[names(rankPro[unique(EntrezReplace)]),whichEntrez]
  }else{
    ranksX.uniqueGene <- ranksX
  }
  
  ranksX.uniqueGene <- ranksX.uniqueGene[which(names(ranksX.uniqueGene)!="")]
  
  return(ranksX.uniqueGene)
}

### function for fgsea pathway enrichment analysis
### input: all_gene_sets -- reference genesets; ranks1 -- gene list ranked by p values derived from differential test
### output: list, fgseaPath-- fgsea object of pathway enrichment test
getEnrichPath <- function(all_gene_sets,ranks){
  
  ### define which scoreType to use
  if(sign(min(ranks))*sign(max(ranks)) == -1){whichType="std"
  }else if(min(ranks)>0){whichType = "pos"
  }else{whichType = "neg"}
  
  # fgseaPath <- fgsea(pathways = all_gene_sets, stats=ranks, scoreType = whichType,eps = 0.0,minSize=10, maxSize=500)
  fgseaPath <- fgsea::fgsea(pathways = all_gene_sets, stats=ranks,scoreType = whichType, eps = 0.0) ### we will not limit the size (minSize and maxSize) here, just output everything, but would select the proper size for further plotting. 
  fgseaPath2 <- fgseaPath[order(fgseaPath$pval),]
  
  return(fgseaPath2)
}