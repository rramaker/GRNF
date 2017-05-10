#'Find Reactome pathways enriched from proximal binding
#' @param distanceFrame A dataframe resulting from the output of the GeneToPeakDist function. Should have two columns with the first having ENSEMBL gene IDs and the second the distance to the nearest ChIP peak.
#' @param minPathwayGenes A numeric indicating the minimum number of genes in a pathway required for inclusion in pathway analysis.
#' @param numCores A numeric indicating the number of cores the function can use via R's parallel package. Defaults to 4.
#' @export
#' @import parallel
#' @return Returns a vector of one-sided Wilcoxon p-values with the alternate hypothesis of pathway genes enriched for shorter distances to the nearest ChIP-peak relative to all non-pathway genes.
#' @note
#' Reactome
#'
#' Terms of Use: http://www.reactome.org/pages/about/license-agreement/
#'


FindPathwayEnrichment<-function(distanceFrame, minPathwayGenes=5, numCores=4){
  data("ENSEMBL_REACT_Paths")
  Pathways<-Pathways[which(row.names(Pathways)%in%distanceFrame[,1]),]
  Pathway_Filt<-as.list(Pathways[,which(colSums(Pathways,na.rm=T)>=minPathwayGenes)])
  GeneNames<-row.names(Pathways[,which(colSums(Pathways,na.rm=T)>=minPathwayGenes)])
  return(unlist(mclapply(Pathway_Filt,function(Path) wilcox.test(distanceFrame[which(distanceFrame[,1]%in%GeneNames[which(Path==1)]),2], distanceFrame[which(!distanceFrame[,1]%in%GeneNames[which(Path==1)]),2], alternative="less")$p.value,mc.cores = numCores)))
}
