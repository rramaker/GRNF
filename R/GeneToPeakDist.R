#'Find the nearest ChIP-seq binding site to each gene in a GTF file with chromatin weighting
#' @param ChIP A dataframe containing ChIP-seq peak information in bed format. At least three columns indicating chromosome, start position, and stop position for each binding site.
#' @param GTF A dataframe containing gene location information in Gene Transfer Format (GTF).
#' Function expects nine columns with chromosome in column 1, gene start and stop postion in columns 4 and 5, strand information in column 6, and gene ID information in column 9.
#' Ensure only rows corresponding to unique genes are included and chromosome column is formatted identically to ChIP-seq file.
#' @param Genes A vector of gene IDs to which the package will calculate the nearest ChiP-seq binding distance. If NULL, all genes present in the provided GTF file will be used. Defaults to NULL.
#' @param TAD A dataframe containing TAD boundaries with at least 3 columns containing chromosome, start position, and stop position of each TAD. Defaults to NULL.
#' @param TAD_Penalty A numeric indicating the distance penalty to put on binding sites that fall outside a genes TAD. Only used if a TAD boundary File is provided. Defaults to 100
#' @param PCHC A dataframe containing promoter looping information. The function expects a 5 column data frame with the bait gene ENSEMBL ID, the capture chromsome, capture start position, capture stop position, and an interaction frequencing metric.
#' Ensure the capture chromosome column is formatted identically to the ChIP file. Defaults to NULL
#' @param PCHC_Bonus A numeric indicating the reward to provide to ChIP peaks falling in regions that loop into a gene's promoter. Only used if a PCHC file is provided. Defaults to 100.
#' @param numCores A numeric indicating the number of cores the function can use via R's parallel package. Defaults to 4.
#' @export
#' @import parallel
#' @import GenomicRanges
#' @import IRanges
#' @return Returns a dataframe with two columns. The first is the ENSEMBL gene ID and the second is the distance to the nearest ChIP-seq binding site.
#' @examples
#'data("GM12878_BATF_ChIP")
#'data("Homo_sapiens.GRCh37.82.chr.gtf")
#'DistanceFrame<-GeneToPeakDist(ChIP = ChIP, Genes = c("ENSG00000186092", "ENSG00000237683", "ENSG00000235249"))
#'
#' @note
#' #Citations:
#'

GeneToPeakDist<-function(ChIP, GTF, Genes=NULL, TAD=NULL, TAD_Penalty=100, PCHC=NULL, PCHC_Bonus=100, numCores=4){

  #Parse ENSEMBL IDs from GTF
  GTF[,"Gene"]<-substr(gsub("^.+ENSG","ENSG",GTF[,9]),1,15)

  if(is.null(Genes)){
    Genes<-GTF$Gene
  }

  #Find all genes of interest in GTF file
  allGeneLoc<-GTF[which(GTF$Gene%in%Genes),]

  #Create empty list to store all gene-peak distances
  allGeneList<-list()

  #Loop through each chromosome of chip bed file
  for(chromosome in unique(ChIP[,1])){

    #Find all positive strand genes
    allGeneLocPos<-allGeneLoc[which(allGeneLoc[,7]=="+"),]

    #Find all genes on current chromosome
    currentChrGenes<-allGeneLocPos[which(allGeneLocPos$V1==chromosome),]

    if(sum(allGeneLocPos$V1==chromosome)>0){

      #Find all peaks on current chromosome
      currentChrChip<-ChIP[which(ChIP[,1]==chromosome),]

      if(!is.null(TAD)){

        #Find TADS on current chromosome
        currentChrTAD<-TAD[which(TAD[,1]==chromosome),]

        #Calculate distance from gene to nearest peak (start and stop) for each gene on chromosome
        currentGeneFrame<-data.frame("LeftChip"=unlist(parallel::mclapply(currentChrGenes[,4], function(x) min((currentChrChip[,2]-x)[which((currentChrChip[,2]-x)>(currentChrChip[,2]-currentChrChip[,3]))]),mc.cores = numCores)),
                                     "LeftTAD"=unlist(parallel::mclapply(currentChrGenes[,4], function(x) min((c(currentChrTAD[,3],currentChrTAD[,2])-x)[which((c(currentChrTAD[,3],currentChrTAD[,2])-x)>0)]),mc.cores = numCores)),
                                     "RightChip"=unlist(parallel::mclapply(currentChrGenes[,4], function(x) max((currentChrChip[,3]-x)[which((currentChrChip[,3]-x)<0)]),mc.cores = numCores)),
                                     "RightTAD"=unlist(parallel::mclapply(currentChrGenes[,4], function(x) max((c(currentChrTAD[,3],currentChrTAD[,2])-x)[which((c(currentChrTAD[,3],currentChrTAD[,2])-x)<0)]),mc.cores = numCores)))
        #Convert hanging dists to NA
        for(i in 1:4){currentGeneFrame[grep("Inf|-Inf",currentGeneFrame[,i]),i]<-NA}
        #Convert overlapping TSSs to zero
        currentGeneFrame[which(currentGeneFrame[,1]<0),]<-0

        #Add TAD boundary if applicable
        currentGeneFrame[which(!is.na(currentGeneFrame[,1])&currentGeneFrame[,1]>currentGeneFrame[,2]),1]<-currentGeneFrame[which(!is.na(currentGeneFrame[,1])&currentGeneFrame[,1]>currentGeneFrame[,2]),1]*TAD_Penalty
        currentGeneFrame[which(!is.na(currentGeneFrame[,3])&currentGeneFrame[,3]<currentGeneFrame[,4]),3]<-currentGeneFrame[which(!is.na(currentGeneFrame[,3])&currentGeneFrame[,3]<currentGeneFrame[,4]),3]*TAD_Penalty
        #Find closest binding site post-TAD
        posGeneFrame<-cbind(currentChrGenes[,"Gene"],apply(currentGeneFrame[,c(1,3)],1,function(x) min(abs(x),na.rm=T)))
      }
      if(is.null(TAD)){
        #Calculate distance from gene to nearest peak (start and stop) for each gene on chromosome
        currentGeneFrame<-data.frame("LeftChip"=unlist(parallel::mclapply(currentChrGenes[,4], function(x) min((currentChrChip[,2]-x)[which((currentChrChip[,2]-x)>(currentChrChip[,2]-currentChrChip[,3]))]),mc.cores = numCores)),
                                     "RightChip"=unlist(parallel::mclapply(currentChrGenes[,4], function(x) max((currentChrChip[,3]-x)[which((currentChrChip[,3]-x)<0)]),mc.cores = numCores)))
        #Convert hanging dists to NA
        for(i in 1:2){currentGeneFrame[grep("Inf|-Inf",currentGeneFrame[,i]),i]<-NA}
        #Convert overlapping TSSs to zero
        currentGeneFrame[which(currentGeneFrame[,1]<0),]<-0
        #Convert peak distances to gene frame and find minimum distance (peak start or stop)
        posGeneFrame<-cbind(currentChrGenes[,"Gene"],apply(currentGeneFrame[,c(1,2)],1,function(x) min(abs(x),na.rm=T)))
      }
      #Add to master distance list
      colnames(posGeneFrame)<-c(1:2)

      #Add to master distance list
      allGeneList[[paste0(chromosome,"+")]]<-posGeneFrame
    }

    #Find all negative strand genes
    allGeneLocNeg<-allGeneLoc[which(allGeneLoc$V7=="-"),]

    #Find all genes on current chromosome
    currentChrGenes<-allGeneLocNeg[which(allGeneLocNeg$V1==chromosome),]
    if(sum(allGeneLocNeg$V1==chromosome)>0){

      #Find all peaks on current chromosome
      currentChrChip<-ChIP[which(ChIP[,1]==chromosome),]
      if(!is.null(TAD)){
        #Find TADS on current chromosome
        currentChrTAD<-TAD[which(TAD[,1]==chromosome),]

        #Calculate distance from gene to nearest peak (start and stop) for each gene on chromosome
        currentGeneFrame<-data.frame("LeftChip"=unlist(parallel::mclapply(currentChrGenes[,5], function(x) min((currentChrChip[,2]-x)[which((currentChrChip[,2]-x)>(currentChrChip[,2]-currentChrChip[,3]))]),mc.cores = numCores)),
                                     "LeftTAD"=unlist(parallel::mclapply(currentChrGenes[,5], function(x) min((c(currentChrTAD[,3],currentChrTAD[,2])-x)[which((c(currentChrTAD[,3],currentChrTAD[,2])-x)>0)]),mc.cores = numCores)),
                                     "RightChip"=unlist(parallel::mclapply(currentChrGenes[,5], function(x) max((currentChrChip[,3]-x)[which((currentChrChip[,3]-x)<0)]),mc.cores = numCores)),
                                     "RightTAD"=unlist(parallel::mclapply(currentChrGenes[,5], function(x) max((c(currentChrTAD[,3],currentChrTAD[,2])-x)[which((c(currentChrTAD[,3],currentChrTAD[,2])-x)<0)]),mc.cores = numCores)))
        #Convert hanging dists to NA
        for(i in 1:4){currentGeneFrame[grep("Inf|-Inf",currentGeneFrame[,i]),i]<-NA}
        #Convert overlapping TSSs to zero
        currentGeneFrame[which(currentGeneFrame[,1]<0),]<-0

        #Add TAD boundary if applicable
        currentGeneFrame[which(!is.na(currentGeneFrame[,1])&currentGeneFrame[,1]>currentGeneFrame[,2]),1]<-currentGeneFrame[which(!is.na(currentGeneFrame[,1])&currentGeneFrame[,1]>currentGeneFrame[,2]),1]*TAD_Penalty
        currentGeneFrame[which(!is.na(currentGeneFrame[,3])&currentGeneFrame[,3]<currentGeneFrame[,4]),3]<-currentGeneFrame[which(!is.na(currentGeneFrame[,3])&currentGeneFrame[,3]<currentGeneFrame[,4]),3]*TAD_Penalty
        #Find closest binding site post-TAD
        negGeneFrame<-cbind(currentChrGenes[,"Gene"],apply(currentGeneFrame[,c(1,3)],1,function(x) min(abs(x),na.rm=T)))
      }
      if(is.null(TAD)){
        #Calculate distance from gene to nearest peak (start and stop) for each gene on chromosome
        #Calculate distance from gene to nearest peak (start and stop) for each gene on chromosome
        currentGeneFrame<-data.frame("LeftChip"=unlist(parallel::mclapply(currentChrGenes[,5], function(x) min((currentChrChip[,2]-x)[which((currentChrChip[,2]-x)>(currentChrChip[,2]-currentChrChip[,3]))]),mc.cores = numCores)),
                                     "RightChip"=unlist(parallel::mclapply(currentChrGenes[,5], function(x) max((currentChrChip[,3]-x)[which((currentChrChip[,3]-x)<0)]),mc.cores = numCores)))
        #Convert hanging dists to NA
        for(i in 1:2){currentGeneFrame[grep("Inf|-Inf",currentGeneFrame[,i]),i]<-NA}
        #Convert overlapping TSSs to zero
        currentGeneFrame[which(currentGeneFrame[,1]<0),]<-0
        #Find closest binding site post-TAD
        negGeneFrame<-cbind(currentChrGenes[,"Gene"],apply(currentGeneFrame[,c(1,2)],1,function(x) min(abs(x),na.rm=T)))
      }
      #Add to master distance list
      colnames(negGeneFrame)<-c(1:2)

      #Add to master distance list
      allGeneList[[paste0(chromosome,"-")]]<-negGeneFrame
    }
  }

  #Bind distance into data frame add column names and return
  allGeneFrame<-data.frame(do.call(rbind, allGeneList))
  colnames(allGeneFrame)<-c("GeneID","Distance (bp)")
  allGeneFrame[,2]<-as.numeric(as.character(allGeneFrame[,2]))

  if(!is.null(PCHC)){
    #Find chip sites that fall in a PCHC region
    PCHC_bed <- GenomicRanges::GRanges(PCHC[,1], IRanges::IRanges(PCHC[,2], PCHC[,3]), strand=rep("+",nrow(PCHC)), PCHC[,4], id=PCHC[,5])
    ChIP_bed <- GenomicRanges::GRanges(ChIP[,1], IRanges::IRanges(ChIP[,2], ChIP[,3]), strand=rep("+",nrow(ChIP)), id=row.names(ChIP))
    bedIntersect <- GenomicRanges::findOverlaps(PCHC_bed,ChIP_bed)
    overlapFrame<-cbind(ChIP[bedIntersect@subjectHits,1:3],PCHC[bedIntersect@queryHits,])
    #Compute loop bonus distance
    overlapFrame[,8]<-as.character(overlapFrame[,8])
    overlapFrame<-overlapFrame[which(overlapFrame[,8]%in%allGeneFrame[,1]),]
    row.names(GTF)<-GTF$Gene
    GTF_PCHC_Filt<-GTF[overlapFrame[,8],]
    overlapFrame[,9]<-rep(NA,nrow(overlapFrame))
    overlapFrame[which(GTF_PCHC_Filt[,7]=="+"),9]<-apply(cbind(overlapFrame[which(GTF_PCHC_Filt[,7]=="+"),2]-GTF_PCHC_Filt[which(GTF_PCHC_Filt[,7]=="+"),4],
                                                               overlapFrame[which(GTF_PCHC_Filt[,7]=="+"),3]-GTF_PCHC_Filt[which(GTF_PCHC_Filt[,7]=="+"),4]),1,function(x) min(abs(x)))

    overlapFrame[which(GTF_PCHC_Filt[,7]=="-"),9]<-apply(cbind(overlapFrame[which(GTF_PCHC_Filt[,7]=="-"),2]-GTF_PCHC_Filt[which(GTF_PCHC_Filt[,7]=="-"),5],
                                                               overlapFrame[which(GTF_PCHC_Filt[,7]=="-"),3]-GTF_PCHC_Filt[which(GTF_PCHC_Filt[,7]=="-"),5]),1,function(x) min(abs(x)))
    overlapFrame[,9]<-overlapFrame[,9]/(1+overlapFrame[,7]*PCHC_Bonus)
    #Find shortest distance with loop bonus
    overlapFrame<-overlapFrame[order(overlapFrame[,9]),]
    overlapFrame<-overlapFrame[which(!duplicated(overlapFrame[,8])),]
    #Convert previously computed distance to loop bonus distance if smaller
    row.names(allGeneFrame)<-allGeneFrame[,1]
    allGeneFrame_PCHC_Filt<-allGeneFrame[overlapFrame[,8],]
    allGeneFrame_PCHC_Filt[which(allGeneFrame_PCHC_Filt[,2]>overlapFrame[,9]),2]<-overlapFrame[which(allGeneFrame_PCHC_Filt[,2]>overlapFrame[,9]),9]
    allGeneFrame<-rbind(allGeneFrame[which(!allGeneFrame[,1]%in%overlapFrame[,8]),],allGeneFrame_PCHC_Filt)
    row.names(allGeneFrame)<-1:nrow(allGeneFrame)
  }
  return(allGeneFrame)
}
