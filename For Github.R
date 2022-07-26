library(dplyr)
library(tidyverse)
library(data.table)
library(gplots)
library(GenomicRanges)
library(GenomicFeatures)
library(BSgenome.Celegans.UCSC.ce11)
library(Biostrings)
library(Biobase) 

#Read in Data from Tintori et al. Developmental Cell 2016
Goldstein_data <- read.table("TableS2_RPKMs_Goldstein_singleCell2016.txt",header=T,row.names=1,stringsAsFactors = FALSE)
Goldstein_data_1 <- Goldstein_data %>%
  select(-contains("tossed"))
colnames(Goldstein_data_1) <- sub("_.*", "", colnames(Goldstein_data_1))

#convert Goldstein gene names to ce11

ce10 = read.table("c_elegans.WS220.geneIDs.txt",sep=",",stringsAsFactors = FALSE)
WS260 = read.table("c_elegans.PRJNA13758.WS260.geneIDs.txt",sep=",",stringsAsFactors = FALSE)

ce10[,4] <-  ce10[,2]
ce10[!nzchar(as.character(ce10[,2])),4] <- ce10[!nzchar(as.character(ce10[,2])),3]
ce10 <- ce10[ce10[,4]%in%rownames(Goldstein_data_1),]
rownames(ce10) <-  ce10[,4]

GoldsteinTNames <- rownames(Goldstein_data_1)
names(GoldsteinTNames) <- GoldsteinTNames
GoldsteinTNames[intersect(GoldsteinTNames,rownames(ce10))] <- ce10[intersect(GoldsteinTNames,rownames(ce10)),1]
names(GoldsteinTNames) <- GoldsteinTNames

rownames(WS260) <- WS260[,2]
WS260[,5] <- WS260[,3]
WS260[!nzchar(WS260[,3]),5] <- WS260[!nzchar(WS260[,3]),4]

GoldsteinTNames[intersect(GoldsteinTNames,rownames(WS260))] <- WS260[intersect(GoldsteinTNames,rownames(WS260)),5]
GoldsteinTNames[c("1-Apr","1-Jun","6-Mar","1-Oct","2-Oct","1-Sep","T05H4.6a","wars-2")] <- c("apr-1","jun-1","mar-6","oct-1","oct-2","sep-1","erfa-1","prx-10")
GoldsteinTNames[GoldsteinTNames == ""] <- names(GoldsteinTNames[GoldsteinTNames == ""])
rownames(Goldstein_data_1) <- GoldsteinTNames

#Calculate median TPMs.
cellMedians <- t(apply(Goldstein_data_1, 1, function(x) tapply(x, colnames(Goldstein_data_1), median)))

#Get Parent Function 
GetParent <- function(x){
  if(x=="P0") return("P0")
  if(is.element(substr(x,nchar(x),nchar(x)),c("a","p","d","v","l","r","x"))) return (substr(x,1,nchar(x)-1))
  if(x=="MSx1" || x=="MSx2") return("MS")
  if(x=="Cx1" || x=="Cx2") return("C")
  if(x=="AB" || x=="P1") return("P0")
  if(x=="EMS" || x=="P2") return("P1")
  if(x=="E" || x=="MS") return("EMS")
  if(x=="C" || x=="P3") return("P2")
  if(x=="D" || x=="P4") return("P3")
  if(x=="Z2" || x=="Z3") return("P4")
  
}

cellMedians_2 <- cellMedians[,!is.na(sapply(colnames(cellMedians),GetParent))]

#Calculate Fold Change (FC) using Pseudocount (PC) of 10
PC = 10 
FC <-  (cellMedians_2+PC) / (cellMedians[,sapply(colnames(cellMedians_2),GetParent)]+PC)

#Calculate daughter minus parent (DP Change)
DP_change <- cellMedians_2 -cellMedians[,sapply(colnames(cellMedians_2),GetParent)]

#Calculate Absolute Change (AC) reading in volumes and adjusting to 8000000 transcripts/embryo 
#Convert to dataframe and rename columns 

AC <- Change_vol* EmbryoTranscripts / 1000000 
AC_DF <- as.data.frame(AC) 
setDT(AC_DF, keep.rownames = TRUE)
setnames(AC_DF, "rn", "Genes")

FC_DF <- as.data.frame(FC) 
setDT(FC_DF, keep.rownames = TRUE)
setnames(FC_DF, "rn", "Genes")

AC_2 <- AC_DF %>% pivot_longer (cols = -Genes)
setnames(AC_2, "value", "AC")
FC_2 <- FC_DF %>% pivot_longer (cols = -Genes)
setnames(FC_2, "value", "FC")
FCvalue <- dplyr::select(FC_2 , cols = "FC")
pdc <- add_column(AC_2, FCvalue, before = NULL, after = NULL)

setnames(pdc, "cols", "FC")
setnames(pdc, "name", "Cell")

# Read in information about Parents and join 
cell_relations_2 <- readxl::read_xlsx("Cell_relations_2.xlsx") %>% as.data.table
pdc_2 <- pdc %>% left_join(cell_relations_2)

#remove maternal
P0_medians <- cellMedians[,"P0"]
P0_amt <- 8*(P0_medians)
maternal <- P0_amt %>% as.data.frame() %>% dplyr::filter(. > 100)
setDT(maternal, keep.rownames = TRUE)
maternal_genes <- maternal$rn %>% unique()

'%ni%' <- Negate('%in%')

#Set thresholds for Low, Medium and High rate genes
#Low genes are never high in any cell

pdc_2 %>% dplyr::filter(Genes %ni% maternal_genes) %>%
  dplyr::mutate(GeneType = case_when(AC> 200 & FC >10 ~ "V. High",
                                     between(AC, 50, 200) & between(FC, 5, 10) ~ "High",
                                     AC> 200 & between(FC, 5, 10) ~ "High",
                                     between(AC, 50, 200) & FC > 10 ~ "High", 
                                     between(AC, 1, 50) & between(FC, 1, 5) ~ "Low",
                                     AC> 50 & between(FC, 1, 5) ~ "Low",
                                     between(AC, 1, 50) & FC > 5 ~ "Low",
                                     TRUE ~ "Other")) %>%
  group_by(Genes) %>%
  mutate(EverHigh = any(GeneType %in% c("V. High","High")),
         GeneType = case_when(GeneType %in% c("V. High","High") ~ GeneType,
                              GeneType == "Low" & EverHigh == FALSE ~ "Low",
                              TRUE ~ "Other" )) %>%
  ungroup() 

#For Transcription Start Sites 

#Collate TSS from Chen et al., Saito et. al and Kruesi et al.
Chen <-fread("Supp_TableS2_Chen.txt", header = T, stringsAsFactors = FALSE, na.strings=c("","NA")) 
setnames(Chen, "Gene", "Genes")
Chen[Chen$Genes==""] <- NA
Chen <- Chen %>% dplyr::select(Genes, strand, mode_position) %>% drop_na(Genes) %>% 
  group_by(Genes) %>%
  dplyr::slice(if(all(strand == '+')) which.max(mode_position) 
               else which.min(mode_position)) 
Chen_2 <- Chen %>% dplyr::select(Genes, mode_position) 

#Read Saito and pick closest TSS
Saito <-fread("Saito_2.txt", header = T, stringsAsFactors = FALSE, na.strings=c("","NA")) 
setnames(Saito, "common gene name", "Genes")
Saito[Saito$Genes==""] <- NA
Saito <- Saito %>% dplyr::select(Genes, strand, `start position`) %>% drop_na(Genes) %>%
  group_by(Genes) %>%
  dplyr::slice(if(all(strand == '+')) which.max(`start position`) 
               else which.min(`start position`)) 
Saito_2 <- Saito %>% dplyr::select("Genes", "start position") 

#Read Kruesi and pick closest TSS
Kruesi <-fread("Kruesi supp T2.txt", header = T, stringsAsFactors = FALSE, na.strings=c("","NA")) 
setnames(Kruesi, "Gene Name", "Genes")
Kruesi[Kruesi$Genes==""] <- NA
Kruesi <- Kruesi %>% dplyr::select(Genes, Strand, `WT Embryo TSS`) %>% drop_na(Genes) %>%
  group_by(Genes) %>%
  dplyr::slice(if(all(Strand == '1')) which.max(`WT Embryo TSS`) 
               else which.min(`WT Embryo TSS`)) 
Kruesi_2 <- Kruesi %>% dplyr::select("Genes", "WT Embryo TSS") 

#Pivot_longer
Chen_3 <- Chen_2 %>% `names<-`(c("Genes","Chen")) %>% pivot_longer(-Genes, names_to = "Source", values_to = "Position")
Saito_3 <- Saito_2 %>% `names<-`(c("Genes","Saito")) %>% pivot_longer(-Genes, names_to = "Source", values_to = "Position")
Kruesi_3 <- Kruesi_2 %>% `names<-`(c("Genes","Kruesi")) %>% pivot_longer(-Genes, names_to = "Source", values_to = "Position")

txdb = loadDb("CE_WS235_Genes.sqlite")
ce_genes=sort(genes(txdb))
geneSymbols <- mapIds(org.Ce.eg.db, keys=names(ce_genes), column="SYMBOL", keytype="ENTREZID", multiVals="first")
geneNumbers = names(geneSymbols)
names(geneNumbers) = geneSymbols

chroms = as.vector(chrom(ce_genes[geneNumbers]))
starts = as.vector(start(ce_genes[geneNumbers]))
ends = as.vector(end(ce_genes[geneNumbers]))

names(chroms)=geneSymbols
names(starts)=geneSymbols
names(ends)=geneSymbols

#for reference TSS
pdc_3 <- pdc_2 %>% mutate(Ref_start = starts[Genes]) %>% mutate(end = ends[Genes]) %>% 
  mutate(chr = chroms[Genes]) %>% dplyr::select(Genes, Ref_start, end, chr) %>% distinct(Genes, .keep_all = TRUE)
pdc_4 <- pdc_3 %>% dplyr::select(Genes, Ref_start) %>%
  `names<-`(c("Genes","Ref_start")) %>% pivot_longer(-Genes, names_to = "Source", values_to = "Position")

all_TSS <- bind_rows(pdc_4,Chen_3, Kruesi_3, Saito_3)
#bind all
all_TSS_2 <- all_TSS %>% pivot_wider(id_cols = Genes, names_from = Source, values_from = Position)

myCDS <- cds(txdb, c("GENEID", "CDSID", "CDSSTART", "CDSSTRAND", "CDSEND", "CDSCHROM"), use.names=TRUE)
CDSIDtoENTREZID <- mapIds(txdb, keys = keys(txdb,keytype = "CDSID"), column = "GENEID", keytype = "CDSID")
ENTREZtoSYMBOL <-  mapIds(org.Ce.eg.db, keys = mapIds(txdb, keys = keys(txdb,keytype = "CDSID"), column = "GENEID", keytype = "CDSID"), column = "SYMBOL", keytype = "ENTREZID")
CDSIDtoSymbol <- setNames(names(CDSIDtoENTREZID),ENTREZtoSYMBOL)
names(myCDS) <- names(CDSIDtoSymbol)

myCDS_df <- myCDS %>%
  as_tibble() %>%
  mutate(seqname = myCDS@ranges@NAMES) %>%
  unnest(GENEID) %>%
  group_by(GENEID) %>%
  mutate(Good = (CDSSTRAND == "+" & CDSSTART == min(CDSSTART)) | (CDSSTRAND == "-" & CDSSTART == max(CDSSTART) )) %>%
  filter(Good) %>% dplyr::select("seqname", "CDSSTART", "CDSEND", "CDSCHROM")    
setnames(myCDS_df, "seqname", "Genes")

all_TSS_2 %>% left_join(myCDS_df, by = "Genes") 

TSS_selected <- all_TSS_2 %>% left_join(myCDS_df, by = "Genes") %>% as.data.frame() %>%
  gather(key, val, "Ref_start":"Saito") %>% 
  group_by(Genes) %>% mutate(Diff = (val-CDSSTART)+1,
                             Diff = replace(Diff, Diff <=0, NA)) %>%
  dplyr::slice(which.min(Diff)) %>%
  dplyr::select(Genes, GENEID, CDSSTART, val, CDSEND, CDSCHROM) %>%
  right_join(all_TSS_2) %>%
  dplyr::select(names(all_TSS_2), everything()) %>% dplyr::select(Genes, GENEID,val, CDSEND, CDSCHROM)

setnames(TSS_selected, "CDSCHROM", "seqnames")
setnames(TSS_selected, "val", "start")
setnames(TSS_selected, "GENEID", "gene_id")
setnames(TSS_selected, "CDSEND", "end")

TSS_selected <- TSS_selected %>% drop_na() %>% dplyr::filter(end > start) %>% distinct(Genes, .keep_all = TRUE)
TSS_selected <- column_to_rownames(TSS_selected, var = "gene_id")


TSS_GRange <- makeGRangesFromDataFrame(TSS_selected, keep.extra.columns = TRUE)
mcols(TSS_GRange, level="within")$gene_id <- names(TSS_GRange)

upstream = 500
downstream = 50
ce_promoters=trim(promoters(TSS_GRange,upstream=upstream,downstream=downstream))

# Add gene_id to Granges object
mcols(ce_promoters, level = "within")$gene_id <- geneSymbols[match(names(ce_promoters),names(geneSymbols))]



















