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

#Figure_5_Gene Structure

#Fig.5_Gene Length

pdc_2 %>%
  dplyr::filter(Genes %ni% maternal$Genes) %>%
  dplyr::mutate(GeneType = case_when(AC> 200 & FC >10 ~ "V. High",
                                     between(AC, 50, 200) & between(FC, 5, 10) ~ "High",
                                     AC> 200 & between(FC, 5, 10) ~ "High",
                                     between(AC, 50, 200) & FC > 10 ~ "High", 
                                     between(AC, 1, 50) & between(FC, 1, 5) ~ "Low",
                                     AC> 50 & between(FC, 1, 5) ~ "Low",
                                     between(AC, 1, 50) & FC > 5 ~ "Low",
                                     TRUE ~ "Other"), Length = gene_width[Genes]) %>%
  group_by(Genes) %>%
  mutate(EverHigh = any(GeneType %in% c("V. High","High")),
         GeneType = case_when(GeneType %in% c("V. High","High") ~ GeneType,
                              GeneType == "Low" & EverHigh == FALSE ~ "Low",
                              TRUE ~ "Other" )) %>%
  ungroup() %>%
  filter(GeneType != "Other") %>%
  ggplot(aes(x=factor(GeneType, levels = c("Low", "High", "V. High")), y=Length, fill = factor(GeneType))) + 
  geom_boxplot(width = 0.5,  outlier.shape = NA) + ylim(0,20000) + theme_bw() + theme(aspect.ratio = 1/1) +
  labs(x = "Gene Type", y = "Gene Length", fill = "Gene Type") +
  scale_fill_brewer(palette="Pastel1", guide = FALSE) + 
  stat_compare_means(aes(label=..p.adj..), comparisons = my_comparisons, label.y = c(11000, 13000, 11000))#+ 
stat_compare_means(label = "p.format", method = "wilcox.test")

#Fig.5_Average Intron Length 

myIntrons <- intronsByTranscript(txdb,use.names = TRUE)
names(myIntrons) <- names(TXIDtoSymbol)[match(names(myIntrons),TXIDtoSymbol)]

avg_intron_lentgh = sapply(myIntrons[],function(x)mean(width(x)))

pdc_2 %>%
  dplyr::filter(Genes %ni% maternal$Genes) %>%
  dplyr::mutate(GeneType = case_when(AC> 200 & FC >10 ~ "V. High",
                                     between(AC, 50, 200) & between(FC, 5, 10) ~ "High",
                                     AC> 200 & between(FC, 5, 10) ~ "High",
                                     between(AC, 50, 200) & FC > 10 ~ "High", 
                                     between(AC, 1, 50) & between(FC, 1, 5) ~ "Low",
                                     AC> 50 & between(FC, 1, 5) ~ "Low",
                                     between(AC, 1, 50) & FC > 5 ~ "Low",
                                     TRUE ~ "Other"), IntronLength = avg_intron_lentgh[Genes]) %>%
  group_by(Genes) %>%
  mutate(EverHigh = any(GeneType %in% c("V. High","High")),
         GeneType = case_when(GeneType %in% c("V. High","High") ~ GeneType,
                              GeneType == "Low" & EverHigh == FALSE ~ "Low",
                              TRUE ~ "Other" )) %>%
  ungroup() %>%
  filter(GeneType != "Other") %>% 
  ggplot(aes(x=factor(GeneType, levels = c("Low", "High", "V. High")), y=IntronLength, fill = factor(GeneType))) + 
  geom_boxplot(width = 0.5,  outlier.shape = NA) + ylim(0,1500) + theme_bw() + theme(aspect.ratio = 1/1) +
  labs(x = "Gene Type", y = "Avg. Inron Length", fill = "Gene Type") +
  scale_fill_brewer(palette="Pastel1", guide = FALSE) + 
  stat_compare_means(aes(label=..p.adj..), comparisons = my_comparisons, label.y = c(1000, 1200, 1000)) #+ 
stat_compare_means(label = "p.format", method = "wilcox.test")

#Fig. 5_Outron vs. not

Saito_outron <- fread("Saito_for_outron.txt",header=T,stringsAsFactors = FALSE)
setnames(Saito_outron, "common gene name", "Genes")
Saito_2 <- Saito_outron %>% dplyr::select("Genes", "TSS type")

pdc_2 %>% dplyr::filter(Genes %ni% maternal$Genes) %>%
  left_join(Saito_2, by = "Genes") %>%
  dplyr::filter(`TSS type` == "outron" |`TSS type` == "exon") %>%
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
  ungroup() %>%
  filter(GeneType != "Other") %>%
  dplyr::select("Genes","TSS type","GeneType") %>% group_by(GeneType, `TSS type`) %>% 
  summarise(each_type = n()) %>% mutate(type_total = sum(each_type)) %>% mutate(proportion = each_type/type_total) %>% ungroup() %>%
  ggplot(aes(x=factor(`TSS type`), y=proportion, fill=`TSS type`)) + ylim(0, 1) + 
  geom_bar(stat='identity') +  labs(x = "Rate Category", fill = "Type") +
  facet_wrap(.~factor(GeneType, levels = c("Low","High","V. High")),nrow=1, strip.position = "bottom") + 
  scale_fill_brewer(palette="Pastel1") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.placement = "outside")

## Use Fisher's exact test for stats

#To run for Isoforms

NumIntrons <- sapply(intronsByTranscript(txdb),NROW)
names(NumIntrons) <- names(intronsByTranscript(txdb,use.names = TRUE))
TXIDtoENTREZID <- mapIds(txdb, keys = keys(txdb,keytype = "TXNAME"), column = "GENEID", keytype = "TXNAME")
ENTREZtoSYMBOL <-  mapIds(org.Ce.eg.db, keys = mapIds(txdb, keys = keys(txdb,keytype = "TXNAME"), column = "GENEID", keytype = "TXNAME"), column = "SYMBOL", keytype = "ENTREZID")
TXIDtoSymbol <- setNames(names(TXIDtoENTREZID),ENTREZtoSYMBOL)
names(NumIntrons) <- names(TXIDtoSymbol[match(names(NumIntrons),TXIDtoSymbol)]) 

#Fig. 5_isoforms

pdc_2 %>%
  dplyr::filter(Genes %ni% maternal$Genes) %>%
  mutate(GeneType = case_when(AC> 200 & FC >10 ~ "V. High",
                              between(AC, 50, 200) & between(FC, 5, 10) ~ "High",
                              between(AC, 2, 50) & between(FC, 2, 5) ~ "Low",
                              TRUE ~ "Other"), Introns = NumIntrons[Genes]) %>%
  group_by(Genes) %>%
  mutate(EverHigh = any(GeneType %in% c("V. High","High")),
         GeneType = case_when(GeneType %in% c("V. High","High") ~ GeneType,
                              GeneType == "Low" & EverHigh == FALSE ~ "Low",
                              TRUE ~ "Other" )) %>%
  ungroup() %>%
  mutate(GeneType = factor(GeneType, levels = c("Low", "High", "V. High", "Other"),
                           labels = c("Low", "High", "Very High", "Other"))) %>%
  arrange(GeneType) %>%         
  filter(GeneType != "Other") %T>% {
    group_by(.,GeneType) %>%
      summarise(Introns = mean(Introns, na.rm = TRUE)) ->> newIntronMeans
  } %>%
  ggplot(aes(x = Introns, fill = GeneType)) +
  scale_x_continuous(limits = c(-0.5,15.5)) +
  geom_vline(data = newIntronMeans, 
             aes(color = GeneType, xintercept = newIntronMeans$Introns[..PANEL..]),
             lwd = 1, lty = 2) +  
  geom_bar(aes(y=..prop..),
           alpha = 0.75, width = 1) + scale_fill_manual(name="GeneType",values = c("#F5B7B1", "#AED6F1", "#A3E4D7")) +
  scale_color_manual(name="Mean",values = c("#F5B7B1", "#AED6F1", "#A3E4D7")) +
  facet_wrap(~GeneType,scales = "fixed" ) + coord_flip() +
  labs(x = "Number of Introns", y = "Density") + theme_bw()

myTranscripts <- transcripts(txdb, c("TXNAME", "GENEID"), use.names=TRUE)
names(myTranscripts) <- names(transcripts(txdb, c("TXNAME", "GENEID"), use.names=TRUE))
names(myTranscripts) <- names(TXIDtoSymbol[match(names(myTranscripts),TXIDtoSymbol)]) 

myTranscripts_df <- as.data.frame(myTranscripts@elementMetadata) %>% mutate(GENEID = unlist(GENEID)) %>% 
  group_by(GENEID) %>% summarise(count = n())

id_geneSymbols <- as.data.frame(geneSymbols)
setDT(id_geneSymbols, keep.rownames = TRUE)
setnames(id_geneSymbols, "rn", "GENEID")
myTranscripts_df <- myTranscripts_df %>% left_join(id_geneSymbols, by = "GENEID")
setnames(myTranscripts_df, "geneSymbols", "Genes")

pdc_2 %>%
  dplyr::filter(Genes %ni% maternal$Genes) %>%
  left_join(myTranscripts_df, by = "Genes") %>%
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
  ungroup() %>%
  mutate(GeneType = factor(GeneType, levels = c("Low", "High", "V. High", "Other"),
                           labels = c("Low", "High", "Very High", "Other"))) %>%
  arrange(GeneType) %>%         
  filter(GeneType != "Other") %T>% {
    group_by(.,GeneType) %>%
      summarise(Isoforms = mean(count, na.rm = TRUE)) ->> newIsoformMeans
  } %>%
  ggplot(aes(x = count, fill = GeneType)) +
  scale_x_continuous(limits = c(-0.5,8.5)) +
  geom_vline(data = newIsoformMeans, 
             aes(color = GeneType, xintercept = newIsoformMeans$Isoforms[..PANEL..]),
             lwd = 1, lty = 2) +  
  geom_bar(aes(y=..prop..),
           alpha = 0.75, width = 1) + scale_fill_manual(name="GeneType",values = c("#F5B7B1", "#AED6F1", "#A3E4D7")) +
  scale_color_manual(name="Mean",values = c("#F5B7B1", "#AED6F1", "#A3E4D7")) +
  facet_wrap(~GeneType,scales = "fixed" ) + coord_flip() +
  labs(x = "Number of Isoforms", y = "Density") + theme_bw()

#Fig.5_chromosomes 

pdc_2 %>% dplyr::filter(Genes %ni% maternal_genes) %>%
  mutate(Chromosome = chroms[Genes]) %>%
  dplyr::filter(Chromosome != "NA") %>%
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
  ungroup() %>%
  filter(GeneType != "Other") %>%
  dplyr::select("Genes","Chromosome","GeneType") %>% group_by(GeneType, Chromosome) %>% 
  dplyr::filter(!duplicated(Genes)) %>% 
  summarise(each_ch = n()) %>% mutate(chr_total = sum(each_ch)) %>% mutate(proportion = each_ch/chr_total) #%>%
  ggplot(aes(x=factor(Chromosome), y=proportion, fill=Chromosome))+
  geom_bar(stat='identity') +  labs(x = "Rate Category", fill = "Chromosome") +
  facet_wrap(.~factor(GeneType, levels = c("Low", "High", "V. High")),nrow=1, strip.position = "bottom") + 
  scale_fill_brewer(palette="Pastel1") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.placement = "outside")
  
#Fig.5_Karyogram
  
data.frame(lTip = c(527,306,494,720,643,572),
             lArm = c(3858,4879,3722,3896,5897,6137),
             center= c(11040,12020,10340,12970,16550,12480),
             rArm = c(14875, 14609, 13217, 16712, 20337, 16417),
             rTip = c(15072,15279,13784,17494,20920,17719),
             chr = c("chrI","chrII","chrIII","chrIV","chrV","chrX")) -> Boundaries

my_cell <- "MSx1"
my_max <- 600
my_label_threshold <- 25

tibble(Tintori) %>%
  filter(Cell == my_cell &
           !is.na(chr)) %>%
  {ggplot(data = ., aes(x=starts, y=AC, color=GeneType, fill = GeneType)) +
      geom_bar(stat="identity",position=position_dodge(),width=1) +
      geom_segment(data = Boundaries, inherit.aes = FALSE, aes(x=-100, y = 0, yend = 0, xend = rTip*1000), color = "black") +
      scale_fill_manual(values = c(High = "#8dd3c7", Medium = "#fb8072"))+
      scale_color_manual(values = c(High = "#8dd3c7", Medium = "#fb8072")) + 
      ylim(0,my_max) +
      xlim(-100,max(Boundaries$rTip*1000))+
      theme_minimal()+
      facet_wrap(facets = "chr", ncol=1)+
      labs(title=paste(my_cell, "expression increase by chromosomal position relative to", Parents[my_cell]),
           color = "Rate Category", fill = "Rate Category")+
      ylab("Absolute change")+
      xlab("Coordinate") +
      geom_text_repel(data=filter(., AC > my_label_threshold), aes(label=Genes), size = 1,
                      vjust=-0.25, nudge_y = 100, nudge_x = 100, segment.size=0.5, max.overlaps = 20, segment.colour="gray")}
  
  
#Fig 6_Example Motif Logo
  
library("ggseqlogo")
  
motif_row_names <- inr_motif2[,1]
  
inr_motif2 <- (data.matrix(fread("inr from list 5.txt", header = T)))
rownames(inr_motif2) <- motif_row_names 
inr_motif2 <- inr_motif2[,-1]

ggseqlogo(inr_motif2)

#Half-lives

est_rates <- fread("halflives_06.18.23.txt", header = T, stringsAsFactors = FALSE) 

pdc_2 %>% mutate(rate = AC/`ccTime (mins)`) %>% left_join(est_rates, by = c("Genes", "Cell")) %>% 
  dplyr::filter(Genes %ni% maternal_genes) %>%
  dplyr::mutate(GeneType = case_when(AC> 200 & FC >10 ~ "V. High",
                                     between(AC, 50, 200) & between(FC, 5, 10) ~ "High",
                                     AC> 200 & between(FC, 5, 10) ~ "High",
                                     between(AC, 50, 200) & FC > 10 ~ "High", 
                                     between(AC, 1, 50) & between(FC, 1, 5) ~ "Low",
                                     AC> 50 & between(FC, 1, 5) ~ "Low",
                                     between(AC, 1, 50) & FC > 5 ~ "Low",
                                     TRUE ~ "Other")) %>%
  dplyr::filter(GeneType != "Other") %>% 
  dplyr::filter(GeneType == "High" | GeneType == "V. High") %>% 
  dplyr::filter(Founder == "AB" | Founder == "C" | Founder == "E" | Founder == "MS" | Founder == "EMS") %>% 
  #filter(n()==1)
  ggplot(aes(x=Cell,y=`Synthesis rate`, fill = as.factor(Cell))) + 
  geom_boxplot(fill = "gray", alpha = 0.3, outlier.shape = NA) + 
  geom_jitter(aes(color = GeneType), size = 0.5) + scale_color_manual(values = c("#fb8072","#8dd3c7")) +
  ylim(0,50) + theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 20)) 


pdc_2 %>% mutate(rate = AC/`ccTime (mins)`) %>% left_join(est_rates, by = c("Genes", "Cell")) %>% 
  mutate(FoldChange = `Synthesis rate`/rate) %>%
  dplyr::filter(Genes %ni% maternal_genes) %>%
  dplyr::mutate(GeneType = case_when(AC> 200 & FC >10 ~ "V. High",
                                     between(AC, 50, 200) & between(FC, 5, 10) ~ "High",
                                     AC> 200 & between(FC, 5, 10) ~ "High",
                                     between(AC, 50, 200) & FC > 10 ~ "High", 
                                     between(AC, 1, 50) & between(FC, 1, 5) ~ "Low",
                                     AC> 50 & between(FC, 1, 5) ~ "Low",
                                     between(AC, 1, 50) & FC > 5 ~ "Low",
                                     TRUE ~ "Other")) %>%
  dplyr::filter(GeneType != "Other") %>% 
  dplyr::filter(GeneType == "High" | GeneType == "V. High") %>% 
  dplyr::filter(Founder == "AB" | Founder == "C" | Founder == "E" | Founder == "MS" | Founder == "EMS") %>% 
  ggplot(aes(x=Cell,y=FoldChange, fill = as.factor(Cell))) + 
  geom_boxplot(fill = "gray", alpha = 0.3, outlier.shape = NA) + 
  geom_jitter(aes(color = GeneType), size = 0.5) + scale_color_manual(values = c("#fb8072","#8dd3c7")) +
  ylim(0,5) + theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 20)) 















