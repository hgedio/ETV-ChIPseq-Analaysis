library(rtracklayer)
etv1 <- import("raw_data/ETV1/GSE110655_RAW/GSM3223705_VcaP_shCt_ERG_peaks.narrowPeak.gz")

library(rtracklayer)
etv4 <- import("raw_data/ETV4/GSE129803_RAW/GSM4260329_siETV4_72hr_ER.narrowPeak")

overlap_etv1_etv4 <- findOverlaps(etv1, etv4)

shared_peaks_etv1 <- etv1[queryHits(overlap_etv1_etv4)]
shared_peaks_etv4 <- etv4[subjectHits(overlap_etv1_etv4)]

etv1_unique <- etv1[-queryHits(overlap_etv1_etv4)]

length(shared_peaks_etv1)
length(etv1_unique)

plot(width(etv1_unique), main="Width of Unique ETV1 Peaks", xlab="Peak Index", ylab="Width (bp)")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("ChIPpeakAnno", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))

library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoData <- toGRanges(txdb)

etv1_annotated <- addGeneIDs(etv1_annotated,
                             orgAnn = org.Hs.eg.db,
                             feature_id_type = "entrez_id",
                             IDs2Add = c("symbol"))

anno_df <- as.data.frame(etv1_annotated)
top_genes <- sort(table(anno_df$symbol), decreasing = TRUE)[1:10]

write.csv(as.data.frame(top_genes), "etv1_top10_genes.csv")

View(read.csv("etv1_top10_genes.csv"))

barplot(top_genes,
        main = "Top 10 Genes Near ETV1-Specific Peaks",
        col = "steelblue",
        las = 2,
        cex.names = 0.8,
        ylab = "Peak Count")
