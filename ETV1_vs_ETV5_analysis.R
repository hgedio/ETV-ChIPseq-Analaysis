library(rtracklayer)
library(GenomicRanges)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

etv1 <- import("raw_data/ETV1/GSE110655_RAW/GSM3223705_VCaP_shCt_ERG_peaks.narrowPeak.gz")

etv5_df <- read.table(
  "raw_data/ETV5/GSE232262_IP_2_vs_Input_2_peaks.txt.gz",
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)

etv5 <- GRanges(
  seqnames = etv5_df$chr,
  ranges = IRanges(start = etv5_df$start, end = etv5_df$end),
  strand = "*"
)

etv5 <- keepSeqlevels(etv5, intersect(seqlevels(etv1), seqlevels(etv5)), pruning.mode = "coarse")
seqlevelsStyle(etv5) <- seqlevelsStyle(etv1)

overlap <- findOverlaps(etv1, etv5)
etv1_unique <- etv1[-queryHits(overlap)]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
anno_data <- toGRanges(txdb)

etv1_annotated <- annotatePeakInBatch(
  etv1_unique,
  AnnotationData = anno_data,
  output = "nearestLocation",
  bindingRegion = c(-5000, 5000),
  use.names = TRUE
)

etv1_annotated <- addGeneIDs(
  etv1_annotated,
  orgAnn = org.Hs.eg.db,
  feature_id_type = "entrez_id",
  IDs2Add = "symbol"
)

head(etv1_annotated)
unique(mcols(etv1_annotated)$feature)


View(etv1_annotated)
write.csv(as.data.frame(etv1_annotated), "results/etv1_unique_vs_etv5_annotated.csv", row.names = FALSE)


peakAnno_etv1 <- annotatePeak(etv1, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Hs.eg.db")

peakAnno_etv5 <- annotatePeak(etv5, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = "org.Hs.eg.db")

plotDistToTSS(list(ETV1 = peakAnno_etv1, ETV5 = peakAnno_etv5),
              title = "Peak Distance to TSS: ETV1 vs ETV5")

