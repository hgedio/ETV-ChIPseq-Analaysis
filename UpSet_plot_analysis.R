library(ChIPpeakAnno)
library(GenomicRanges)
library(rtracklayer)
library(UpSetR)


etv1 <- toGRanges(
  "raw_data/ETV1/GSE110655_RAW/GSM3223708_VCaP_shCT_BAF155_peaks.narrowPeak",
  format = "narrowPeak"
)

etv4 <- toGRanges(
  "raw_data/ETV4/GSE129803_RAW/GSM4260329_siETV4_72hr_ER.narrowPeak",
  format = "narrowPeak"
)

etv5_df <- read.table(
  "raw_data/ETV5/GSE232262_IP_2_vs_Input_2_peaks.txt.gz",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)


etv5 <- GRanges(
  seqnames = etv5_df$chr,    # adjust if your column name is different
  ranges = IRanges(start = etv5_df$start,
                   end = etv5_df$end)
)

erg <- toGRanges(
  "raw_data/ETV1/GSE110655_RAW/GSM3223706_VCaP_shERG_ERG_peaks.narrowPeak",
  format = "narrowPeak"
)

peak_list <- list(
  ETV1 = etv1,
  ETV4 = etv4,
  ETV5 = etv5,
  ERG  = erg
)

overlaps <- findOverlapsOfPeaks(peak_list)

overlapping_peaks <- overlaps$peaklist

all_peaks <- reduce(unlist(GRangesList(peak_list)))

binary_matrix <- matrix(0, nrow = length(all_peaks), ncol = length(peak_list))
colnames(binary_matrix) <- names(peak_list)

for (i in seq_along(peak_list)) {
  binary_matrix[, i] <- countOverlaps(all_peaks, peak_list[[i]]) > 0
}

upset_df <- as.data.frame(binary_matrix)

upset(
  upset_df,
  sets = names(peak_list),
  main.bar.color = "steelblue",
  sets.bar.color = "forestgreen",
  order.by = "freq",
  sets.x.label = "Number of peaks",
  text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.2)
)
