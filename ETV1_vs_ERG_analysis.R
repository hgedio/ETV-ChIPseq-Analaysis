library(GenomicRanges)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)
library(ggplot2)
library(dplyr)

etv1_peaks <- toGRanges("raw_data/ETV1/GSE110655_RAW/GSM3223708_VCaP_shCT_BAF155_peaks.narrowPeak", format = "narrowPeak")
erg_peaks <- toGRanges("raw_data/ETV1/GSE110655_RAW/GSM3223706_VCaP_shERG_ERG_peaks.narrowPeak", format = "narrowPeak")

etv1_peaks
erg_peaks

overlap <- findOverlaps(etv1_peaks, erg_peaks)

etv1_specific <- etv1_peaks[-queryHits(overlap)]
erg_specific <- erg_peaks[-subjectHits(overlap)]


library(tibble)
library(ggplot2)

summary_df <- tibble(
  category = c("ETV1 Total", "ERG Total", "Overlapping", "ETV1 Specific", "ERG Specific"),
  count = c(length(etv1_peaks), length(erg_peaks), length(overlap), length(etv1_specific), length(erg_specific))
)

ggplot(summary_df, aes(x = category, y = count, fill = category)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 14) +
  labs(title = "ETV1 vs ERG Peak Comparison", x = "", y = "Number of Peaks") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

etv1_widths <- width(etv1_peaks)
erg_widths <- width(erg_peaks)

library(GenomicRanges)
library(ggplot2)
library(dplyr)

peak_widths_df <- bind_rows(
  data.frame(
    width = width(etv1_peaks),
    group = "ETV1"
  ),
  data.frame(
    width = width(erg_peaks),
    group = "ERG"
  )
)

ggplot(peak_widths_df, aes(x = width, fill = group)) +
  geom_histogram(alpha = 0.5) +
  theme_minimal(base_size =14) +
  labs(
    title = "Peak WIdth Distribution",
    x = "Peak WIdth (bp)",
    y = "Density"
  ) + 
  scale_fill_manual(values = c("ETV1" = "#1f77b4", "ERG" = "#ff7f0e"))
