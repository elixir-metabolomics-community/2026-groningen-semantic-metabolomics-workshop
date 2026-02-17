library("Spectra")
library("MsExperiment")
library("MsQuality")
library("ggplot2")
library("tidyr")

str_split <- function(x) {
    stringr::str_split(x, '\\.')[[1]][1]
}

drop_na_cols <- function(df) {
    as.data.frame(df[, colSums(is.na(df)) == 0])
}

download_url_files <- function(sd_df, col_name = "url") {
  urls <- sd_df[[col_name]]
  tmp_dir <- tempdir()
  local_paths <- file.path(tmp_dir, basename(urls))

  for (i in seq_along(urls)) {
    utils::download.file(urls[i], local_paths[i], mode = "wb", quiet = TRUE)
  }

  local_paths
}

sd <- read.csv2("sampleData.tsv", sep = "\t", header=TRUE)
paths <- download_url_files(sd, "url")

spectra_obj <- Spectra(paths, backend = MsBackendMzR())

lmse <- MsExperiment()
sampleData(lmse) <- DataFrame(sd)
experimentFiles(lmse) <- MsExperimentFiles(mzML_files = paths)
spectra(lmse) <- spectra_obj
lmse <- linkSampleData(
    lmse,
    with = "experimentFiles.mzML_file",
    withIndex = seq_len(21)
)

metrics <- calculateMetrics(
    object = lmse,
    metrics = qualityMetrics(lmse),
    filterEmptySpectra = FALSE,
    msLevel = 1L
)

md <- drop_na_cols(metrics)
md$sample_name <- basename(rownames(md)) |> sapply(str_split)


# Plot all metrics by sample ID, colored by dilution level
plot_data <- merge(x=md, y=sd, by="sample_name", all=TRUE)

# Get numeric columns (excluding sample_name, dilution, injection)
numeric_cols <- names(plot_data)[sapply(plot_data, is.numeric)]

# Create a long format dataframe for plotting
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = all_of(numeric_cols),
    names_to = "metric",
    values_to = "value"
  )

# Create faceted plots for all metrics
p <- ggplot(plot_data_long, aes(x = sample_name, y = value, color = dilution)) +
  geom_point(size = 2) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    strip.text = element_text(size = 8)
  ) +
  labs(
    title = "QC Metrics by Sample",
    x = "Sample Name",
    y = "Value",
    color = "Dilution Level"
  )

# Save the plot to file
ggsave("all_metrics_plot.png", plot = p, width = 16, height = 12)

mzqc_metrics <- calculateMetricsFromSpectra(
    spectra_obj,
    metrics = qualityMetrics(spectra_obj),
    format="mzQC"
)

rmzqc::writeMZQC("metrics.mzqc", mzqc_metrics)