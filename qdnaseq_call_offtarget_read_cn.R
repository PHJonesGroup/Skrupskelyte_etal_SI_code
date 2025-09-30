# This pipeline is suitable for targeted sequencing

### Load libraries ####
suppressPackageStartupMessages({
  library(QDNAseq)
  library(Biobase)
  library(png)
  library(grid)
  library(gridExtra)
  library(gtools)
  library(ggplot2)
  library(VariantAnnotation)
  library(logging)
})
###############################################################################
### Read in command line arguments
args <- commandArgs(TRUE)

samplename <- args[1]
tumourbam <- args[2]
normalbam <- args[3]
targetregions <- args[4]
run_name <- args[9]

# Check if sample sex is provided
if (length(args) > 4) {
  sex <- args[5]
} else {
  sex <- NULL
}

### CNA calling run parameters
# Whether to correct for replication timing bias
correctReplication <- FALSE

# Higher gamma means less segments
segmentation_gamma <- as.numeric(args[6])
# Minimum number of bins per segment
segmentation_kmin <- as.numeric(args[7])
# minimum logr has to deviate from 0 to be considered for significance testing
logr_minimum_deviation <- 0.025

# maximum number of bins a segment can contain before downsampling
# for large numbers of bins a very small deviation will always be significant
max_num_bins_test <- 500

# Statistical testing parameters
test_significance_threshold <- 0.05
segm_summary_method <- 2 # 1 = mean, 2 = median

# Load precalculated bins
bins_path <- args[8]
load(bins_path)

# Log all variable
loginfo(paste0("Sample name: ", samplename))
loginfo(paste0("Tumour bam: ", tumourbam))
loginfo(paste0("Normal bam: ", normalbam))
loginfo(paste0("Target regions: ", targetregions))
loginfo(paste0("Sex as specified: ", sex))
loginfo(paste0("Segmentation gamma: ", segmentation_gamma))
loginfo(paste0("Segmentation kmin: ", segmentation_kmin))
loginfo(paste0("Logr minimum deviation: ", logr_minimum_deviation))
loginfo(paste0("Max number of bins for testing: ", max_num_bins_test))
loginfo(paste0("Test significance threshold: ", test_significance_threshold))
loginfo(paste0("Segmentation summary method: ", segm_summary_method))
loginfo(paste0("Bins path: ", bins_path))
loginfo(paste0("Run name: ", run_name))

# Output directory variable
output_dir <- paste0(run_name)

# create output directory
dir.create(output_dir)

###############################################################################
# Read in data for off-target regions from bam files

# Read in target regions BED file
target <- read.table(targetregions, header = TRUE,
                     stringsAsFactors = FALSE, sep = "\t")

# Add column names
colnames(target)[1:3] <- c("chromosome", "start", "end")

# Remove 'chr' from chromosome names if present
if (any(grepl("chr", target$chromosome))) {target$chromosome <- gsub("chr", "", target$chromosome)}

# Extract Start and End positions
target$start <- target$start
target$end <- target$end

# Convert to GRanges and find overlaps with bins
target <- makeGRangesFromDataFrame(target)
bins_gr <- makeGRangesFromDataFrame(bins)
overlap <- findOverlaps(target, bins_gr)
selection <- which(!(1:nrow(bins)) %in% subjectHits(overlap))
bins <- bins[selection, ]
# Clean up
rm(target, bins_gr, overlap)

# Bin the read counts for tumour and normal using the off-target bins
readCounts_tumour <- binReadCounts(bins, bamfiles = tumourbam)
readCounts_normal <- binReadCounts(bins, bamfiles = normalbam)

# If Female, remove Y chromosome bins in both tumour and normal
if (sex == "female") {
  readCounts_tumour <-  new('QDNAseqReadCounts', bins = bins[bins$chromosome != "Y",,drop = F], counts = assayDataElement(readCounts_tumour, "counts")[bins$chromosome != "Y",,drop = F], phenodata = phenoData(readCounts_tumour))
	readCounts_normal <- new('QDNAseqReadCounts', bins = bins[bins$chromosome != "Y",,drop = F], counts = assayDataElement(readCounts_normal, "counts")[bins$chromosome != "Y",,drop = F], phenodata = phenoData(readCounts_normal))
}
# Normalise tumour read counts by normal read counts
normaliseReadCounts <- function(readCounts_tumour, readCounts_normal, sample_index = 1, minReadsThreshold = 10, correctReplication = F) {
  tumour_counts <- assayDataElement(readCounts_tumour, "counts")[, sample_index]
  normal_counts <- assayDataElement(readCounts_normal, "counts")[, sample_index]
  # Set bins with very low counts to NA in the normal sample
  normal_counts[normal_counts < minReadsThreshold] <- NA
  # Calculate ratios of tumour to normal
  tumour_r <- tumour_counts / normal_counts
  # Set infinite values to NA
  tumour_r[!is.finite(tumour_r)] <- NA
  # Log2 ratio transformation
  tumour_logr <- log2(tumour_r / mean(tumour_r, na.rm = TRUE))
  # Set infinite values to NA
  tumour_logr[!is.finite(tumour_logr)] <- NA
  # Convert to matrix
  tumour_logr <- matrix(tumour_logr, ncol = 1)
  # Set dimnames
  dimnames(tumour_logr) <- dimnames(assayDataElement(readCounts_tumour, "counts"))
  # Create a QDNAseq object
  res <- new("QDNAseqCopyNumbers", bins = featureData(readCounts_tumour), copynumber = tumour_logr, phenodata = phenoData(readCounts_tumour))

  # Store logR values in copynumber assayData element
  counts <- tumour_logr

# Determine the sample sex
if (is.null(sex)) {
  y_counts <- assayDataElement(readCounts_normal, "counts")[bins$chromosome == "Y"]
  if (median(y_counts) > 1) {
    sex <- "male"
    loginfo(paste0("Sex from calculation: ", sex))
  } else {
    sex <- "female"
    loginfo(paste0("Sex from calculation: ", sex))
  }
}

  # Get GC content and mappability
  gc <- round(fData(res)$gc)
  mappability <- round(fData(res)$mappability)

  # Fit a loess model to correct for GC content, mappability and replication timing
  if (correctReplication) {
    replication <- round(fData(res)$replication)
    condition <- QDNAseq:::binsToUse(res) & !is.na(gc) & !is.na(mappability) & is.finite(replication)
    smoothT <- loess(counts[, sample_index] ~ gc * mappability * replication)
  } else {
    condition <- QDNAseq:::binsToUse(res) & !is.na(gc) & !is.na(mappability)
    smoothT <- loess(counts[, sample_index] ~ gc * mappability)
  }
  # Get residuals from the loess fit
  fit <- matrix(NA, ncol = 1, nrow = nrow(res))
  dimnames(fit) <- dimnames(counts)
  # Fill in the fitted values
  fit[names(smoothT$residuals), sample_index] <- 2^smoothT$residuals
  assayDataElement(res, "copynumber") <- fit
  # Update binsToUse to exclude bins with NA values after fitting
  QDNAseq:::binsToUse(res) <- QDNAseq:::binsToUse(res) & !is.na(rowMeans(fit))
  return(res)
}

# Function to get median logr for each segment
get_median_logr <- function(segmentation, logr_chrom) {
  logr_segm_chrom <- rep(NA, length(logr_chrom))

  segs <- rle(segmentation$yhat)
  for (i in 1:length(segs$lengths)) {
    end <- cumsum(segs$lengths[1:i])
    end <- end[length(end)]
    start <- (end-segs$lengths[i]) + 1 # segs$lengths contains end points
    logr_segm_chrom[start:end] <- median(logr_chrom[start:end])
  }
  return(logr_segm_chrom)
}

# Custom segmentation function using PCF from Battenberg package
# with option to summarise segments by mean or median
segmentBinsPCF <- function(copyNumbersSmooth, segmentation_gamma, segmentation_kmin, segm_summary_method = 1) {
  condition <- QDNAseq:::binsToUse(copyNumbersSmooth)
  # Get logR values for bins to use
  logr <- log2(assayDataElement(copyNumbersSmooth, "copynumber")[condition, 1])
  # Get unique chromosomes
  chroms <- unique(fData(copyNumbersSmooth)$chromosome)
  # Initialise vector to hold segmented logR values
  logr_segm <- c()
  # Loop over chromosomes and segment each separately
  for (chrom in chroms) {
    # Build mask of bins on this chromosome that are also in 'condition'
    chrom_mask_condition <- (fData(copyNumbersSmooth)$chromosome == chrom)[condition]
    # Subset logr values for this chromosome
    logr_chrom <- logr[chrom_mask_condition]
    # Remove non-finite values
    logr_chrom <- logr_chrom[is.finite(logr_chrom)]
    # If no values, skip
    if (length(logr_chrom) == 0) {
      next
    }
    # If too few bins to run PCF reliably, fallback to a simple constant segment
    if (length(logr_chrom) < segmentation_kmin) {
      if (segm_summary_method == 1) {
        logr_segm <- c(logr_segm, rep(median(logr_chrom, na.rm = TRUE), length(logr_chrom)))
      } else if (segm_summary_method == 2) {
        logr_segm <- c(logr_segm, rep(median(logr_chrom, na.rm = TRUE), length(logr_chrom)))
      }
      next
    }
    # Perform segmentation using PCF using Battenberg functions, with safety checks
    sdev <- Battenberg:::getMad(logr_chrom, k = 25)
    if (!is.finite(sdev) || sdev == 0) sdev <- 1e-6
    chrom_segm <- tryCatch(
      Battenberg:::selectFastPcf(logr_chrom, segmentation_kmin, segmentation_gamma * sdev, TRUE),
      error = function(e) {
        warning(sprintf("selectFastPcf failed on chrom %s: %s -- falling back to median segments", chrom, e$message))
        list(yhat = rep(median(logr_chrom, na.rm = TRUE), length(logr_chrom)))
      }
    )
    # Summarise segments by mean or median
    if (segm_summary_method == 1) { # take the mean (yhat already provided by PCF)
      logr_segm <- c(logr_segm, chrom_segm$yhat)
    } else if (segm_summary_method == 2) { # take the median per-segment
      logr_segm <- c(logr_segm, get_median_logr(chrom_segm, logr_chrom))
    } else {
      warning("Unknown segmentation segment summary method supplied; using PCF yhat")
      logr_segm <- c(logr_segm, chrom_segm$yhat)
    }
  }
  # Store segmented logR values in the segmented assayData element
  temp <- assayDataElement(copyNumbersSmooth, "copynumber")
  temp[,] <- NA
  temp[QDNAseq:::binsToUse(copyNumbersSmooth), 1] <- 2^logr_segm
  assayDataElement(copyNumbersSmooth, "segmented") <- temp
  return(copyNumbersSmooth)
}

# Normalise read counts
normalisedCounts <- normaliseReadCounts(readCounts_tumour, readCounts_normal)
# Smooth outlier bins
copyNumbersSmooth <- smoothOutlierBins(normalisedCounts)
# Segment bins using PCF
copyNumbersSegmented <- segmentBinsPCF(copyNumbersSmooth, segmentation_gamma=segmentation_gamma, segmentation_kmin=segmentation_kmin, segm_summary_method=segm_summary_method)
# Normalize segmented bins
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
# Call copy numbers
copyNumbersCalled <- callBins(copyNumbersSegmented)


################
# Write count summary to file
# Write summarised readcount per bin file for tumour and normal, and add normalised counts

# Extract counts
tumour_counts_vec <- as.numeric(assayDataElement(readCounts_tumour, "counts"))
normal_counts_vec <- as.numeric(assayDataElement(readCounts_normal, "counts"))

# Extract normalised log2 ratio (logR)
normalised_logr_vec <- as.numeric(assayDataElement(normalisedCounts, "copynumber"))

# Combine into a single data frame
combined_counts_df <- data.frame(
  chromosome = bins$chromosome,
  start = bins$start,
  end = bins$end,
  tumour_readcount = tumour_counts_vec,
  normal_readcount = normal_counts_vec,
  normalised_log2ratio = normalised_logr_vec
)

write.table(
  combined_counts_df,
  file = file.path(output_dir, paste0(samplename, "_combined_readcounts.txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)


###############################################################################################
# Make results figures

# Make main figure
png(paste0(output_dir, "/", samplename, "_figure.png"), width = 1000, height = 400)
plot(copyNumbersCalled)
dev.off()

# Make raw data figures
processRawDataDefault <- function(readCounts, correctReplication) {
	copyNumbers <- correctBins(readCounts, correctReplication=correctReplication, method="median")
	copyNumbersNormalized <- normalizeBins(copyNumbers)
	copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
	copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
	copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
	return(copyNumbersSegmented)
}

# Run the Process raw funciton above for tumour and normal
tumour_raw <- processRawDataDefault(readCounts_tumour, correctReplication)
normal_raw <- processRawDataDefault(readCounts_normal, correctReplication)

png(paste0(output_dir, "/", samplename,  "_rawdata.png"), width = 1000, height = 800)
par(mfrow = c(2, 1))
plot(tumour_raw)
plot(normal_raw)
dev.off()

# make posthoc plot and segmentation output
get_segments <- function(x) {
  ns <- asNamespace("CGHbase")
  # get the .makeSegments function from CGHbase
  .makeSegments <- get(".makeSegments", envir = ns, mode = "function")
  # extract data from QDNAseq object
  condition <- QDNAseq:::binsToUse(x)	
  all.chrom <- QDNAseq:::chromosomes(x)
  # Get chromosome names as character vector
  if (is.integer(all.chrom)) # when x is a cghRaw, cghSeg, or cghCall object
    all.chrom <- as.character(all.chrom)
  chrom <- all.chrom[condition]
  # Extract copynumber and segmented data for bins to use
  condition <- QDNAseq:::binsToUse(x)
  copynumber <- as.data.frame(Biobase::assayDataElement(x, "copynumber")[condition, , drop=FALSE])

  # assuming a single sample here, so only one column
  segmented <- Biobase::assayDataElement(x, "segmented")[condition, 1]
  # convert to log2 space
  segmented <- QDNAseq:::log2adhoc(segmented)
  # make segments
  segment <- as.data.frame(.makeSegments(segmented, chrom))

  # add segment values to copynumber data frame
  copynumber$segment_values <- NA

  # Loop over segments and fill in segment values
  for (i in seq_len(nrow(segment))) {
    copynumber$segment_values[segment$start[i]:segment$end[i]] <- segment$values[i]
  }
  # Add chromosome, start and end positions
  colnames(copynumber)[1] <- "values"
  # Add logR values, chrom, start and end positions
  copynumber$logr_values <- QDNAseq:::log2adhoc(copynumber$values)
  # Extract chrom, start and end from rownames
  copynumber$chrom <- unlist(lapply(rownames(copynumber), function(x) unlist(stringr::str_split(x, ":"))[1]))
  # Extract start and end positions
  copynumber$start <- as.numeric(unlist(lapply(rownames(copynumber), function(x) unlist(stringr::str_split(unlist(stringr::str_split(x, ":"))[2], "-"))[1])))
  copynumber$end <- as.numeric(unlist(lapply(rownames(copynumber), function(x) unlist(stringr::str_split(unlist(stringr::str_split(x, ":"))[2], "-"))[2])))
  # order by chrom and start position
  copynumber$chrom <- factor(copynumber$chrom, levels=gtools::mixedsort(unique(copynumber$chrom)))

  return(copynumber)
}

# Extract segmented copy-number table for plotting and downstream analysis
copynumber <- get_segments(copyNumbersCalled)

# Build a segmentation data frame with statistical tests for gain/loss
get_segmentation_df <- function(copynumber, max_num_bins_test, logr_minimum_deviation, test_significance_threshold) {
  segments <- rle(copynumber$segment_values)
  segmentation <- data.frame(stringsAsFactors = FALSE)

  for (i in seq_along(segments$lengths)) {
    # Determine current segment start and end indices
    if (i == 1) {
      curr_segment_start_index <- 1
    } else {
      segment_end_index <- cumsum(segments$lengths[1:(i - 1)])
      curr_segment_start_index <- segment_end_index[length(segment_end_index)] + 1
    }

    if (i == length(segments$lengths)) {
      curr_segment_end_index <- nrow(copynumber)
    } else {
      segment_end_index <- cumsum(segments$lengths[1:i])
      curr_segment_end_index <- segment_end_index[length(segment_end_index)]
    }

    # Extract log-ratio values for the current segment
    logr_values <- log2(copynumber$values[curr_segment_start_index:curr_segment_end_index])

    # Subsample if segment is large to limit computational cost
    num_snps_sample <- ifelse(length(logr_values) > max_num_bins_test, max_num_bins_test, length(logr_values))

    # Test for gains and losses separately (comparing to a normal distribution shifted by the minimum deviation)
    p_value_gain <- t.test(
      sample(logr_values, num_snps_sample),
      rnorm(n = num_snps_sample, mean = 0 + logr_minimum_deviation, sd = sd(logr_values)),
      alternative = "greater"
    )$p.value

    p_value_loss <- t.test(
      sample(logr_values, num_snps_sample),
      rnorm(n = num_snps_sample, mean = 0 - logr_minimum_deviation, sd = sd(logr_values)),
      alternative = "less"
    )$p.value

    # Append row for this segment
    segmentation <- rbind(
      segmentation,
      data.frame(
        chrom = copynumber$chrom[curr_segment_start_index],
        start = copynumber$start[curr_segment_start_index],
        end = copynumber$end[curr_segment_end_index],
        value = copynumber$segment_values[curr_segment_start_index],
        p_value_gain = p_value_gain,
        p_value_loss = p_value_loss,
        stringsAsFactors = FALSE
      )
    )
  }

  # Multiple testing correction
  segmentation$p_value_gain_adj <- p.adjust(segmentation$p_value_gain, method = "bonferroni")
  segmentation$p_value_loss_adj <- p.adjust(segmentation$p_value_loss, method = "bonferroni")

  # Classification of segments
  segmentation$class <- "normal"
  is_gain <- segmentation$p_value_gain_adj < test_significance_threshold
  is_loss <- segmentation$p_value_loss_adj < test_significance_threshold
  segmentation$class[is_gain] <- "gain"
  segmentation$class[is_loss] <- "loss"

  # Sanity check: a segment shouldn't be both gain and loss
  if (any(is_gain & is_loss)) {
    segmentation$class[is_gain & is_loss] <- "error"
    warning("Found segment(s) classified as both gain and loss:")
    print(segmentation[is_gain & is_loss, ])
  }

  segmentation$len <- segmentation$end - segmentation$start
  return(segmentation)
}

# Generate segmentation df
segmentation <- get_segmentation_df(
  copynumber,
  max_num_bins_test = max_num_bins_test,
  logr_minimum_deviation = logr_minimum_deviation,
  test_significance_threshold = test_significance_threshold
)

# Estimate tumor purity per segment and overall; pass copynumber explicitly
get_purity_estimates <- function(segmentation, copynumber, sex) {
  segmentation$purity <- NA

  # Compute chromosome lengths from copynumber table
  chrom_lengths <- data.frame(stringsAsFactors = FALSE)
  for (chrom in unique(copynumber$chrom)) {
    pos_start <- min(copynumber[copynumber$chrom == chrom, "start"])
    pos_end <- max(copynumber[copynumber$chrom == chrom, "end"])
    chrom_lengths <- rbind(chrom_lengths, data.frame(chrom = chrom, start = pos_start, end = pos_end))
  }
  chrom_lengths$len <- chrom_lengths$end - chrom_lengths$start
  chrom_lengths$frac <- (chrom_lengths$len / 1000) / sum(chrom_lengths$len / 1000)

  # Estimate a normal/ploidy baseline (placeholder; here we set to 2)
  normal_ploidy <- 2

  # Helper transform and error functions used to fit purity
  mytransform <- function(x, rho, phi) (phi * 2^x - 2 * (1 - rho)) / rho

  geterrors <- function(rho, phi, meansSeg, weights, sds) {
    signal <- mytransform(meansSeg, rho, phi)
    mean(((round(signal) - signal) / sds)^2 * weights / 1000, na.rm = TRUE)
  }

  get_purity <- function(meansSeg, weights, normal_ploidy) {
    purs <- seq(0.05, 1, 0.01)
    ploidies <- normal_ploidy
    errs <- matrix(NA, length(purs), length(ploidies))
    rownames(errs) <- purs
    colnames(errs) <- ploidies
    for (pp in seq_along(purs)) {
      for (pl in seq_along(ploidies)) {
        errs[pp, pl] <- geterrors(rho = purs[pp], phi = ploidies[pl], meansSeg, weights, (pl / pp) * 0 + 1)
      }
    }
    mins <- arrayInd(which.min(errs), dim(errs))
    return(purs[mins[1]])
  }

  # Compute purity per segment
  for (i in seq_len(nrow(segmentation))) {
    segmentation$purity[i] <- get_purity(segmentation$value[i], 1, normal_ploidy)
  }

  # Overall purity estimate using segment lengths as weights
  segmentation$overall_purity <- get_purity(segmentation$value, segmentation$len, normal_ploidy)
  segmentation$assumed_ploidy <- normal_ploidy
  segmentation$len <- segmentation$len / 1e6  # convert length to Mb
  return(segmentation)
}

# Run purity estimation and write segmentation table to disk
segmentation <- get_purity_estimates(segmentation, copynumber, sex)
write.table(segmentation, file = paste0(output_dir, "/", samplename, "_segmentation.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Determine plotting scale based on observed max log-ratio
max_value <- max(abs(copynumber$logr_values), na.rm = TRUE)
if (max_value < 2) {
  max_y <- 2
  significance_bar_height <- 0.05
} else if (max_value < 3) {
  max_y <- 3
  significance_bar_height <- 0.05 * 1.5
} else {
  max_y <- 5
  significance_bar_height <- 0.05 * 2.5
}

# Custom plotting function for QDNAseq-style summary with significance bars
make_custom_plot <- function(copynumber, segmentation, max_y, plot_title = NULL, plot_subtitle = NULL, significance_bar_height = 0.01) {
  background <- data.frame(y = seq(max_y * -1, max_y, 1))
  p <- ggplot(copynumber) +
    geom_hline(data = background, mapping = aes(yintercept = y), colour = "black", alpha = 0.3) +
    geom_point(mapping = aes(x = start, y = logr_values), alpha = 0.5, size = 0.9, colour = "black") +
    geom_point(mapping = aes(x = start, y = segment_values), alpha = 1, size = 0.2, colour = "red") +
    facet_grid(~chrom, scales = "free_x", space = "free_x") +
    scale_x_continuous(expand = c(0, 0)) +
    ylim(-1 * max_y, max_y) +
    ylab("log2 ratio") +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(colour = "black", size = 18, face = "plain"),
      axis.title.y = element_text(colour = "black", size = 20, face = "plain"),
      strip.text.x = element_text(colour = "black", size = 16, face = "plain"),
      plot.title = element_text(colour = "black", size = 36, face = "plain", hjust = 0.5)
    )

  # Add significance bars for gains and losses
  if (any(segmentation$class == "loss")) {
    p <- p + geom_rect(
      data = segmentation[segmentation$class == "loss", ],
      mapping = aes(xmin = start, xmax = end, ymin = -(max_y - significance_bar_height), ymax = -max_y),
      fill = "blue"
    )
  }
  if (any(segmentation$class == "gain")) {
    p <- p + geom_rect(
      data = segmentation[segmentation$class == "gain", ],
      mapping = aes(xmin = start, xmax = end, ymin = (max_y - significance_bar_height), ymax = max_y),
      fill = "red"
    )
  }

  # Add title/subtitle if provided
  if (!is.null(plot_title) && !is.null(plot_subtitle)) {
    p <- p + ggtitle(bquote(atop(.(plot_title), atop(.(plot_subtitle), ""))))
  } else if (!is.null(plot_title)) {
    p <- p + ggtitle(plot_title)
  }
  return(p)
}

# Prepare plot title/subtitle and render PNG
plot_title <- dimnames(Biobase::assayDataElement(readCounts_tumour, "counts"))[[2]]
plot_subtitle <- paste0("Estimated purity: ", round(segmentation$overall_purity[1], 2),
                        "  Assumed ploidy: ", round(segmentation$assumed_ploidy, 2))
p <- make_custom_plot(copynumber, segmentation, max_y = max_y, plot_title = plot_title, plot_subtitle = plot_subtitle, significance_bar_height = significance_bar_height)

png(paste0(output_dir, "/", samplename, "_qdnaseq.png"), height = 500, width = 2000)
print(p)
dev.off()

# Prepare side/top summary plots for combined figure
segmentation$chrom <- factor(segmentation$chrom, levels = gtools::mixedsort(unique(segmentation$chrom)))
segmentation$frac <- segmentation$len / sum(segmentation$len)

breaks <- seq(0, 1, 0.02)
segmentation$purity_group <- cut(segmentation$purity, breaks)
purity_axis_breaks <- levels(segmentation$purity_group)[(breaks * 100) %% 20 == 0]
purity_axis_labels <- breaks[(breaks * 100) %% 20 == 0]

p_side_top <- ggplot(segmentation) +
  aes(x = purity_group, y = frac) +
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete(expand = c(0, 0.1), drop = FALSE, breaks = purity_axis_breaks, labels = purity_axis_labels) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(colour = "black", size = 14, face = "plain", angle = 90, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 14, face = "plain"),
    axis.title.y = element_blank(),
    strip.text.x = element_text(colour = "black", size = 16, face = "plain"),
    plot.title = element_text(colour = "black", size = 36, face = "plain", hjust = 0.5),
    plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.07, unit = "cm")
  )

offset <- 0.01
background <- data.frame(y = seq(0, 1, 0.1))
p_side_bottom <- ggplot(segmentation) +
  geom_hline(data = background, mapping = aes(yintercept = y), colour = "black", alpha = 0.3) +
  geom_rect(mapping = aes(xmin = start, xmax = end, ymin = purity - offset, ymax = purity + offset),
            alpha = 0.5, size = 0.9, colour = "black", fill = "black") +
  facet_grid(rows = vars(chrom), scales = "free_y", space = "free_y", switch = "y") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1 + offset), breaks = purity_axis_labels, labels = purity_axis_labels) +
  ylab("Fraction of cells") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(colour = "black", size = 18, face = "plain"),
    axis.title.x = element_text(colour = "black", size = 20, face = "plain"),
    strip.text = element_text(colour = "black", size = 14, face = "plain"),
    plot.title = element_text(colour = "black", size = 36, face = "plain", hjust = 0.5),
    plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.1, unit = "cm")
  ) +
  coord_flip()

# Combine raw data image and summary plots into one combined image
rawdata <- readPNG(paste0(output_dir, "/", samplename, "_rawdata.png"), native = FALSE)
png(paste0(output_dir, "/", samplename, "_combined.png"), height = 1600, width = 2000)
grid.arrange(
  arrangeGrob(
    rasterGrob(rawdata, interpolate = FALSE),
    arrangeGrob(p_side_top, p_side_bottom, ncol = 1, heights = c(1/5, 4/5)),
    ncol = 2, widths = c(4/5, 1/5)
  ),
  p,
  ncol = 1, nrow = 2, heights = c(0.7, 0.30)
)
dev.off()

# Save workspace for debugging / reproducibility
save.image(file = paste0(output_dir, "/", samplename, "_output.RData"))
