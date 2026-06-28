
**PROJECT-GROUP 2**

**Phase 0: Environment Setup**

During environment setup, storing the necessary files on a drive was necessary, so everyone can work on the Colab with the same files and libraries, overcoming some of the major problems when working on ephemeral instances like Google Colab, and so avoiding the need to reinstall every library from scratch when the runtime is restarted.

In detail, the analysis data are stored in a shared Google Drive folder and downloaded into the folder then set as baseDir, while the libraries are stored in a compressed file that is downloaded and decompressed, allowing to reduce the time necessary for the setup from 40+ minutes to less than a hundred seconds.

Additionally, a fixed seed (`6767`) was set to ensure reproducibility across different runs, especially for procedures involving random number generation, such as `preprocessSWAN`.
"""

# Setup R library + data
# precompiled R library
setup_persistent_r_libs <- function(file_id_targz = "1jxRXgic6slClamJrgKFRrJtALZBRzFfT") {
  lib_path <- "/content/R_libs_minfi"
  targz_path <- "/content/R_libs_minfi.tar.gz"


  if (!dir.exists(lib_path) || length(list.files(lib_path)) < 50) {
    message("library download ..")
    cmd <- paste0("gdown ", file_id_targz, " -O ", targz_path)
    system(cmd)


    message("Decompression...")
    dir.create(lib_path, showWarnings = FALSE, recursive = TRUE)
    system(paste0("tar -xzf ", targz_path, " -C /content"))
  } else {
    message("skip download")
  }


  .libPaths(c(lib_path, .libPaths()))
  library(minfi)
  library(qqman)
  library(gplots)
  library(ggrepel)

  set.seed(6767)



  message("Library ready")
}

# Download shared data folder
setup_colab_drive_folder <- function(folder_url) {
  folder_id <- sub(".*/folders/([^?]+).*", "\\1", folder_url)
  cmd <- paste0("gdown --folder https://drive.google.com/drive/folders/", folder_id)

  message("Data Folder download...")
  output <- system(cmd, intern = TRUE)
  cat(output, sep = "\n")

  all_dirs <- list.dirs(full.names = FALSE, recursive = FALSE)
  return(all_dirs[1])
}


setup_persistent_r_libs()

my_folder_url <- "https://drive.google.com/drive/folders/1KumBGpRTrO_EAoWVIwdYvo33fE5Cu0Bd?usp=drive_link"
data <- setup_colab_drive_folder(my_folder_url)
baseDir <- "./Input_Data"
print(list.files(baseDir))

"""# *Analysis Pipeline*

**Phase 1 & 2: loading data and dataframe creation**

The intensity files and the sample metadata sheet were loaded.
From the files, the RGChannelSet object was stored as a binary RGset.RData file.
From this, two dataframes (Red and Green) were created to store the red and green light intensities.
"""

#1 Loading raw data
load("/content/Input_Data/Illumina450Manifest_clean.RData")
SampleSheet <- read.csv(file.path(baseDir, "SampleSheet_Report_II.csv"),
               stringsAsFactors = FALSE)

targets <- read.metharray.sheet(baseDir)
RGset <- read.metharray.exp(targets = targets)
save(RGset, file = "RGset.RData")

# 2 Dataset

Red <- data.frame(getRed(RGset))
Green <- data.frame(getGreen(RGset))

"""**Phase 3: Fluorescence extraction and probe inspection**

Targeted inspection of the group's assigned probe address (10804411) was performed to verify signal behaviour.


Based on the generated probe summary table, the target probe is an Illumina Infinium Type II (consequently receiving the signal from both Red and Green channels).


Across almost all samples, the Green fluorescence is consistently higher than the Red fluorescence.
Compared to the others, sample 5 exhibits low intensities in both channels. This could be due to:
Lower DNA yield
Technical Issues during microarray analysis (highlighting the need for Quality Checks (Phase 5))


Sample 6 shows close values for the two channels, suggesting a roughly balanced mix of methylated and unmethylated alleles in that specific subject.

"""

# 3 Ispection of a specific probe address
target_address <- "10804411" #Our probe adress



probe_type <- Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID == target_address |
                                          Illumina450Manifest_clean$AddressB_ID == target_address, "Infinium_Design_Type"][1]
probe_color <- Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID == target_address |
                                           Illumina450Manifest_clean$AddressB_ID == target_address, "Color_Channel"][1]


probe_data <- data.frame(
  Sample = colnames(Red),
  Red_fluor = as.numeric(Red[target_address, ]),
  Green_fluor = as.numeric(Green[target_address, ]),
  Type = probe_type,
  Color = ifelse(probe_type == "I", probe_color, "Both")
)
print("--- Point 3 & 4: Probe Table ---")
print(probe_data)

"""**Phase 4 & 5: Quality check**

**QCplot**

With the plotQC() function the median log2 values for Methylated and Unmethylated channels across samples was computed and plotted against a threshold line, where samples falling below this line are technically compromised.

Only 3 samples are classified as "good" (black circles). The remaining 5 samples are explicitly flagged as "bad" (red circles with indices). Notably, samples 2 and 5 are extreme outliers, highlighting sample-wide signal failure (low DNA yield or hybridization failure).


**Negative Control Intensity Verification**

The strip plot visualizes the behavior of internal, engineered negative control probes, showing that their intensities are low and sable.
This means that variance seen in the QC plot is caused by sample quality decay or binding failures, rather than a malfunctioning laser scanner or chip chemical contamination.

**Decection of pValues**

We use the assigned threshold of p > 0.01 to flag failed probes, meaning probes whose signal is not significantly louder than the background noise floor.
For each sample, the output table reports probes that have a detection p-value higher than the threshold.
Particularly, we observe that:
Sample 5 holds 4,555 failed probes, representing the highest failure rate. This is consistent with previous results (the low fluorescence value observed in Phase 3)
"""

# 4
MSet.raw <- preprocessRaw(RGset)

# 5 Quality control
options(warn = -1)
qc <- getQC(MSet.raw)
plotQC(qc)

controlStripPlot(RGset, controls = "NEGATIVE")

detP <- detectionP(RGset)
failed <- detP > 0.01 #Our threshold
summary(failed)
head(failed)

output <- apply(failed, 2, table)
output
options(warn = 0)

"""**Phase 6: Group subdivision**

Intensity signals were converted into Beta values and the total dataset was split according to the clinical classifications stored in the sample metadata sheet, yielding two distinct subsets: `beta_df_control` (CTRL) and `beta_df_disease` (DIS).
There’s no difference in the dimensionality of the two groups, meaning that the same number of individual CpG loci was recognised between the two subgroups, each of 4 individuals.

"""

# 6 Row data extraction e group suddivision
beta <- getBeta(MSet.raw)
M <- getM(MSet.raw)

# division based on the group (control vs unhealthy)
beta_df <- data.frame(getBeta(MSet.raw))

# Control group (CTRL)
beta_df_control <- beta_df[, SampleSheet$Group == 'CTRL']
dim(beta_df_control)

# unhealthy group (DIS)
beta_df_disease <- beta_df[, SampleSheet$Group == 'DIS']
dim(beta_df_disease)

"""
To visualize methylation profiles, the mean Beta value per probe was calculated across the 4 individuals within each clinical group using the `apply` function. Kernel density estimations were then computed to generate a comparative distribution plot.

The plot presents **bimodal distribution pattern**, with two peaks: around **0.06** (representing heavily unmethylated CpG sites) and a broader peak near **0.85** (representing highly methylated genomic areas like gene bodies).

The global distributions show highly overlapping geometry if the CTRl (blue) and DIs (red) groups, however subtle variations are visible. These difference may be be attributed to true pathology or reflect technical artefacts (dye bias, overall sample quality).
"""

mean_of_beta_c <- apply(beta_df_control, 1, mean, na.rm = TRUE)
mean_of_beta_d <- apply(beta_df_disease, 1, mean, na.rm = TRUE)

# Calcolo della densità e Visualizzazione dei profili di metilazione

d_mean_of_beta_c <- density(mean_of_beta_c, na.rm = TRUE)
d_mean_of_beta_d <- density(mean_of_beta_d, na.rm = TRUE)

# Generazione del grafico comparativo
plot(d_mean_of_beta_c, col = 'blue', main = "Density plot of Mean Beta Values")
lines(d_mean_of_beta_d, col = 'red')

"""**Phase 7: SWAN Normalization**

Raw methylation data were normalized using the **Subset-quantile Within Array Normalization (SWAN)** method implemented in the **minfi** package. SWAN is specifically designed for Illumina HumanMethylation450 BeadChip arrays and corrects the technical bias introduced by the different chemistries of **Type I** and **Type II** probes while preserving biological variation. After normalization, Beta values were extracted from the normalized methylation set and compared with the raw data.
"""

MSet.norm <- preprocessSWAN(RGset)
beta_norm <- getBeta(MSet.norm)

beta_raw <- getBeta(MSet.raw)
annot <- getAnnotation(MSet.raw)
typeI_idx  <- rownames(beta_raw) %in% rownames(annot)[annot$Type == "I"]
typeII_idx <- rownames(beta_raw) %in% rownames(annot)[annot$Type == "II"]


short_names <- sub("_.*", "", colnames(beta_raw))

par(mfrow=c(2,3), cex.main=1.2)

plot(density(apply(beta_raw[typeI_idx, ], 1, mean, na.rm=TRUE), na.rm=TRUE), col='darkgreen', main="Raw: Mean Beta by Type")
lines(density(apply(beta_raw[typeII_idx, ], 1, mean, na.rm=TRUE), na.rm=TRUE), col='purple')

plot(density(apply(beta_raw[typeI_idx, ], 1, sd, na.rm=TRUE), na.rm=TRUE), col='darkgreen', main="Raw: SD Beta by Type")
lines(density(apply(beta_raw[typeII_idx, ], 1, sd, na.rm=TRUE), na.rm=TRUE), col='purple')

boxplot(beta_raw, names=short_names, col=ifelse(targets$Group == 'CTRL', 'blue', 'red'), main="Raw: Boxplot", las=2, cex.axis=0.7)

plot(density(apply(beta_norm[typeI_idx, ], 1, mean, na.rm=TRUE), na.rm=TRUE), col='darkgreen', main="Norm: Mean Beta by Type")
lines(density(apply(beta_norm[typeII_idx, ], 1, mean, na.rm=TRUE), na.rm=TRUE), col='purple')

plot(density(apply(beta_norm[typeI_idx, ], 1, sd, na.rm=TRUE), na.rm=TRUE), col='darkgreen', main="Norm: SD Beta by Type")
lines(density(apply(beta_norm[typeII_idx, ], 1, sd, na.rm=TRUE), na.rm=TRUE), col='purple')

boxplot(beta_norm, names=short_names, col=ifelse(targets$Group == 'CTRL', 'blue', 'red'), main="Norm: Boxplot", las=2, cex.axis=0.7)

par(mfrow=c(1,1))

"""To evaluate the effect of normalization, six comparison control plots were generated. For both the raw and normalized datasets, the distribution of the mean Beta values and the standard deviation of Beta values were compared separately for Type I and Type II probes. In addition, boxplots of Beta values were produced for each sample, with CTRL samples shown in **blue** and DIS samples in **red**.


The density plots show that the characteristic bimodal methylation pattern is preserved after SWAN normalization, indicating that the overall biological structure of the dataset was maintained. The distributions of Type I and Type II probes remain comparable before and after normalization, while slight adjustments in their shapes suggest a reduction of technical differences associated with probe chemistry.


The density distributions of the Beta standard deviations remain highly comparable before and after normalization, indicating that the overall variability of methylation measurements was preserved. Likewise, the boxplots display very similar distributions across samples, with only minor changes in the median and interquartile ranges after normalization. No evident global differences between the CTRL and DIS groups are observed, suggesting that SWAN normalization preserved the overall methylation profiles while reducing probe-specific technical bias.

**Phase 8: Principal Component Analysis (PCA)**

Principal Component Analysis (PCA) was performed on the normalized Beta value matrix to explore the main sources of variation within the dataset after SWAN normalization. PCA is a dimensionality reduction technique that summarizes the variability of high-dimensional data into a small number of principal components, allowing the visualization of potential clustering according to biological or technical variables. In this analysis, the first two principal components explained **41.0%** and **20.8%** of the total variance, respectively.
"""

# 8 PCA on Normalized Data
pca_res <- prcomp(t(na.omit(beta_norm)))

pca_data <- data.frame(
  PC1   = pca_res$x[, 1],
  PC2   = pca_res$x[, 2],
  Group = targets$Group,
  Sex   = targets$Sex,
  Batch = as.factor(targets$Slide)
)

var_explained <- summary(pca_res)$importance[2, 1:2] * 100

# Plot 1: Group
plot(pca_data$PC1, pca_data$PC2,
     col  = ifelse(pca_data$Group == "CTRL", "blue", "red"),
     pch  = 19,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA by Group")
legend("topright", legend = c("CTRL", "DIS"),
       col = c("blue", "red"), pch = 19)

# Plot 2: Sex
plot(pca_data$PC1, pca_data$PC2,
     col  = ifelse(pca_data$Sex == "Male", "darkgreen", "purple"),
     pch  = 19,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA by Sex")
legend("bottomright", legend = c("M", "F"),
       col = c("darkgreen", "purple"), pch = 19)

# Plot 3: Batch
batches <- levels(pca_data$Batch)
batch_colors <- rainbow(length(batches))
plot(pca_data$PC1, pca_data$PC2,
     col  = batch_colors[as.numeric(pca_data$Batch)],
     pch  = 19,
     xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
     ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
     main = "PCA by Batch")
legend("topright", legend = batches,
       col = batch_colors, pch = 19, cex = 0.7)

"""Three PCA plots were generated by colouring the samples according to **clinical group (CTRL/DIS)**, **sex**, and **processing batch**.

The PCA coloured by clinical group does not show a clear separation between CTRL and DIS samples. Although some disease samples are positioned farther from the main cluster, the two groups partially overlap, indicating that disease status is not the dominant source of variation within the dataset.

The PCA coloured by sex reveals a more evident clustering pattern. Male and female samples tend to occupy different regions of the PCA space, suggesting that sex contributes substantially to the observed methylation variability and represents an important biological factor influencing the dataset.

The PCA coloured by processing batch does not reveal a clear grouping of samples according to the Illumina BeadChip slide. Samples processed on different slides are broadly intermingled, suggesting that no evident batch effect is present after normalization.

Overall, the PCA indicates that the normalized methylation data retain the biological variability of the samples while showing no evident clustering driven by technical batch effects. Among the variables examined, **sex appears to explain a greater proportion of the observed variability than disease status**.

**Phase 9: Differential Methylation Analysis (Mann-Whitney Test)**

To identify differentially methylated positions (DMPs) between the CTRL and DIS groups, we used the Mann-Whitney U test (Wilcoxon rank-sum test), as assigned to Group 2.
This non-parametric approach was chosen because beta values tend to pile up near 0 and 1, which is a bimodal shape that breaks the normality assumption required by tests like the t-test.
Since the Mann-Whitney works on ranked data rather than raw values, it is a more appropriate and reliable choice for this type of data. The test was run probe-wise across the entire normalized beta matrix (beta_norm from Step 7), meaning one independent test was performed per CpG probe, each time comparing the 4 CTRL beta values against the 4 DIS beta values.
"""

# 9 - Differential Methylation Analysis (Mann-Whitney test)
ctrl_idx <- targets$Group == 'CTRL'
dis_idx  <- targets$Group == 'DIS'

My_mannwhitney_function <- function(x) {
  wilcox <- wilcox.test(x[ctrl_idx], x[dis_idx],
                        exact = FALSE, correct = FALSE)
  return(wilcox$p.value)
}

p_values <- apply(beta_norm, 1, My_mannwhitney_function)

# Dataframe
dmp_results <- data.frame(
  probe   = names(p_values),
  p_value = p_values,
  row.names = names(p_values)
)
dmp_results <- dmp_results[order(dmp_results$p_value), ]

print("--- Point 9: Top Differentially Methylated Probes ---")
print(head(dmp_results, n = 10))

"""Looking at the results, the lowest p-value we obtained was approximately 0.021, and, a large number of probes share this exact value. Despite this might look questionable at first, it is actually expected.
With only 4 samples per group, the Mann-Whitney test is severely limited in the range of p-values it can produce. Since the test works by ranking all 8 observations together and checking whether one group tends to rank higher than the other, the number of possible rank arrangements with groups of this size is very small (meaning many probes inevitably map to the same p-value). More critically, the minimum p-value the test can ever reach is mathematically fixed around 0.021, regardless of how biologically different CTRL and DIS actually are.
This is already a red flag for the correction step ahead, since surviving genome-wide multiple testing correction would require p-values on the order of 10⁻⁷ or lower across ~485,000 probes.

**Phase 10: Multiple Testing Correction**

Since we tested roughly 485,000 probes at once, running that many statistical tests simultaneously means that a large number of probes could appear significant purely by chance (even with no real biology behind them). To address this, we applied two multiple testing correction methods to the nominal p-values from Step 9: Bonferroni and Benjamini-Hochberg (BH). A significance threshold of 0.05 was used for all comparisons.

Bonferroni is the strictest of the two: it multiplies every p-value by the total number of tests (~485,000), making sure the probability of getting even a single false positive across all tests stays below 0.05.

BH takes a more tolerant approach: instead of trying to eliminate all false positives, it controls the proportion of false positives among the probes you actually call significant (the False Discovery Rate), which makes it a more practical choice for large-scale genomic studies where some tolerance for error is acceptable.
"""

# 10 - Multiple test correction
corrected_pValues_BH   <- p.adjust(dmp_results$p_value, "BH")
corrected_pValues_Bonf <- p.adjust(dmp_results$p_value, "bonferroni")
dmp_results_corrected  <- data.frame(dmp_results, corrected_pValues_BH, corrected_pValues_Bonf)
head(dmp_results_corrected)


boxplot(dmp_results_corrected[, c("p_value", "corrected_pValues_BH", "corrected_pValues_Bonf")],
        names = c("Nominal_pValue", "BH_correction", "Bonf_correction"),
        main = "Distribution of p-values")


sig_nominal    <- dim(dmp_results_corrected[dmp_results_corrected$p_value             <= 0.05, ])[1]
sig_bh         <- dim(dmp_results_corrected[dmp_results_corrected$corrected_pValues_BH   <= 0.05, ])[1]
sig_bonferroni <- dim(dmp_results_corrected[dmp_results_corrected$corrected_pValues_Bonf <= 0.05, ])[1]

cat("Significant probes (threshold 0.05):\n")
cat("Nominal pValues:     ", sig_nominal,    "\n")
cat("Bonferroni adjusted: ", sig_bonferroni, "\n")
cat("BH (FDR) adjusted:   ", sig_bh,         "\n")

"""Looking at the boxplot, the three distributions paint a clear picture. The nominal p-values (left box) are spread broadly across the 0–1 range with a median around 0.57, which is exactly what weexpect when most probes are not truly differentially methylated.

The BH-corrected values (center box) shift upward and compress, clustering around a median of 0.80, which is still far above 0.05, so no probe survives even this more lenient correction.

The Bonferroni-corrected values (right) barely form a proper box. They collapse into a flat line sitting at exactly 1.0, because multiplying every p-value by ~485,000 pushes them all to the ceiling with zero variation left.

The final numbers confirm what the boxplot already suggested: at the nominal threshold, 21,950 probes appear significant, but after both Bonferroni and BH correction, 0 probes survive. As we already anticipated at the end of Step 9, this outcome is completely expected. For a probe to pass Bonferroni correction it would need a nominal p-value below roughly 1.03 × 10⁻⁷ (0.05 / ~485,000) is a value that is simply out of reach with only 4 samples per group. The 21,950 nominally significant probes are therefore most likely a product of statistical noise and low power rather than genuine differential methylation between CTRL and DIS.

## Phase 11: Volcano Plot and Manhattan Plot
To visualize the results of the differential methylation analysis, a volcano plot and a Manhattan plot were produced from the nominal p-values obtained with the Mann-Whitney test.



### Volcano Plot
The volcano plot displays the methylation difference between groups for each probe on the x axis against statistical significance on the y axis. With only 8 samples (4 CTRL + 4 DIS) the Mann-Whitney test can produce just a few p-values, so the probes align on discrete horizontal bands. Only the two highest bands stay above the significance line, corresponding to the two smallest p-values the Mann-Whitney test can produce with this sample size that is a direct visual illustration of the discreteness that limits the test's resolution.
"""

# 11 - Volcano Plot

beta_ctrl <- beta_norm[, ctrl_idx]
beta_dis  <- beta_norm[, dis_idx]

mean_beta_ctrl <- apply(beta_ctrl, 1, mean, na.rm = TRUE)
mean_beta_dis  <- apply(beta_dis,  1, mean, na.rm = TRUE)

delta <- mean_beta_dis - mean_beta_ctrl


toVolcPlot <- data.frame(
  delta    = delta[dmp_results_corrected$probe],
  neg_log10_p = -log10(dmp_results_corrected$p_value),
  row.names = dmp_results_corrected$probe
)

head(toVolcPlot)

plot(toVolcPlot[, 1], toVolcPlot[, 2], pch = 16, cex = 0.5,
     xlab = "Delta Beta (DIS - CTRL)",
     ylab = "-log10(p-value)",
     main = "Volcano Plot")
abline(h = -log10(0.05), col = "red")

# Evidenzia i punti con delta assoluto > 0.1 e p < 0.05
toHighlight <- toVolcPlot[abs(toVolcPlot[, 1]) > 0.1 & toVolcPlot[, 2] > (-log10(0.05)), ]
points(toHighlight[, 1], toHighlight[, 2], pch = 16, cex = 0.7, col = "red")

"""### Manhattan Plot
The Manhattan plot maps the nominal p-values agaist their genomic position coloured by chromosome. All points sit on a few low horizontal levels and the plot appears flat, with no peaks standing out. No probe approaches genome-wide significance and there is no regional enrichment across chromosomes, consistent with the absence of robust differential methylation at this sample size.
"""

# 11 - Manhattan Plot

library(qqman)


dmp_results_corrected <- data.frame(
  IlmnID = rownames(dmp_results_corrected),
  dmp_results_corrected
)


dmp_annotated <- merge(dmp_results_corrected, Illumina450Manifest_clean, by = "IlmnID")
dim(dmp_annotated)
head(dmp_annotated)

input_Manhattan <- dmp_annotated[, colnames(dmp_annotated) %in% c("IlmnID", "CHR", "MAPINFO", "p_value")]
head(input_Manhattan)


order_chr <- c("1","2","3","4","5","6","7","8","9","10","11",
               "12","13","14","15","16","17","18","19","20","21","22","X","Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR, levels = order_chr)
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)
table(input_Manhattan$CHR)

manhattan(input_Manhattan,
          snp = "IlmnID", chr = "CHR", bp = "MAPINFO", p = "p_value",
          genomewideline = FALSE,
          suggestiveline = FALSE,
          col = rainbow(24))

"""## Phase 12: Heatmap
We display the top 100 probes ranked by nominal p-value as a heatmap with hierarchical clustering on both rows "probes" and columns "samples". The colour bar above the columns marks the group, CTRL = blue and DIS = red. Significant probes drops to zero after correction, for this reason the ranking is based on nominal p-values. The samples cluster cleanly into the two groups even if there is lack of formal significance, with one block of hypermethylated and one of hypomethylated probes, showing that the nominal top-ranked probes still capture a clear CTRL/DIS separation.




"""

# 12 - Heatmap Top 100
library(gplots)

# Top 100 for nominal P-value
input_heatmap <- as.matrix(beta_norm[dmp_results_corrected$IlmnID[1:100], ])

colorbar <- ifelse(targets$Group == "CTRL", "blue", "red")

col_heatmap <- colorRampPalette(c("green", "black", "red"))(100)

short_names_heatmap <- sub("_.*", "", colnames(input_heatmap))

heatmap.2(input_heatmap,
          col          = col_heatmap,
          Rowv         = TRUE, Colv = TRUE,
          dendrogram   = "both",
          key          = TRUE,
          ColSideColors = colorbar,
          labCol        = short_names_heatmap,
          density.info = "none",
          trace        = "none",
          scale        = "none",
          symm         = FALSE,
          margins       = c(10, 8),
          main         = "Top 100 DMPs - Complete linkage")

"""## Comparative analysis (extra)
Here we ran the analysis again with a t-test just to see how the choice of test changes the results, even though Mann-Whitney was our assigned test. The clearest difference is in the volcano plot: the Mann-Whitney one shows discrete horizontal bands and tops out around −log10(p) ≈ 1.7, while the t-test works on the actual beta values and gives a continuous spread reaching much higher values up to ≈ 5.4. The t-test assumes roughly normal data, which doesn't hold well for beta values and with only 4 samples per group a few outliers can easily inflate significance. Still, as with Mann-Whitney, no probe survives multiple-testing correction.
"""

# ============================================================
# Comparative Analysis - t-test (only for visual confrontation)
# ============================================================

# Step 9 - t-test
My_ttest_function <- function(x) {
  t_test <- t.test(x[ctrl_idx], x[dis_idx])
  return(t_test$p.value)
}

p_values_ttest <- apply(beta_norm, 1, My_ttest_function)

dmp_results_ttest <- data.frame(
  probe   = names(p_values_ttest),
  p_value = p_values_ttest,
  row.names = names(p_values_ttest)
)
dmp_results_ttest <- dmp_results_ttest[order(dmp_results_ttest$p_value), ]

print(head(dmp_results_ttest, n = 10))

# Step 10 - Correzione
corrected_pValues_BH_ttest   <- p.adjust(dmp_results_ttest$p_value, "BH")
corrected_pValues_Bonf_ttest <- p.adjust(dmp_results_ttest$p_value, "bonferroni")
dmp_results_corrected_ttest  <- data.frame(
  dmp_results_ttest,
  corrected_pValues_BH   = corrected_pValues_BH_ttest,
  corrected_pValues_Bonf = corrected_pValues_Bonf_ttest
)

cat("Significant probes (threshold 0.05):\n")
cat("Nominal pValues:     ", sum(dmp_results_corrected_ttest$p_value                < 0.05), "\n")
cat("Bonferroni adjusted: ", sum(dmp_results_corrected_ttest$corrected_pValues_Bonf < 0.05), "\n")
cat("BH (FDR) adjusted:   ", sum(dmp_results_corrected_ttest$corrected_pValues_BH   < 0.05), "\n")

# Volcano Plot
mean_beta_ctrl <- apply(beta_norm[, ctrl_idx], 1, mean, na.rm = TRUE)
mean_beta_dis  <- apply(beta_norm[, dis_idx],  1, mean, na.rm = TRUE)
delta <- mean_beta_dis - mean_beta_ctrl

toVolcPlot_ttest <- data.frame(
  delta       = delta[dmp_results_corrected_ttest$probe],
  neg_log10_p = -log10(dmp_results_corrected_ttest$p_value),
  row.names   = dmp_results_corrected_ttest$probe
)

plot(toVolcPlot_ttest$delta, toVolcPlot_ttest$neg_log10_p,
     pch = 16, cex = 0.5,
     xlab = "Delta Beta (DIS - CTRL)",
     ylab = "-log10(p-value)",
     main = "Volcano Plot (t-test)")
abline(h = -log10(0.05), col = "red")
toHighlight <- toVolcPlot_ttest[abs(toVolcPlot_ttest$delta) > 0.1 &
                                  toVolcPlot_ttest$neg_log10_p > (-log10(0.05)), ]
points(toHighlight$delta, toHighlight$neg_log10_p, pch = 16, cex = 0.7, col = "red")

# Manhattan Plot
library(qqman)
dmp_results_corrected_ttest$IlmnID <- rownames(dmp_results_corrected_ttest)
dmp_annotated_ttest <- merge(dmp_results_corrected_ttest, Illumina450Manifest_clean, by = "IlmnID")

input_Manhattan_ttest <- dmp_annotated_ttest[, colnames(dmp_annotated_ttest) %in%
                                               c("IlmnID", "CHR", "MAPINFO", "p_value")]
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11",
               "12","13","14","15","16","17","18","19","20","21","22","X","Y")
input_Manhattan_ttest$CHR <- factor(input_Manhattan_ttest$CHR, levels = order_chr)
input_Manhattan_ttest$CHR <- as.numeric(input_Manhattan_ttest$CHR)

manhattan(input_Manhattan_ttest,
          snp = "IlmnID", chr = "CHR", bp = "MAPINFO", p = "p_value",
          genomewideline = FALSE,
          suggestiveline = FALSE,
          annotatePval = 0.00001, col = rainbow(24),
          main = "Manhattan Plot (t-test)")

# Heatmap Top 100
library(gplots)
input_heatmap_ttest <- as.matrix(beta_norm[dmp_results_corrected_ttest$probe[1:100], ])
colorbar <- ifelse(targets$Group == "CTRL", "blue", "red")
col_heatmap <- colorRampPalette(c("green", "black", "red"))(100)

short_names_ttest <- sub("_.*", "", colnames(input_heatmap_ttest))

heatmap.2(input_heatmap_ttest,
          col           = col_heatmap,
          Rowv          = TRUE, Colv = TRUE,
          dendrogram    = "both",
          key           = TRUE,
          ColSideColors = colorbar,
          labCol        = short_names_ttest,
          density.info  = "none",
          trace         = "none",
          scale         = "none",
          symm          = FALSE,
          margins       = c(10, 8), #add margin
          main          = "Top 100 DMPs - t-test")