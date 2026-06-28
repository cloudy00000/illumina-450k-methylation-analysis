# Setup R library + data

# Precompiled R library
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


mean_of_beta_c <- apply(beta_df_control, 1, mean, na.rm = TRUE)
mean_of_beta_d <- apply(beta_df_disease, 1, mean, na.rm = TRUE)

# Calcolo della densità e Visualizzazione dei profili di metilazione

d_mean_of_beta_c <- density(mean_of_beta_c, na.rm = TRUE)
d_mean_of_beta_d <- density(mean_of_beta_d, na.rm = TRUE)

# Generazione del grafico comparativo
plot(d_mean_of_beta_c, col = 'blue', main = "Density plot of Mean Beta Values")
lines(d_mean_of_beta_d, col = 'red')


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
          margins       = c(10, 8),
          main          = "Top 100 DMPs - t-test")


