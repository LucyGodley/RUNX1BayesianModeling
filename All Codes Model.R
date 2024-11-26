# Load required libraries
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("The 'readxl' package is required. Install it using install.packages('readxl').")
}

library(readxl)

# File paths
input_file <- "data/RUNX1_raw_data.xlsx"
output_file <- "output/Allcodes.csv"

# Ensure output directory exists
dir.create("output", showWarnings = FALSE)

# Check if input file exists
if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

# Load data
raw_data <- read_excel(input_file)

# Data preprocessing: Convert "." to NA for numerical columns
numeric_columns <- c("gnomAD_maf", "revel_score", "p_curations_dif_residue", 
                     "lp_curations_dif_residue", "nt_num", "phyloP_score", 
                     "spliceAI_score", "transactivation_percent", 
                     "pathogenic_secondary_assay_count", 
                     "benign_secondary_assay_count", "proband_count", 
                     "nonsegregation_count", "segregation_count", 
                     "denovo_assumed_count", "denovo_proven_count", 
                     "gnomAD_homozygotes_count")
raw_data[numeric_columns] <- lapply(raw_data[numeric_columns], function(x) {
  as.numeric(ifelse(x == ".", NA, x))
})

# Helper function to apply conditions
apply_code <- function(condition, true_value = 1, false_value = 0) {
  ifelse(condition, true_value, false_value)
}

# Compute Codes
# Population codes
raw_data$BA1 <- apply_code(raw_data$gnomAD_maf > 0.0015)
raw_data$BS1 <- apply_code(raw_data$gnomAD_maf >= 0.00015 & raw_data$gnomAD_maf <= 0.0015)
raw_data$PM2_supporting <- apply_code(!is.na(raw_data$gnomAD_maf) & raw_data$gnomAD_maf == 0)

# Missense codes
raw_data$PP3 <- apply_code(raw_data$revel_score >= 0.88)
raw_data$BP4 <- apply_code(raw_data$revel_score <= 0.5)
raw_data$PM5_strong <- apply_code(raw_data$p_curations_dif_residue >= 2)
raw_data$PM5 <- apply_code(raw_data$p_curations_dif_residue == 1)
raw_data$PM5_supporting <- apply_code(
  raw_data$lp_curations_dif_residue >= 1 | (!is.na(raw_data$aa_alt) & raw_data$aa_alt == "Ter" & raw_data$nt_num > 98)
)

# PS1 codes
raw_data$PS1 <- apply_code(raw_data$p_curations_same_residue >= 1)
raw_data$PS1_moderate <- apply_code(raw_data$lp_curations_same_residue >= 1)

# Null variant codes
raw_data$PVS1 <- apply_code(raw_data$aa_alt == "Ter" & raw_data$nt_num >= 98 & raw_data$nt_num <= 916)
raw_data$PVS1_strong <- apply_code(raw_data$aa_alt == "Ter" & raw_data$nt_num >= 917 & raw_data$nt_num <= 1440)

# Synonymous and intronic variant codes
raw_data$BP7 <- apply_code(
  raw_data$aa_ref == raw_data$aa_alt & raw_data$phyloP_score <= 2 & raw_data$spliceAI_score <= 0.2
)

# Experimental codes
pm1_residues <- c("R107", "K110", "A134", "R162", "R166", "S167", "R169", "G170", "K194", "T196", "D198", "R201", "R204")
raw_data$PM1 <- apply_code(paste0(raw_data$aa_ref, raw_data$aa_num) %in% pm1_residues)
raw_data$PM1_supporting <- apply_code(!raw_data$PM1 & raw_data$aa_num >= 89 & raw_data$aa_num <= 204)

# Transactivation percent-based conditions
raw_data$PS3 <- apply_code(raw_data$transactivation_percent < 20 & raw_data$pathogenic_secondary_assay_count >= 1 & raw_data$PVS1 == 0)
raw_data$PS3_Moderate <- apply_code(raw_data$PS3 == 0 & raw_data$pathogenic_secondary_assay_count >= 2)
raw_data$BS3 <- apply_code(raw_data$transactivation_percent >= 80 & raw_data$transactivation_percent <= 115 & raw_data$benign_secondary_assay_count >= 1)

# Proband and segregation codes
raw_data$PS4 <- apply_code(raw_data$proband_count >= 4)
raw_data$PP1_Strong <- apply_code(raw_data$segregation_count >= 7)

# Classification
raw_data$PVS_count <- rowSums(raw_data[, c("PVS1")], na.rm = TRUE)
raw_data$PS_count <- rowSums(raw_data[, c("PS1", "PM5_strong", "PVS1_strong", "PS3", "PS4", "PP1_Strong")], na.rm = TRUE)
raw_data$PM_count <- rowSums(raw_data[, c("PS1_moderate", "PM5", "PM1")], na.rm = TRUE)
raw_data$PP_count <- rowSums(raw_data[, c("PM2_supporting", "PP3", "PM5_supporting", "PM1_supporting")], na.rm = TRUE)
raw_data$BS_count <- rowSums(raw_data[, c("BA1", "BS1", "BS3")], na.rm = TRUE)
raw_data$BP_count <- rowSums(raw_data[, c("BP4", "BP7")], na.rm = TRUE)

raw_data$OP <- 350^((raw_data$PP_count / 8) + (raw_data$PM_count / 4) + (raw_data$PS_count / 2) + raw_data$PVS_count - ((raw_data$BP_count / 8) + raw_data$BS_count / 2))
raw_data$Post_P <- (raw_data$OP * 0.1) / (((raw_data$OP - 1) * 0.1) + 1)

raw_data$Classification <- ifelse(
  raw_data$Code_Count < 2 & raw_data$BP_count == 1 & raw_data$BS_count == 0 & raw_data$BA1 == 0,
  "VUS",
  ifelse(raw_data$Post_P > 0.99, "Pathogenic",
         ifelse(raw_data$Post_P <= 0.99 & raw_data$Post_P > 0.9, "Likely Pathogenic",
                ifelse(raw_data$Post_P >= 0.10 & raw_data$Post_P <= 0.9, "VUS",
                       ifelse(raw_data$Post_P > 0.001 & raw_data$Post_P < 0.10, "Likely Benign", "Benign"))))
)

# Save output
write.csv(raw_data, file = output_file, row.names = FALSE)
cat("Processing complete. Output saved to", output_file, "\n")
