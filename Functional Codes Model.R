# Load required libraries
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("The 'readxl' package is required. Install it using install.packages('readxl').")
}

library(readxl)

# File paths
input_file <- "data/RUNX1_raw_data.xlsx"
output_file <- "output/Functional.csv"

# Ensure output directory exists
dir.create("output", showWarnings = FALSE)

# Check if input file exists
if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}

# Load data
raw_data <- read_excel(input_file)

# Compute Codes
# Population codes
raw_data$BA1 <- ifelse(raw_data$gnomAD_maf != "." & raw_data$gnomAD_maf > 0.0015, 1, 0)
raw_data$BS1 <- ifelse(raw_data$gnomAD_maf != "." & raw_data$gnomAD_maf >= 0.00015 & raw_data$gnomAD_maf <= 0.0015, 1, 0)
raw_data$PM2_supporting <- ifelse(!is.na(raw_data$gnomAD_maf) & raw_data$gnomAD_maf == 0, 1, 0)

# Missense Codes
raw_data$PP3 <- ifelse(raw_data$revel_score != "." & raw_data$revel_score >= 0.88, 1, 0)
raw_data$BP4 <- ifelse(raw_data$revel_score != "." & raw_data$revel_score <= 0.5, 1, 0)
raw_data$PM5_strong <- ifelse(raw_data$p_curations_dif_residue != "." & raw_data$p_curations_dif_residue >= 2, 1, 0)
raw_data$PM5 <- ifelse(raw_data$p_curations_dif_residue != "." & raw_data$p_curations_dif_residue == 1, 1, 0)
raw_data$PM5_supporting <- ifelse((raw_data$lp_curations_dif_residue != "." & raw_data$lp_curations_dif_residue >= 1) | 
                                    (raw_data$aa_alt == "Ter" & raw_data$nt_num != "." & raw_data$nt_num > 98), 1, 0)
raw_data$PS1 <- ifelse(raw_data$p_curations_same_residue != "." & raw_data$p_curations_same_residue >= 1, 1, 0)
raw_data$PS1_moderate <- ifelse(raw_data$lp_curations_same_residue != "." & raw_data$lp_curations_same_residue >= 1, 1, 0)

# Null Variant Codes
raw_data$PVS1 <- ifelse(raw_data$aa_alt == "Ter" & raw_data$nt_num != "." & raw_data$nt_num >= 98 & raw_data$nt_num <= 916, 1, 0)
raw_data$PVS1_strong <- ifelse(raw_data$aa_alt == "Ter" & raw_data$nt_num != "." & raw_data$nt_num >= 917 & raw_data$nt_num <= 1440, 1, 0)

# Synonymous and Intronic Variant Codes
raw_data$BP7 <- ifelse(raw_data$aa_ref != "." & raw_data$aa_alt != "." & raw_data$aa_ref == raw_data$aa_alt & 
                         raw_data$phyloP_score != "." & raw_data$phyloP_score <= 2 & 
                         raw_data$spliceAI_score != "." & raw_data$spliceAI_score <= 0.2, 1, 0)

# Experimental Codes
pm1_residues <- c("R107", "K110", "A134", "R162", "R166", "S167", "R169", "G170", "K194", "T196", "D198", "R201", "R204")
raw_data$PM1 <- ifelse(paste0(raw_data$aa_ref, raw_data$aa_num) %in% pm1_residues, 1, 0)
raw_data$PM1_supporting <- ifelse(!raw_data$PM1 & raw_data$aa_num != "." & raw_data$aa_num >= 89 & raw_data$aa_num <= 204, 1, 0)
raw_data$BS3 <- ifelse(raw_data$AM_score != "." & raw_data$AM_score <= 0.425, 1, 0)
raw_data$PS3 <- ifelse(raw_data$AM_score != "." & raw_data$AM_score >= 0.90, 1, 0)

# Applied Code Counts
raw_data$PVS_count <- rowSums(raw_data[, c("PVS1")], na.rm = TRUE)
raw_data$PS_count <- rowSums(raw_data[, c("PS1", "PM5_strong", "PVS1_strong", "PS3")], na.rm = TRUE)
raw_data$PM_count <- rowSums(raw_data[, c("PS1_moderate", "PM5", "PM1")], na.rm = TRUE)
raw_data$PP_count <- rowSums(raw_data[, c("PM2_supporting", "PP3", "PM5_supporting", "PM1_supporting")], na.rm = TRUE)
raw_data$BS_count <- rowSums(raw_data[, c("BA1", "BS1", "BS3")], na.rm = TRUE)
raw_data$BP_count <- rowSums(raw_data[, c("BP4", "BP7")], na.rm = TRUE)
raw_data$Code_Count <- raw_data$PVS_count + raw_data$PS_count + raw_data$PM_count + raw_data$PP_count + raw_data$BS_count + raw_data$BP_count

# Compute OP and Post_P
raw_data$OP <- 350^((raw_data$PP_count/8) + (raw_data$PM_count/4) + (raw_data$PS_count/2) + raw_data$PVS_count - ((raw_data$BP_count/8) + raw_data$BS_count/2))
raw_data$Post_P <- (raw_data$OP * 0.1) / (((raw_data$OP - 1) * 0.1) + 1)

# Classify variants
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
