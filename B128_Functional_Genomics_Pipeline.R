R
# ==============================================================================
# Functional Genomics of Limosilactobacillus fermentum B1 28 (ATCC 14931)
# Analysis Pipeline
# Author: Mandana AmeliMojarad
# ==============================================================================

library(tidyverse)
library(readxl)
library(ggplot2)
library(stringr)
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(GO.db)
library(AnnotationDbi)

# Set working directory (Update this to your local path)
setwd("/Users/apple/Desktop/LF/ncbi_dataset/ncbi_dataset/data/GCF_005341425.1")
dir.create("figures", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)

# ==============================================================================
# 1. GENOME LOADING & CDS EXTRACTION
# ==============================================================================
cat("Step 1: Extracting CDS sequences...\n")

genome <- readDNAStringSet("GCF_005341425.1_ASM534142v1_genomic.fna")
names(genome) <- "B128_chromosome"

gff <- import("genomic.gff")
gff_df <- as.data.frame(gff)

# Extract CDS respecting strand
cds_gr <- gff_df %>%
  filter(type == "CDS") %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)
seqlevels(cds_gr) <- names(genome)

cds_list <- vector("list", length(cds_gr))
strand_vec <- as.character(strand(cds_gr))

for (i in seq_along(cds_gr)) {
  seq_range <- start(cds_gr)[i]:end(cds_gr)[i]
  if (strand_vec[i] == "+") {
    cds_list[[i]] <- genome[[1]][seq_range]
  } else {
    cds_list[[i]] <- reverseComplement(genome[[1]][seq_range])
  }
}

cds_seqs <- DNAStringSet(cds_list)
names(cds_seqs) <- cds_gr$ID

# Translate & clean
protein_seqs <- translate(cds_seqs)
protein_seqs <- sub("\\*$", "", protein_seqs)
internal_stops <- grepl("\\*", protein_seqs)
clean_proteins <- AAStringSet(protein_seqs[!internal_stops])

writeXStringSet(clean_proteins, "B128_proteins.faa")
cat("Total clean proteins:", length(clean_proteins), "\n")

# ==============================================================================
# 2. ENZYME CLASSIFICATION
# ==============================================================================
cat("Step 2: Classifying enzymes...\n")

cds_annotations <- gff_df %>%
  filter(type == "CDS") %>%
  dplyr::select(protein_id, locus_tag, gene, product) %>%
  distinct()

enzyme_classified <- cds_annotations %>%
  mutate(
    enzyme_class = case_when(
      str_detect(product, regex("dehydrogenase|oxidase|reductase|peroxidase|catalase|superoxide|oxidoreductase", ignore_case = TRUE)) ~ "Oxidoreductase",
      str_detect(product, regex("transferase|kinase|phosphorylase", ignore_case = TRUE)) ~ "Transferase",
      str_detect(product, regex("hydrolase|peptidase|protease|nuclease|phosphatase|glycosidase|esterase|thioesterase", ignore_case = TRUE)) ~ "Hydrolase",
      str_detect(product, regex("lyase|decarboxylase|synthase", ignore_case = TRUE)) ~ "Lyase",
      str_detect(product, regex("isomerase|mutase", ignore_case = TRUE)) ~ "Isomerase",
      str_detect(product, regex("ligase|synthetase", ignore_case = TRUE)) ~ "Ligase",
      TRUE ~ "Non-enzyme"
    )
  )

enzyme_distribution <- enzyme_classified %>%
  count(enzyme_class, sort = TRUE) %>%
  mutate(percentage = round(n / sum(n) * 100, 1))

write_csv(enzyme_distribution, "tables/Table2_enzyme_distribution.csv")

# Figure 2
p2 <- enzyme_distribution %>%
  filter(enzyme_class != "Non-enzyme") %>%
  ggplot(aes(x = reorder(enzyme_class, -n), y = n, fill = enzyme_class)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = paste0(n, "\n(", percentage, "%)")), vjust = -0.3) +
  theme_minimal() + theme(legend.position = "none") +
  labs(title = "Enzyme Class Distribution", x = "", y = "Count")
ggsave("figures/Fig2_enzyme_class.png", p2, width = 8, height = 5)

# ==============================================================================
# 3. ANTICANCER CANDIDATE SCREENING
# ==============================================================================
cat("Step 3: Screening anticancer candidates...\n")

anticancer_patterns <- c(
  "ATPase", "dehydrogenase", "oxidase", "reductase", "peroxidase", "catalase", "superoxide", 
  "glutathione", "nuclease", "endonuclease", "exonuclease", "ribonuclease", "peptidase", 
  "protease", "arginine deiminase", "glycosidase", "beta-glucosidase", "beta-galactosidase", 
  "phosphatase", "kinase", "transferase", "methyltransferase"
)

tier1_anticancer <- enzyme_classified %>%
  filter(enzyme_class != "Non-enzyme") %>%
  filter(str_detect(product, regex(paste(anticancer_patterns, collapse = "|"), ignore_case = TRUE))) %>%
  mutate(
    anticancer_mechanism = case_when(
      str_detect(product, regex("superoxide dismutase|catalase|peroxidase|glutathione|thioredoxin", ignore_case = TRUE)) ~ "ROS Reduction",
      str_detect(product, regex("dehydrogenase|oxidase|reductase", ignore_case = TRUE)) ~ "ROS Generation",
      str_detect(product, regex("endonuclease|exonuclease|ribonuclease|nuclease", ignore_case = TRUE)) ~ "DNA/RNA Damage",
      str_detect(product, regex("arginine deiminase|asparaginase", ignore_case = TRUE)) ~ "Amino Acid Depletion",
      str_detect(product, regex("beta-glucosidase|beta-galactosidase|glycosidase", ignore_case = TRUE)) ~ "Glycoside Hydrolysis",
      str_detect(product, regex("peptidase|protease", ignore_case = TRUE)) ~ "Protein Degradation",
      str_detect(product, regex("ATPase", ignore_case = TRUE)) ~ "Energy Metabolism",
      str_detect(product, regex("phosphatase|kinase", ignore_case = TRUE)) ~ "Signal Disruption",
      str_detect(product, regex("transferase|synthase|synthetase", ignore_case = TRUE)) ~ "Metabolic Modulation",
      TRUE ~ "Other"
    )
  )

write_csv(tier1_anticancer, "tables/Anticancer_Candidates.csv")

# Figure 3A (Mechanisms)
mechanism_dist <- tier1_anticancer %>% count(anticancer_mechanism, sort = TRUE) %>% mutate(pct = round(n/sum(n)*100,1))
p3a <- ggplot(mechanism_dist, aes(x = reorder(anticancer_mechanism, n), y = n, fill = anticancer_mechanism)) +
  geom_bar(stat = "identity") + coord_flip() + theme_minimal() + theme(legend.position = "none") +
  labs(title = "Anticancer Mechanisms", x = "", y = "Count")
ggsave("figures/Fig3A_anticancer_mechanisms.png", p3a, width = 9, height = 5)

# Figure 3B (Enzyme classes of anticancer)
p3b <- tier1_anticancer %>% count(enzyme_class) %>%
  ggplot(aes(x = reorder(enzyme_class, -n), y = n, fill = enzyme_class)) +
  geom_bar(stat = "identity") + theme_minimal() + theme(legend.position = "none") +
  labs(title = "Anticancer Enzyme Classes", x = "", y = "Count")
ggsave("figures/Fig3B_anticancer_classes.png", p3b, width = 7, height = 5)

# ==============================================================================
# 4. GO, KEGG, & COG ANNOTATION (EggNOG-mapper)
# ==============================================================================
cat("Step 4: Parsing EggNOG annotation...\n")

eggnog <- read_excel("emapper_annotations.xlsx")

# --- 4A. GO ANNOTATION ---
go_long <- eggnog %>%
  filter(!is.na(GOs), GOs != "-") %>%
  separate_rows(GOs, sep = ",") %>%
  mutate(GOs = str_trim(GOs)) %>%
  filter(GOs != "", GOs != "-")

go_counts <- go_long %>% count(GOs, sort = TRUE) %>% top_n(20, n)

go_info <- suppressMessages(AnnotationDbi::select(GO.db, keys = go_counts$GOs, columns = c("TERM", "ONTOLOGY"), keytype = "GOID"))
go_annotated <- go_counts %>% left_join(go_info, by = c("GOs" = "GOID")) %>% filter(!is.na(ONTOLOGY))
write_csv(go_annotated, "tables/Table_GO_enrichment.csv")

# --- 4B. KEGG ANNOTATION ---
kegg_long <- eggnog %>%
  filter(!is.na(KEGG_Pathway), KEGG_Pathway != "-") %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  mutate(KEGG_Pathway = str_trim(KEGG_Pathway)) %>%
  filter(KEGG_Pathway != "", KEGG_Pathway != "-")

kegg_long$pathway_group <- "Other metabolic"
kegg_long$pathway_group[str_sub(kegg_long$KEGG_Pathway, 1, 3) %in% c("000","001","002","003","004","005","006","007","008","009")] <- "Global maps"
kegg_long$pathway_group[str_sub(kegg_long$KEGG_Pathway, 1, 3) %in% c("010","011","012")] <- "Carbohydrate metabolism"
kegg_long$pathway_group[str_sub(kegg_long$KEGG_Pathway, 1, 3) %in% c("020","021","022")] <- "Amino acid metabolism"
kegg_long$pathway_group[str_sub(kegg_long$KEGG_Pathway, 1, 3) %in% c("030","031","032")] <- "Nucleotide metabolism"

kegg_distribution <- kegg_long %>% count(pathway_group, sort = TRUE) %>% mutate(pct = round(n/sum(n)*100,1))
write_csv(kegg_distribution, "tables/Table_KEGG_subcategories.csv")

# --- 4C. COG ANNOTATION ---
cog_long <- eggnog %>%
  filter(!is.na(COG_category), COG_category != "-") %>%
  separate_rows(COG_category, sep = "") %>%
  filter(COG_category != "") %>%
  count(COG_category, sort = TRUE)
write_csv(cog_long, "tables/Table_COG_categories.csv")

cog_pie <- cog_long %>% mutate(pct = round(n/sum(n)*100, 1))
p8 <- ggplot(cog_pie, aes(x = "", y = n, fill = COG_category)) +
  geom_bar(stat = "identity", width = 1, color = "white") + coord_polar("y", start = 0) +
  theme_void() + labs(title = "COG Categories (1,642 proteins)")
ggsave("figures/Fig8_COG_pie.png", p8, width = 10, height = 8)

# ==============================================================================
# 5. SAVE WORKSPACE
# ==============================================================================
save.image("B128_Full_Analysis_Workspace.RData")
cat("Pipeline complete. All tables, figures, and workspace saved.\n")
