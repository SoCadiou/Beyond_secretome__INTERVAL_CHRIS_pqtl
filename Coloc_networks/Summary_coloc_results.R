rm(list=ls(all=TRUE))

time.start <- Sys.time()

library(readr)
library(dplyr)
library(scales)
library(ggrepel)
library(data.table)
library(readxl)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann <- read_excel("data/supplementary_table_2.xlsx", sheet = "ST2", skip = 1)
colocalization_results <- fread("data/combined_colocalization_results_mock_dataset.csv")

dim(colocalization_results)


colocalization_results <- colocalization_results %>%
  mutate(
    trait_a_new = pmin(trait_a, trait_b),
    trait_b_new = pmax(trait_a, trait_b),
    locus_a_new = if_else(trait_a == trait_a_new, locus_a, locus_b),
    locus_b_new = if_else(trait_a == trait_a_new, locus_b, locus_a),
    target_a_new = if_else(trait_a == trait_a_new, target_a, target_b),
    target_b_new = if_else(trait_a == trait_a_new, target_b, target_a),
    top_cond_a_new =if_else(trait_a == trait_a_new, top_cond_a, top_cond_b),
    top_cond_b_new =if_else(trait_a == trait_a_new, top_cond_b, top_cond_a),
    top_freq_a_new =if_else(trait_a == trait_a_new, top_freq_a, top_freq_b),
    top_freq_b_new =if_else(trait_a == trait_a_new, top_freq_b, top_freq_a),
    top_freq_geno_a_new =if_else(trait_a == trait_a_new, top_freq_geno_a, top_freq_geno_b),
    top_freq_geno_b_new =if_else(trait_a == trait_a_new, top_freq_geno_b, top_freq_geno_a),
    top_mlog10pC_a_new =if_else(trait_a == trait_a_new, top_mlog10pC_a, top_mlog10pC_b),
    top_mlog10pC_b_new =if_else(trait_a == trait_a_new, top_mlog10pC_b, top_mlog10pC_a)
  ) %>%
  select(-trait_a,-trait_b,-locus_a,-locus_b,-target_a,
         -target_b, -top_cond_a, -top_cond_b, -top_freq_a,- top_freq_b,
         -top_freq_geno_a, -top_freq_geno_b, -top_mlog10pC_a, -top_mlog10pC_b
         ) %>% 
  rename(trait_a = trait_a_new, trait_b = trait_b_new,
         locus_a=locus_a_new, locus_b = locus_b_new,
         target_a = target_a_new, target_b = target_b_new, 
         top_cond_a= top_cond_a_new, top_cond_b= top_cond_b_new, top_freq_a=  top_freq_a_new,
         top_freq_b=top_freq_b_new, top_freq_geno_a=top_freq_geno_a_new, top_freq_geno_b=top_freq_geno_b_new,
         top_mlog10pC_a=top_mlog10pC_a_new,  top_mlog10pC_b=top_mlog10pC_b_new)

pair_tested <- colocalization_results %>%
  distinct(trait_a, trait_b, locus_a, locus_b, .keep_all = TRUE)

colocalization_results_filtered <- colocalization_results %>%
  filter(PP.H4.abf > 0.8) 

pair_significant <- colocalization_results_filtered %>%
  distinct(trait_a, trait_b, locus_a, locus_b, .keep_all = TRUE)

dim(pair_significant)

colocalization_results_filtered <- colocalization_results_filtered %>%
  mutate(
    chr_locus_a = sub("chr(\\d+)_.*", "\\1", locus_a),
    start_locus_a = sub("chr\\d+_(\\d+)_.*", "\\1", locus_a),
    end_locus_a = sub("chr\\d+_\\d+_(\\d+)", "\\1", locus_a)
  )

colocalization_results_filtered <- colocalization_results_filtered %>%
  mutate(
    chr_locus_b = sub("chr(\\d+)_.*", "\\1", locus_b),
    start_locus_b = sub("chr\\d+_(\\d+)_.*", "\\1", locus_b),
    end_locus_b = sub("chr\\d+_\\d+_(\\d+)", "\\1", locus_b)
  )


colnames(mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann)
hotspots_df <- mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann[,c("CHR","locus_START_END_37","SeqID", "cis_or_trans","UniProt_ID", "hotspot", "full_hotspot_gene_window", "uniprot_new_in_somascan7k_vs5k" , "uniprot_match")]
colnames(hotspots_df)<-c("chr", "locus_START_END_37", "phenotype_id", "cis_or_trans","UniProt_ID", "hotspot", "full_hotspot_gene_window", "new_somamer", "uniprot_match")
hotspots_df$start<- sub("^[^_]*_([^_]*)_.*$", "\\1",hotspots_df$locus_START_END_37)
hotspots_df$end  <- sub("^.*_", "", hotspots_df$locus_START_END_37)
hotspots_df$chr <- as.character(hotspots_df$chr)
hotspots_df$start <- as.character(hotspots_df$start)
hotspots_df$end <- as.character(hotspots_df$end)
hotspots_df <- hotspots_df[,c("chr", "start","end", "phenotype_id", "cis_or_trans","UniProt_ID", "hotspot", "full_hotspot_gene_window", "new_somamer", "uniprot_match")]
dim_dataset <- dim(colocalization_results_filtered)[2]
colocalization_results_filtered<-as.data.frame(colocalization_results_filtered)

colocalization_results_filtered <- merge(colocalization_results_filtered, hotspots_df, by.x = c("chr_locus_a", "start_locus_a", "end_locus_a", "trait_a"), by.y = c("chr", "start","end", "phenotype_id"))
table(colocalization_results_filtered$hotspot)
colnames(colocalization_results_filtered)[dim_dataset + 1:6] <- c("cis_or_trans_locus_a","UniProt_ID_locus_a","hotspot_locus_a","full_hotspot_gene_window_locus_a", "new_somamer_a","uniprot_match_a" )
dim_dataset <- dim(colocalization_results_filtered)[2]
colocalization_results_filtered <- merge(colocalization_results_filtered, hotspots_df, by.x = c("chr_locus_b", "start_locus_b", "end_locus_b", "trait_b"), by.y = c("chr", "start","end", "phenotype_id"))
colnames(colocalization_results_filtered)[dim_dataset + 1:6] <- c("cis_or_trans_locus_b","UniProt_ID_locus_b","hotspot_locus_b","full_hotspot_gene_window_locus_b","new_somamer_b","uniprot_match_b")

dir.create("results", showWarnings = FALSE, recursive = TRUE)

###########################################
# SUMMARY OF THE COLOCALIZATION RESULTS
###########################################

pair_significant <- colocalization_results_filtered %>%
  distinct(trait_a, trait_b, locus_a, locus_b, .keep_all = TRUE)

table(colocalization_results_filtered$hotspot_locus_a, colocalization_results_filtered$hotspot_locus_b)
# The discordant pair refers to one locus that overlap the end of a hotspot region and the other not 

count_df <- colocalization_results_filtered %>%
  group_by(UniProt_ID_locus_a,UniProt_ID_locus_b) %>%
  filter(UniProt_ID_locus_a != UniProt_ID_locus_b)   %>%
  summarise(count = n())
dim(count_df)
 
table(colocalization_results_filtered$cis_or_trans_locus_a, colocalization_results_filtered$cis_or_trans_locus_b)

cis_trans_coloc_df <- colocalization_results_filtered %>%
  filter(cis_or_trans_locus_a == "cis" & cis_or_trans_locus_b == "trans" | cis_or_trans_locus_a == "trans" & cis_or_trans_locus_b == "cis" )

cis_trans_df <- colocalization_results_filtered %>%
  filter(cis_or_trans_locus_a == "cis" ) %>%
  select(locus_a,trait_a)%>%
  rename(locus = locus_a, trait = trait_a)

cis_trans_df_2 <- colocalization_results_filtered %>%
  filter(cis_or_trans_locus_b == "cis" ) %>%
  select(locus_b,trait_b)%>%
  rename(locus = locus_b, trait = trait_b)

cis_trans_df <- rbind(cis_trans_df, cis_trans_df_2)
cis_trans_pairs <- cis_trans_df %>%
  distinct(trait,locus, .keep_all = TRUE)

cis_df <- colocalization_results_filtered %>%
   filter(cis_or_trans_locus_a == "cis" & cis_or_trans_locus_b == "cis")
 
cis_df <- cis_df %>%
  group_by(UniProt_ID_locus_a, 
           # symbol_a, Protein.names_a, 
           UniProt_ID_locus_b #, 
           # symbol_b, Protein.names_b
           ) %>%
  filter(UniProt_ID_locus_a != UniProt_ID_locus_b)%>%
  summarise(count = n())
dim(cis_df)
 
try_df <- colocalization_results_filtered %>%
  filter(cis_or_trans_locus_a == "cis" & cis_or_trans_locus_b == "cis") %>%
  filter(UniProt_ID_locus_a != UniProt_ID_locus_b)

cis_df$hotspot_locus_a <- FALSE
cis_df$hotspot_locus_b <- FALSE
for(i in 1:nrow(cis_df)){
    temp <- try_df %>%
      filter(UniProt_ID_locus_a == cis_df$UniProt_ID_locus_a[i] & UniProt_ID_locus_b==cis_df$UniProt_ID_locus_b[i])
    if (TRUE %in% temp$hotspot_locus_a){
      cis_df$hotspot_locus_a[i] <- TRUE
    }
    if (TRUE %in% temp$hotspot_locus_b){
      cis_df$hotspot_locus_b[i] <- TRUE
    }
  }
table(cis_df$hotspot_locus_a, cis_df$hotspot_locus_b)

write.csv(cis_df, "results/Colocalizing_cis.csv")

colocalization_results_temp <- colocalization_results_filtered %>%
  filter(hotspot_locus_a == TRUE & hotspot_locus_b == TRUE)

colocalization_results_temp <- colocalization_results_temp %>%
  filter((cis_or_trans_locus_a == "cis" & new_somamer_a == TRUE) |(cis_or_trans_locus_b == "cis" & new_somamer_b == TRUE))

count_df <- colocalization_results_temp %>%
  group_by(full_hotspot_gene_window_locus_a, chr_locus_a) %>%
  summarise(
    colocalization_count = n(),
    new_somamer_uniprot = paste(unique(c(
      UniProt_ID_locus_a[cis_or_trans_locus_a == "cis" & new_somamer_a],
      UniProt_ID_locus_b[cis_or_trans_locus_b == "cis" & new_somamer_b]
    )), collapse = "|"),
    .groups = "drop"
  )

write.csv(count_df, "results/Colocalizing_cis_new_somamer.csv")

colocalization_results_temp <- colocalization_results_filtered %>%
  filter(hotspot_locus_a == TRUE & hotspot_locus_b == TRUE)

colocalization_results_temp <- colocalization_results_temp %>%
  filter((cis_or_trans_locus_a == "cis" & uniprot_match_a == "NO") |(cis_or_trans_locus_b == "cis" &  uniprot_match_b == "NO"))

count_df <- colocalization_results_temp %>%
  group_by(full_hotspot_gene_window_locus_a, chr_locus_a) %>%
  summarise(
    colocalization_count = n(),
    new_cis_uniprot = paste(unique(c(
      UniProt_ID_locus_a[cis_or_trans_locus_a == "cis" & uniprot_match_a == "NO"],
      UniProt_ID_locus_b[cis_or_trans_locus_b == "cis" & uniprot_match_b == "NO"]
    )), collapse = "|"),
    .groups = "drop"
  )

write.csv(count_df, "results/Colocalizing_cis_new_uniprot.csv")

#############################################################

colocalization_results_filtered <- colocalization_results_filtered %>%
  filter(hotspot_locus_a == TRUE | hotspot_locus_b == TRUE)
fwrite(colocalization_results_filtered, "results/Colocalization_results_filtered.csv", sep=",")


count_df <- colocalization_results_filtered %>%
  group_by(full_hotspot_gene_window_locus_a, chr_locus_a) %>%
  summarise(colocalization_count = n())%>%
  filter(full_hotspot_gene_window_locus_a != "[]")
write.csv(count_df, "results/Hotspots_df.csv")

LB_in_hotsposts <- mapped_LB_gp_ann_va_ann_bl_ann_collapsed_hf_ann %>%
  filter(hotspot == TRUE)
write.csv(LB_in_hotsposts, "results/LB_in_hotsposts.csv")

cat("\nSummary_coloc_results.R completed in ", round(as.numeric(difftime(Sys.time(), time.start, units = "secs")), 2), " seconds.\n", sep = "")

