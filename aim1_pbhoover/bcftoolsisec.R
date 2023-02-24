library(tidyverse)
library(reshape2)
set.seed(4)

figs_dir <- "/home/ebishop/thesis_projects/figures/bcftools_norm_isec_het/"
start_dir <- "/home/ebishop/thesis_projects/aim1_pbhoover/arrow_comparison_resource_use/bcftools_norm_isec_hetero"
subdirs <- list.dirs(start_dir)[2:33]  # Exclude parent and input vcf directories

###########################################
#   COMPARE PBHOOVER AND ARROW VARIANTS   #
###########################################

# Parsing bcftools isec output

# 0000.vcf = pbhoover-only variants
# 0001.vcf = arrow-only variants
# 0002.vcf = pbhoover variants shared by arrow
# Confirmed that 0000 and 0002 add up to all the non-HETERO pbhoover variants


parse_vcf_isec <- function(vcf_path, skip, prefix) {
  vcf_raw <- read_tsv(vcf_path, skip = skip, show_col_types = FALSE)
  # Totals
  tot_vars <- nrow(vcf_raw)
  tot_ins <- vcf_raw %>% filter(nchar(REF) < nchar(ALT)) %>% nrow()
  tot_del <- vcf_raw %>% filter(nchar(REF) > nchar(ALT)) %>% nrow()
  tot_snp <- vcf_raw %>% filter(nchar(REF) == nchar(ALT)) %>% nrow()
  # Proportions
  ins <- tot_ins / tot_vars
  del <- tot_del / tot_vars
  snp <- tot_snp / tot_vars
  # Put together
  df <- data.frame(
    total = c(tot_vars),
    ins = c(ins),
    del = c(del),
    snp = c(snp)
  )
  colnames(df) <- paste0(prefix, colnames(df))
  return(df)
}

out_df <- NULL
for (d in subdirs) {
  isolate <- basename(d)
  # Load pbhoover-only, arrow-only, and shared variants
  pbhvcf <- parse_vcf_isec(paste0(d, "/0000.vcf"), skip = 12, "pb_")
  arrvcf <- parse_vcf_isec(paste0(d, "/0001.vcf"), skip = 14, "ar_")
  shavcf <- parse_vcf_isec(paste0(d, "/0002.vcf"), skip = 12, "sh_")
  # Put together
  df <- bind_cols(pbhvcf, arrvcf, shavcf) %>%
    mutate(isolate = c(isolate), .before = pb_total) %>%
    mutate(pbh_private = pb_total / (pb_total + sh_total)) %>%
    mutate(pbh_shared = sh_total / (pb_total + sh_total))
  out_df <- bind_rows(out_df, df)
}

# Save object in current working directory
saveRDS(out_df, "out_df.rds")

###############
##  FIGURES  ##
###############

# Boxplot showing percent of PBHoover variants that are shared with arrow

sharedfig <- ggplot(out_df, aes(
  x = NA, y = pbh_shared)) +
  theme_classic() +
  geom_boxplot(width = 0.4) +
  geom_jitter(aes(x = NA, y = pbh_shared), width = 0.1) +
  stat_boxplot(geom = 'errorbar', width = 0.1) +
  ylab('Proportion') +
  ylim(0, 1) +
  ggtitle('Proportion of non-HETERO PBHoover \nVariants Shared with Arrow') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 13)
  )
sharedfig

# Save
ggsave(
  sharedfig,
  filename = paste0(figs_dir, 'propshared.png'),
  device = 'png',
  width = 5,
  height = 4,
  dpi = 300
)

summary(out_df$pbh_shared)

# Boxplot showing what proportion of shared variants are INS/DEL/SNP

longshared <- out_df %>%
  select(isolate, sh_snp, sh_ins, sh_del) %>%
  melt(id = 'isolate') %>%
  mutate(variable = factor(variable, levels = c("sh_del", "sh_ins", "sh_snp")))

sharedprop <-
  ggplot(longshared, aes(x = factor(variable), y = value)) +
  theme_classic() +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  xlab('Variant') +
  scale_x_discrete(labels = c(
    "sh_del" = "Deletion",
    "sh_ins" = "Insertion",
    "sh_snp" = "SNP"
  )) +
  ylab('Proportion') +
  ylim(0, 1) +
  ggtitle('Variants Shared Between PBHoover and Arrow') +
  theme(plot.title = element_text(hjust = 0.5))
sharedprop

# Save
ggsave(
  sharedprop,
  filename = paste0(figs_dir, 'sharedtypes.png'),
  device = 'png',
  width = 5,
  height = 4,
  dpi = 300
)

longshared %>% 
  filter(variable == "sh_snp") %>%
  summary('value')

# Boxplot showing what proportion of non-shared variants are INS/DEL/SNP for pbh and arrow

longprivate <- out_df %>%
  select(isolate, pb_snp, pb_ins, pb_del, ar_snp, ar_ins, ar_del) %>%
  melt(id = 'isolate') %>%
  mutate(varcaller = case_when(
    grepl("pb_", variable) ~ "PBHoover",
    TRUE ~ "Arrow"
  )) %>%
  mutate(type = substr(variable, 4, 6))

privateprop <-
  ggplot(longprivate, aes(x = type, y = value, color = varcaller)) +
  theme_classic() +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1)) +
  xlab('Variant') +
  scale_x_discrete(labels = c(
    "snp" = "SNP",
    "ins" = "Insertion",
    "del" = "Deletion"
  )) +
  ylab('Proportion') +
  ylim(0, 1) +
  ggtitle('Proportions of Variants Unique to PBHoover or Arrow') +
  labs(colour = "Unique to") +
  theme(plot.title = element_text(hjust = 0.5))
privateprop

# Save
ggsave(
  privateprop,
  filename = paste0(figs_dir, 'privateprop.png'),
  device = 'png',
  width = 5,
  height = 4,
  dpi = 300
)

longprivate %>%
  filter(variable == "pb_del") %>%
  summary('value')

