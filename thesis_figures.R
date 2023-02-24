library(tidyverse)
library(RColorBrewer)
library(reshape2)
options(scipen = 10)
set.seed(2)

# PBHoover
# pbhoover_vcf_dir = '/home/ebishop/thesis_projects/aim1_pbhoover/arrow_comparison_resource_use/vcfs_pbhoover/'
# arrow_vcf_dir = '/home/ebishop/thesis_projects/aim1_pbhoover/arrow_comparison_resource_use/vcfs_arrow/'
# pbhoover_mem_dir = '/home/ebishop/thesis_projects/aim1_pbhoover/arrow_comparison_resource_use/memory_use/'
# runtime_f = '/home/ebishop/thesis_projects/aim1_pbhoover/arrow_comparison_resource_use/stdout_runtimes/runtimes.txt'
# blast_dir = '/home/ebishop/thesis_projects/aim2_hybran/homology/'  # Hybran
# figs_dir = '/home/ebishop/thesis_projects/figures/'  # Output


# Emma's mac
arrow_vcf_dir = '/Users/emmabishop/workspace/thesis/aim1_phoover/arrow_comparison_resource_use/vcfs_arrow/'
pbhoover_mem_dir = '/Users/emmabishop/workspace/thesis/aim1_phoover/arrow_comparison_resource_use/mem/'
runtime_f = '/Users/emmabishop/workspace/thesis/aim1_phoover/arrow_comparison_resource_use/runtimes.txt'
blast_dir = '/Users/emmabishop/workspace/thesis/aim2_hybran/homology/'
figs_dir = '/Users/emmabishop/workspace/thesis/figures/'


##########################################################
#                      PBHOOVER                          #
##########################################################

######################
#    MEMORY USAGE    #
######################
parse_mem <- function(mem_path) {
  isolate <- read_delim(mem_path, 
                        col_names = FALSE, 
                        n_max = 1)$X14
  fname <- basename(mem_path)
  mem_raw <- read_delim(
    mem_path,
    skip = 1,
    col_names = FALSE,
    show_col_types = FALSE
  )
  peak_gb <- max(mem_raw$X2) / 1000
  df <- data.frame(peak_mem = c(peak_gb),
                   isolate = c(isolate),
                   file = c(fname))
  return(df)
}
mem_files <- list.files(
  path = pbhoover_mem_dir,
  pattern = '*.dat',
  full.names = TRUE,
  recursive = FALSE
)

all_mem <- NA
# Loop over memory profiler files
for (f in mem_files) {
  all_mem <- rbind(all_mem, parse_mem(f))
}

# All peak memories
all_mem <- drop_na(all_mem)
summary(all_mem)

# Plot and save
memfig <- ggplot(all_mem, aes(x = NA, y = peak_mem)) +
  theme_classic() +
  geom_boxplot(width = 0.4) +
  geom_jitter(aes(x = NA, y = peak_mem), width = 0.1) +
  stat_boxplot(geom = 'errorbar', width = 0.1) +
  ylab('Peak Memory Usage (GiB)') +
  ylim(0, 25) +
  ggtitle('Peak Memory Usage Measured \nfor 27 PBHoover Runs') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 13)
  )
memfig
ggsave(
  memfig,
  filename = paste0(figs_dir, 'pbh_peak_mem.png'),
  device = 'png',
  width = 4,
  height = 4,
  dpi = 300
)

##################
#    RUNTIMES    #
##################
run <-
  read_csv(runtime_f, col_names = FALSE, show_col_types = FALSE)
summary(run)

# Plot and save
runfig <- ggplot(run, aes(x = NA, y = X1)) +
  theme_classic() +
  geom_boxplot(width = 0.4) +
  geom_jitter(aes(x = NA, y = X1), width = 0.1) +
  stat_boxplot(geom = 'errorbar', width = 0.1) +
  ylab('Wall Time (minutes)') +
  ylim(0, 250) +
  ggtitle('PBHoover Runtimes for 32 Isolates') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 13)
  )
runfig
# Save
ggsave(
  runfig,
  filename = paste0(figs_dir, 'pbh_runtimes.png'),
  device = 'png',
  width = 4,
  height = 4,
  dpi = 300
)

##############
#    VCFS    #
##############
parse_vcf <- function(vcf_path, skip) {
  # Load and get isolate name
  vcf_raw <- read_tsv(vcf_path, skip = skip, show_col_types = FALSE)
  fname <- basename(vcf_path)[1]
  isolate <- str_split(fname, '\\.')[[1]][1]
  # Totals
  tot_vars <- nrow(vcf_raw)
  ins <- vcf_raw %>% filter(nchar(REF) < nchar(ALT)) %>% nrow()
  del <- vcf_raw %>% filter(nchar(REF) > nchar(ALT)) %>% nrow()
  snp <- vcf_raw %>% filter(nchar(REF) == nchar(ALT)) %>% nrow()
  # Proportions
  prop_ins <- ins / tot_vars
  prop_del <- del / tot_vars
  prop_snp <- snp / tot_vars
  # Number of hetero calls
  numhet <- length(which(vcf_raw$FILTER == "HETERO"))
  # Make into df to rbind
  df1 <- data.frame(
    isolate = c(isolate),
    tot_vars = c(tot_vars),
    snp = c(snp),
    ins = c(ins),
    del = c(del),
    prop_snp = c(prop_snp),
    prop_ins = c(prop_ins),
    prop_del = c(prop_del),
    numhet = c(numhet)
  )
  return(df1)
}

make_vcf_df <- function(in_dir, variant_caller) {
  vcf_files <- list.files(
    path = in_dir,
    pattern = '*.vcf',
    full.names = TRUE,
    recursive = FALSE
  )
  # Empty df we'll add to as we loop over files
  all_vars <- NULL
  
  if (variant_caller == 'PBHoover') {
    skip = 6
  } else {
    skip = 9
  }
  for (f in vcf_files) {
    all_vars <- rbind(all_vars, parse_vcf(f, skip))
  }
  all_vars <- drop_na(all_vars)
  all_vars$tot_vars <- as.numeric(all_vars$tot_vars)
  all_vars$variant_caller <- variant_caller
  
  return(all_vars)
}

# Get variant stats for each set of vcfs
pbh_df <- make_vcf_df(pbhoover_vcf_dir, 'PBHoover')
summary(pbh_df)

arr_df <- make_vcf_df(arrow_vcf_dir, 'Arrow')
summary(arr_df)

total_vars <- rbind(pbh_df, arr_df)

# Sanity check same number arrow and pbhoover rows
table(total_vars$variant_caller)

# Plot total vars
totvars <-
  ggplot(total_vars, aes(x = variant_caller, y = tot_vars)) +
  theme_classic() +
  geom_boxplot(width = 0.4) +
  geom_jitter(aes(x = variant_caller, y = tot_vars), width = 0.1) +
  #geom_point(position = position_jitter(seed = 10)) +
  stat_boxplot(geom = 'errorbar', width = 0.1) +
  xlab('Variant Caller') +
  ylab('Total Variants') +
  ylim(0, 70000) +
  ggtitle('Total Variants Across 32 Isolates') +
  theme(plot.title = element_text(hjust = 0.5))
totvars
# Save
ggsave(
  totvars,
  filename = paste0(figs_dir, 'pbh_total_vars.png'),
  device = 'png',
  width = 5,
  height = 4,
  dpi = 300
)


# Plot hetero calls
hetdf <- total_vars %>% filter(variant_caller == "PBHoover")
summary(hetdf)

hetcalls <- ggplot(total_vars, aes(x = NA, y = numhet)) +
  theme_classic() +
  geom_boxplot(width = 0.4) +
  geom_jitter(aes(x = NA, y = numhet), width = 0.1) +
  stat_boxplot(geom = 'errorbar', width = 0.1) +
  ylab('Total Variants') +
  ylim(0, 12000) +
  ggtitle('Numer of "HETERO" Calls \nAcross 32 PBHoover Runs') +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 13)
  )
hetcalls
# Save
ggsave(
  filename = paste0(figs_dir, 'pbh_hetero.png'),
  device = 'png',
  width = 4,
  height = 4,
  dpi = 300
)


# Plot all proportions
prop_only <- total_vars %>% select(prop_snp,
                                   prop_del,
                                   prop_ins,
                                   variant_caller)
vars_long <- melt(prop_only, id = 'variant_caller')

allprop <-
  ggplot(vars_long, aes(x = variable, y = value, color = variant_caller)) +
  theme_classic() +
  geom_boxplot(width = 0.4) +
  #geom_jitter(aes(x = variable, y = value), width = 0.1) +
  #geom_point(position = position_jitter(seed = 10))
  stat_boxplot(geom = 'errorbar', width = 0.4) +
  xlab('Variant') +
  scale_x_discrete(labels = c(
    "prop_snp" = "SNP",
    "prop_ins" = "Insertion",
    "prop_del" = "Deletion"
  )) +
  ylab('Proportion') +
  ylim(0, 1) +
  ggtitle('Proportions of Variants Called Across 32 Isolates') +
  labs(colour = "Variant Caller") +
  theme(plot.title = element_text(hjust = 0.5))
allprop
# Save
ggsave(
  allprop,
  filename = paste0(figs_dir, 'pbh_all_prop.png'),
  device = 'png',
  width = 5,
  height = 4,
  dpi = 300
)

#######################
# Print summary stats #
#######################
pbh_df <- total_vars %>% filter(variant_caller == 'PBHoover')
mean(pbh_df$tot_vars)
mean(pbh_df$prop_snp)
mean(pbh_df$prop_ins)
mean(pbh_df$prop_del)

arr_df <- total_vars %>% filter(variant_caller == 'Arrow')
mean(arr_df$tot_vars)
mean(arr_df$prop_snp)
mean(arr_df$prop_ins)
mean(arr_df$prop_del)

##########################################################
#                       HYBRAN                           #
##########################################################

###############
#    BLAST    #
###############
for_pie_chart <-
  read_csv(paste0(blast_dir, 'blastx_homology_status.csv'),
           show_col_types = FALSE)
total = sum(for_pie_chart$count)

my_colors <- brewer.pal(4, "Blues")
my_colors <- rev(my_colors)

piechart <-
  ggplot(for_pie_chart, aes(x = "", y = count, fill = result_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(
    values = my_colors,
    name = 'BLASTX Top Hit',
    labels = c(
      'Known gene (n = 364)',
      'Homologous gene (n = 92)',
      'Unknown gene (n = 2)',
      'No hit/Incomplete (n = 45)'
    ),
    breaks = c('known', 'homologous', 'unknown', 'no_hit_incomplete')
  ) +
  theme_void() +
  ggtitle('Novel ORF Homology') +
  theme(plot.title = element_text(hjust = 0.5))
piechart
# Save
ggsave(
  piechart,
  filename = paste0(figs_dir, 'hyb_homology_pie.png'),
  device = 'png',
  dpi = 300
)
