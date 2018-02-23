library(tidyverse)
library(stringr)
library(viridis)
library(cowplot)
library(ggExtra)
library(modelr)
library(lazyeval)
library(splines)
library(broom)

cbPalette7 <- c('#440154FF', '#39568CFF', '#287D8EFF', '#20A387FF', '#73D055FF',
                '#B8DE29FF', '#FDE725FF')

#Written for transient library analyses of indexed reads corresponding to:
#bc_JD02: RNA induced rep. 1
#bc_JD03: RNA induced rep. 2
#bc_JD04: RNA uninduced rep. 1
#bc_JD05: RNA uninduced rep. 2
#bc_JD11: DNA pool


#Load index files---------------------------------------------------------------

bc_R25A <- read_tsv('bc_JD02.txt')
bc_R25B <- read_tsv('bc_JD03.txt')
bc_R0A <- read_tsv('bc_JD04.txt')
bc_R0B <- read_tsv('bc_JD05.txt')
bc_DNA <- read_tsv('bc_JD11.txt')


#Load barcode mapping table, remember sequences are rcomp

barcode_map <- read_tsv('../../BCMap/uniqueSP2345.txt', 
                        col_names = c(
                          'fluff', 'barcode', 'name', 'most_common'), 
                        skip = 1) %>%
  select(-fluff) %>%
  mutate(subpool = ifelse(startsWith(name, 'subpool'), 
                          substr(name, 1, 8), 
                          'control')) 


#Join reads---------------------------------------------------------------------

#Join BC reads to barcode map per index pool, only keeping barcodes that appear 
#in the barcode map, and replacing NA with 0. Then normalizing to reads per 
#million of total reads per index

bc_map_join_bc <- function(df1, df2) {
  df2 <- df2 %>%
    mutate(normalized = as.numeric((num_reads * 1000000) / (sum(num_reads))))
  keep_bc <- left_join(df1, df2, by = 'barcode') %>%
    mutate(normalized = if_else(is.na(normalized), 
                               0, 
                               normalized)) %>%
    mutate(num_reads = if_else(is.na(num_reads), 
                                as.integer(0), 
                                num_reads))
  return(keep_bc)
}

bc_join_R25A <- bc_map_join_bc(barcode_map, bc_R25A)
bc_join_R25B <- bc_map_join_bc(barcode_map, bc_R25B)
bc_join_R0A <- bc_map_join_bc(barcode_map, bc_R0A)
bc_join_R0B <- bc_map_join_bc(barcode_map, bc_R0B)
bc_join_DNA <- bc_map_join_bc(barcode_map, bc_DNA)


#Determine variant counts by summing--------------------------------------------

#sum unique barcodes and normalized bc reads across barcodes per variant. Output 
#is barcodes and sum normalized reads per variant. Set a BC minimum of 8 per
#variant, still need to check how this affects loss of variants per subpool
#after joining all dfs

var_sum_bc_num <- function(df1) {
  bc_count <- df1 %>%
    filter(df1$num_reads > 0) %>%
    group_by(subpool, name, most_common) %>%
    summarise(barcodes = n())
  variant_sum <- df1 %>%
    group_by(subpool, name, most_common) %>%
    count(name, wt = normalized) %>%
    rename(sum = n)
  bc_sum <- inner_join(variant_sum, bc_count, 
                       by = c("name", "subpool", "most_common")) %>%
    ungroup()
  return(bc_sum)
}

variant_counts_R25A <- var_sum_bc_num(bc_join_R25A)
variant_counts_R25B <- var_sum_bc_num(bc_join_R25B)
variant_counts_R0A <- var_sum_bc_num(bc_join_R0A)
variant_counts_R0B <- var_sum_bc_num(bc_join_R0B)
variant_counts_DNA <- var_sum_bc_num(bc_join_DNA) %>%
  filter(barcodes > 7)


#Normalizing RNA reads to DNA---------------------------------------------------

#combine DNA and RNA cumm. BC counts, only keeping instances in both sets and 
#determining RNA/DNA per variant. Ratio is summed normalized reads of RNA over 
#summed normalized reads DNA

var_expression <- function(df1, df2) {
  RNA_DNA <- inner_join(df1, df2, 
                        by = c("name", "subpool", "most_common"), 
                        suffix = c("_RNA", "_DNA")) %>%
    mutate(ratio = sum_RNA / sum_DNA)
  print('x defined as RNA, y defined as DNA in var_expression(x,y)')
  return(RNA_DNA)
}

RNA_DNA_R25A <- var_expression(variant_counts_R25A, variant_counts_DNA)
RNA_DNA_R25B <- var_expression(variant_counts_R25B, variant_counts_DNA)
RNA_DNA_R0A <- var_expression(variant_counts_R0A, variant_counts_DNA)
RNA_DNA_R0B <- var_expression(variant_counts_R0B, variant_counts_DNA)

#Combine biological replicates

var_rep <- function(df0A, df0B, df25A, df25B) {
  join_0 <- inner_join(df0A, df0B, 
                       by = c("name", "subpool", "most_common", "sum_DNA", 
                              "barcodes_DNA"), 
                       suffix = c("_0A", "_0B"))
  join_25 <- inner_join(df25A, df25B, 
                        by = c("name", "subpool", "most_common", "sum_DNA", 
                               "barcodes_DNA"), 
                        suffix = c("_25A", "_25B"))
  join_0_25 <- inner_join(join_0, join_25, 
                          by = c("name", "subpool", "most_common", "sum_DNA", 
                                 "barcodes_DNA"))
  print('processed dfs in order: df0A, df0B, df25A, df25B')
  return(join_0_25)
}

rep_1_2 <- var_rep(RNA_DNA_R0A, RNA_DNA_R0B, RNA_DNA_R25A, RNA_DNA_R25B)

#Normalize to background

back_norm_1_ind_2 <- function(df1) {
  gsub_df1 <- df1 %>%
    ungroup () %>%
    filter(subpool != 'control') %>%
    mutate(
      name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
      name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name),
      name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
    ) %>%
    mutate(background = name) %>%
    mutate(background = str_sub(background, 
                                nchar(background)-5, 
                                nchar(background)))
  backgrounds <- gsub_df1 %>%
    filter(startsWith(name,
                      'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')
    ) %>%
    select(background, ratio_0A, ratio_0B, ratio_25A, ratio_25B) %>%
    rename(ratio_0A_back = ratio_0A) %>%
    rename(ratio_0B_back = ratio_0B) %>%
    rename(ratio_25A_back = ratio_25A) %>%
    rename(ratio_25B_back = ratio_25B)
  back_join_norm <- left_join(gsub_df1, backgrounds, by = 'background') %>%
    mutate(ratio_0A_norm = ratio_0A/ratio_0A_back) %>%
    mutate(ratio_0B_norm = ratio_0B/ratio_0B_back) %>%
    mutate(ratio_25A_norm = ratio_25A/ratio_25A_back) %>%
    mutate(ratio_25B_norm = ratio_25B/ratio_25B_back) %>%
    mutate(ave_ratio_0_norm = (ratio_0A_norm + ratio_0B_norm)/2) %>%
    mutate(ave_ratio_25_norm = (ratio_25A_norm + ratio_25B_norm)/2) %>%
    mutate(induction = ave_ratio_25_norm/ave_ratio_0_norm)
}

rep_1_2_back_norm <- back_norm_1_ind_2(rep_1_2)


#Log transform df then average barcodes in each biological replicate

var_log2 <- function(df) {
  log_ratio_df <- df %>% 
    mutate_if(is.double, funs(log2(.)))
  return(log_ratio_df)
}

var_log10 <- function(df) {
  log_ratio_df <- df %>% 
    mutate_if(is.double, funs(log10(.)))
  return(log_ratio_df)
}

log2_rep_1_2_back_norm <- var_log2(rep_1_2_back_norm) %>%
  mutate(ave_barcodes_0 = (barcodes_RNA_0A + barcodes_RNA_0B)/2) %>%
  mutate(ave_barcodes_25 = (barcodes_RNA_25A + barcodes_RNA_25B)/2)

log10_rep_1_2_back_norm <- var_log10(rep_1_2_back_norm) %>%
  mutate(ave_barcodes_0 = (barcodes_RNA_0A + barcodes_RNA_0B)/2) %>%
  mutate(ave_barcodes_25 = (barcodes_RNA_25A + barcodes_RNA_25B)/2)
  

#Subpool separation-------------------------------------------------------------

#Separate string qualifiers per subpool

#Subpool 2 contains either a consensus site (TGACGTCA) or consensus surrounded 
#by 2 bp flanks on either side (ATTGACGTCAGC) as is used in the rest of the 
#subpools. Each site is placed on the background starting closest to minP and 
#are then moved along the backgrounds at 1 bp increments. Separation lists the 
#type of site, the distance (start of the consensus site) and the background 
#used. Added 2 bp to consensusflank so that the start of the binding site is 
#represented instead of the start of the flank

subpool2 <- 
  filter(rep_1_2_back_norm, subpool == "subpool2") %>%
  ungroup() %>%
  select(-subpool) %>%
  separate(name, into = c("fluff1", "site", "fluff2", "dist", "fluff3"),
           sep = "_", convert = TRUE) %>% 
  select(-fluff1, -fluff2, -fluff3) %>%
  mutate(dist = ifelse(startsWith(site, 'consensusflank'), dist + 2, dist)) %>%
  group_by(background, site) %>%
  mutate(med_ratio_25A_norm = median(ratio_25A_norm)) %>%
  mutate(med_ratio_25B_norm = median(ratio_25B_norm)) %>%
  ungroup() %>%
  mutate(ratio_25A_norm_med = ratio_25A_norm/med_ratio_25A_norm) %>%
  mutate(ratio_25B_norm_med = ratio_25B_norm/med_ratio_25B_norm) %>%
  mutate(ave_ratio_25_norm_med = (ratio_25A_norm_med + ratio_25B_norm_med)/2)


#Subpool 3 contains 2 consensus binding sites with flanks (ATTGACGTCAGC) that 
#vary in distance from one another by 0 (no inner flanks), 5, 10, 15, 20 and 70 
#bp (all but 0 appear as -4 bp spacing). Each site distance combination is then 
#moved along the backgrounds at 1 bp increments starting from closest to the 
#minP. Separation lists the spacing between sites, distance (start of consensus)
#and the background. Added 2 to all distances to measure the start of BS and not 
#to flank. Added 4 to all spacing but 0 to include flank spaces.

subpool3 <- 
  filter(rep_1_2_back_norm, subpool == "subpool3") %>%
  ungroup () %>%
  select(-subpool) %>%
  mutate(name = gsub('2BS ', '', name), 
         name = gsub(' bp spacing ', '_', name)) %>%
  separate(name, 
           into = c("subpool", "spacing", "fluff2", "fluff3", "dist", "fluff4"),
           sep = "_", convert = TRUE) %>%
  select(-subpool, -fluff2, -fluff3, -fluff4) %>%
  mutate(dist = as.integer(dist + 2)) %>%
  mutate(spacing = ifelse(spacing != as.integer(0), 
                          as.integer(spacing + 4), 
                          as.integer(spacing))) 

#Subpool 4 contains 2 binding sites with flanks that vary in site type from 
#consensus (ATTGACGTCAGC) moderate (ATTGACGTCTGC) weak (ATTGAAGTCAGC) and 
#half(NNNGCCGTCATA). These sites vary in spacing from 0, 5, 10, 15, and 20 bp. 
#Both sites are moved along the promoter at 5 bp increments starting closest to 
#minP. The order of sites in the name is from further to closer to the minP and 
#distance is from the end of the oligo to the site closer to the proximal 
#promoter. 

subpool4 <- 
  filter(rep_1_2_back_norm, subpool == "subpool4") %>%
  ungroup () %>%
  select(-subpool) %>%
  mutate(name = gsub('consensus0 ', 'consensus_0', name),
         name = gsub('weak0 ', 'weak_0', name),
         name = gsub('moderate0 ', 'moderate_0', name),
         name = gsub('half0 ', 'half_0', name),
         name = gsub('bp spacing', '', name)) %>%
  separate(name, 
           into = c("subpool", "site1", "site2", "spacing", "fluff1", "dist", 
                    "fluff2"), 
           sep = "_", convert = TRUE) %>%
  select(-subpool, -fluff1, -fluff2) %>%
  mutate(consensus = str_detect(site1, "consensus") + 
           str_detect(site2, "consensus")) %>%
  mutate(moderate = str_detect(site1, "moderate") + 
           str_detect(site2, "moderate")) %>%
  mutate(weak = str_detect(site1, "weak") +
           str_detect(site2, "weak")) %>%
  mutate(half = str_detect(site1, "half") +
           str_detect(site2, "half")) %>%
  mutate(dist = dist + 2) %>%
  mutate(spacing = ifelse(spacing != as.integer(0), spacing + 4, spacing))

#Subpool 5 contains 6 equally spaced sites spaced 13 bp apart and starting from 
#furthest to the minP. These sites are filled with sites of either the consensus
#site, a weak site or no site. Both the weak and consensus sites have the same 
#flanking sequence.

subpool5 <- 
  filter(rep_1_2_back_norm, subpool == "subpool5") %>%
  ungroup () %>%
  select(-subpool) %>%
  mutate(name = gsub('no_site', 'nosite', name)) %>%
  separate(name, into = c(
    "subpool", "site1", "site2", "site3", "site4", "site5", "site6", "fluff"), 
    sep = "_") %>%
  select(-subpool, -fluff) %>%
  mutate(consensus = str_detect(site1, "consensus") + 
           str_detect(site2, "consensus") + 
           str_detect(site3, "consensus") + 
           str_detect(site4, "consensus") + 
           str_detect(site5, "consensus") + 
           str_detect(site6, "consensus")) %>%
  mutate(weak = str_detect(site1, "weak") +
           str_detect(site2, "weak") +
           str_detect(site3, "weak") +
           str_detect(site4, "weak") +
           str_detect(site5, "weak") +
           str_detect(site6, "weak")) %>%
  mutate(nosite = str_detect(site1, "nosite") +
           str_detect(site2, "nosite") +
           str_detect(site3, "nosite") +
           str_detect(site4, "nosite") +
           str_detect(site5, "nosite") +
           str_detect(site6, "nosite")) %>%
  mutate(total_sites = consensus + weak) %>%
  mutate(site_combo = ifelse(weak == 0 & consensus > 0, 
                             'consensus', 
                             'mixed')) %>%
  mutate(site_type = ifelse(consensus == 0 & weak > 0, 
                  'weak', 
                  site_combo))

#Controls

controls <- 
  filter(rep_1_2, subpool == "control") %>%
  ungroup() %>%
  mutate(ave_ratio_0 = (ratio_0A + ratio_0B)/2) %>%
  mutate(ave_ratio_25 = (ratio_25A + ratio_25B)/2) %>%
  mutate(induction = ave_ratio_25/ave_ratio_0) %>%
  mutate(ave_barcodes_0 = (barcodes_RNA_0A + barcodes_RNA_0B)/2) %>%
  mutate(ave_barcodes_25 = (barcodes_RNA_25A + barcodes_RNA_25B)/2)
  

#Overall library analysis plots-------------------------------------------------

#replicate plots using ratio

log10_rep_1_2 <- var_log10(rep_1_2)

#No background

p_rep_ratio_0 <- ggplot(data = NULL, aes(ratio_0A, ratio_0B)) +
  facet_grid(. ~ subpool) +
  geom_point(data = filter(log10_rep_1_2, subpool == 'subpool5'),
             color = '#440154FF', alpha = 0.3) +
  geom_point(data = filter(log10_rep_1_2, subpool == 'subpool3'),
             color = '#33638DFF', alpha = 0.3) +
  geom_point(data = filter(log10_rep_1_2, subpool == 'subpool4'),
             color = '#29AF7FFF', alpha = 0.3) +
  geom_point(data = filter(log10_rep_1_2, subpool == 'subpool2'),
             color = '#DCE319FF', alpha = 0.3) +
  geom_density2d(data = filter(log10_rep_1_2, subpool != 'control'), 
                 color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log10_rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("log10 sum\nRNA/DNA BR 1") +
  ylab("log10 sum\nRNA/DNA BR 2") +
  scale_x_continuous(breaks = c(-1:2), limits = c(-1, 2)) + 
  scale_y_continuous(breaks = c(-1:2), limits = c(-1, 2)) +
  background_grid(major = 'xy', minor = 'none') + 
  panel_border() +
  annotate("text", x = 0, y = 1.8, color = '#440154FF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2, 
                                          subpool == 'subpool5')$ratio_0A,
                                   filter(log10_rep_1_2, 
                                          subpool == 'subpool5')$ratio_0B,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = 0, y = 1.5, color = '#33638DFF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2, 
                                          subpool == 'subpool3')$ratio_0A,
                                   filter(log10_rep_1_2, 
                                          subpool == 'subpool3')$ratio_0B,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = 0, y = 1.2, color = '#29AF7FFF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2, 
                                          subpool == 'subpool4')$ratio_0A,
                                   filter(log10_rep_1_2, 
                                          subpool == 'subpool4')$ratio_0B,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = 0, y = 0.9, color = '#DCE319FF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2, 
                                          subpool == 'subpool2')$ratio_0A,
                                   filter(log10_rep_1_2, 
                                          subpool == 'subpool2')$ratio_0B,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2)))

p_rep_ratio_25 <- ggplot(data = NULL, aes(ratio_25A, ratio_25B)) +
  facet_grid(. ~ subpool) +
  geom_point(data = filter(log10_rep_1_2, subpool == 'subpool5'),
             color = '#440154FF', alpha = 0.3) +
  geom_point(data = filter(log10_rep_1_2, subpool == 'subpool3'),
             color = '#33638DFF', alpha = 0.3) +
  geom_point(data = filter(log10_rep_1_2, subpool == 'subpool4'),
             color = '#29AF7FFF', alpha = 0.3) +
  geom_point(data = filter(log10_rep_1_2, subpool == 'subpool2'),
             color = '#DCE319FF', alpha = 0.3) +
  geom_density2d(data = filter(log10_rep_1_2, subpool != 'control'), 
                 color = 'black', size = 0.2, bins = 10) +
  geom_point(data = filter(log10_rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = '#B8DE29FF', alpha = 0.7) +
  annotation_logticks(scaled = TRUE) +
  xlab("log10 sum\nRNA/DNA BR 1") +
  ylab("log10 sum\nRNA/DNA BR 2") +
  scale_x_continuous(breaks = c(-1:2), limits = c(-1, 2)) + 
  scale_y_continuous(breaks = c(-1:2), limits = c(-1, 2)) +
  background_grid(major = 'xy', minor = 'none') + 
  panel_border() +
  annotate("text", x = 0, y = 1.8, color = '#440154FF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2, 
                                          subpool == 'subpool5')$ratio_25A,
                                   filter(log10_rep_1_2, 
                                          subpool == 'subpool5')$ratio_25B,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = 0, y = 1.5, color = '#33638DFF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2, 
                                          subpool == 'subpool3')$ratio_25A,
                                   filter(log10_rep_1_2, 
                                          subpool == 'subpool3')$ratio_25B,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = 0, y = 1.2, color = '#29AF7FFF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2, 
                                          subpool == 'subpool4')$ratio_25A,
                                   filter(log10_rep_1_2, 
                                          subpool == 'subpool4')$ratio_25B,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = 0, y = 0.9, color = '#DCE319FF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2, 
                                          subpool == 'subpool2')$ratio_25A,
                                   filter(log10_rep_1_2, 
                                          subpool == 'subpool2')$ratio_25B,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2)))

p_rep_ratio_0_25 <- plot_grid(p_rep_ratio_0, p_rep_ratio_25, 
                              nrow = 2, labels = c(' 0 µM', '25 µM'), 
                              align = 'v', hjust = -3, vjust = -1, scale = 0.9)

save_plot('plots/p_rep_ratio_0_25.png', p_rep_ratio_0_25,
          base_width = 10.5, base_height = 7)

#Background-normalized

p_rep_backnorm_ratio_0 <- ggplot(data = NULL, 
                                 aes(ratio_0A_norm, ratio_0B_norm)) +
  facet_grid(. ~ subpool) +
  geom_point(data = filter(log10_rep_1_2_back_norm, subpool == 'subpool5'),
             color = '#440154FF', alpha = 0.3) +
  geom_point(data = filter(log10_rep_1_2_back_norm, subpool == 'subpool3'),
             color = '#33638DFF', alpha = 0.3) +
  geom_point(data = filter(log10_rep_1_2_back_norm, subpool == 'subpool4'),
             color = '#29AF7FFF', alpha = 0.3) +
  geom_point(data = filter(log10_rep_1_2_back_norm, subpool == 'subpool2'),
             color = '#DCE319FF', alpha = 0.3) +
  geom_density2d(data = log10_rep_1_2_back_norm, 
                 color = 'black', size = 0.2, bins = 10) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  annotation_logticks(scaled = TRUE) +
  xlab("log10 background-norm.\nsum RNA/DNA BR 1") +
  ylab("log10 background-norm.\nsum RNA/DNA BR 2") +
  scale_x_continuous(breaks = c(0:2), limits = c(-0.5, 2)) + 
  scale_y_continuous(breaks = c(0:2), limits = c(-0.5, 2)) +
  background_grid(major = 'xy', minor = 'none') + 
  panel_border() +
  annotate("text", x = 0, y = 1.8, color = '#440154FF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool5')$ratio_0A_norm,
                                   filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool5')$ratio_0B_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = 0, y = 1.5, color = '#33638DFF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool3')$ratio_0A_norm,
                                   filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool3')$ratio_0B_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = 0, y = 1.2, color = '#29AF7FFF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool4')$ratio_0A_norm,
                                   filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool4')$ratio_0B_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = 0, y = 0.9, color = '#DCE319FF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool2')$ratio_0A_norm,
                                   filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool2')$ratio_0B_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2)))

p_rep_backnorm_ratio_25 <- ggplot(data = NULL, 
                                  aes(ratio_25A_norm, ratio_25B_norm)) +
  facet_grid(. ~ subpool) +
  geom_point(data = filter(log10_rep_1_2_back_norm, subpool == 'subpool5'),
             color = '#440154FF', alpha = 0.3) +
  geom_point(data = filter(log10_rep_1_2_back_norm, subpool == 'subpool3'),
             color = '#33638DFF', alpha = 0.3) +
  geom_point(data = filter(log10_rep_1_2_back_norm, subpool == 'subpool4'),
             color = '#29AF7FFF', alpha = 0.3) +
  geom_point(data = filter(log10_rep_1_2_back_norm, subpool == 'subpool2'),
             color = '#DCE319FF', alpha = 0.3) +
  geom_density2d(data = log10_rep_1_2_back_norm, 
                 color = 'black', size = 0.2, bins = 10) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  annotation_logticks(scaled = TRUE) +
  xlab("log10 background-norm.\nsum RNA/DNA BR 1") +
  ylab("log10 background-norm.\nsum RNA/DNA BR 2") +
  scale_x_continuous(breaks = c(0:2), limits = c(-0.5, 2)) + 
  scale_y_continuous(breaks = c(0:2), limits = c(-0.5, 2)) +
  background_grid(major = 'xy', minor = 'none') + 
  panel_border() +
  annotate("text", x = 0, y = 1.8, color = '#440154FF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool5')$ratio_25A_norm,
                                   filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool5')$ratio_25B_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = 0, y = 1.5, color = '#33638DFF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool3')$ratio_25A_norm,
                                   filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool3')$ratio_25B_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = 0, y = 1.2, color = '#29AF7FFF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool4')$ratio_25A_norm,
                                   filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool4')$ratio_25B_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = 0, y = 0.9, color = '#DCE319FF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool2')$ratio_25A_norm,
                                   filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool2')$ratio_25B_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2)))

p_rep_backnorm_ratio_0_25 <- plot_grid(p_rep_backnorm_ratio_0, 
                                       p_rep_backnorm_ratio_25, 
                                       nrow = 2, labels = c(' 0 µM', '25 µM'), 
                                       align = 'v', hjust = -3, vjust = -1, 
                                       scale = 0.9)

save_plot('plots/p_rep_backnorm_ratio_0_25.png', p_rep_backnorm_ratio_0_25,
          base_width = 10.5, base_height = 7)

#Replicate plots for sp5 for CMB presentation

p_rep_backnorm_ratio_0_sp5 <- ggplot(data = filter(log10_rep_1_2_back_norm, 
                                                 subpool == 'subpool5'),
                                   aes(ratio_0A_norm, ratio_0B_norm)) +
  geom_point(color = '#440154FF', alpha = 0.3) +
  geom_density2d(data = log10_rep_1_2_back_norm, 
                 color = 'black', size = 0.2, bins = 10) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  annotation_logticks(scaled = TRUE) +
  xlab("log10 background-norm.\nsum RNA/DNA BR 1") +
  ylab("log10 background-norm.\nsum RNA/DNA BR 2") +
  scale_x_continuous(breaks = c(0:2), limits = c(-0.5, 2)) + 
  scale_y_continuous(breaks = c(0:2), limits = c(-0.5, 2)) +
  background_grid(major = 'xy', minor = 'none') + 
  annotate("text", x = 0.5, y = 1.5, color = '#440154FF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool5')$ratio_0A_norm,
                                   filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool5')$ratio_0B_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) 

p_rep_backnorm_ratio_25_sp5 <- ggplot(data = filter(log10_rep_1_2_back_norm, 
                                                   subpool == 'subpool5'),
                                     aes(ratio_25A_norm, ratio_25B_norm)) +
  geom_point(color = '#440154FF', alpha = 0.3) +
  geom_density2d(data = log10_rep_1_2_back_norm, 
                 color = 'black', size = 0.2, bins = 10) +
  geom_hline(yintercept = 0, alpha = 0.5) +
  geom_vline(xintercept = 0, alpha = 0.5) +
  annotation_logticks(scaled = TRUE) +
  xlab("log10 background-norm.\nsum RNA/DNA BR 1") +
  ylab("log10 background-norm.\nsum RNA/DNA BR 2") +
  scale_x_continuous(breaks = c(0:2), limits = c(-0.5, 2)) + 
  scale_y_continuous(breaks = c(0:2), limits = c(-0.5, 2)) +
  background_grid(major = 'xy', minor = 'none') + 
  annotate("text", x = 0.5, y = 1.5, color = '#440154FF', 
           label = paste('r =', 
                         round(cor(filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool5')$ratio_25A_norm,
                                   filter(log10_rep_1_2_back_norm, 
                                          subpool == 'subpool5')$ratio_25B_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2)))

p_rep_backnorm_ratio_0_25_sp5 <- plot_grid(p_rep_backnorm_ratio_0_sp5, 
                                       p_rep_backnorm_ratio_25_sp5, 
                                       nrow = 2, scale = 0.9)

save_plot('plots/p_rep_backnorm_ratio_0_25_sp5.png', 
          p_rep_backnorm_ratio_0_25_sp5, base_width = 3.5, base_height = 7)

#Showing what library looks like before and after normalizing to background

induction_noback <- rep_1_2 %>%
  mutate(ave_sum_RNA_0 = (sum_RNA_0A + sum_RNA_0B)/2) %>%
  mutate(ave_sum_RNA_25 = (sum_RNA_25A + sum_RNA_25B)/2) %>%
  mutate(induction = ave_sum_RNA_25/ave_sum_RNA_0) %>%
  var_log2() %>%
  filter(subpool != 'control')

p_induction_noback <- ggplot(induction_noback, aes(induction)) +
  facet_grid(. ~ subpool) +
  geom_density(kernel = 'gaussian') +
  geom_vline(xintercept = -1.172789, alpha = 0.5) +
  annotation_logticks(scaled = TRUE, sides = 'b') +
  xlab("Log2 induction of sum RNA reads") +
  panel_border() +
  scale_x_continuous(limits = c(-2, 4), breaks = c(-2:4))

p_induction_back <- ggplot(log2_rep_1_2_back_norm, aes(induction)) +
  facet_grid(. ~ subpool) +
  geom_density(kernel = 'gaussian') +
  geom_vline(xintercept = 0, alpha = 0.5) +
  annotation_logticks(scaled = TRUE, sides = 'b') +
  xlab("Log2 induction of background-norm. sum RNA reads") +
  panel_border() +
  scale_x_continuous(limits = c(-2, 4), breaks = c(-2:4))

p_induction_back_compare <- plot_grid(p_induction_noback, p_induction_back,
                                      nrow = 2)

save_plot('plots/p_induction_back_compare.png', p_induction_back_compare,
          base_width = 9, base_height = 5)


#Comparison to integrated-------------------------------------------------------

int_back_norm_rep_1_2 <- read_tsv('../20171129_intLib/int_back_norm_rep_1_2.txt')

int_trans <- left_join(int_back_norm_rep_1_2, trans_back_norm_rep_0_22, 
                       by = c('subpool', 'name', 'most_common', 'background'))

#Analyzing 1-site expression patterns in subpool5

subpool5_int <- 
  filter(int_back_norm_rep_1_2, subpool == "subpool5") %>%
  ungroup () %>%
  select(-subpool) %>%
  mutate(name = gsub('no_site', 'nosite', name)) %>%
  separate(name, into = c(
    "subpool", "site1", "site2", "site3", "site4", "site5", "site6", "fluff"), 
    sep = "_") %>%
  select(-subpool, -fluff) %>%
  mutate(consensus = str_detect(site1, "consensus") + 
           str_detect(site2, "consensus") + 
           str_detect(site3, "consensus") + 
           str_detect(site4, "consensus") + 
           str_detect(site5, "consensus") + 
           str_detect(site6, "consensus")) %>%
  mutate(weak = str_detect(site1, "weak") +
           str_detect(site2, "weak") +
           str_detect(site3, "weak") +
           str_detect(site4, "weak") +
           str_detect(site5, "weak") +
           str_detect(site6, "weak")) %>%
  mutate(nosite = str_detect(site1, "nosite") +
           str_detect(site2, "nosite") +
           str_detect(site3, "nosite") +
           str_detect(site4, "nosite") +
           str_detect(site5, "nosite") +
           str_detect(site6, "nosite")) %>%
  mutate(total_sites = consensus + weak) %>%
  mutate(site_combo = ifelse(weak == 0 & consensus > 0, 
                             'consensus', 
                             'mixed')) %>%
  mutate(site_type = ifelse(consensus == 0 & weak > 0, 
                            'weak', 
                            site_combo))

s5_single_site_exp_int <- subpool5_int %>%
  filter(weak == 0 & consensus == 1) %>%
  mutate(site1 = str_detect(site1, "consensus") * 1) %>%
  mutate(site2 = str_detect(site2, "consensus") * 2) %>%
  mutate(site3 = str_detect(site3, "consensus") * 3) %>%
  mutate(site4 = str_detect(site4, "consensus") * 4) %>%
  mutate(site5 = str_detect(site5, "consensus") * 5) %>%
  mutate(site6 = str_detect(site6, "consensus") * 6) %>%
  mutate(site = site1 + site2 + site3 + site4 + site5 + site6)

p_s5_single_site_exp_int <- ggplot(s5_single_site_exp_int, 
                               aes(as.factor(site), ave_med_ratio_norm)) +
  geom_bar(stat = 'identity', fill = '#440154FF') +
  facet_grid(. ~ background) +
  xlab('Site position') + 
  panel_border() +
  ylab('Average background-normalized\nmedian RNA/DNA at 25 µM') +
  geom_hline(yintercept = 1, color = 'gray')

save_plot('plots/p_s5_single_site_exp.png', p_s5_single_site_exp,
          scale = 1.2, base_width = 4.5, base_height = 3)


#Barcode analysis---------------------------------------------------------------

#Plot reads per BC
sample_DNA <- filter(bc_join_DNA, num_reads > 0)
sample_R0A <- filter(bc_join_R0A, num_reads > 0)
sample_R0B <- filter(bc_join_R0B, num_reads > 0)
sample_R25A <- filter(bc_join_R25A, num_reads > 0)
sample_R25B <- filter(bc_join_R25B, num_reads > 0)

p_BC_num_reads_viol_full <- ggplot(NULL, 
                                   aes(x = "", 
                                       y = num_reads, 
                                       color = subpool
                                       )
                                   ) +
  geom_violin(data = sample_DNA, aes(x = "DNA"), 
              position = position_dodge(0.75)) +
  geom_violin(data = sample_R0A, aes(x = "R0A"), 
              position = position_dodge(0.75)) + 
  geom_violin(data = sample_R0B, aes(x = "R0B"), 
              position = position_dodge(0.75)) +
  geom_violin(data = sample_R25A, aes(x = "R25A"), 
              position = position_dodge(0.75)) + 
  geom_violin(data = sample_R25B, aes(x = "R25B"), 
              position = position_dodge(0.75)) + 
  scale_color_manual(values = cbPalette7) +
  xlab("") +
  ylab("Reads per barcode")

save_plot('plots/BC_num_reads_viol_full.png', p_BC_num_reads_viol_full)

p_BC_num_reads_box_zoom <- ggplot(NULL, 
                                  aes(x = "", 
                                      y = num_reads, 
                                      color = subpool
                                      )
                                  ) +
  geom_boxplot(data = sample_DNA, aes(x = "DNA"), 
               position = position_dodge(0.9)) +
  geom_boxplot(data = sample_R0A, aes(x = "R0A"), 
               position = position_dodge(0.9)) + 
  geom_boxplot(data = sample_R0B, aes(x = "R0B"), 
               position = position_dodge(0.9)) +
  geom_boxplot(data = sample_R25A, aes(x = "R25A"), 
               position = position_dodge(0.9)) + 
  geom_boxplot(data = sample_R25B, aes(x = "R25B"), 
               position = position_dodge(0.9)) + 
  scale_color_manual(values = cbPalette7) +
  xlab("") +
  ylab("Reads per barcode") + 
  ylim(0, 25)

save_plot('plots/BC_num_reads_box_zoom.png', p_BC_num_reads_box_zoom)

p_BC_per_variant <- ggplot(rep_1_2, 
                           aes(x = "", 
                               y = num_reads, 
                               color = subpool
                               )
                           ) +
  geom_boxplot(aes(x = "DNA", y = barcodes_DNA), 
               position = position_dodge(0.9)) +
  geom_boxplot(aes(x = "R0A", y = barcodes_RNA_0A), 
               position = position_dodge(0.9)) + 
  geom_boxplot(aes(x = "R0B", y = barcodes_RNA_0B), 
               position = position_dodge(0.9)) + 
  geom_boxplot(aes(x = "R25A", y = barcodes_RNA_25A), 
               position = position_dodge(0.9)) + 
  geom_boxplot(aes(x = "R25B", y = barcodes_RNA_25B), 
               position = position_dodge(0.9)) + 
  scale_color_manual(values = cbPalette7) +
  xlab("") +
  ylab("Barcodes per variant") +
  background_grid(major = c('y'), minor = c('y')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot('plots/BC_per_variant.png', p_BC_per_variant)


#Control analysis---------------------------------------------------------------

p_contr_25_induction <- ggplot(filter(controls, ave_ratio_25 < 20), 
                                aes(ave_ratio_0, ave_ratio_25)) +
  geom_point(alpha = 0.7) + xlim(0,2) + ylim(0,5) +
  ylab(expression(
    paste("Average Expression at 25 ", mu, "M Forskolin")
  )) +
  xlab(expression(
    paste("Average Expression at 0 ", mu, "M Forskolin")))

save_plot('plots/contr_25_induction.png',
          p_contr_25_induction)


#Comparing to integrated luciferase assays
luc_comp_seq <- filter(ave_variant, 
                       name %in% c('pGL4.29 Promega 1-63 + 1-87', 
                                   'subpool3_2BS 16 bp spacing consensus+flank x2_dist_79_scramble pGL4.29 Promega 1-63 + 1-87',
                                   'subpool3_2BS 6 bp spacing consensus+flank x2_dist_9_scramble pGL4.29 Promega 1-63 + 1-87',
                                   'subpool3_2BS 11 bp spacing consensus+flank x2_dist_14_scramble pGL4.29 Promega 1-63 + 1-87',
                                   'subpool5_no_site_consensus_weak_consensus_no_site_weak_Smith R. Vista chr9:83712599-83712766')) %>%
  select(subpool, name, ave_ratio_25, ave_ratio_0, ave_induction)

luc_comp_luc <- tibble('name' = c('pGL4.29 Promega 1-63 + 1-87',
                                  'subpool3_2BS 11 bp spacing consensus+flank x2_dist_14_scramble pGL4.29 Promega 1-63 + 1-87',
                                  'subpool3_2BS 16 bp spacing consensus+flank x2_dist_79_scramble pGL4.29 Promega 1-63 + 1-87',
                                  'subpool3_2BS 6 bp spacing consensus+flank x2_dist_9_scramble pGL4.29 Promega 1-63 + 1-87',
                                  'subpool5_no_site_consensus_weak_consensus_no_site_weak_Smith R. Vista chr9:83712599-83712766'),
                       'luc_0' = c(481, 1923, 172, 1114, 452),
                       'luc_25' = c(33383, 7398, 669, 4645, 2559))

luc_comp_both <- left_join(luc_comp_seq, luc_comp_luc, by = 'name') %>%
  mutate(luc_induction = luc_25/luc_0)

p_luc_comp_both_0 <- ggplot(luc_comp_both, 
                                aes(luc_0, ave_ratio_0, 
                                    color = subpool)) +
  geom_point() +
  labs(x = 'Integrated Variant Luminescence (AU)', 
       y = 'Transient Variant Expression (Seq)')

p_luc_comp_both_25 <- ggplot(luc_comp_both, 
                              aes(luc_25, ave_ratio_25, 
                                  color = subpool)) +
  geom_point() +
  labs(x = 'Integrated Variant Luminescence (AU)', 
       y = 'Transient Variant Expression (Seq)')

p_luc_comp_both_induction <- ggplot(luc_comp_both, 
                                    aes(luc_induction, ave_induction, 
                                        color = subpool)) +
  geom_point() +
  labs(x = 'Integrated Variant Induction Ratio', 
       y = 'Transient Variant Induction Ratio')

#Positive control
p_control_3_5_0_25 <- ggplot(NULL, 
                                  aes(ave_ratio_0, ave_ratio_25)) +
  geom_point(data = sep_5, color = '#999999', 
             alpha = 0.3) +
  geom_point(data = sep_3, color = '#56B4E9', 
             alpha = 0.3) +
  geom_point(data = filter(
    controls, name == 'pGL4.29 Promega 1-63 + 1-87'
  ), color = 'black', size = 2) +
  labs(x = 'Average Expression at 0 µM', 
       y = 'Average Expression at 25 µM')

save_plot('plots/m_control_3_5_0_25.png',
          p_control_3_5_0_25, scale = 0.3, 
          base_width = 18.5, base_height = 16)


#Subpool plots from summing---------------------------------------------

#Subpool 2

p_subpool2_dist_0_25 <- ggplot(filter(subpool2, site == 'consensusflank')) + 
  facet_grid(~ background) + 
  geom_point(aes(dist, ratio_25A_norm_med), 
             alpha = 0.5, color = '#287D8EFF') +
  geom_point(aes(dist, ratio_25B_norm_med), 
             alpha = 0.5, color = '#73D055FF') +
  geom_hline(yintercept = 1) +
  geom_smooth(aes(dist, ave_ratio_25_norm_med), span = 0.1, size = 0.4,
              se = FALSE, color = '#440154FF') +
  scale_x_continuous("Distance from Proximal Promoter End (bp)", 
                     breaks = seq(from = 0, to = 150, by = 10)) +
  panel_border() + ylab('Expression at\n25 µM norm. to median') +
  background_grid(major = 'xy', minor = 'none')

save_plot('plots/subpool2_dist_0_25.png',
          p_subpool2_dist_0_25, base_width = 46, base_height = 10,
          scale = 0.3)

#Subpool 3

p_subpool3_induction <- ggplot(sep_3) + 
  geom_point(aes(dist, induction_1), alpha = 0.3, size = 1.5) +
  geom_point(aes(dist, induction_2), alpha = 0.3, size = 1.5) +
  facet_grid(spacing ~ background) + 
  geom_smooth(aes(dist, ave_induction), span = 0.1, size = 0.7, 
              se = FALSE) +
  ylab('Average Induction Ratio') + 
  panel_border() +
  background_grid(major = 'xy', minor = 'none') +
  scale_x_continuous(
    "Distance from First Site to Proximal Promoter End (bp)", 
                     breaks = seq(from = 0, to = 150, by = 10))

save_plot('plots/subpool3_induction.png',
          p_subpool3_induction, base_width = 46, base_height = 17,
          scale = 0.35)


p_subpool3_chr9_0_25 <- ggplot(filter(sep_3, 
                                           background == 'vista chr9')) + 
  geom_point(aes(dist, ratio_0_1), alpha = 0.5, size = 1.2,
             color = '#999999') +
  geom_point(aes(dist, ratio_0_2), alpha = 0.5, size = 1.2,
             color = '#999999') +
  geom_point(aes(dist, ratio_25_1), alpha = 0.5, size = 1.2,
             color = '#56B4E9') +
  geom_point(aes(dist, ratio_25_2), alpha = 0.5, size = 1.2,
             color = '#56B4E9') +
  facet_grid(spacing ~ .) + 
  geom_smooth(aes(dist, ave_ratio_0), span = 0.1, size = 0.4,
              se = FALSE, color = '#999999') +
  geom_smooth(aes(dist, ave_ratio_25), span = 0.1, size = 0.4,
              se = FALSE, color = '#56B4E9') +
  panel_border() +
  scale_x_continuous("Distance from First Site to Proximal Promoter End (bp)", 
                     breaks = seq(from = 0, to = 150, by = 10),
                     limits = c(0, 135))

save_plot('plots/subpool3_chr9_0_25.png',
          p_subpool3_chr9_0_25, base_width = 32, base_height = 20,
          scale = 0.3)

p_subpool3_chr9_5_0_25 <- ggplot(filter(sep_3, background == 'vista chr9')) + 
  geom_vline(xintercept = 4, color = 'black') +
  geom_vline(xintercept = 14, color = 'black') +
  geom_vline(xintercept = 24, color = 'black') +
  geom_point(aes(dist, ratio_0_1), alpha = 0.5, size = 1.2,
             color = '#999999') +
  geom_point(aes(dist, ratio_0_2), alpha = 0.5, size = 1.2,
             color = '#999999') +
  geom_point(aes(dist, ratio_25_1), alpha = 0.5, size = 1.2,
             color = '#56B4E9') +
  geom_point(aes(dist, ratio_25_2), alpha = 0.5, size = 1.2,
             color = '#56B4E9') +
  facet_grid(spacing ~ .) + 
  geom_smooth(aes(dist, ave_ratio_0), span = 0.1, size = 0.4,
              se = FALSE, color = '#999999') +
  geom_smooth(aes(dist, ave_ratio_25), span = 0.1, size = 0.4,
              se = FALSE, color = '#56B4E9') +
  panel_border() +
  scale_y_continuous("Expression") +
  scale_x_continuous("Distance from First Site to Proximal Promoter End (bp)", 
                     breaks = seq(from = 0, to = 150, by = 10),
                     limits = c(0, 135))

save_plot('plots/subpool3_chr9_5_0_25.png',
          p_subpool3_chr9_5_0_25, base_width = 32, base_height = 20,
          scale = 0.3)

p_subpool3_chr9_10_0_25 <- ggplot(filter(sep_3, background == 'vista chr9')) + 
  geom_vline(xintercept = 7, color = 'black') +
  geom_vline(xintercept = 17, color = 'black') +
  geom_vline(xintercept = 27, color = 'black') +
  geom_point(aes(dist, ratio_0_1), alpha = 0.5, size = 1.2,
             color = '#999999') +
  geom_point(aes(dist, ratio_0_2), alpha = 0.5, size = 1.2,
             color = '#999999') +
  geom_point(aes(dist, ratio_25_1), alpha = 0.5, size = 1.2,
             color = '#56B4E9') +
  geom_point(aes(dist, ratio_25_2), alpha = 0.5, size = 1.2,
             color = '#56B4E9') +
  facet_grid(spacing ~ .) + 
  geom_smooth(aes(dist, ave_ratio_0), span = 0.1, size = 0.4,
              se = FALSE, color = '#999999') +
  geom_smooth(aes(dist, ave_ratio_25), span = 0.1, size = 0.4,
              se = FALSE, color = '#56B4E9') +
  panel_border() +
  scale_y_continuous("Expression") +
  scale_x_continuous("Distance from First Site to Proximal Promoter End (bp)", 
                     breaks = seq(from = 0, to = 150, by = 10),
                     limits = c(0, 135))

save_plot('plots/subpool3_chr9_10_0_25.png',
          p_subpool3_chr9_10_0_25, base_width = 32, base_height = 20,
          scale = 0.3)


#Subpool 5

#plot combinations of consensus and weak vs. their induced expression as a 
#function of the total number of sites filled

p_num_sites_num_weak_c <- ggplot(filter(subpool5, weak == 0),
                               aes(x = as.factor(total_sites), 
                                   y = ave_ratio_25_norm)) +
  geom_boxplot(aes(color = as.factor(weak)), show.legend = FALSE,
               position = position_dodge(1)) +
  facet_grid(background ~ .) + 
  panel_border() +
  scale_color_manual(values = cbPalette7) +
  annotation_logticks(sides = 'lr') +
  background_grid(major = 'y', minor = 'none') + 
  geom_vline(xintercept = c(1.5:6.5), alpha = 0.5) +
  ylab('log10 background-normalized\nexpression at 25 µM') +
  xlab('Number of binding sites')

save_plot('plots/p_num_sites_num_weak_c.png', p_num_sites_num_weak_c, scale = 1.3,
          base_width = 6, base_height = 4)

p_num_sites_num_weak_cw1 <- ggplot(filter(subpool5, weak <= 1),
                                 aes(x = as.factor(total_sites), 
                                     y = ave_ratio_25_norm)) +
  geom_boxplot(aes(color = as.factor(weak)), show.legend = FALSE,
               position = position_dodge(1)) +
  scale_y_log10() +
  facet_grid(background ~ .) + 
  panel_border() +
  scale_color_manual(values = cbPalette7) +
  annotation_logticks(sides = 'lr') +
  background_grid(major = 'y', minor = 'none') + 
  geom_vline(xintercept = c(1.5:6.5), alpha = 0.5) +
  ylab('log10 background-normalized\nexpression at 25 µM') +
  xlab('Number of binding sites')

save_plot('plots/p_num_sites_num_weak_cw1.png', p_num_sites_num_weak_cw1, scale = 1.3,
          base_width = 6, base_height = 4)

p_num_sites_num_weak_cw2 <- ggplot(filter(subpool5, weak <= 2),
                                  aes(x = as.factor(total_sites), 
                                      y = ave_ratio_25_norm)) +
  geom_boxplot(aes(color = as.factor(weak)), show.legend = FALSE,
               position = position_dodge(1)) +
  scale_y_log10() +
  facet_grid(background ~ .) + 
  panel_border() +
  scale_color_manual(values = cbPalette7) +
  annotation_logticks(sides = 'lr') +
  background_grid(major = 'y', minor = 'none') + 
  geom_vline(xintercept = c(1.5:6.5), alpha = 0.5) +
  ylab('log10 background-normalized\nexpression at 25 µM') +
  xlab('Number of binding sites')

save_plot('plots/p_num_sites_num_weak_cw2.png', p_num_sites_num_weak_cw2, scale = 1.3,
          base_width = 6, base_height = 4)

p_num_sites_num_weak_cw6 <- ggplot(subpool5,
                                  aes(x = as.factor(total_sites), 
                                      y = ave_ratio_25_norm)) +
  geom_boxplot(aes(color = as.factor(weak)), show.legend = FALSE,
               position = position_dodge(1)) +
  scale_y_log10() +
  facet_grid(background ~ .) + 
  panel_border() +
  scale_color_manual(values = cbPalette7) +
  annotation_logticks(sides = 'lr') +
  background_grid(major = 'y', minor = 'none') + 
  geom_vline(xintercept = c(1.5:6.5), alpha = 0.5) +
  ylab('log10 background-normalized\nexpression at 25 µM') +
  xlab('Number of binding sites')

save_plot('plots/p_num_sites_num_weak_cw6.png', p_num_sites_num_weak_cw6, scale = 1.3,
          base_width = 6, base_height = 4)


#Filter out sites that have higher induced expression than their 6 consensus 
#site-counterpart

exp_25_greater_c6 <- function(df) {
  c6 <- df %>%
    filter(consensus == 6) %>%
    select(background, induction, ave_ratio_0_norm, ave_ratio_25_norm) %>%
    rename(induction_c6 = induction) %>%
    rename(ave_ratio_0_norm_c6 = ave_ratio_0_norm) %>%
    rename(ave_ratio_25_norm_c6 = ave_ratio_25_norm)
  c6_exp_25_join <- left_join(df, c6, by = 'background') %>%
    group_by(background) %>%
    filter(ave_ratio_25_norm > ave_ratio_25_norm_c6) %>%
    ungroup()
  return(c6_exp_25_join) 
}

subpool5_exp_25_greater_c6 <- exp_25_greater_c6(subpool5)

subpool5_exp_25_greater_c6 %>%
  group_by(background) %>%
  summarize(sequences = n())

#Count site type per site location

site_loc_type_count <- function(df) {
  site1 <- df %>%
    group_by(background) %>%
    count(site1, site1 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site1) %>%
    select(-3) %>%
    mutate(site = 1)
  site2 <- df %>%
    group_by(background) %>%
    count(site2, site2 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site2) %>%
    select(-3) %>%
    mutate(site = 2)
  site3 <- df %>%
    group_by(background) %>%
    count(site3, site3 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site3) %>%
    select(-3) %>%
    mutate(site = 3)
  site4 <- df %>%
    group_by(background) %>%
    count(site4, site4 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site4) %>%
    select(-3) %>%
    mutate(site = 4)
  site5 <- df %>%
    group_by(background) %>%
    count(site5, site5 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site5) %>%
    select(-3) %>%
    mutate(site = 5)
  site6 <- df %>%
    group_by(background) %>%
    count(site6, site6 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site6) %>%
    select(-3) %>%
    mutate(site = 6)
  site_join <- bind_rows(site1, site2, site3, site4, site5, site6)
  return(site_join)
}

subpool5_exp_25_greater_c6_sites <- site_loc_type_count(subpool5_exp_25_greater_c6)

p_site_exp_25_greater_c6 <- ggplot(subpool5_exp_25_greater_c6_sites,
                                   aes(as.factor(site), counts, fill = type)) +
  facet_grid(. ~ background) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_viridis(discrete = TRUE) + 
  xlab('Site position') + 
  panel_border()

save_plot('plots/p_site_exp_25_greater_c6.png', p_site_exp_25_greater_c6, 
          base_width = 4, base_height = 2, scale = 1.5)

compare_2_5 <- inner_join(subpool2, subpool5, 
                          by = c('most_common', 'background', 'ratio_0A_back',
                                 'ratio_0B_back', 'ratio_25A_back', 
                                 'ratio_25B_back'),
                          suffix = c('_sp2', '_sp5')) %>%
  ungroup()

s5_single_site_exp <- compare_2_5 %>%
  filter(weak == 0 & consensus == 1) %>%
  mutate(site1 = str_detect(site1, "consensus") * 1) %>%
  mutate(site2 = str_detect(site2, "consensus") * 2) %>%
  mutate(site3 = str_detect(site3, "consensus") * 3) %>%
  mutate(site4 = str_detect(site4, "consensus") * 4) %>%
  mutate(site5 = str_detect(site5, "consensus") * 5) %>%
  mutate(site6 = str_detect(site6, "consensus") * 6) %>%
  mutate(site = site1 + site2 + site3 + site4 + site5 + site6)

p_s5_single_site_exp <- ggplot(s5_single_site_exp, 
                               aes(x = as.factor(site))) +
  geom_bar(aes(y = ave_ratio_25_norm_med),
           stat = 'identity', fill = '#440154FF') +
  geom_point(aes(y = ratio_25A_norm_med), 
             fill = '#20A387FF', shape = 21, size = 2) +
  geom_point(aes(y = ratio_25B_norm_med), 
             fill = '#20A387FF', shape = 21, size = 2) +
  facet_grid(. ~ background) +
  xlab('Site position') + 
  panel_border() +
  ylab('Expression at 25 µM normalized\nto median 1-site expression') +
  geom_hline(yintercept = 1, color = 'gray')

save_plot('plots/p_s5_single_site_exp.png', p_s5_single_site_exp,
          scale = 1.2, base_width = 4.5, base_height = 3)


#Do weak sites operate as no_sites or do they still contribute to expression?

p_num_cons_num_weak <- ggplot(subpool5, 
                              aes(x = as.factor(consensus), 
                                  y = ave_ratio_25_norm)) +
  geom_boxplot(aes(color = as.factor(weak)), show.legend = TRUE,
               position = position_dodge(1)) +
  scale_y_log10() +
  facet_grid(background ~ .) + 
  panel_border() +
  annotation_logticks(sides = 'lr') +
  scale_color_manual(values = cbPalette7) +
  background_grid(major = 'y', minor = 'none') + 
  geom_vline(xintercept = c(1.5:6.5), alpha = 0.5) +
  ylab('log10 expression at 25 µM') +
  xlab('Number of consensus sites')

save_plot('plots/p_num_cons_num_weak.png', p_num_cons_num_weak, scale = 1.3,
          base_width = 6, base_height = 4)


#Look at site architecture changes per background per n% population bins that 
#span induced expression

s5_binned_exp_25_back <- subpool5 %>%
  group_by(background) %>%
  mutate(bin = ntile(ave_ratio_25_norm, 10))

site_loc_type_count_back_bin <- function(df) {
  site1 <- df %>%
    group_by(background, bin) %>%
    count(site1, site1 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site1) %>%
    select(-4) %>%
    mutate(site = 1)
  site2 <- df %>%
    group_by(background, bin) %>%
    count(site2, site2 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site2) %>%
    select(-4) %>%
    mutate(site = 2)
  site3 <- df %>%
    group_by(background, bin) %>%
    count(site3, site3 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site3) %>%
    select(-4) %>%
    mutate(site = 3)
  site4 <- df %>%
    group_by(background, bin) %>%
    count(site4, site4 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site4) %>%
    select(-4) %>%
    mutate(site = 4)
  site5 <- df %>%
    group_by(background, bin) %>%
    count(site5, site5 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site5) %>%
    select(-4) %>%
    mutate(site = 5)
  site6 <- df %>%
    group_by(background, bin) %>%
    count(site6, site6 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site6) %>%
    select(-4) %>%
    mutate(site = 6)
  site_join <- bind_rows(site1, site2, site3, site4, site5, site6)
  return(site_join)
}

s5_binned_exp_25_back_sites <- site_loc_type_count_back_bin(s5_binned_exp_25_back)

p_s5_binned_exp_25_back_sites <- ggplot(s5_binned_exp_25_back_sites,
                                   aes(as.factor(site), counts, fill = type)) +
  facet_grid(background ~ bin) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_viridis(discrete = TRUE) + 
  xlab('Site position')

#Binning by expression alone

s5_binned_exp_25 <- subpool5 %>%
  mutate(bin = ntile(ave_ratio_25_norm, 5))

ggplot(s5_binned_exp_25, aes(ratio_25A_norm, ratio_25B_norm, fill = bin)) +
  facet_grid(background ~ bin) + 
  scale_fill_viridis() +
  panel_border() +
  scale_x_log10() + scale_y_log10() +
  geom_point(alpha = 0.3, shape = 21)

site_loc_type_count_bin <- function(df) {
  site1 <- df %>%
    group_by(bin) %>%
    count(site1, site1 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site1) %>%
    select(-3) %>%
    mutate(site = 1)
  site2 <- df %>%
    group_by(bin) %>%
    count(site2, site2 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site2) %>%
    select(-3) %>%
    mutate(site = 2)
  site3 <- df %>%
    group_by(bin) %>%
    count(site3, site3 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site3) %>%
    select(-3) %>%
    mutate(site = 3)
  site4 <- df %>%
    group_by(bin) %>%
    count(site4, site4 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site4) %>%
    select(-3) %>%
    mutate(site = 4)
  site5 <- df %>%
    group_by(bin) %>%
    count(site5, site5 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site5) %>%
    select(-3) %>%
    mutate(site = 5)
  site6 <- df %>%
    group_by(bin) %>%
    count(site6, site6 == 'consensus') %>%
    rename(counts = n) %>%
    rename(type = site6) %>%
    select(-3) %>%
    mutate(site = 6)
  site_join <- bind_rows(site1, site2, site3, site4, site5, site6)
  return(site_join)
}

s5_binned_exp_25_sites <- site_loc_type_count_bin(s5_binned_exp_25)

p_s5_binned_exp_25_sites <- ggplot(filter(s5_binned_exp_25_sites, 
                                          type == 'consensus'),
                                   aes(as.factor(site), counts, fill = type)) +
  facet_grid(. ~ bin) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_viridis(discrete = TRUE) + 
  xlab('Site position')

ggplot(data = NULL, aes(x = "", y = ave_ratio_25_norm)) +
  facet_grid(background ~ .) +
  geom_boxplot(data = filter(subpool5, site1 == 'consensus' & site_combo == 'consensus'), 
               aes(x = 'site 1')) +
  geom_boxplot(data = filter(subpool5, site2 == 'consensus' & site_combo == 'consensus'),
               aes(x = 'site 2')) +
  geom_boxplot(data = filter(subpool5, site3 == 'consensus' & site_combo == 'consensus'), 
               aes(x = 'site 3')) +
  geom_boxplot(data = filter(subpool5, site4 == 'consensus' & site_combo == 'consensus'),
               aes(x = 'site 4')) +
  geom_boxplot(data = filter(subpool5, site5 == 'consensus' & site_combo == 'consensus'), 
               aes(x = 'site 5')) +
  geom_boxplot(data = filter(subpool5, site6 == 'consensus' & site_combo == 'consensus'), 
               aes(x = 'site 6')) +
  scale_y_log10() + annotation_logticks(sides = 'l') +
  ylab("log10 average background-norm.\nreads RNA/DNA at 25 µM") + 
  xlab('Site position') +
  panel_border()


#There is almost a linear relationship between induced expression and induction
#with high-expression uninduced variants as outliers

p_25_ind_s5 <- ggplot(subpool5, aes(ave_ratio_25_norm, induction)) +
  facet_grid(. ~ background) +
  panel_border() +
  geom_point(aes(fill = log10(ave_ratio_0_norm)), shape = 21) +
  scale_fill_viridis() +
  xlab('Log10 normalized expression at 25 µM') +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks(sides = 'bl') +
  ylab('Log10 induction')

save_plot('plots/p_25_ind_s5.png', p_25_ind_s5, scale = 1.5, base_width = 4,
          base_height = 4)

#There is a minor population of very active sequences at 0 µM

p_exp_distrib_0 <- ggplot(subpool5, aes(ave_ratio_0_norm)) +
  facet_grid(. ~ background) +
  scale_x_log10() +
  geom_density(kernel = 'gaussian') +
  annotation_logticks(scaled = TRUE, sides = 'b') +
  xlab("Log10 average background-norm.\nreads RNA/DNA") +
  geom_vline(xintercept = 4) +
  panel_border()

save_plot('plots/p_exp_distrib_0.png', p_exp_distrib_0, 
          base_width = 3, base_height = 2, scale = 1.5)

s5_highexp_0 <-subpool5 %>%
  filter(ave_ratio_0_norm >= 4)

s5_lowexp_0 <-subpool5 %>%
  filter(ave_ratio_0_norm < 4)

p_pop_unind_exp25_vs_ind <- ggplot(data = NULL, aes(ave_ratio_25_norm, induction)) +
  facet_grid(background ~ .) +
  panel_border() +
  geom_point(data = s5_lowexp_0, fill = '#482677FF', shape = 21) +
  geom_point(data = s5_highexp_0, fill = '#B8DE29FF', shape = 21) +
  xlab('Log10 normalized expression at 25 µM') +
  scale_x_log10() + scale_y_log10() +
  annotation_logticks(sides = 'bl') +
  ylab('Log10 induction')

save_plot('plots/p_pop_unind_exp25_vs_ind.png', p_pop_unind_exp25_vs_ind, 
          scale = 1.5, base_width = 2.5, base_height = 4)

#Count site type per location in these high uninduced-expression variants

s5_highexp_0_sites <- site_loc_type_count(s5_highexp_0)
s5_lowexp_0_sites <- site_loc_type_count(s5_lowexp_0)

p_site_highexp_0 <- ggplot(s5_highexp_0_sites,
                           aes(as.factor(site), counts, fill = type)) +
  facet_grid(. ~ background) +
  geom_bar(stat = 'identity', position = 'fill') +
  scale_fill_viridis(discrete = TRUE) + 
  xlab('Site position')

save_plot('plots/p_site_highexp_0.png', p_site_highexp_0, 
          scale = 2, base_width = 3, base_height = 1.5)


#Fitting a model to the 6-site library------------------------------------------

#subset to only fit to consensus sites

s5_onlycons_log2 <- subpool5 %>%
  filter(weak == 0) %>%
  var_log2()

bin_site_s5 <- s5_onlycons_log2 %>%
  mutate(site1 = str_detect(site1, 'consensus') + 0) %>%
  mutate(site2 = str_detect(site2, 'consensus') + 0) %>%
  mutate(site3 = str_detect(site3, 'consensus') + 0) %>%
  mutate(site4 = str_detect(site4, 'consensus') + 0) %>%
  mutate(site5 = str_detect(site5, 'consensus') + 0) %>%
  mutate(site6 = str_detect(site6, 'consensus') + 0)

pred_resid <- function(df1, x) {
  df2 <- df1 %>%
    add_predictions(x)
  df3 <- df2 %>%
    add_residuals(x)
  return(df3)
  print('processed pre_res_trans_int(df1, df2) in order of (data, model)')
}

#Linear models

#Just use total sites, independent background variable

totsite_ind_back <- function(df) {
  model <- lm(ave_ratio_25_norm ~ total_sites + background, data = df)
}

totsite_ind_back_fit <- totsite_ind_back(bin_site_s5)
summary(totsite_ind_back_fit)
anova(totsite_ind_back_fit)
totsite_ind_back_p_r <- pred_resid(bin_site_s5, totsite_ind_back_fit)

ggplot(totsite_ind_back_p_r, aes(x = as.factor(total_sites))) +
  facet_grid(. ~ background) +
  geom_boxplot(aes(y = ave_ratio_25_norm)) +
  geom_boxplot(aes(y = pred), color = 'red')

p_totsite_ind_back_r <- ggplot(totsite_ind_back_p_r, aes(ave_ratio_25_norm, pred,
                                                         fill = total_sites)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis() +
  annotation_logticks(sides = 'bl') +
  xlab('log2 observed expression') + ylab('log2 predicted expression') +
  annotate("text", x = 1, y = 5, 
           label = paste('r =', 
                         round(cor(totsite_ind_back_p_r$pred,
                                   totsite_ind_back_p_r$ave_ratio_25_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))


#All independent sites, no background variable

ind_site <- function(df) {
  model <- lm(ave_ratio_25_norm ~ 
                site1 + site2 + site3 + site4 + site5 + site6, data = df)
}

ind_site_fit <- ind_site(bin_site_s5)
summary(ind_site_fit)
anova(ind_site_fit)
ind_site_p_r <- pred_resid(bin_site_s5, ind_site_fit)

ggplot(ind_site_p_r, aes(x = as.factor(total_sites))) +
  facet_grid(. ~ background) +
  geom_boxplot(aes(y = ave_ratio_25_norm)) +
  geom_boxplot(aes(y = pred), color = 'red')

ggplot(ind_site_p_r, aes(ave_ratio_25_norm, pred,
                         fill = total_sites)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis() +
  annotation_logticks(sides = 'bl') +
  xlab('log2 observed expression') + ylab('log2 predicted expression') +
  annotate("text", x = 1.6, y = 4, 
           label = paste('r =', 
                         round(cor(ind_site_p_r$pred,
                                   ind_site_p_r$ave_ratio_25_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

#All independent sites, with independent background variable

ind_site_ind_back <- function(df) {
  model <- lm(ave_ratio_25_norm ~ 
                site1 + site2 + site3 + site4 + site5 + site6 + background, 
              data = df)
}

ind_site_ind_back_fit <- ind_site_ind_back(bin_site_s5)
summary(ind_site_ind_back_fit)
anova(ind_site_ind_back_fit)
ind_site_ind_back_p_r <- pred_resid(bin_site_s5, ind_site_ind_back_fit)

ggplot(ind_site_ind_back_p_r, aes(x = as.factor(total_sites))) +
  facet_grid(. ~ background) +
  geom_boxplot(aes(y = ave_ratio_25_norm)) +
  geom_boxplot(aes(y = pred), color = 'red')

p_ind_site_ind_back <- ggplot(ind_site_ind_back_p_r, aes(ave_ratio_25_norm, pred,
                                                         fill = total_sites)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis() +
  annotation_logticks(sides = 'bl') +
  xlab('log2 observed expression') + ylab('log2 predicted expression') +
  scale_y_continuous(limits = c(-1, 7.5)) +
  annotate("text", x = 1, y = 5, 
           label = paste('r =', 
                         round(cor(ind_site_ind_back_p_r$pred,
                                   ind_site_ind_back_p_r$ave_ratio_25_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

p_tsite_indsite <- plot_grid(p_totsite_ind_back_r, 
                             p_ind_site_ind_back, nrow = 2)

save_plot('plots/p_tsite_indsite.png', p_tsite_indsite, 
          base_width = 5, base_height = 7, scale = 1.2)

#Including weak sites, all independent, independent background 

subpool5_log2 <- var_log2(subpool5)

ind_site_ind_back_fit_nwc <- ind_site_ind_back(subpool5_log2)
summary(ind_site_ind_back_fit_nwc)
summary_nwc <- tidy(ind_site_ind_back_fit_nwc) %>%
  filter(str_detect(term, '^site')) %>%
  mutate(term = gsub('nosite', '_nosite', term)) %>%
  mutate(term = gsub('weak', '_weak', term)) %>%
  separate(term, into = c('site', 'type'), sep = "_")

anova_nwc <- tidy(anova(ind_site_ind_back_fit_nwc))

ind_site_ind_back_nwc_p_r <- pred_resid(subpool5_log2, ind_site_ind_back_fit_nwc)

ggplot(ind_site_ind_back_nwc_p_r, aes(x = as.factor(total_sites))) +
  facet_grid(. ~ background) +
  geom_boxplot(aes(y = ave_ratio_25_norm)) +
  geom_boxplot(aes(y = pred), color = 'red')

p_ind_site_ind_back_nwc <- ggplot(ind_site_ind_back_nwc_p_r, 
                                  aes(ave_ratio_25_norm, pred, 
                                      fill = consensus)) +
  geom_point(shape = 21, alpha = 0.5) +
  scale_fill_viridis() +
  annotation_logticks(sides = 'bl') +
  xlab('log2 observed expression') + ylab('log2 predicted\nexpression') +
  scale_y_continuous(limits = c(-1.5, 8.5)) +
  annotate("text", x = 1, y = 5, 
           label = paste('r =', 
                         round(cor(ind_site_ind_back_nwc_p_r$pred,
                                   ind_site_ind_back_nwc_p_r$ave_ratio_25_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

p_ind_site_ind_back_nwc_sum <- ggplot(summary_nwc, aes(site, estimate, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_x_discrete(position = 'top') + 
  scale_fill_viridis(discrete = TRUE) + 
  ylab('log2 weight relative\nto consensus')

p_ind_site_ind_back_nwc_anova <- anova_nwc %>% 
  mutate(term_fctr = factor(term, levels = term)) %>% 
  ggplot(aes(term_fctr, sumsq)) + 
    geom_bar(stat = 'identity') + 
  ylab('Sum of squares') +
  xlab('Model term') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_ind_sit_ind_back_nwc_grid <- plot_grid(p_ind_site_ind_back_nwc, 
          p_ind_site_ind_back_nwc_sum, 
          p_ind_site_ind_back_nwc_anova,
          nrow = 3)

save_plot('plots/p_ind_sit_ind_back_nwc_grid.png', p_ind_sit_ind_back_nwc_grid,
          base_width = 5, base_height = 8)

#All independent sites, dependent background, not sure how useful this one is,
#again, too many parameters to fit to

ind_site_dep_back <- function(df) {
  model <- lm(ave_ratio_25_norm ~ 
                (site1 + site2 + site3 + site4 + site5 + site6) * background, 
              data = df)
}

ind_site_dep_back_fit <- ind_site_dep_back(bin_site_s5)
summary(ind_site_dep_back_fit)
anova(ind_site_dep_back_fit)
ind_site_dep_back_p_r <- pred_resid(bin_site_s5, ind_site_dep_back_fit)

ggplot(ind_site_dep_back_p_r, aes(x = as.factor(total_sites))) +
  facet_grid(. ~ background) +
  geom_boxplot(aes(y = ave_ratio_25_norm)) +
  geom_boxplot(aes(y = pred), color = 'red')

ggplot(ind_site_dep_back_p_r, aes(ave_ratio_25_norm, pred,
                                  fill = total_sites)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis() +
  annotation_logticks(sides = 'bl') +
  xlab('log2 observed expression') + ylab('log2 predicted expression') +
  annotate("text", x = 1.6, y = 4, 
           label = paste('r =', 
                         round(cor(ind_site_dep_back_p_r$pred,
                                   ind_site_dep_back_p_r$ave_ratio_25_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

#All dependent sites, independent background, this essentially just over-fits 
#the data, looks good though

dep_site_ind_back <- function(df) {
  model <- lm(ave_ratio_25_norm ~ 
                site1 * site2 * site3 * site4 * site5 * site6 + background, 
              data = df)
}

mm <- model_matrix(bin_site_s5, ave_ratio_25_norm ~ 
                     site1 * site2 * site3 * site4 * site5 * site6 + background)

dep_site_ind_back_fit <- dep_site_ind_back(bin_site_s5)
test <- summary(dep_site_ind_back_fit)
anova(dep_site_ind_back_fit)
dep_site_ind_back_p_r <- pred_resid(bin_site_s5, dep_site_ind_back_fit)

ggplot(dep_site_ind_back_p_r, aes(x = as.factor(total_sites))) +
  facet_grid(. ~ background) +
  geom_boxplot(aes(y = ave_ratio_25_norm)) +
  geom_boxplot(aes(y = pred), color = 'red')

ggplot(dep_site_ind_back_p_r, aes(ave_ratio_25_norm, pred,
                                  fill = total_sites)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis() +
  annotation_logticks(sides = 'bl') +
  xlab('log2 observed expression') + ylab('log2 predicted expression') +
  annotate("text", x = 1.6, y = 4, 
           label = paste('r =', 
                         round(cor(dep_site_ind_back_p_r$pred,
                                   dep_site_ind_back_p_r$ave_ratio_25_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

#Non-linear models

hill_like_model <- function(df) {
  n_init <- 1
  site_half_max_init <- 1
  max_ave_ratio_25_norm_init <- 60
  hill_nls <- nls(
    ave_ratio_25_norm ~ I(max_ave_ratio_25_norm * (total_sites)^n) / I(site_half_max + (total_sites)^n), 
    data = df, start = c(site_half_max = site_half_max_init, n = n_init,
                         max_ave_ratio_25_norm = max_ave_ratio_25_norm_init))
  return(hill_nls)
}

hill_total_site_ind_back_fit <- hill_like_model(bin_site_s5)
summary(hill_total_site_ind_back_fit)
hill_total_site_ind_back_p_r <- pred_resid(bin_site_s5, 
                                           hill_total_site_ind_back_fit)

ggplot(hill_total_site_ind_back_p_r, aes(ave_ratio_25_norm, pred, 
                                         fill = total_sites)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis() +
  annotation_logticks(sides = 'bl') +
  xlab('log2 observed expression') + ylab('log2 predicted expression') +
  annotate("text", x = 1.6, y = 4, 
           label = paste('r =', 
                         round(cor(hill_total_site_ind_back_p_r$pred,
                                   hill_total_site_ind_back_p_r$ave_ratio_25_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

log_curve_model <- function(df) {
  n_init <- 1
  site_half_max_init <- 1
  max_ave_ratio_25_norm_init <- 60
  log_curve_nls <- nls(
    ave_ratio_25_norm ~ max_ave_ratio_25_norm/(1 + exp(-n * (total_sites - site_half_max))),
    data = df, start = c(n = n_init, site_half_max = site_half_max_init, 
                         max_ave_ratio_25_norm = max_ave_ratio_25_norm_init))
  return(log_curve_nls)
}

log_curve_fit <- log_curve_model(bin_site_s5)
summary(log_curve_fit)
log_curve_p_r <- pred_resid(bin_site_s5, log_curve_fit)

ggplot(log_curve_p_r , aes(x = as.factor(total_sites))) +
  facet_grid(. ~ background) +
  geom_boxplot(aes(y = ave_ratio_25_norm)) +
  geom_boxplot(aes(y = pred), color = 'red')

p_log_curve <- ggplot(log_curve_p_r, aes(ave_ratio_25_norm, pred, fill = total_sites)) +
  geom_point(shape = 21, alpha = 0.7) +
  scale_fill_viridis() +
  annotation_logticks(sides = 'bl') +
  xlab('log2 observed expression') + ylab('log2 predicted expression') +
  scale_y_continuous(limits = c(-1, 7.5)) +
  annotate("text", x = 1, y = 6, 
           label = paste('r =', 
                         round(cor(log_curve_p_r$pred,
                                   log_curve_p_r$ave_ratio_25_norm,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"), 2)))

save_plot('plots/p_log_curve.png', p_log_curve, 
          base_width = 4, base_height = 3, scale = 1.2)


#Median analysis of expression--------------------------------------------------

#Filter DNA reads to set a minimum of 3 reads per BC join with RNA and determine
#RNA/DNA

bc_dna_join_rna <- function(df1, df2) {
  filter_reads <- filter(df1, num_reads > 2)
  DNA_RNA_join <- left_join(filter_reads, df2,
                            by = c("barcode", "name", "subpool", 
                                   "most_common"), 
                            suffix = c('_DNA', '_RNA')) %>%
    mutate(ratio = normalized_RNA/normalized_DNA)
  print('processed dfs in order of (DNA, RNA) in bc_dna_join_rna(df1, df2)')
  return(DNA_RNA_join)
}

RNA_DNA_bc_R0A <- bc_dna_join_rna(bc_join_DNA, bc_join_R0A)
RNA_DNA_bc_R0B <- bc_dna_join_rna(bc_join_DNA, bc_join_R0B)
RNA_DNA_bc_R25A <- bc_dna_join_rna(bc_join_DNA, bc_join_R25A)
RNA_DNA_bc_R25B <- bc_dna_join_rna(bc_join_DNA, bc_join_R25B)

#Count barcodes per variant per DNA and RNA, set minimum of 7 BC's per variant 
#in DNA sample, take median RNA/DNA per variant, find absolute deviation for
#each BC per variant then per variant determine the median absolute deviation.

ratio_bc_med_var <- function(df1) {
  bc_count_DNA <- df1 %>%
    group_by(subpool, name, most_common) %>%
    summarize(barcodes_DNA = n()) %>%
    filter(barcodes_DNA > 7)
  bc_count_RNA <- df1 %>%
    group_by(subpool, name, most_common) %>%
    filter(num_reads_RNA != 0) %>%
    summarize(barcodes_RNA = n())
  bc_DNA_RNA <- inner_join(bc_count_DNA, bc_count_RNA, 
                           by = c('subpool', 'name', 'most_common'))
  med_ratio <- df1 %>%
    group_by(subpool, name, most_common) %>%
    summarize(med_ratio = median(ratio))
  mad_ratio <- inner_join(df1, med_ratio, 
                          by = c('subpool', 'name', 'most_common')) %>%
    mutate(absdev = abs(ratio - med_ratio)) %>%
    group_by(subpool, name, most_common) %>%
    summarize(mad = median(absdev))
  med_mad <- inner_join(med_ratio, mad_ratio, 
                        by = c('subpool', 'name', 'most_common')) %>%
    mutate(mad_over_med = as.double(mad/med_ratio)) %>%
    mutate(mad_over_med = if_else(
      is.na(mad_over_med),
      as.double(0), 
      mad_over_med))
  bc_med <- inner_join(med_mad, bc_DNA_RNA, 
                       by = c('subpool', 'name', 'most_common')) %>%
    ungroup()
  return(bc_med)
}

med_RNA_DNA_R0A <- ratio_bc_med_var(RNA_DNA_bc_R0A)
med_RNA_DNA_R0B <- ratio_bc_med_var(RNA_DNA_bc_R0B)
med_RNA_DNA_R25A <- ratio_bc_med_var(RNA_DNA_bc_R25A)
med_RNA_DNA_R25B <- ratio_bc_med_var(RNA_DNA_bc_R25B)


#Combine all samples

med_var_rep <- function(df0A, df0B, df25A, df25B) {
  join_0 <- inner_join(df0A, df0B, 
                       by = c("name", "subpool", "most_common", 'barcodes_DNA'), 
                       suffix = c("_0A", "_0B"))
  join_25 <- inner_join(df25A, df25B, 
                        by = c("name", "subpool", "most_common", 
                               'barcodes_DNA'), 
                        suffix = c("_25A", "_25B"))
  join_0_25 <- inner_join(join_0, join_25, 
                          by = c("name", "subpool", "most_common", 
                                 'barcodes_DNA'))
  print('processed dfs in order: df0A, df0B, df25A, df25B')
  return(join_0_25)
}

med_rep_1_2 <- med_var_rep(med_RNA_DNA_R0A, med_RNA_DNA_R0B, 
                           med_RNA_DNA_R25A, med_RNA_DNA_R25B)

med_rep_1_2_0corr <- med_rep_1_2 %>%
  mutate(med_ratio_0A = if_else(
    med_ratio_0A  == as.double(0),
    as.double(0.01), 
    med_ratio_0A )) %>%
  mutate(med_ratio_0B = if_else(
    med_ratio_0B  == as.double(0),
    as.double(0.01), 
    med_ratio_0B )) %>%
  mutate(med_ratio_25A = if_else(
    med_ratio_25A  == as.double(0),
    as.double(0.01), 
    med_ratio_25A )) %>%
  mutate(med_ratio_25B = if_else(
    med_ratio_25B  == as.double(0),
    as.double(0.01), 
    med_ratio_25B ))
 

#Following mad/med per med

med_mad_gather <- function(df1) {
  med <- df1 %>%
    select(subpool, name, med_ratio_0A, med_ratio_0B, med_ratio_25A, 
           med_ratio_25B) %>%
    gather(med_ratio_0A, med_ratio_0B, med_ratio_25A, med_ratio_25B,
           key = condition, value = med) %>%
    separate(condition, into = c("fluff1", "fluff2", "condition"),
             sep = "_", convert = TRUE) %>% 
    select(-fluff1, -fluff2)
  mad_over_med <- df1 %>%
    select(subpool, name, mad_over_med_0A, mad_over_med_0B, mad_over_med_25A, 
           mad_over_med_25B) %>%
    gather(mad_over_med_0A, mad_over_med_0B, mad_over_med_25A, mad_over_med_25B,
           key = condition, value = mad_over_med) %>%
    separate(condition, into = c("fluff1", "fluff2", "fluff3", "condition"),
             sep = "_", convert = TRUE) %>% 
    select(-fluff1, -fluff2, -fluff3)
  med_mad <- inner_join(med, mad_over_med, 
                        by = c('subpool', 'name', 'condition'))
  return(med_mad)
}

med_mad_rep_1_2 <- med_mad_gather(med_rep_1_2)

p_mad_over_med_distr <- ggplot(med_mad_rep_1_2, aes(mad_over_med)) +
  geom_density(kernel = 'gaussian') +
  facet_grid(subpool ~ condition) + 
  panel_border() +
  xlab('MAD/Med')

save_plot('plots/p_mad_over_med_distr_BCcutoff.png', p_mad_over_med_distr,
          base_width = 7, base_height = 5)

sum_rep_1_2 <- rep_1_2 %>%
  select(subpool, name, ratio_0A, ratio_0B, ratio_25A, ratio_25B) %>%
  gather(ratio_0A, ratio_0B, ratio_25A, ratio_25B,
         key = condition, value = sum) %>%
  separate(condition, into = c("fluff1", "condition"),
           sep = "_", convert = TRUE) %>% 
  select(-fluff1)

sum_med_var_cond <- inner_join(sum_rep_1_2, med_mad_rep_1_2,
                          by = c('subpool', 'name', 'condition')) %>%
  mutate(med = if_else(
    med == as.double(0),
    as.double(0.01), 
    med))


#Compare med replicates and mad to sum, focusing on mad/median = 1 contribution 
#to replicability

sum_med_var <- inner_join(rep_1_2, med_rep_1_2_0corr,
                          by = c('subpool', 'name', 'most_common'),
                          suffix = c('sum', 'med')) %>%
  mutate(ave_sum_0 = (ratio_0A + ratio_0B)/2) %>%
  mutate(ave_sum_25 = (ratio_25A + ratio_25B)/2) %>%
  mutate(ave_med_0 = (med_ratio_0A + med_ratio_0B)/2) %>%
  mutate(ave_med_25 = (med_ratio_25A + med_ratio_25B)/2) %>%
  var_log10()

p_med_vs_sum_0_rep <- ggplot(data = NULL, aes(ave_sum_0, ave_med_0)) +
  facet_grid(. ~ subpool) +
  geom_point(data = filter(sum_med_var, subpool == 'subpool5'),
             color = '#440154FF', alpha = 0.3) +
  geom_point(data = filter(sum_med_var, subpool == 'subpool3'),
             color = '#33638DFF', alpha = 0.3) +
  geom_point(data = filter(sum_med_var, subpool == 'subpool4'),
             color = '#29AF7FFF', alpha = 0.3) +
  geom_point(data = filter(sum_med_var, subpool == 'subpool2'),
             color = '#DCE319FF', alpha = 0.3) +
  geom_point(data = filter(sum_med_var, 
                           mad_over_med_0A == as.double(0) & subpool != 'control'),
             color = 'red', alpha = 0.7) +
  geom_point(data = filter(sum_med_var, 
                           mad_over_med_0B == as.double(0) & subpool != 'control'),
             color = 'red', alpha = 0.7) +
  xlab("log10 Ave Sum RNA/DNA 0A") +
  ylab("log10 Ave Med RNA/DNA 0B") +
  scale_x_continuous(limits = c(-2, 2)) + 
  scale_y_continuous(limits = c(-2, 2)) +
  annotation_logticks(sides = 'bl') +
  background_grid(major = 'xy', minor = 'none') + 
  panel_border() +
  annotate("text", x =  -1, y = 1.75, color = '#440154FF', 
           label = paste('r =', 
                         round(cor(filter(sum_med_var, 
                                          subpool == 'subpool5')$ave_med_0,
                                   filter(sum_med_var, 
                                          subpool == 'subpool5')$ave_sum_0,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = -1, y = 1.25, color = '#33638DFF', 
           label = paste('r =', 
                         round(cor(filter(sum_med_var, 
                                          subpool == 'subpool3')$ave_med_0,
                                   filter(sum_med_var, 
                                          subpool == 'subpool3')$ave_sum_0,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = -1, y = 0.75, color = '#29AF7FFF', 
           label = paste('r =', 
                         round(cor(filter(sum_med_var, 
                                          subpool == 'subpool4')$ave_med_0,
                                   filter(sum_med_var, 
                                          subpool == 'subpool4')$ave_sum_0,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = -1, y = 0.25, color = '#DCE319FF', 
           label = paste('r =', 
                         round(cor(filter(sum_med_var, 
                                          subpool == 'subpool2')$ave_med_0,
                                   filter(sum_med_var, 
                                          subpool == 'subpool2')$ave_sum_0,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2)))

p_med_vs_sum_25_rep <- ggplot(data = NULL, aes(ave_sum_25, ave_med_25)) +
  facet_grid(. ~ subpool) +
  geom_point(data = filter(sum_med_var, subpool == 'subpool5'),
             color = '#440154FF', alpha = 0.3) +
  geom_point(data = filter(sum_med_var, subpool == 'subpool3'),
             color = '#33638DFF', alpha = 0.3) +
  geom_point(data = filter(sum_med_var, subpool == 'subpool4'),
             color = '#29AF7FFF', alpha = 0.3) +
  geom_point(data = filter(sum_med_var, subpool == 'subpool2'),
             color = '#DCE319FF', alpha = 0.3) +
  geom_point(data = filter(sum_med_var, 
                           mad_over_med_25A == as.double(0) & subpool != 'control'),
             color = 'red', alpha = 0.7) +
  geom_point(data = filter(sum_med_var, 
                           mad_over_med_25B == as.double(0) & subpool != 'control'),
             color = 'red', alpha = 0.7) +
  xlab("log10 Ave Sum RNA/DNA 25A") +
  ylab("log10 Ave Med RNA/DNA 25B") +
  scale_x_continuous(limits = c(-2, 2)) + 
  scale_y_continuous(limits = c(-2, 2)) +
  annotation_logticks(sides = 'bl') +
  background_grid(major = 'xy', minor = 'none') + 
  panel_border() +
  annotate("text", x = -1, y = 1.75, color = '#440154FF', 
           label = paste('r =', 
                         round(cor(filter(sum_med_var, 
                                          subpool == 'subpool5')$ave_med_25,
                                   filter(sum_med_var, 
                                          subpool == 'subpool5')$ave_sum_25,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = -1, y = 1.25, color = '#33638DFF', 
           label = paste('r =', 
                         round(cor(filter(sum_med_var, 
                                          subpool == 'subpool3')$ave_med_25,
                                   filter(sum_med_var, 
                                          subpool == 'subpool3')$ave_sum_25,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = -1, y = 0.75, color = '#29AF7FFF', 
           label = paste('r =', 
                         round(cor(filter(sum_med_var, 
                                          subpool == 'subpool4')$ave_med_25,
                                   filter(sum_med_var, 
                                          subpool == 'subpool4')$ave_sum_25,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2))) +
  annotate("text", x = -1, y = 0.25, color = '#DCE319FF', 
           label = paste('r =', 
                         round(cor(filter(sum_med_var, 
                                          subpool == 'subpool2')$ave_med_25,
                                   filter(sum_med_var, 
                                          subpool == 'subpool2')$ave_sum_25,
                                   use = "pairwise.complete.obs", 
                                   method = "pearson"),
                               2)))

p_med_vs_sum_0_25_rep <- plot_grid(p_med_vs_sum_0_rep, p_med_vs_sum_25_rep, 
                                   nrow = 2, scale = 0.9, 
                                   labels = c(' 0 µM', '25 µM'),
                                   align = 'v', hjust = -3, vjust = -1)

save_plot('plots/p_med_vs_sum_0_25_rep_read3_nomed0.png', p_med_vs_sum_0_25_rep,
          base_width = 10.5, base_height = 7)


p_med_rep_0 <- ggplot(data = NULL, aes(med_ratio_0A, med_ratio_0B)) +
  facet_grid(. ~ subpool) +
  geom_point(data = filter(med_rep_1_2_0corr, subpool == 'subpool5'),
             color = '#440154FF', alpha = 0.3) +
  geom_point(data = filter(med_rep_1_2_0corr, subpool == 'subpool3'),
             color = '#33638DFF', alpha = 0.3) +
  geom_point(data = filter(med_rep_1_2_0corr, subpool == 'subpool4'),
             color = '#29AF7FFF', alpha = 0.3) +
  geom_point(data = filter(med_rep_1_2_0corr, subpool == 'subpool2'),
             color = '#DCE319FF', alpha = 0.3) +
  geom_point(data = filter(med_rep_1_2_0corr, 
                           mad_over_med_0A == as.double(1) & subpool != 'control'),
             color = 'red', alpha = 0.7) +
  geom_point(data = filter(med_rep_1_2_0corr, 
                           mad_over_med_0B == as.double(1) & subpool != 'control'),
             color = 'red', alpha = 0.7) +
  xlab("Median RNA/DNA 0A") +
  ylab("Median RNA/DNA 0B") +
  scale_x_log10(limits = c(0.01, 32)) + 
  scale_y_log10(limits = c(0.01, 32)) +
  annotation_logticks(sides = 'bl') +
  background_grid(major = 'xy', minor = 'none') + 
  panel_border() +
  annotate("text", x = 0.05, y = 10.1, color = '#440154FF',
           label = paste('r =', 
                         round(
                           cor(
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool5')$med_ratio_0A),
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool5')$med_ratio_0B),
                             use = "pairwise.complete.obs", method = "pearson"), 
                           2))) +
  annotate("text", x = 0.05, y = 5, color = '#33638DFF',
           label = paste('r =', 
                         round(
                           cor(
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool3')$med_ratio_0A),
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool3')$med_ratio_0B),
                             use = "pairwise.complete.obs", method = "pearson"), 
                           2))) +
  annotate("text", x = 0.05, y = 2.5, color = '#29AF7FFF',
           label = paste('r =', 
                         round(
                           cor(
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool4')$med_ratio_0A),
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool4')$med_ratio_0B),
                             use = "pairwise.complete.obs", method = "pearson"), 
                           2))) +
  annotate("text", x = 0.05, y = 1.25, color = '#DCE319FF',
           label = paste('r =', 
                         round(
                           cor(
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool2')$med_ratio_0A),
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool2')$med_ratio_0B),
                             use = "pairwise.complete.obs", method = "pearson"), 
                           2)))

p_med_rep_25 <- ggplot(filter(med_rep_1_2_0corr, subpool != 'control'), 
                      aes(med_ratio_25A, med_ratio_25B)) +
  facet_grid(. ~ subpool) +
  geom_point(data = filter(med_rep_1_2_0corr, subpool == 'subpool5'),
             color = '#440154FF', alpha = 0.3) +
  geom_point(data = filter(med_rep_1_2_0corr, subpool == 'subpool3'),
             color = '#33638DFF', alpha = 0.3) +
  geom_point(data = filter(med_rep_1_2_0corr, subpool == 'subpool4'),
             color = '#29AF7FFF', alpha = 0.3) +
  geom_point(data = filter(med_rep_1_2_0corr, subpool == 'subpool2'),
             color = '#DCE319FF', alpha = 0.3) +
  geom_point(data = filter(med_rep_1_2_0corr, 
                           mad_over_med_25A == as.double(1) & subpool != 'control'),
             color = 'red', alpha = 0.7) +
  geom_point(data = filter(med_rep_1_2_0corr, 
                           mad_over_med_25B == as.double(1) & subpool != 'control'),
             color = 'red', alpha = 0.7) +
  xlab("Median RNA/DNA 25A") +
  ylab("Median RNA/DNA 25B") +
  scale_x_log10(limits = c(0.01, 32)) + 
  scale_y_log10(limits = c(0.01, 32)) +
  annotation_logticks(sides = 'bl') +
  background_grid(major = 'xy', minor = 'none') + 
  panel_border() +
  annotate("text", x = 0.05, y = 10.1, color = '#440154FF',
           label = paste('r =', 
                         round(
                           cor(
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool5')$med_ratio_25A),
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool5')$med_ratio_25B),
                             use = "pairwise.complete.obs", method = "pearson"), 
                           2))) +
  annotate("text", x = 0.05, y = 5, color = '#33638DFF',
           label = paste('r =', 
                         round(
                           cor(
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool3')$med_ratio_25A),
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool3')$med_ratio_25B),
                             use = "pairwise.complete.obs", method = "pearson"), 
                           2))) +
  annotate("text", x = 0.05, y = 2.5, color = '#29AF7FFF',
           label = paste('r =', 
                         round(
                           cor(
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool4')$med_ratio_25A),
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool4')$med_ratio_25B),
                             use = "pairwise.complete.obs", method = "pearson"), 
                           2))) +
  annotate("text", x = 0.05, y = 1.25, color = '#DCE319FF',
           label = paste('r =', 
                         round(
                           cor(
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool2')$med_ratio_25A),
                             log10(
                               filter(med_rep_1_2_0corr, 
                                      subpool == 'subpool2')$med_ratio_25B),
                             use = "pairwise.complete.obs", method = "pearson"), 
                           2)))

p_rep_med_0_25 <- plot_grid(p_med_rep_0, p_med_rep_25, nrow = 2, 
                                  scale = 0.9, labels = c(' 0 µM', '25 µM'),
                                  align = 'v', hjust = -3, vjust = -1)

save_plot('plots/p_rep_med_0_25_read3.png', p_rep_med_0_25,
          base_width = 10.5, base_height = 7)


#Normalize to background

med_back_norm_1_ind_2 <- function(df1) {
  gsub_1_2 <- df1 %>%
    ungroup () %>%
    filter(subpool != 'control') %>%
    mutate(
      name = gsub('Smith R. Vista chr9:83712599-83712766', 'v chr9', name),
      name = gsub('Vista Chr5:88673410-88674494', 'v chr5', name),
      name = gsub('scramble pGL4.29 Promega 1-63 \\+ 1-87', 's pGl4', name)
    ) %>%
    mutate(background = name) %>%
    mutate(background = str_sub(background, 
                                nchar(background)-5, 
                                nchar(background)))
  backgrounds <- gsub_1_2 %>%
    filter(startsWith(name, 
                      'subpool5_no_site_no_site_no_site_no_site_no_site_no_site')) %>%
    select(background, med_ratio_0A, med_ratio_0B, 
           med_ratio_25A, med_ratio_25B) %>%
    rename(med_ratio_0A_back = med_ratio_0A) %>%
    rename(med_ratio_0B_back = med_ratio_0B) %>%
    rename(med_ratio_25A_back = med_ratio_25A) %>%
    rename(med_ratio_25B_back = med_ratio_25B)
  back_join_norm <- left_join(gsub_1_2, backgrounds, by = 'background') %>%
    mutate(med_ratio_0A_norm = med_ratio_0A/med_ratio_0A_back) %>%
    mutate(med_ratio_0B_norm = med_ratio_0B/med_ratio_0B_back) %>%
    mutate(med_ratio_25A_norm = med_ratio_25A/med_ratio_25A_back) %>%
    mutate(med_ratio_25B_norm = med_ratio_25B/med_ratio_25B_back) %>%
    mutate(ave_med_ratio_0_norm = (med_ratio_0A_norm + med_ratio_0B_norm)/2) %>%
    mutate(ave_med_ratio_25_norm = (med_ratio_25A_norm + med_ratio_25B_norm)/2) %>%
    mutate(induction = ave_med_ratio_25_norm/ave_med_ratio_0_norm)
}

med_rep_1_2_back_norm <- med_back_norm_1_ind_2(med_rep_1_2)

log2_med_rep_1_2_back_norm <- var_log2(med_rep_1_2_back_norm)

log10_med_rep_1_2_back_norm <- var_log10(med_rep_1_2_back_norm)


#Determining subpool bias-------------------------------------------------------

#Compare expression of sequences shared between subpool 2 and 5 at both 0 and 25
#µM forksolin. This analysis was done on data normalized to background
#expression levels to show the disparity between relative subpool expression
#that may result from the DNA ratios not representing biological ratios between
#the subpools. For this use data normalized to background and not log
#transformed.

compare_2_5_induction <- inner_join(subpool2, subpool5, 
                          by = c('most_common', 'background', 'ratio_0A_back', 
                                 'ratio_0B_back', 'ratio_25A_back', 
                                 'ratio_25B_back'),
                          suffix = c('_sp2', '_sp5'))

p_compare_back_norm_2_5_induction <- ggplot(compare_2_5_induction, 
                                    aes(induction_sp2, induction_sp5,
                                        fill = background)) +
  geom_point(shape = 21) +
  geom_abline(slope = 1, intercept = 0) +
  geom_hline(yintercept = 1, alpha = 0.5) +
  geom_vline(xintercept = 1, alpha = 0.5) +
  scale_fill_viridis(discrete = TRUE) +
  xlab("Subpool2 induction of\nbackground-normalized reads") +
  ylab("Subpool5 induction of\nbackground-normalized reads")

save_plot('plots/p_compare_back_norm_2_5_induction.png', 
          p_compare_back_norm_2_5_induction, base_width = 4.3, base_height = 3,
          scale = 1.2)

#compare_2_5_ratio

compare_2_5_ratio <- inner_join(subpool2, subpool5, 
                                by = c('most_common', 'background', 
                                       'ratio_0A_back', 'ratio_0B_back', 
                                       'ratio_25A_back', 'ratio_25B_back'),
                                suffix = c('_sp2', '_sp5'))

p_compare_back_norm_2_5_0 <- ggplot(compare_2_5_ratio, 
                                    aes(ave_ratio_0_norm_sp2, 
                                        ave_ratio_0_norm_sp5,
                                        fill = background)) +
  geom_point(shape = 21) +
  scale_fill_viridis(discrete = TRUE) +
  geom_hline(yintercept = 1, alpha = 0.5) +
  geom_vline(xintercept = 1, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  annotation_logticks(sides = 'bl') +
  xlab("Subpool2 average\nbackground-norm. sum RNA/DNA") +
  ylab("Subpool5 average\nbackground-norm. sum RNA/DNA") +
  scale_x_continuous(breaks = c(0:2), limits = c(0, 2.5)) + 
  scale_y_continuous(breaks = c(0:2), limits = c(0, 2.5))

p_compare_back_norm_2_5_25 <- ggplot(compare_2_5_ratio, 
                                     aes(ave_ratio_25_norm_sp2, 
                                         ave_ratio_25_norm_sp5,
                                         fill = background)) +
  geom_point(shape = 21) +
  scale_fill_viridis(discrete = TRUE) +
  geom_hline(yintercept = 1, alpha = 0.5) +
  geom_vline(xintercept = 1, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0) +
  annotation_logticks(sides = 'bl') +
  xlab("Subpool2 average\nbackground-norm. sum RNA/DNA") +
  ylab("Subpool5 average\nbackground-norm. sum RNA/DNA") +
  scale_x_continuous(breaks = c(0:2), limits = c(0, 2.5)) + 
  scale_y_continuous(breaks = c(0:2), limits = c(0, 2.5))

p_compare_back_norm_2_5 <- plot_grid(p_compare_back_norm_2_5_0, 
                                     p_compare_back_norm_2_5_25, 
                                     nrow = 2, labels = c(' 0 µM', '25 µM'),
                                     align = 'v', hjust = -2.5)

save_plot('plots/p_compare_back_norm_2_5.png', p_compare_back_norm_2_5,
          scale = 1.1, base_height = 7, base_width = 4.8)

data_0 <- tibble(ave_ratio_0_norm_sp2 = 0:3, ave_ratio_0_norm_sp5 = 0:3)
data_25 <- tibble(ave_ratio_25_norm_sp2 = 0:3, ave_ratio_25_norm_sp5 = 0:3)

mod_0 <- lm(ave_ratio_0_norm_sp5 ~ ave_ratio_0_norm_sp2, data_0)
mod_25 <- lm(ave_ratio_25_norm_sp5 ~ ave_ratio_25_norm_sp2, data_25)

compare_2_5_res <- compare_2_5 %>%
  add_residuals(mod_0) %>%
  rename(resid_0 = resid) %>%
  add_residuals(mod_25) %>%
  rename(resid_25 = resid)

p_compare_back_norm_2_5_0_res <- ggplot(compare_2_5_res, aes(x = most_common)) +
  geom_point(aes(y = resid_0), color = 'black') +
  geom_point(aes(y = resid_25), color = 'red')


#Selecting variants for luciferase assay--------------------------------------

sep_5_log <- sep_5 %>%
  mutate(ave_log_ratio_25 = (log10(ratio_25_1) + log10(ratio_25_2))/2) %>%
  mutate(ave_log_ratio_0 = (log10(ratio_0_2) + log10(ratio_0_2))/2)

sep_5_luc <- sep_5_log %>%
  arrange(ave_log_ratio_25) %>%
  slice(c(1, 78, 100, 500, 650, 800, 1000, 1100, 1200, 1380, 1450, 1576, 1800, 
          1880, 1965, 2010, 2075, 2145, 2178, 2187))

sep_5_back <- sep_5_log %>%
  filter(nosite == 6)

p_luc_var_pick <- ggplot(NULL, aes(ave_log_ratio_0, ave_log_ratio_25)) +
  geom_point(data = sep_5_log, color = '#999999') +
  geom_point(data = sep_5_luc, color = 'black') +
  geom_point(data = sep_5_back, color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant Average\nlog10 Uninduced Expression") +
  ylab("Variant Average\nlog10 Induced Expression")  +
  scale_x_continuous(breaks = c(0, 1), limits = c(-0.2, 1.2)) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(0, 2.5))

save_plot('plots/luc_var_pick.png',
          p_luc_var_pick)

bc_median_list <- function(df1, df2) {
  selection_bc <- left_join(df1, df2)
  med <- df2 %>%
    group_by(name, subpool, most_common) %>%
    summarize(med_reads_25 = median(num_reads))
  selection_bc_med <- left_join(selection_bc, med)
  sum_selection_bc_med <- filter(
    selection_bc_med, 
    as.numeric(
      num_reads
    ) >= med_reads_25 - 3 & as.numeric(
      num_reads
    ) <= med_reads_25 + 3
  ) %>%
    ungroup() %>%
    select(-subpool, -barcodes_RNA_25_1, -barcodes_RNA_25_2, 
           -barcodes_RNA_0_1, -barcodes_RNA_0_2, -ratio_25_1, 
           -ratio_25_2, -ratio_0_1, -ratio_0_2, -ave_ratio_25, 
           -ave_ratio_0, -normalized)
  return(sum_selection_bc_med)
}

sep_5_back_bc <- bc_median_list(sep_5_back, bc_join_R25A)
sep_5_luc_bc <- bc_median_list(sep_5_luc, bc_join_R25A)

write_csv(sep_5_back_bc, 
          "sep_5_bc_list_med_background.csv")
write_csv(sep_5_luc_bc, 
          "sep_5_bc_list_med_luc.csv")


#Testing validity of library manipulations-------------------------------

#Counting normalized barcodes
test <- bc_join_R25A %>%
  filter(num_reads > 0) %>%
  filter(name == "control_kheradpour_10" |
           name == "control_kheradpour_15") %>%
  head(10) %>%
  group_by(subpool, name, most_common) %>%
  count(name, wt = normalized)

#Median BC summarization
test <- inner_join(bc_join_R25A, bc_join_DNA,
                   by = c("name", "subpool", "most_common", "barcode"),
                   suffix = c("_RNA", "_DNA")) %>%
  filter(name == "control_kheradpour_10" |
           name == "control_kheradpour_15") %>%
  head(7) %>%
  mutate(BC_expression = normalized_RNA / normalized_DNA) %>%
  group_by(name, subpool, most_common) %>%
  summarize(med_BC_expression = median(BC_expression))


#Random work--------------------------------------------------------------------

#Selecting high uninduced high induction for Rishi
arrange(ave_variant, desc(ave_induction)) %>%
  head(1) %>%
  ungroup %>%
  select(name, most_common, ave_ratio_0, ave_induction) %>%
  write.table(
    "highunind_highinduction.txt", 
    sep = '\t', row.names = FALSE)

reportcomp <- filter(log_rep_1_2, 
                     name == 'subpool5_weak_consensus_consensus_consensus_consensus_consensus_Vista Chr5:88673410-88674494' | name == 'pGL4.29 Promega 1-63 + 1-87' | name == 'pGL4.29 Promega 1-87')


#Current unfinished code------------------------------------------------

fourier <- sep_2 %>%
  filter(site == "consensus", 
         background == "scramble pGL4.29 Promega 1-63 + 1-87") %>%
  arrange(dist) %>%
  fft(fourier$ave_induction)

#Subpool 3, tested but did not find interesting displaying as a
#function of 2nd site or midpoint between sites

#determine distance according to second site, measuring from start of background,
#12 bp of the first site, spacing and 2 bp of flanks, for 0 spacing, just add 10
#to account for outer flank and first site
sep_3_2ndsite <-subpool3 %>% 
  mutate(name = gsub('2BS ', '', name), 
         name = gsub(' bp spacing ', '_', name)) %>%
  separate(name, into = c(
    "subpool", "spacing", "fluff2", "fluff3", "dist", "background"
  ), sep = "_", convert = TRUE
  ) %>%
  select(-fluff2, -fluff3) %>%
  mutate(background = gsub('Smith R. Vista chr9:83712599-83712766',
                           'vista chr9', background),
         background = gsub('Vista Chr5:88673410-88674494', 'vista chr5',
                           background),
         background = str_sub(background, 1, 13)) %>%
  mutate(dist = ifelse(
    spacing != as.integer(0), dist + 12 + spacing + 2, dist + 10)
  ) %>%
  mutate(spacing = 
           ifelse(spacing != as.integer(0), 
                  spacing + 4, spacing))

#determine distance as a function of midpoint between the two spaces
sep_3_mid <-subpool3 %>% 
  mutate(name = gsub('2BS ', '', name), 
         name = gsub(' bp spacing ', '_', name)) %>%
  separate(name, into = c(
    "subpool", "spacing", "fluff2", "fluff3", "dist", "background"
  ), sep = "_", convert = TRUE
  ) %>%
  select(-fluff2, -fluff3) %>%
  mutate(background = gsub('Smith R. Vista chr9:83712599-83712766',
                           'vista chr9', background),
         background = gsub('Vista Chr5:88673410-88674494', 'vista chr5',
                           background),
         background = str_sub(background, 1, 13)) %>%
  mutate(dist = ifelse(
    spacing != as.integer(0), dist + 12 + (spacing/2), dist + 10)
  ) %>%
  mutate(spacing = 
           ifelse(spacing != as.integer(0), 
                  spacing + 4, spacing))

#Subpool 4 - plot number of sites, might take easy route and just add
#extra columns of the consensus consensus expression

test_set_4 <- sep_4 %>%
  filter(background == 'vista chr5', dist < 3)

#try to spread instead

#normalize_to_consensus <- function(df) {
#if (site1 == 'consensus' & site2 == 'consensus') {
#consensus_ref = ave_induction}
#mutate(ave_induction_over_consensus = ave_induction/consensus_ref)
#return(df)
#}

test_set_4 <- test_set_4 %>%
  group_by(background, spacing, dist) %>%
  do(normalize_to_consensus(.)) %>%
  ungroup

test_set_4 <- sep_4 %>%
  group_by(background, dist, spacing) %>%
  mutate(c_c_ave_induction = )

ggplot(NULL, aes(dist, ave_induction)) +
  geom_smooth(data = filter(sep_3, spacing != 70), 
              span = 0.1, size = 0.7, color = '#0072B2',
              show.legend = TRUE) +
  geom_point(data = filter(sep_4, site1 == 'consensus', 
                           site2 == 'consensus'), 
             alpha = 0.6, size = 2, color = '#D55E00', 
             show.legend = TRUE) +
  geom_point(data = filter(sep_4, site1 == 'moderate', 
                           site2 == 'consensus'), 
             alpha = 0.6, size = 2, color = 'black', 
             show.legend = TRUE) +
  geom_point(data = filter(sep_4, site1 == 'moderate', 
                           site2 == 'moderate'), 
             alpha = 0.6, size = 2, color = 'red', 
             show.legend = TRUE) +
  facet_grid(spacing ~ background) +
  ylab('Average Induction Ratio') + 
  ggtitle('Effect of Distance of Two Consensus Sites along 3 Backgrounds upon Induction') +
  scale_x_continuous("Distance from First Site to Proximal Promoter End (bp)", 
                     breaks = seq(from = 0, to = 150, by = 10))

#Subpool5
#play around with this fit model
fit5 <- lm(ave_induction ~ site1 + site2 + site3 + site4 + site5 + site6 +
             background, data=sep_5)
plot(fit) +
  layout(matrix(c(1,2,3,4),2,2)) # optional 4 graphs/page 
anova(fit)

#plot per singular site, but this is subpool2 as well
consensus_1_sp5 <- sep_5 %>%
  filter(consensus == 1, weak == 0) %>%
  mutate(site_fill = ifelse(site1 = consensus, 1, 0))
ggplot(consensus_1_sp5, aes())

#Plot change in induction from consensus to weak or no_site starting at
#all 6 sites being filled with consensus. Do this per background
library('reshape2')

sep_5_pos_type <- sep_5 %>%
  filter(consensus >= 5) %>%
  mutate(position = 0) %>%
  mutate(type = 'consensus') %>%
  mutate(position = ifelse(site1 != 'consensus', 1, position)) %>%
  mutate(position = ifelse(site2 != 'consensus', 2, position)) %>%
  mutate(position = ifelse(site3 != 'consensus', 3, position)) %>%
  mutate(position = ifelse(site4 != 'consensus', 4, position)) %>%
  mutate(position = ifelse(site5 != 'consensus', 5, position)) %>%
  mutate(position = ifelse(site6 != 'consensus', 6, position)) %>%
  mutate(type = ifelse(weak == 1, 'weak', type)) %>%
  mutate(type = ifelse(nosite == 1, 'none', type))

sep_5_pos_type_scr <- sep_5_pos_type %>%
  filter(background == 'scramble pGL4') %>%
  dcast(position ~ type, value.var = 'ave_induction') %>%
  mutate(weak_diff = weak - first(consensus)) %>%
  mutate(nosite_diff = none - first(consensus)) %>%
  select(position, weak_diff, nosite_diff) %>%
  mutate(background = 'scramble pGL4')

p_cons_diff_25_scr4 <- ggplot(filter(sep_5_pos_type_scr, position != 0)) +
  geom_point(aes(position, weak_diff), color = '#56B4E9') +
  geom_point(aes(position, nosite_diff), color = 'gray') +
  geom_hline(yintercept = 0) + 
  ylab('∆ Induction\nfrom site change') +
  ggtitle('Scramble pGL4')

sep_5_pos_type_chr5 <- sep_5_pos_type %>%
  filter(background == 'vista chr5') %>%
  dcast(position ~ type, value.var = 'ave_induction') %>%
  mutate(weak_diff = weak - first(consensus)) %>%
  mutate(nosite_diff = none - first(consensus)) %>%
  select(position, weak_diff, nosite_diff) %>%
  mutate(background = 'vista chr5')

p_cons_diff_25_chr5 <- ggplot(filter(sep_5_pos_type_chr5, position != 0)) +
  geom_point(aes(position, weak_diff), color = '#56B4E9') +
  geom_point(aes(position, nosite_diff), color = 'gray') +
  geom_hline(yintercept = 0) + 
  ylab('∆ Induction\nfrom site change') +
  ggtitle('Vista Chr5')

sep_5_pos_type_chr9 <- sep_5_pos_type %>%
  filter(background == 'vista chr9') %>%
  dcast(position ~ type, value.var = 'ave_induction') %>%
  mutate(weak_diff = weak - first(consensus)) %>%
  mutate(nosite_diff = none - first(consensus)) %>%
  select(position, weak_diff, nosite_diff) %>%
  mutate(background = 'vista chr9')

p_cons_diff_25_chr9 <- ggplot(filter(sep_5_pos_type_chr9, position != 0)) +
  geom_point(aes(position, weak_diff), color = '#56B4E9') +
  geom_point(aes(position, nosite_diff), color = 'gray') +
  geom_hline(yintercept = 0) + 
  ylab('∆ Induction\nfrom site change') +
  ggtitle('Vista Chr9')

p_cond_diff_25_all <- plot_grid(p_cons_diff_25_chr5, 
                                p_cons_diff_25_scr4, ncol = 2)

save_plot('plots/cond_diff_25_all.png',
          p_cond_diff_25_all, scale = 0.3, 
          base_width = 28, base_height = 10)

#sample a certain amount to make plotting easier
p_inducedBCrep <- bc_rep_1_2 %>%
  sample_n(10000) %>%
  ggplot(aes(ratio_25_1, ratio_25_2)) + 
  geom_point(alpha = 0.1) + 
  scale_x_log10(limits = c(0.02,1000)) + 
  scale_y_log10(limits = c(0.02,1000)) +
  annotation_logticks() +
  xlab("Barcode Expression Rep. 1") +
  ylab("Barcode Expression Rep. 2") +
  annotate("text", x = 0.1, y = 100,
           label = paste('r =', round(bc_rep_1_2_ratio_25_R2, 2)))
