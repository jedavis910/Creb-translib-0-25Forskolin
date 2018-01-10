library(tidyverse)
library(stringr)
library(viridis)
library(cowplot)
library(ggExtra)
library(modelr)
library(lazyeval)
library(splines)

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
                          'fluff', 'barcode', 'name', 'most_common'
                          ), 
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
  keep_bc <- left_join(df1, df2, by = 'barcode') %>%
    mutate(num_reads = if_else(is.na(num_reads), 
                               as.integer(0), 
                               num_reads)) %>%
    mutate(normalized = as.numeric((num_reads * 1000000) / (sum(num_reads))))
  return(keep_bc)
}

bc_join_R25A <- bc_map_join_bc(barcode_map, bc_R25A)
bc_join_R25B <- bc_map_join_bc(barcode_map, bc_R25B)
bc_join_R0A <- bc_map_join_bc(barcode_map, bc_R0A)
bc_join_R0B <- bc_map_join_bc(barcode_map, bc_R0B)
bc_join_DNA <- bc_map_join_bc(barcode_map, bc_DNA)


#Determine variant counts by summing--------------------------------------------

#sum unique barcodes and normalized bc reads across barcodes per variant output 
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
                       by = c("name", "subpool", "most_common")
                       ) %>%
    ungroup() %>%
    filter(barcodes > 7)
  return(bc_sum)
}

variant_counts_R25A <- var_sum_bc_num(bc_join_R25A)
variant_counts_R25B <- var_sum_bc_num(bc_join_R25B)
variant_counts_R0A <- var_sum_bc_num(bc_join_R0A)
variant_counts_R0B <- var_sum_bc_num(bc_join_R0B)
variant_counts_DNA <- var_sum_bc_num(bc_join_DNA)


#Join RNA to DNA and determine expression from summing--------------------------

#combine DNA and RNA cumm. BC counts, only keeping instances in both sets and 
#determining RNA/DNA per variant. Ratio is summed normalized reads of RNA over 
#summed normalized reads DNA

var_expression <- function(df1, df2) {
  RNA_DNA <- inner_join(df1, df2, 
                        by = c("name", "subpool", "most_common"), 
                        suffix = c("_RNA", "_DNA")
                        ) %>%
    mutate(ratio = sum_RNA / sum_DNA)
  print('x defined as RNA, y defined as DNA in var_expression(x,y)')
  return(RNA_DNA)
}

RNA_DNA_R25A <- var_expression(variant_counts_R25A, variant_counts_DNA)
RNA_DNA_R25B <- var_expression(variant_counts_R25B, variant_counts_DNA)
RNA_DNA_R0A <- var_expression(variant_counts_R0A, variant_counts_DNA)
RNA_DNA_R0B <- var_expression(variant_counts_R0B, variant_counts_DNA)


#combine biological replicates--------------------------------------------------

var_rep <- function(df0A, df0B, df25A, df25B) {
  join_0 <- inner_join(df0A, df0B, 
                       by = c(
                         "name", "subpool", "most_common", "sum_DNA", 
                         "barcodes_DNA"), 
                       suffix = c("_0A", "_0B")
                       )
  join_25 <- inner_join(df25A, df25B, 
                       by = c(
                         "name", "subpool", "most_common", "sum_DNA", 
                         "barcodes_DNA"), 
                       suffix = c("_25A", "_25B")
                       )
  join_0_25 <- inner_join(join_0, join_25,
                        by = c("name", "subpool", "most_common", "sum_DNA", 
                               "barcodes_DNA")
                        )
  print('processed dfs in order: df0A, df0B, df25A, df25B')
  return(join_0_25)
}

rep_1_2 <- var_rep(RNA_DNA_R0A, RNA_DNA_R0B, RNA_DNA_R25A, RNA_DNA_R25B)

#Normalize variants to background-----------------------------------------------

#After combining, rename backgrounds to simplified names, make background column 
#(excluding controls), separate out background values in each dataset and left 
#join to original dataset. Normalize expression of each variant to its 
#background in that biological replicate. Determine average normalized 
#expression across biological replicates.

#Add text here





var_log <- function(df) {
  log_ratio <- df %>%
    mutate(log_ratio_25_1 = log10(ratio_25_1)) %>%
    mutate(log_ratio_25_2 = log10(ratio_25_2)) %>%
    mutate(log_ratio_0_1 = log10(ratio_0_1)) %>%
    mutate(log_ratio_0_2 = log10(ratio_0_2))
  return(log_ratio)
}

log_rep_1_2 <- var_log(rep_1_2)
  

#Subpool manipulations using expression from sum-----------------------

#Separate average variant values into subpools
subpool2 <- 
  filter(rep_1_2, subpool == "subpool2") %>%
  ungroup () %>%
  select(-subpool)
subpool3 <- 
  filter(rep_1_2, subpool == "subpool3") %>%
  ungroup () %>%
  select(-subpool)
subpool4 <- 
  filter(rep_1_2, subpool == "subpool4") %>%
  ungroup () %>%
  select(-subpool)
subpool5 <- 
  filter(rep_1_2, subpool == "subpool5") %>%
  ungroup () %>%
  select(-subpool)
controls <- 
  filter(rep_1_2, subpool == "control") %>%
  ungroup ()

#Separate string qualifiers per subpool

#Subpool 2 contains either a consensus site (TGACGTCA) or consensus 
#surrounded by 2 bp flanks on either side (ATTGACGTCAGC) as is used in the
#rest of the subpools. Each site is placed on the background starting 
#closest to minP and are then moved along the backgrounds at 1 bp 
#increments. Separation lists the type of site, the distance (start of 
#consensus or consensus and flanks) and the background used. Added 2 bp to
#consensusflank so that the start of the binding site is represented
#instead of the start of the flank
sep_2 <- separate(
  subpool2, name, into = c(
    "subpool", "site", "fluff2", "dist", "background"
    ), sep = "_", convert = TRUE
  ) %>% 
  select(-fluff2) %>%
  mutate(background = gsub('Smith R. Vista chr9:83712599-83712766',
                           'vista chr9', background),
         background = gsub('Vista Chr5:88673410-88674494', 'vista chr5',
                           background),
         background = str_sub(background, 1, 13)) %>%
  mutate(dist = 
           ifelse(startsWith(site, 'consensusflank'), 
                  dist + 2, dist))

#Subpool 3 contains 2 consensus binding sites with flanks (ATTGACGTCAGC) 
#that vary in distance from one another by 0 (no inner flanks), 5, 10, 15,
#20 and 70 bp (all but 0 appear as -4 bp spacing). Each site distance 
#combination is then moved along the backgrounds at 1 bp increments 
#starting from closest to the minP. Separation lists the spacing between 
#sites, distance (start of consensus and flanks) and the background. Added
#2 to all distances to measure to start of BS and not to flank. Added 4 
#to all spacing but 0 to measure difference between start of sites.
sep_3 <-subpool3 %>% 
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
  mutate(dist = dist + 2) %>%
  mutate(spacing = 
           ifelse(spacing != as.integer(0), 
                  spacing + 4, spacing))


#Subpool 4 contains 2 binding sites with flanks that vary in site type 
#from consensus (ATTGACGTCAGC) moderate (ATTGACGTCTGC) weak (ATTGAAGTCAGC) 
#and half(NNNGCCGTCATA). These sites vary in spacing from 0, 5, 10, 15, 
#and 20 bp. Both sites are moved along the promoter at 5 bp increments 
#starting closest to minP. The order of sites in the name is from further
#to closer to the minP and distance is from the closer site to 0. 
#Separation lists the first and second site, the spacing between sites, 
#distance (start of site) and the background.
sep_4 <- subpool4 %>%
  mutate(name = gsub('consensus0 ', 'consensus_0', name),
         name = gsub('weak0 ', 'weak_0', name),
         name = gsub('moderate0 ', 'moderate_0', name),
         name = gsub('half0 ', 'half_0', name),
         name = gsub('bp spacing', '', name)) %>%
  separate(name, into = c(
    "subpool", "site1", "site2", "spacing", "fluff2", "dist", "background"
  ), sep = "_", convert = TRUE
  ) %>%
  select(-fluff2) %>%
  mutate(consensus = str_detect(site1, "consensus") + 
           str_detect(site2, "consensus")) %>%
  mutate(moderate = str_detect(site1, "moderate") + 
           str_detect(site2, "moderate")) %>%
  mutate(weak = str_detect(site1, "weak") + 
     str_detect(site2, "weak")) %>%
  mutate(half = str_detect(site1, "half") + 
           str_detect(site2, "half")) %>%
  mutate(background = gsub('Smith R. Vista chr9:83712599-83712766',
                           'vista chr9', background),
         background = gsub('Vista Chr5:88673410-88674494', 'vista chr5',
                           background),
         background = str_sub(background, 1, 13)) %>%
  mutate(dist = dist + 2) %>%
  mutate(spacing = 
           ifelse(spacing != as.integer(0), 
                  spacing + 4, spacing))

#Subpool 5 contains 6 equally spaced sites spaced 13 bp apart and starting
#from furthest to the minP. These sites are filled with sites of either the
#consensus site, a weak site or no site. Both the weak and consensus sites
#are flanked by the same flanking sequence.
sep_5 <- subpool5 %>%
  mutate(name = gsub('no_site', 'nosite', name)) %>%
  separate(name, into = c(
    "subpool", "site1", "site2", "site3", "site4", "site5", "site6", 
    "background"), sep = "_"
  ) %>%
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
  mutate(background = gsub('Smith R. Vista chr9:83712599-83712766',
                           'vista chr9', background),
         background = gsub('Vista Chr5:88673410-88674494', 'vista chr5',
                           background),
         background = str_sub(background, 1, 13))

ave_site_5 <- sep_5 %>%
  group_by(background, consensus, weak) %>%
  summarize(ave_site_induction = mean(ave_induction), 
            ave_site_25 = mean(ave_ratio_25),
            ave_site_0 = mean(ave_ratio_0))

#model fitting to subpool 5
sep_5_chr9 <- sep_5 %>%
  filter(background == 'vista chr9')
  

#Replicate plots from summing-------------------------------------------

#Determine pearson and plot a lot of replicate graphs
rep_1_2_ratio_0_R2 <- cor(log_rep_1_2$log_ratio_0_1 ,
                              log_rep_1_2$log_ratio_0_2,
                              use = "pairwise.complete.obs",
                              method = "pearson")

rep_1_2_ratio_25_R2 <- cor(log_rep_1_2$log_ratio_25_1 ,
                            log_rep_1_2$log_ratio_25_2,
                            use = "pairwise.complete.obs",
                            method = "pearson")

rep_1_2_25_BCnum_R2 <- cor(rep_1_2$barcodes_RNA_25_1, 
                         rep_1_2$barcodes_RNA_25_2,
                         use = "pairwise.complete.obs",
                         method = "pearson")

rep_1_2_0_BCnum_R2 <- cor(rep_1_2$barcodes_RNA_0_1, 
                           rep_1_2$barcodes_RNA_0_2,
                           use = "pairwise.complete.obs",
                           method = "pearson")

p_0_var_rep <- ggplot(NULL, aes(log_ratio_0_1, 
                                           log_ratio_0_2)) +
  geom_point(data = log_rep_1_2, alpha = 0.3) +
  geom_point(data = filter(log_rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  annotation_logticks(scaled = TRUE) +
  xlab("Variant log10 Expression Rep. 1") +
  ylab("Variant log10 Expression Rep. 2") +
  scale_x_continuous(breaks = c(0, 1), limits = c(-0.6, 2.4)) + 
  scale_y_continuous(breaks = c(0, 1), limits = c(-0.6, 2.4)) + 
  annotate("text", x = 0, y = 1,
           label = paste('r =', round(rep_1_2_ratio_0_R2, 2)))

save_plot('plots/0_var_rep_neg.png',
          base_aspect_ratio = 1, p_0_var_rep)

p_25_var_rep <- ggplot(NULL, 
                        aes(log_ratio_25_1, log_ratio_25_2)) + 
  geom_point(data = log_rep_1_2, alpha = 0.3) +
  geom_point(data = filter(log_rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(-0.6, 2.4)) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-0.6, 2.4)) +
  annotation_logticks() +
  xlab("Variant log10 Expression Rep. 1") + 
  ylab("Variant log10 Expression Rep. 2") +
  annotate("text", x = 0, y = 1,
           label = paste('r =', round(rep_1_2_ratio_25_R2, 2)))

save_plot('plots/ind_var_rep_neg.png',
          base_aspect_ratio = 1, p_25_var_rep)

p_0_bc_var <- ggplot(ave_variant, aes(ave_barcodes_0, 
                                          ave_ratio_0)) +
  geom_point(alpha = 0.1) + 
  xlab("Average Number of Barcodes per Variant") +
  ylab("Average Expression per Variant") + xlim(0, 702) +
  facet_wrap( ~ subpool, nrow = 5) + panel_border()

save_plot('plots/0_bc_var.png',
          p_0_bc_var, scale = 1.7)

p_25_bc_var <- ggplot(ave_variant, aes(ave_barcodes_25, 
                                        ave_ratio_25)) +
  geom_point(alpha = 0.1) +
  xlab("Average Number of Barcodes per Variant") +
  ylab("Average Expression per Variant") + xlim(0, 702) +
  facet_wrap( ~ subpool, nrow = 5) + panel_border()

save_plot('plots/ind_bc_var.png',
          p_25_bc_var, scale = 1.7)

p_bc_rep_0 <- ggplot(
  rep_1_2, aes(barcodes_RNA_0_1, barcodes_RNA_0_2)
  ) +
  geom_point(alpha = 0.1) +
  geom_density2d() +
  scale_x_continuous("Number of Barcodes per Variant Rep. 1", 
                     breaks = seq(from = 0, to = 700, by = 100)) +
  scale_y_continuous("Number of Barcodes per Variant Rep. 2", 
                     breaks = seq(from = 0, to = 700, by = 100)) +
  annotate("text", x = 100, y = 400,
           label = paste('r =', round(rep_1_2_0_BC_R2, 2)))

save_plot('plots/bc_rep_0.png',
          p_bc_rep_0)

p_bc_rep_25 <- ggplot(
  rep_1_2, aes(barcodes_RNA_25_1, barcodes_RNA_25_2)
  ) +
  geom_point(alpha = 0.1) +
  geom_density2d() +
  scale_x_continuous("Number of Barcodes per Variant Rep. 1", 
                     breaks = seq(from = 0, to = 700, by = 100), 
                     limits = c(0, 700)) +
  scale_y_continuous("Number of Barcodes per Variant Rep. 2",
                     breaks = seq(from = 0, to = 700, by = 100), 
                     limits = c(0, 700)) +
  annotate("text", x = 100, y = 400,
           label = paste('r =', round(rep_1_2_25_BC_R2, 2)))

save_plot('plots/bc_rep_25.png',
          p_bc_rep_25)

#Plot reads per BC, filtering out 0 reads that skew plot
sample_11 <- filter(bc_join_DNA, num_reads > 0)
sample_02 <- filter(bc_join_R25A, num_reads > 0)
sample_03 <- filter(bc_join_R25B, num_reads > 0)
sample_04 <- filter(bc_join_R0A, num_reads > 0)
sample_05 <- filter(bc_join_R0B, num_reads > 0)

p_BC_num_reads_viol_full <- ggplot(NULL, aes(x = "", y = num_reads)) +
  geom_violin(data = sample_11, aes(x = "DNA")) +
  geom_violin(data = sample_02, aes(x = "Ind_1")) + 
  geom_violin(data = sample_03, aes(x = "Ind_2")) +
  geom_violin(data = sample_04, aes(x = "0_1")) + 
  geom_violin(data = sample_05, aes(x = "0_2")) + 
  xlab("") +
  ylab("Reads")

save_plot('plots/BC_num_reads_viol_full.png',
          p_BC_num_reads_viol_full)

p_BC_num_reads_box_zoom <- ggplot(NULL, aes(x = "", y = num_reads)) +
  geom_boxplot(data = sample_11, aes(x = "DNA"), alpha = 0.1) +
  geom_boxplot(data = sample_02, aes(x = "Ind_1"), alpha = 0.1) + 
  geom_boxplot(data = sample_03, aes(x = "Ind_2"), alpha = 0.1) +
  geom_boxplot(data = sample_04, aes(x = "0_1"), alpha = 0.1) + 
  geom_boxplot(data = sample_05, aes(x = "0_2"), alpha = 0.1) + 
  xlab("") +
  ylab("Reads") + ylim(0, 25)

save_plot('plots/BC_num_reads_box_zoom.png',
          p_BC_num_reads_box_zoom)

#checking negatives relative to library expression
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
p_subpool2_dist_0_25 <- ggplot(filter(sep_2, 
                                           site == 'consensusflank')) + 
  facet_grid(~ background) + 
  geom_point(aes(dist, ratio_0_1), alpha = 0.5, size = 1.2,
             color = '#999999') +
  geom_point(aes(dist, ratio_0_2), alpha = 0.5, size = 1.2,
             color = '#999999') +
  geom_point(aes(dist, ratio_25_1), alpha = 0.5, size = 1.2,
             color = '#56B4E9') +
  geom_point(aes(dist, ratio_25_2), alpha = 0.5, size = 1.2,
             color = '#56B4E9') +
  geom_smooth(aes(dist, ave_ratio_0), span = 0.1, size = 0.4,
              se = FALSE, color = '#999999') +
  geom_smooth(aes(dist, ave_ratio_25), span = 0.1, size = 0.4,
              se = FALSE, color = '#56B4E9') +
  scale_x_continuous("Distance from Proximal Promoter End (bp)", 
                     breaks = seq(from = 0, to = 150, by = 10)
                     ) +
  panel_border() + ylab('Expression') +
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
p_subpool5_box_induction <- ggplot(sep_5, aes(as.character(consensus), 
                                       ave_induction)) + 
  facet_wrap(~ background) + 
  geom_boxplot() +
  panel_border() +
  xlab("Number of Consensus Binding Sites") +
  ylab("Average Induction Ratio")

save_plot('plots/subpool5_box_induction.png',
          p_subpool5_box_induction, scale = 1.2)

#Plot site number vs. induction (counting weak and consensus) and 
#overlay # of consensus sites
p_subpool5_cons_weak_sitenum_25 <- ggplot(NULL) + 
  facet_wrap(~ background) + 
  geom_point(data = sep_5,
               aes(total_sites, ave_ratio_25), 
             alpha = 0.3, color = '#999999') +
  geom_point(data = filter(sep_5, weak == 0 & consensus >0),
             aes(consensus, ave_ratio_25), 
             alpha = 1, color = 'black') +
  geom_point(data = filter(sep_5, consensus == 0 & weak >0),
             aes(weak, ave_ratio_25), 
             alpha = 1, color = '#56B4E9') +
  panel_border() +
  xlab("Total Number of Binding Sites") +
  scale_y_continuous("Expression at 25 µM Forskolin", limits = c(0,250))

save_plot('plots/subpool5_cons_weak_sitenum_25.png',
          p_subpool5_cons_weak_sitenum_25, scale = 1.2)

p_subpool5_cons_weak_sitenum_0 <- ggplot(NULL) + 
  facet_wrap(~ background) + 
  geom_point(data = sep_5,
             aes(total_sites, ave_ratio_0), 
             alpha = 0.3, color = '#999999') +
  geom_point(data = filter(sep_5, weak == 0 & consensus >0),
             aes(consensus, ave_ratio_0), 
             alpha = 1, color = 'black') +
  geom_point(data = filter(sep_5, consensus == 0 & weak >0),
             aes(weak, ave_ratio_0), 
             alpha = 1, color = '#56B4E9') +
  panel_border() +
  xlab("Total Number of Binding Sites") +
  scale_y_continuous("Expression")

save_plot('plots/subpool5_cons_weak_sitenum_0.png',
          p_subpool5_cons_weak_sitenum_0, scale = 1.2)

p_subpool5_weak_sitenum_induction <- ggplot(
  filter(sep_5, consensus == 0), aes(as.character(total_sites), 
                                     ave_induction)) + 
  geom_point(color = '#56B4E9') +
  facet_wrap(~ background) + 
  panel_border() +
  xlab("Number of Weak Sites Only") +
  ylab("Average Induction Ratio")

save_plot('plots/subpool5_weak_sitenum_induction.png',
          p_subpool5_weak_sitenum_induction, scale = 1.2)

p_subpool5_box_0_25 <- ggplot(sep_5) + 
  facet_wrap(~ background) + 
  geom_boxplot(aes(as.character(consensus), ave_ratio_0), 
               color = '#999999', show.legend = TRUE) +
  geom_boxplot(aes(as.character(consensus), ave_ratio_25), 
               color = '#56B4E9', show.legend = TRUE) +
  panel_border() + 
  xlab("Number of Consensus Binding Sites") +
  ylab("Average Expression")

save_plot('plots/subpool5_box_0_25.png',
          p_subpool5_box_0_25, scale = 1.2)

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

#Median calculations of Expression-----------------------------------------

#Re-do BC normalization to account for pseudocount of 1
pseudo_bc_join_R25A <- 
  left_join(barcode_map, bc_R25A, by = "barcode") %>%
  mutate(num_reads = 
           if_else(is.na(num_reads),
                   as.integer(1), as.integer(num_reads + 1))) %>%
  mutate(normalized = 
           as.numeric((num_reads * 1000000 * 5.2) / (sum(
             num_reads)))) %>%
  mutate(subpool = 
           ifelse(startsWith(name, 'subpool'), 
                  substr(name, 1, 8), 'control'))

pseudo_bc_join_R25B <- 
  left_join(barcode_map, bc_R25B, by = "barcode") %>%
  mutate(num_reads = 
           if_else(is.na(num_reads),
                   as.integer(1), as.integer(num_reads + 1))) %>%
  mutate(normalized = 
           as.numeric((num_reads * 1000000 * 5.2) / (sum(
             num_reads)))) %>%
  mutate(subpool = 
           ifelse(startsWith(name, 'subpool'), 
                  substr(name, 1, 8), 'control'))

pseudo_bc_join_R0A <- 
  left_join(barcode_map, bc_R0A, by = "barcode") %>%
  mutate(num_reads = 
           if_else(is.na(num_reads),
                   as.integer(1), as.integer(num_reads + 1))) %>%
  mutate(normalized = 
           as.numeric((num_reads * 1000000) / (sum(
             num_reads)))) %>%
  mutate(subpool = 
           ifelse(startsWith(name, 'subpool'),
                  substr(name, 1, 8), 'control'))

pseudo_bc_join_R0B <- 
  left_join(barcode_map, bc_R0B, by = "barcode") %>%
  mutate(num_reads = 
           if_else(is.na(num_reads),
                   as.integer(1), as.integer(num_reads + 1))) %>%
  mutate(normalized = 
           as.numeric((num_reads * 1000000) / (sum(
             num_reads)))) %>%
  mutate(subpool = 
           ifelse(startsWith(name, 'subpool'), 
                  substr(name, 1, 8), 'control'))

pseudo_bc_join_DNA <- 
  left_join(barcode_map, bc_DNA, by = "barcode") %>%
  mutate(num_reads = 
           if_else(is.na(num_reads),
                   as.integer(1), as.integer(num_reads + 1))) %>%
  mutate(normalized = 
           as.numeric((num_reads * 1000000) / (sum(
             num_reads)))) %>%
  mutate(subpool =
           ifelse(startsWith(name, 'subpool'),
                  substr(name, 1, 8), 'control'))

#Determine BC expression RNA/DNA then take median BC expression per
#variant, appending BC numbers
med_RNA_DNA_R25A <-
  inner_join(pseudo_bc_join_R25A, pseudo_bc_join_DNA,
             by = c("name", "subpool", "most_common", "barcode"),
             suffix = c("_RNA", "_DNA")) %>%
  mutate(BC_expression = normalized_RNA / normalized_DNA) %>%
  group_by(name, subpool, most_common) %>%
  summarize(med_BC_expression = median(BC_expression)) %>%
  left_join(bc_num_R25A, 
             by = c("name", "subpool", "most_common"))

med_RNA_DNA_R25B <-
  inner_join(pseudo_bc_join_R25B, pseudo_bc_join_DNA,
             by = c("name", "subpool", "most_common", "barcode"),
             suffix = c("_RNA", "_DNA")) %>%
  mutate(BC_expression = normalized_RNA / normalized_DNA) %>%
  group_by(name, subpool, most_common) %>%
  summarize(med_BC_expression = median(BC_expression)) %>%
  left_join(bc_num_R25B, 
            by = c("name", "subpool", "most_common"))

med_RNA_DNA_R0A <-
  inner_join(pseudo_bc_join_R0A, pseudo_bc_join_DNA,
             by = c("name", "subpool", "most_common", "barcode"),
             suffix = c("_RNA", "_DNA")) %>%
  mutate(BC_expression = normalized_RNA / normalized_DNA) %>%
  group_by(name, subpool, most_common) %>%
  summarize(med_BC_expression = median(BC_expression)) %>%
  left_join(bc_num_R0A, 
            by = c("name", "subpool", "most_common"))

med_RNA_DNA_R0B <-
  inner_join(pseudo_bc_join_R0B, pseudo_bc_join_DNA,
             by = c("name", "subpool", "most_common", "barcode"),
             suffix = c("_RNA", "_DNA")) %>%
  mutate(BC_expression = normalized_RNA / normalized_DNA) %>%
  group_by(name, subpool, most_common) %>%
  summarize(med_BC_expression = median(BC_expression)) %>%
  left_join(bc_num_R0B, 
            by = c("name", "subpool", "most_common"))

#Combine induced and uninduced per replicate
med_25_rep_1 <-
  inner_join(med_RNA_DNA_R25A, med_RNA_DNA_R0A,
             by = c("name", "subpool", "most_common"),
             suffix = c("_25", "_0"))

med_25_rep_2 <-
  inner_join(med_RNA_DNA_R25B, med_RNA_DNA_R0B,
             by = c("name", "subpool", "most_common"),
             suffix = c("_25", "_0"))

#Join both replicates and take average unique BC numbers, determine 
#correlation between replicates at uninduced and induced
med_rep_1_2 <- inner_join(med_25_rep_1, med_25_rep_2,
                          by = c("name", "subpool", "most_common"),
                          suffix = c("_1", "_2")) %>%
  mutate(ave_barcodes_25 = 
           (barcodes_25_1 + barcodes_25_2)/2) %>%
  mutate(ave_barcodes_0 = 
           (barcodes_0_1 + barcodes_0_2)/2) %>%
  left_join(bc_num_DNA, by = c("name", "subpool", "most_common"))

names(med_rep_1_2)[colnames(med_rep_1_2) == 'barcodes'] <- 'barcodes_DNA'
               

log_med_rep_1_2 <- med_rep_1_2 %>%
  mutate(log_BC_expression_25_1 = log10(med_BC_expression_25_1)) %>%
  mutate(log_BC_expression_25_2 = log10(med_BC_expression_25_2)) %>%
  mutate(log_BC_expression_0_1 = log10(med_BC_expression_0_1)) %>%
  mutate(log_BC_expression_0_2 = log10(med_BC_expression_0_2))

log_med_rep_1_2_0_R2 <- cor(log_med_rep_1_2$log_BC_expression_0_1,
                                log_med_rep_1_2$log_BC_expression_0_2,
                                use = "pairwise.complete.obs",
                                method = "pearson")

log_med_rep_1_2_25_R2 <- cor(log_med_rep_1_2$log_BC_expression_25_1,
                              log_med_rep_1_2$log_BC_expression_25_2,
                              use = "pairwise.complete.obs",
                              method = "pearson")

p_log_med_0_var_rep <- ggplot(NULL, aes(log_BC_expression_0_1, 
                                            log_BC_expression_0_2)) +
  geom_point(data = log_med_rep_1_2, aes(color = log2(ave_barcodes_0)),
             alpha = 0.3) +
  geom_point(data = filter(log_med_rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  xlab("Log10 Median BC Expression Rep. 1") +
  ylab("Log10 Median BC Expression Rep. 2") +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(-0.3, 2.4)) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-0.3, 2.4)) +
  annotation_logticks() +
  scale_color_viridis() +
  labs(color = 'log2 Average\nNumber\nBarcodes') +
  annotate("text", x = 0.1, y = 1.1,
           label = paste('r =', round(log_med_rep_1_2_0_R2, 2)))

save_plot('plots/log_med_0_var_rep.png',
          p_log_med_0_var_rep, scale = 0.4, 
          base_width = 12, base_height = 8)

p_log_med_25_var_rep <- ggplot(NULL, aes(log_BC_expression_25_1, 
                                          log_BC_expression_25_2)) +
  geom_point(data = log_med_rep_1_2, aes(color = log2(ave_barcodes_25)),
             alpha = 0.3) +
  geom_point(data = filter(log_med_rep_1_2, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  xlab("Log10 Median BC Expression Rep. 1") +
  ylab("Log10 Median BC Expression Rep. 2") +
  scale_x_continuous(breaks = c(0, 1, 2), limits = c(-0.3, 2.4)) + 
  scale_y_continuous(breaks = c(0, 1, 2), limits = c(-0.3, 2.4)) +
  annotation_logticks() +
  scale_color_viridis() +
  labs(color = 'log2 Average\nNumber\nBarcodes') +
  annotate("text", x = 0.1, y = 1.1,
           label = paste('r =', round(log_med_rep_1_2_25_R2, 2)))

save_plot('plots/log_med_25_var_rep.png',
          p_log_med_25_var_rep, scale = 0.4, 
          base_width = 12, base_height = 8)

#Take the average BC expression at induced and uninduced across reps.
med_ave_variant <- med_rep_1_2 %>%
  select(subpool, name, most_common, med_BC_expression_25_1, 
         med_BC_expression_25_2, med_BC_expression_0_1, 
         med_BC_expression_0_2, barcodes_DNA, ave_barcodes_25, 
         ave_barcodes_0) %>%
  mutate(ave_BC_expression_25 = (med_BC_expression_25_1 + med_BC_expression_25_2)/2) %>%
  mutate(ave_BC_expression_0 = (med_BC_expression_0_1 + med_BC_expression_0_2)/2)

#Join median and summed variant expression tables
med_BC_join_sum_variant <- inner_join(ave_variant, med_ave_variant,
                                      by = c('subpool', 'name', 
                                             'most_common', 'barcodes_DNA',
                                             'ave_barcodes_25',
                                             'ave_barcodes_0'))

med_sum_0_R2 <- cor(med_BC_join_sum_variant$ave_ratio_0,
                        med_BC_join_sum_variant$ave_BC_expression_0,
                        use = "pairwise.complete.obs",
                        method = "pearson")

med_sum_25_R2 <- cor(med_BC_join_sum_variant$ave_ratio_25,
                      med_BC_join_sum_variant$ave_BC_expression_25,
                      use = "pairwise.complete.obs",
                      method = "pearson")

p_med_sum_0 <- ggplot(NULL, aes(ave_ratio_0, 
                                    ave_BC_expression_0)) +
  geom_point(data = med_BC_join_sum_variant, 
             aes(color = log2(ave_barcodes_0)),
             alpha = 0.3) +
  geom_point(data = filter(med_BC_join_sum_variant, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  xlab("Average Summed BC\nRNA/DNA per Variant") +
  ylab("Average Median BC\nExpression per Variant") +
  scale_color_viridis() +
  labs(color = 'log2 Average\nNumber\nBarcodes') +
  annotate("text", x = 2, y = 8,
           label = paste('r =', round(med_sum_0_R2, 2)))

save_plot('plots/med_sum_0_neg.png',
          p_med_sum_0, scale = 0.4, 
          base_width = 12, base_height = 8)

p_med_sum_25 <- ggplot(NULL, aes(ave_ratio_25, 
                                  ave_BC_expression_25)) +
  geom_point(data = med_BC_join_sum_variant, 
             aes(color = log2(ave_barcodes_25)),
             alpha = 0.3) +
  geom_point(data = filter(med_BC_join_sum_variant, 
                           grepl(
                             'subpool5_no_site_no_site_no_site_no_site_no_site_no_site',
                             name)), 
             color = 'red') +
  xlab("Average Summed BC\nRNA/DNA per Variant") +
  ylab("Average Median BC\nExpression per Variant") +
  scale_color_viridis() +
  labs(color = 'log2 Average\nNumber\nBarcodes') +
  annotate("text", x = 25, y = 150,
           label = paste('r =', round(med_sum_25_R2, 2)))

save_plot('plots/med_sum_25_neg.png',
          p_med_sum_25, scale = 0.4, 
          base_width = 12, base_height = 8)


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
