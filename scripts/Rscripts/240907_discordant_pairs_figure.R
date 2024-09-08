library(tidyverse)

# Files
discordant_pairs_file <- "/labs/kingsley/ambenj/TSIV/analysis_STMY_2012_42/assembly/02_spades/STMY12-42_down0.005/remapped/discordant_pairs.txt"
output_fig <-"/labs/kingsley/ambenj/TSIV/analysis_STMY_2012_42/assembly/02_spades/STMY12-42_down0.005/remapped/discordant_pairs_fig.pdf"


# Read sam formatted file with discordant read pairs
discordant_pairs <- read_tsv(discordant_pairs_file, col_names = FALSE) 

col_names <- c(R1_contig="X3", R1_pos="X4", MAPQ="X5", R2_contig="X7", R2_pos="X8")

# Rename columns and filter for high quality mapped reads
discordant_pairs_tidy <- discordant_pairs %>% 
  select(X3, X4, X5, X7, X8) %>% 
  rename(all_of(col_names)) %>% 
  filter(MAPQ >= 30)
  
# Filter for NODE1/contig1 forward reads
contig1_F <- discordant_pairs_tidy %>% 
  filter(R1_contig == "NODE_1_length_79043_cov_150.476925") %>% 
  select(R1_pos, R2_pos) %>% 
  rename(all_of(c(contig1_pos = "R1_pos", contig2_pos = "R2_pos")))

# Filter for NODE2/contig2 forward reads
contig2_F <- discordant_pairs_tidy %>% 
  filter(R1_contig == "NODE_2_length_36080_cov_145.724195") %>% 
  select(R2_pos, R1_pos) %>% 
  rename(all_of(c(contig1_pos = "R2_pos", contig2_pos = "R1_pos")))

# Join tables together
joined_pos <- rbind(contig1_F, contig2_F)

# Make plot
joined_pos %>% 
  ggplot(aes(contig1_pos, contig2_pos)) +
  geom_point(alpha=0.008) +
  xlab("Contig1 read position") +
  ylab("Contig2 paired read position") +
  theme_classic()

ggsave(output_fig, plot = last_plot(), width = 4, height = 3, units = "in")
