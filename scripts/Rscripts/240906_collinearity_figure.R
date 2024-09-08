library(tidyverse)
library(gggenomes)

seq_file <- "seq_table.txt"
gene_file <- "gene_seq_rearranged.tsv"
links_file <- "links_rearranged.tsv"


# Read seq table
s <- read_tsv(seq_file)

# Read feature table with ORF annotations
g <- read_tsv(gene_file) %>% 
  mutate(Order = case_when(Order == "1" ~ "01",
                           Order == "2" ~ "02",
                           Order == "3" ~ "03",
                           Order == "4" ~ "04",
                           Order == "5" ~ "05",
                           Order == "6" ~ "06",
                           Order == "7" ~ "07",
                           Order == "8" ~ "08",
                           Order == "9" ~ "09",
                           TRUE ~ Order))


# Read links table
l <- read_tsv(links_file) %>% 
  mutate(Order = case_when(Order == "1" ~ "01",
                           Order == "2" ~ "02",
                           Order == "3" ~ "03",
                           Order == "4" ~ "04",
                           Order == "5" ~ "05",
                           Order == "6" ~ "06",
                           Order == "7" ~ "07",
                           Order == "8" ~ "08",
                           Order == "9" ~ "09",
                           TRUE ~ Order))



# Create data object for plotting
p <- gggenomes(seqs=s, genes=g, links=l)

p + 
  geom_seq() +         # draw contig/chromosome lines
  geom_seq_label() +    # label each sequence 
  geom_gene(aes(fill = Order)) +         # draw genes as arrow
  geom_link(aes(fill = Order, color=Order))          # draw some connections between syntenic regions


p + 
  geom_link_line(aes(color = Order), show.legend = FALSE)  +        # draw some connections between syntenic regions
  geom_seq() +         # draw contig/chromosome lines
#  geom_seq_label() +    # label each sequence 
  geom_gene(aes(fill = Order), size = 3) +  # draw genes as arrow
  theme(legend.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.2, "in")) +
  labs(fill = "Core gene")

ggsave("collinearity_plot.pdf", plot = last_plot(), width = 7.2, height = 5, units = "in")
