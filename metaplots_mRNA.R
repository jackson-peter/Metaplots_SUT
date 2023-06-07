library(tidyverse)
library(data.table)
library(ggh4x)
library(yaml)
library(ggpubr)

################################################ Config
args = commandArgs(trailingOnly=TRUE)

config_f <- args[1]
#config_f <- "~/Scripts/SeedUTails/metaplots_mRNA_config.yaml"
config=read_yaml(file = config_f)

gene_list_f=config$gene_list_f
gene_list=read_tsv(gene_list_f, col_names = c("AGI"))

sample_list_f=config$sample_list_f
samples_list=read_tsv(sample_list_f, col_names = c("samples"))

repr_gene_model_f=config$repr_gene_model_f
repr_gene_model <- read_tsv(repr_gene_model_f,
                            col_names = c("chr", "start", "end", "feature", "score", "strand")) %>%
  separate(feature, into=c("repr_gene", NA), sep='_') %>%
  filter(repr_gene %in% gene_list$AGI)


tail_result_f=config$tail_result_f
tail_result <- read_tsv(tail_result_f) %>% filter(mRNA %in% gene_list$AGI,
                                                 sample %in% samples_list$samples)
gc()

annot_file=config$annot_file

nb_bin_UTR=as.numeric(config$nb_bin_UTR)
nb_bin_CDS=as.numeric(config$nb_bin_CDS)

outplot=config$outplot

########################################## FUNCTIONS
# Divide a region (start, end) on X number of bins
binRegion <- function(start, end, bins, idDF = NULL, info = NULL, strand_annot = "*") {
  finalColNames <- c("chr", "start", "end", "info", "strand_annot", "id", "binID", "ubinID")
  binSize <- (end - start) / (bins)
  breaks <- round(rep(start, each = (bins + 1))
                  + (0:(bins)) * rep(binSize, each = (bins + 1)))
  
  endpoints <- (bins + 1) * (1:(length(start)))
  startpoints <- 1 + (bins + 1)  * (0:(length(start) - 1))
  # do all regions in the same order
  dt <- data.table(
    start = breaks[-endpoints],
    end = breaks[-startpoints] - 1,
    #end = breaks[-startpoints] ,
    id = rep((seq_along(start)), each = bins),
    binID = 1:bins,
    ubinID = seq_along(breaks[-startpoints]),
    strand_annot = rep(strand_annot, each = bins),
    key = "id"
  )
  
  if (!is.null(idDF)) {
    chr <- rep(idDF, each = bins)
    #strand <- rep(Strand, each=bins)
    info <- rep(info, each = bins)
    
    dt <- dt[, chr := chr]
    dt <- dt[, info := info]
    setcolorder(dt, finalColNames)  # putting chr first, does not copy
  }
  
  return(dt[])
}


#tail_result$sample <- factor(tail_result$sample, levels=c("Col0", "urt1", "heso1", "heso1urt1"))

tail_result_gl <- tail_result %>% 
  inner_join(distinct(repr_gene_model%>%select(c("repr_gene", "strand"))), by=c("mRNA"="repr_gene"))
tail_result_gl$mRNA_start <- as.numeric(tail_result$mRNA_start)
tail_result_gl$mRNA_end <- as.numeric(tail_result$mRNA_end)
tail_result_gl$read_start <- as.numeric(tail_result$read_start)
tail_result_gl$read_end <- as.numeric(tail_result$read_end)
setDT(tail_result_gl)

annot_df <- read_tsv(annot_file, col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")) %>%
  #filter(feature %in% c("CDS", "three_prime_UTR", "five_prime_UTR")) %>%
  filter(feature %in% c("protein", "mRNA")) %>%
  separate(attributes, into = c(NA, "gene"), extra = "drop", sep=c("[;=,-]")) %>%
  filter(gene %in% gene_list$AGI) %>%
  select(-c("score", "frame")) %>%
  #mutate(five_p=case_when(strand=="+" ~ start,
  #                        strand=="-" ~ end),
  #       three_p=case_when(strand=="+" ~ end,
  #                         strand=="-" ~ start)) %>%
  mutate(five_p=start,
         three_p=end) %>%
  select(-c("start", "end")) %>%
  pivot_wider(names_from = feature, values_from = c("five_p", "three_p")) %>%
  relocate(three_p_protein, .after=five_p_protein) 
  #pivot_longer(c("five_p_mRNA", "five_p_protein"), names_to = "region_type", values_to = "val")

annot_df <- setDT(annot_df)

#save.image(file = "~/Scripts/SeedUTails/data.RData")
#load("~/Scripts/SeedUTails/data.RData")

# Get upstream region
upstream <- copy(annot_df)
upstream <- upstream %>%
  mutate(End=case_when(strand=="+" ~five_p_protein,
                       strand == "-" ~five_p_protein),
         Start=case_when(strand=="+" ~five_p_mRNA,
                         strand=="-" ~five_p_mRNA))

upstream <- setDT(upstream)


# Get downstream region
downstream <- copy(annot_df)

downstream<- downstream %>%
  mutate(Start=three_p_protein,
         End=three_p_mRNA)


# data.table "j command" using column names and numberOfBins variable
binnedUpstreamRegionDT <- upstream[, binRegion(Start, End, nb_bin_UTR, seqname, gene, strand)]
binnedUpstreamRegionDT <- binnedUpstreamRegionDT %>%
  mutate(region=case_when(strand_annot=="+" ~ "5pUTR",
                          strand_annot=="-" ~ "3pUTR"))


binnedCDSRegionDT <- annot_df[, binRegion(five_p_protein , three_p_protein , nb_bin_CDS, seqname, gene, strand)]
binnedCDSRegionDT[, binID := binID + nb_bin_UTR]
binnedCDSRegionDT[, region := "GENE"]

binnedCDSRegionDT %>% filter(info=="AT1G12410.1")
binnedCDSRegionDT %>% filter(info=="AT5G66510.2")

binnedDownstreamRegionDT <- downstream[, binRegion(Start, End, nb_bin_UTR, seqname, gene, strand)]
binnedDownstreamRegionDT[, binID := binID + nb_bin_UTR+nb_bin_CDS]
binnedDownstreamRegionDT <- binnedDownstreamRegionDT %>%
  mutate(region=case_when(strand_annot=="+" ~ "3pUTR",
                          strand_annot=="-" ~ "5pUTR"))

res <- rbindlist(list(binnedUpstreamRegionDT, binnedCDSRegionDT, binnedDownstreamRegionDT)) %>%
  rename(start_overlap=start,
         end_overlap=end) %>%
  select(-c(region))

res[strand_annot == "-", c("start_overlap", "end_overlap", "binID") := list(rev(start_overlap), rev(end_overlap), rev(binID)), by = id]


res_ol_5p <- tail_result_gl%>%
  #filter(mRNA!="AT2G31141.1") %>%
  filter(mRNA!="AT1G49700.1") %>%
  select(c("sample", "replicate", "mRNA", "read_start", "read_end", "strand", "condition")) %>%
  left_join(res, by=join_by(mRNA==info, strand==strand_annot, between(read_start, start_overlap, end_overlap))) %>%
  select(-c("read_start", "read_end")) %>%
  group_by(mRNA, sample, replicate, condition, strand,binID) %>%
  summarise(nb_reads_mRNA_sample_rep_cond_bin=n()) %>%
  group_by(mRNA, sample, replicate, condition, strand) %>%
  mutate(nb_reads_mRNA_sample_rep_cond=sum(nb_reads_mRNA_sample_rep_cond_bin),
         pct=100*nb_reads_mRNA_sample_rep_cond_bin/nb_reads_mRNA_sample_rep_cond) %>%
  group_by(binID, condition, replicate, sample) %>%
  summarise(mean_pct=mean(pct),
            median_pct=median(pct)) %>%
  filter(!is.na(binID)) %>%
  mutate(extr_type="5p")



res_ol_3p <- tail_result_gl%>%
  #filter(mRNA!="AT2G31141.1") %>%
  filter(mRNA!="AT1G49700.1") %>%
  select(c("sample", "replicate", "mRNA", "read_start", "read_end", "strand", "condition")) %>%
  left_join(res, by=join_by(mRNA==info, strand==strand_annot, between(read_end, start_overlap, end_overlap))) %>%
  select(-c("read_start", "read_end")) %>%
  group_by(mRNA, sample, replicate, condition, strand,binID) %>%
  summarise(nb_reads_mRNA_sample_rep_cond_bin=n()) %>%
  group_by(mRNA, sample, replicate, condition, strand) %>%
  mutate(nb_reads_mRNA_sample_rep_cond=sum(nb_reads_mRNA_sample_rep_cond_bin),
         pct=100*nb_reads_mRNA_sample_rep_cond_bin/nb_reads_mRNA_sample_rep_cond) %>%
  group_by(binID, condition, replicate, sample) %>%
  summarise(mean_pct=mean(pct),
            median_pct=median(pct)) %>%
  filter(!is.na(binID)) %>%
  mutate(extr_type="3p")




#test %>% filter(mRNA=="AT1G51110.1")

res_ol_extr <- rbind(res_ol_5p, res_ol_3p)


med_plot <- ggplot(res_ol_extr, aes(x=binID, y=median_pct, color=extr_type)) + 
  #geom_line(aes(linetype=extr_type)) +
  geom_line() +
  geom_vline(xintercept = nb_bin_UTR) +
  geom_vline(xintercept = nb_bin_UTR+nb_bin_CDS) +
  theme_bw() +
  #scale_linetype_manual(values=c("dotdash", "longdash"))+
  facet_nested(replicate ~ condition + sample) +
  ggtitle("Metaplots of 5' & 3' extremities (median)")


mean_plot <- ggplot(res_ol_extr, aes(x=binID, y=mean_pct, color=extr_type)) + 
  #geom_line(aes(linetype=extr_type)) +
  geom_line() +
  geom_vline(xintercept = nb_bin_UTR) +
  geom_vline(xintercept = nb_bin_UTR+nb_bin_CDS) +
  theme_bw() +
  #scale_linetype_manual(values=c("dotdash", "longdash"))+
  facet_nested(replicate ~ condition + sample) +
  ggtitle("Metaplots of 5' & 3' extremities (mean)")


figure <- ggarrange(med_plot, mean_plot,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2)
figure
ggexport(figure, filename = outplot)


