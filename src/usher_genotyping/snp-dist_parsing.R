library(reshape2)
library(ggplot2)

# Read in SNP matrix created by SNP-dists
snps <-read.csv('/home/phe.gov.uk/nicholas.ellaby/Documents/vector_born_diseases/dengue/data/20231004/ncbi/ncbi_usher_genotyping/complete_genomic.snp-dists.tsv', sep='\t', header=T, row.names=T)
# Make sure the reference isn't duplicated in the matrix as this might cause problems in naming rows

# Add first column as 
row.names(snps) <-snps[,1]
snps <- snps[-1]

# Melt SNP Matrix
# Sample 1 | Sample 2 | SNPs
snps_melt <- melt(as.matrix(snps))[melt(upper.tri(as.matrix(snps)))$value,]

# Generate graphs
snp_counts <- snps_melt %>% group_by(value) %>%  summarize(counts = n())

# Boxplot for SNP counts
p1<-ggplot(snp_counts, aes(y=value))+geom_boxplot()+scale_x_discrete( )+theme_bw()+labs(y='Number of SNPs')

# Number of Samples with SNP counts
p2<-ggplot(snp_counts, aes(x=value, y=counts))+geom_bar(stat='identity')+theme_bw()+theme(axis.text.x = element_text(angle = 90))+labs(x='Number of SNPs', y='Number of Pairwise Sample Combinations')

p_combined<-grid.arrange(p1,p2, nrow=1)
ggsave('snp_dist_summary_graphs.png',p_combined, heigh=9, width=15)

# Summary Table
row_names<-c('Maximum','Upper Quartile','Median','Lower_ Quartile','min')
values <-c(max(snp_counts$value), quantile(snp_counts$value)[[4]], median(snp_counts$value),  quantile(snp_counts$value)[[2]], min(snp_counts$value))
summary_table <- do.call(rbind, Map(data.frame, "Stat"=row_names, 'SNP Distance'=values))
summary_table <- summary_table[-1]
write.csv(x=summary_table, file="SNP-dist_summary_metrics.csv", row.names=T)

# Identify samples of interest
# Select pairwise combinations that exceed the upper quartile value
up_quartile_snps <-snps_melt[snps_melt$value > quantile(snp_counts$value)[[4]],]
# Combine all unique instances a sample is in a pairwise combination > the upper quartile
up_quartile_snps_ids<- c(up_quartile_snps$Var1,up_quartile_snps$Var2)
# Count number of instances a samples has been above the upper quartile SNP threshold
up_quartile_snps_ids_instances<- as.data.frame(table(up_quartile_snps_ids))
# remove empty values
up_quartile_snps_ids_instances <- up_quartile_snps_ids_instances[up_quartile_snps_ids_instances$Fre > 0, ]
# order by freq
up_quartile_snps_ids_instances <- up_quartile_snps_ids_instances[rev(order(up_quartile_snps_ids_instances$Freq)),]
# rename columns
names(up_quartile_snps_ids_instances)<- c("accession","PW Comb. > UQ")
write.csv(x=up_quartile_snps_ids_instances, file="Instances_sample_was_present_in_pairwise_combinations_greater_than_SNP_upper_quantile_value.csv", row.names=F)
