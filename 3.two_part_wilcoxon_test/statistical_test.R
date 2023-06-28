# Two part wilcoxon test of the taxon or go abundances between Tumor and NAT samples
library(dplyr)
library(plyr)
library(reshape2)
library(stringr)
library(openxlsx)

# Load the files of the taxa(GO) abundance matrix and the sample information
# abundance matrix : "../2.taxon_go_abundance_calculation_normalization/output_taxon_abundance_data/*.noheader" 
#                          or "../2.taxon_go_abundance_calculation_normalization/output_go_abundance_data/*.noheader"
# sample information : "../2.taxon_go_abundance_calculation_normalization/output_taxon_abundance_data/*.header"  
#                         or "../2.taxon_go_abundance_calculation_normalization/output_go_abundance_data/*.header"
# Take the statistical test of CNHPP taxon as an example
abundance_df = read.table(
  "../2.taxon_go_abundance_calculation_normalization/output_taxon_abundance_data/CNHPP.taxon_abundance.sl.NAfiltered.txt.noheader",
  sep = '\t',
  row.names = 1
)
sample_info_df = read.table(
  "../2.taxon_go_abundance_calculation_normalization/output_taxon_abundance_data/CNHPP.taxon_abundance.sl.NAfiltered.txt.header",
  sep = '\t',
  row.names = 1
)

# Data Frame to store the test result
test_result = data.frame(feature = character(),
                         statistic = double(),
                         p.value = double())

# Log2 transform of the feature abundance
abundance_df = log(abundance_df, base=2)

# The function of calculating two-part statistics is referred to the paper of Taylor & Pollard (2009). 
# DOI: 10.2202/1544-6115.1425
source("two_part_test.R")

# Test the abundance of each feature (each row in the taxon/go abundance matrix) separately
for (i in 1:nrow(abundance_df)) {
  data_new = c(t(abundance_df[i, ]))
  feature_new = rownames(abundance_df)[i]
  meta = c(t(sample_info_df['sample_type', ]))
  data_new_nat = data_new[meta=='NAT']
  data_new_tumor = data_new[meta=='Tumor']
  data_new_nat_median = median(2 ^ data_new_nat, na.rm = TRUE)
  data_new_tumor_median = median(2 ^ data_new_tumor, na.rm = TRUE)
  log2fc_new = log2(data_new_tumor_median/data_new_nat_median)
  
  data_new[is.na(data_new)] = 0
  meta[meta=='NAT'] = 0
  meta[meta=='Tumor'] = 1
  result_new = TwoPartTest(data_new, meta, test="wilcoxon", point.mass=0)
  test_result = rbind(test_result,
                      data.frame(
                        feature = feature_new,
                        statistic = result_new$statistic,
                        p.value = result_new$pvalue,
                        tumor_median = data_new_tumor_median,
                        nat_median = data_new_nat_median,
                        log2fc = log2fc_new
                      ))
}

# Benjamini-Hochberg correction for multiple test
test_result$p.adjust = p.adjust(test_result$p.value,
                            method = "BH")
# Decide the direction of change between tumor and NAT samples
test_result$direction = test_result$log2fc > 0
test_result$direction[test_result$direction == T] = 'Up'
test_result$direction[test_result$direction == F] = 'Down'

# Significantly differential feature 
test_result_filtered = test_result %>% filter(p.adjust < 0.05)
table(test_result_filtered$direction)

# Store the test result
write.xlsx(test_result, './test_result/test_result_CNHPP_taxon.xlsx')
