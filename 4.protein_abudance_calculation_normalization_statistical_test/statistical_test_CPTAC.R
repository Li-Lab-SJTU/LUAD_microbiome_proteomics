rm(list = ls())
library(limma)
library(ggplot2)
library(dplyr)
library(openxlsx)

data_df = read.table(
  './03.normalized_protein_abundance/CPTAC.protein_abundance.irs_sl.NAfiltered.txt.noheader',
  header = F,
  sep = '\t',
  row.names = 1
)
sample_info_df = read.table(
  './03.normalized_protein_abundance/CPTAC.protein_abundance.irs_sl.NAfiltered.txt.header',
  sep = '\t',
  row.names = 1
)

data_df[data_df == 0] = NA

##log2-transformation
data_df = log2(data_df)

test_result = data.frame(taxonomy = character(),
                         statistic = double(),
                         p.value = double())

#############two-part test###################
TwoPart <- function(data,
                    group,
                    test = "t.test",
                    point.mass = 0) {
  Index1 <- c(group == 1)
  Group1 <- data[Index1]
  Group0 <- data[!Index1]
  n1 <- length(Group1)
  n2 <- length(Group0)
  obs <- c(n1, n2)
  success <- c(sum(Group1 != point.mass), sum(Group0 != point.mass))
  pointmass <- obs - success
  if (sum(success) == 0) {
    T2 <- 0
    B2 <- 0
  }
  else if ((success[1] == 0) | (success[2] == 0)) {
    T2 <- 0
    B2 <- prop.test(pointmass, obs)$statistic
  }
  else if ((success[1] == 1) | (success[2] == 1)) {
    T2 <- 0
    B2 <- prop.test(pointmass, obs)$statistic
  }
  else {
    uniq1 <- length(unique(Group1[Group1 != point.mass]))
    uniq2 <- length(unique(Group0[Group0 != point.mass]))
    if ((uniq1 < 2) & (uniq2 < 2)) {
      T2 <- 0
      if (sum(pointmass) == 0)
        B2 <- 0
      else
        B2 <- prop.test(pointmass, obs)$statistic
    }
    else if (sum(pointmass) == 0) {
      B2 <- 0
      if (test == "t.test")
        T2 <- t.test(data ~ group)$statistic ^ 2
      if (test == "wilcoxon") {
        W <- wilcox.test(data ~ group, exact = FALSE)$statistic
        mu <- (n1 * n2) / 2
        sigma <- sqrt((n1 * n2 * (n1 + n2 + 1)) / 12)
        T2 <- ((abs(W - mu) - 0.5) / sigma) ^ 2
      }
    }
    else {
      B2 <- prop.test(pointmass, obs)$statistic
      contIndex <- data != point.mass
      cont <- data[contIndex]
      cGroup <- group[contIndex]
      n1c <- sum(cGroup == 1)
      n2c <- sum(cGroup == 0)
      if (test == "t.test")
        T2 <- t.test(cont ~ cGroup)$statistic ^ 2
      if (test == "wilcoxon") {
        W <- wilcox.test(cont ~ cGroup, exact = FALSE)$statistic
        mu <- (n1c * n2c) / 2
        sigma <- sqrt((n1c * n2c * (n1c + n2c + 1)) / 12)
        T2 <- ((abs(W - mu) - 0.5) / sigma) ^ 2
      }
    }
  }
  X2 <- B2 + T2
  if ((T2 == 0) | (B2 == 0)) {
    X2pv <- 1 - pchisq(X2, 1)
  } else {
    X2pv <- 1 - pchisq(X2, 2)
  }
  ans <- list(statistic = as.numeric(X2), pvalue = as.numeric(X2pv))
  return(ans)
}

for (i in 1:nrow(data_df)) {
  data_new = c(t(data_df[i,]))
  taxonomy_new = rownames(data_df)[i]
  meta = c(t(sample_info_df['sample_type',]))
  data_new_nat = data_new[meta == 'NAT']
  data_new_tumor = data_new[meta == 'Tumor']
  data_new_nat_median = median(2 ^ data_new_nat, na.rm = TRUE)
  data_new_tumor_median = median(2 ^ data_new_tumor, na.rm = TRUE)
  log2fc_new = log2(data_new_tumor_median/data_new_nat_median)
  
  data_new[is.na(data_new)] = 0
  meta[meta == 'NAT'] = 0
  meta[meta == 'Tumor'] = 1
  result_new = TwoPart(data_new, meta, test = "wilcoxon", point.mass = 0)
  test_result = rbind(
    test_result,
    data.frame(
      taxonomy = taxonomy_new,
      statistic = result_new$statistic,
      p.value = result_new$pvalue,
      tumor_median = data_new_tumor_median,
      nat_median = data_new_nat_median,
      log2fc = log2fc_new
    )
  )
}

test_result$p.adjust = p.adjust(test_result$p.value,
                                method = "BH")
test_result$direction = test_result$log2fc > 0
test_result$direction[test_result$direction == T] = 'Up'
test_result$direction[test_result$direction == F] = 'Down'
test_result_filtered = test_result %>% filter(p.adjust < 0.05)
test_result_filtered = test_result_filtered %>% filter(abs(log2fc) > 0.5)

table(test_result_filtered$direction)

write.xlsx(test_result, './04.statistical_test_results/test_result_CPTAC_protein.xlsx')
