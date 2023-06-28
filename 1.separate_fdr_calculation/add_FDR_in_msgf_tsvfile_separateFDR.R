# Calculate separate fdr of microbial PSMs in ms-gf+ .tsv output file

#setwd("$work_dir")
# load libraries
library(stringr)
library(dplyr)

# delete the trypsin prefix and modification to get the pure pepetide sequences
add_PurePeptide_col = function(input_frame) {
  input_frame$Pure_Peptide = input_frame$Peptide
  input_frame$Pure_Peptide = gsub("[\\+-]\\d+.?\\d*", "", input_frame$Pure_Peptide, perl = TRUE)
  return(input_frame)
}

# read the ms-gf+ .tsv output file and calculate the separate FDR
handle_tsv = function(msgf_tsv_path) {
  print(paste("Handling", msgf_tsv_path))
  in_frame = read.csv(msgf_tsv_path, sep = '\t')
  in_frame = add_PurePeptide_col(in_frame)
  in_frame$SpecFile_ScanNum = paste(in_frame$X.SpecFile,
                                    in_frame$ScanNum)
  
  in_frame$Decoy = str_detect(in_frame$Protein, 'DECOY_')
  in_frame$Human = str_detect(in_frame$Protein, 'HUMAN')
  
  # micro : microbial
  micro_decoy_num = sum(in_frame$Decoy & !in_frame$Human)
  decoy_num = sum(in_frame$Decoy)
  
  # calculate the ratio of microbial decoy PSMs
  micro_decoy_p = micro_decoy_num / decoy_num
  in_frame$Micro_Target_count = ''
  in_frame$Decoy_count = cumsum(in_frame$Decoy)
  
  in_frame$scan_count = ''
  in_frame$sp_fdr = ''
  
  # filter the microbial PSMs and calculate the separate FDR
  is_micro = !str_detect(in_frame$Protein, "HUMAN")
  in_micro_frame = in_frame[is_micro,]
  
  # uniqueSS : unqiue SpecFile_ScanNum
  # one spectra may match to different peptides
  # especially when amino acid "I" and "L" were present in the peptide sequence
  # we choose one of these PSMs as a representative in the following analysis
  in_micro_frame_uniqueSS = in_micro_frame[!duplicated(in_micro_frame$SpecFile_ScanNum), ]
  
  
  in_micro_frame_uniqueSS$Micro_Target_count = cumsum(!in_micro_frame_uniqueSS$Decoy &
                                                        !in_micro_frame_uniqueSS$Human)
  # sp_fdr  separate FDR
  in_micro_frame_uniqueSS[in_micro_frame_uniqueSS$Micro_Target_count == 0, 'sp_fdr'] = 0
  # calculate separate FDR
  # sp_FDR=(D+ × D_n/D)/(T_n+)
  # D+ : the number of decoy sequences identified with scores above the score threshold.
  # D_n : the number of identified decoy microbial sequences
  # D : the total number of identified decoy sequences
  # T_n+ : the number of microbial sequence identifications in the target database above the score threshold
  in_micro_frame_uniqueSS[in_micro_frame_uniqueSS$Micro_Target_count != 0, 'sp_fdr'] =
    in_micro_frame_uniqueSS[in_micro_frame_uniqueSS$Micro_Target_count != 0, 'Decoy_count'] *
    micro_decoy_p /
    in_micro_frame_uniqueSS[in_micro_frame_uniqueSS$Micro_Target_count != 0, 'Micro_Target_count']
  
  # decide how many microbial PSMs can be accepted under separate FDR control
  max_micro_accept = max(which(in_micro_frame_uniqueSS$sp_fdr < 0.01))
  in_micro_frame_uniqueSS$scan_count = 1:dim(in_micro_frame_uniqueSS)[1]
  
  in_spfdr001_micro_target_frame_uniqueSS = in_micro_frame_uniqueSS %>% filter(scan_count <= max_micro_accept)
  in_spfdr001_micro_target_frame_uniqueSS = in_spfdr001_micro_target_frame_uniqueSS %>% filter(Decoy == FALSE)
  # return the filtered microbial PSMs after separate FDR control
  return(in_spfdr001_micro_target_frame_uniqueSS)
}

# script test in one .tsv file
msgf_tsv_path = "20170207_LC_TMTB1_prot_F15_01.mzML.s3.smalldb.tsv"
out_frame = handle_tsv(msgf_tsv_path)
