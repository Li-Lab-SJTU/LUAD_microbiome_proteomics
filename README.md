# Proteomic Lung Adenocarcinoma Microbiome Analysis
The codes of statistical analysis and visualization of "Proteomic Characterization Reveals Microbiome Diversity in Lung Adenocarcinoma Cohorts"

# Directory contents
## 1.separate_fdr_calculation
__add_FDR_in_msgf_tsvfile_separateFDR_revise.R :__ Calculate the separate fdr of microbial PSMs in the MF-GF+ output files.

__20170207_LC_TMTB1_prot_F15_01.mzML.s3.smalldb.tsv :__ An example file for the separate fdr calculation. MS-GF+ Output file of one fraction of MS experiments in the TW cohort. 

## 2.taxon_go_abundance_calculation_normalization
__MS_experiment_sample_info :__ The sample information of each plex in the two TMT-10 plex analysis cohorts (CPTAC and TW).

__normalization_plots :__ The taxon/go abundance distribution before and after normalization, corresponding to the Supplementary Figure S2 BC of the manuscript.

__output_go/taxon_abundance_data :__ The abundances of microbial taxa or GO terms in samples of each cohort.

-  `*.noheader` : The abundances data.       
-  `*.header` : The header file, contains the sample names and sample types (tumor or NAT) of each cohort.

__peptide_intensity_data :__ The intensity values of identified microbial peptides in each sample in the three cohorts.

__unipept_output_data :__ The taxonomic annotations and GO (Gene Ontology) annotations of the identified microbial peptides in Unipept metaproteomics annotation. 

***_go/abundance_calculation_normalization.R :** The script to calculate and normalize the go terms/taxon abundances in each sample of the three cohorts, combining the peptide_intensity_data and unipept_output_data.

## 3.two_part_wilcoxon_test
__test_result :__ The results of statistical tests of the taxon/go term abundances between the tumor and NAT samples in each cohort.

__statistical_test.R :__ The script of statistical test. 

__two_part_test.R :__ The function of two-part wilcoxon test.

## 4.protein_abudance_calculation_normalization_statistical_test
__01.protein_channel_intensity_calculation :__ Combine the MS-GF+ database search result and Masic peptides quantification result to get the protein intensity in each sample. (In TMT labeled samples of CPTAC and TW cohorts, to get the protein intensity in each channel of TMT labels.)

-  `masic_result_dir_TW` : Put all the output files of Masic in one directory.     
-  `01.get_protein_scan_protein.TW.test.R & 02.TW_combine_scan_protein_intensity.R` : Two steps to combine the MS-GF+ output files and Masic output files to get the protein intensity in each sample. The result will be saved to `../02.protein_channel_intensity`.
-  `20170207_LC_TMTB1_prot_F1_01.mzML.s3.smalldb.tsv.filtered.txt` : The merged processed output file of MS-GF+ (Processed in `../1.separate_fdr_calculation`). Here we take one fraction of TW experiments as an example.
-  `TW.microbiome.proteins_scan.F1_01.txt` : The output file of 01.`get_protein_scan_protein.TW.test.R`.

***_protein_abundance_calculation_normalization.R :** The script of calculating the microbial protein abundance and normalization. The result will be saved in `03.normalized_protein_abundance` and `normalization_plots`.

**statistical_test_*.R :** The script of statistical test of the differences between the abundance of microbial proteins of tumor and NAT samples. The result will be saved in `04.statistical_test_results`. 


## 5.survival_analysis
__Cox_bioForest_plots :__ The bioForest plots of significant microbial taxa in Cox models.

__KM_plots :__ The Kaplan-Meier plots of significant microbial taxa in Cox models, corresponding to the Figure 6 of the manuscript.


__CNHPP_clinical_mmc1.xlsx :__ The clinical information of CNHPP patients, including the survival data (disease free survival and overall survival) of each patient.

__CNHPP_DFS/OS_result.xlsx :__ The results of Univariate Cox proportional hazards model between the abundances of microbial taxa and DFS and OS.

__survival_analysis.R :__ The script of Cox analysis and Kaplan-Meier plots.

## 6.other_visulization
The R scripts of creating figure 1BC, 2B, 3B, 4A.

## ggplot_theme
The code of ggplot themes in creating the figures.
