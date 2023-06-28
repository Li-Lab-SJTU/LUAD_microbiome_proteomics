# Proteomic Lung Adenocarcinoma Microbiome Analysis
The codes of statistical analysis and visualization of "Characterization of tumor-associated microorganisms with the proteomics data from independent lung adenocarcinoma cohorts"

Contact:

Mou Xinyi : xinyimou(at)sjtu(dot)edu(dot)cn

Li jing :  jing.li(at)sjtu(dot)edu(dot)cn

# Directory contents
## 1.separate_fdr_calculation
add_FDR_in_msgf_tsvfile_separateFDR_revise.R : Calculate the separate fdr of microbial PSMs in the MF-GF+ output files.

20170207_LC_TMTB1_prot_F15_01.mzML.s3.smalldb.tsv : An example file for the separate fdr calculation. MS-GF+ Output file of one fraction of MS experiments in the TW cohort. 

## 2.taxon_go_abundance_calculation_normalization
MS_experiment_sample_info : The sample information of each plex in the two TMT-10 plex analysis cohorts (CPTAC and TW).

normalization_plots : The taxon/go abundance distribution before and after normaliztion, corresponding to the Supplementary Figure S1 BC of the manuscript.

output_go/taxon_abundance_data : The abundances of microbal taxa or GO terms in samples of each cohort.
        
    *.noheader : The abundances data.
    *.header : The header file, contains the sample names and sample types (tumor or NAT) of each cohort.

peptide_intensity_data : The intensity values of identified microbial peptides in each sample identified in the three cohorts.

unipept_output_data : The taxonomic annotations and GO (Gene Ontology) annotations of the identified microbial peptides in Unipept meteproteomics annotation. 

*_go/abundance_calculation_normalization.R : The script to calculate and normalize the go terms/taxon abundaces in each sample of the three cohorts, combining the peptide_intensity_data and unipept_output_data.

## 3.two_part_wilcoxon_test
test_result : The results of statistical test of the taxon/go term abundances between the tumor and NAT samples in each cohort.

statistical_test.R : The script of statistical test. 

two_part_test.R : The function of two-part wilcoxon test.
## 4.survival_analysis
Cox_bioForest_plots : The bioForest plots of significant microbial taxa in Cox models.

KM_plots : The Kaplan-Meier plots of significant microbial taxa in Cox models, corresponding to the Figure 7 of the manuscript.


CNHPP_clinical_mmc1.xlsx : The clinical information of CNHPP patients, including the survival data (disease free survival and overall survival) of each patient.

CNHPP_DFS/OS_result.xlsx : The results of Univariate Cox proportional hazards model between the abundances of microbial taxa and DFS and OS.

survival_analysis.R : The script of Cox analysis and Kaplan-Meier plots.
## 5.other_visulization
The R scripts of creating figure 2, 3b, 4b, 5a and 5b.

## ggplot_theme
The code of ggplot themes in creating the figures.