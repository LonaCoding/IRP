# IRP
Independent Research Project: Effects of Valproic Acid and Similarity of Gene Expression Changes Compared to Obesity

## Acknowledgement: 
The owner of this repository would like to thank Dr Marcin Wozniak for writing all of the scripts available in this repository with the exception of analysis_foldAndPathComparison.R and analysis_WCNA.R.

The following files were originally written and supplied by Dr Marcin Wozniak:
  1. GEO_datasetsSearch.ipnyb (with substantial modifications/additions done by the owner of this repository)
  2. GEO_dataExtraction.ipnyb
  3. GEO_rawDataDownload.ipnyb
  4. reading_Affymetrix.R
  5. reading_NimbleGenIllumina.R
  6. analysis_differentialGE.R (with substantial modifications/additions done by the owner of this repository)

analysis_WCNA.R, with substantial modifications/additions done by the owner of this repository, was mostly based on code from "Tutorials for the WGCNA package" by Peter Langfelder and Steve Horvath (https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html), specifically: 
<br>"Tutorial for the WGCNA package for R:</br>
<br>I. Network analysis of liver expression data in female mice</br>
<br>1. Data input and cleaning" </br>
<br>(https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf)</br>
<br>and </br>
<br>"Tutorial for the WGCNA package for R:</br>
<br>I. Network analysis of liver expression data in female mice</br>
<br>2.a Automatic network construction and module detection"</br>
<br>(https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf)</br>


All CSV files used as input in analysis_foldAndPathComparison.R, with the exception of vpa_significant pathway.csv and VPA_gene_set_enrichment.csv, were given by Dr Marcin Wozniak from unpublished materials for use in the project and would not be available via any scripts in this repository.


## Pipeline:

### I. Getting raw and meta data from NCBI GEO:
  1. GEO_datasetsSearch.ipnyb
  2. GEO_dataExtraction.ipnyb
  3. GEO_rawDataDownload.ipnyb

### II. Reading gene expression read files:

  4. reading_Affymetrix.R
  5. reading_NimbleGenIllumina.R

### III. Analyse data:

  6. analysis_differentialGE.R
  7. analysis_WCNA.R
  8. analysis_foldAndPathComparison.R
