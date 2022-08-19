# IRP
Independent Research Project: Effects of Valproic Acid and Similarity of Gene Expression Changes Compared to Obesity

## Acknowledgement: 
The owner of this repository would like to thank Dr Marcin Wozniak for writing all of the scripts available in this repository with the exception of analysis_foldAndPathComparison.R and analysis_WCNA.R as well as certain parts of GEO_datasetsSearch.ipnyb and analysis_differentialGE.R.

## Pipeline:

### I. Getting raw data from NCBI GEO:
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
