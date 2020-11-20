# AgingClock
The repository contains the code to predict the biological age based on gene expression data of C. elegans.


The predict function in src/biological_age_prediction.py takes 2 files as input  
- a CSV file with count-per-million normalized RNA-seq counts with genes as rows, and samples as columns  (e.g. Data/GSE65765_CPM.csv)
- the CSV file Data/Predictor_Genes.csv with the predictor genes and corresponding regression coefficients  

The function will return a Pandas DataFrame with the sample names in the first and the predicted biological age in the second column.



To calculate the second biological age correction described in the paper run the calculate_Bio_Age_correction function in src/biological_age_correction.py.  
It takes a Pandas DataFrame as input with samples as rows and at least one column with the biological age of the samples. 


For more details see:

A transcriptome based aging clock near the theoretical limit of accuracy  
https://www.biorxiv.org/content/10.1101/2020.05.29.123430v1



The example data is downloaded from:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse65765  
Sen, P. et al. H3K36 methylation promotes longevity by enhancing transcriptional fidelity. Genes Dev. 29, 1362–1376 (2015).

