import numpy as np
import pandas as pd



def make_binary(df, filter_genes='WBG'):
    '''
    Takes a Pandas DataFrame with Count-Per-Million normalized read counts and returns a binarized version of it
    :param df: Pandas DataFrame with a row for each Sample.
    Columns contain Count-Per-Million(CPM) read counts, and possibly meta-data,i.e.
    the Strain, Treatment, RNAi, Biological Age, GEO accession number
    :param filter_genes: Filter columns names by filter_genes to compute the binarization only on the genes.
    :return: A binarized copy of the original data without meta-information
    '''
    df_bin = df.copy()
    df_bin = df_bin.filter(regex=filter_genes)  # get the gene columns
    df_bin[df_bin == 0] = np.nan  # mask 0-genes that skew the median
    df_bin['Median'] = df_bin.median(axis=1)  # calculate the median for each row and append a column
    # divide each value, i.e. CPM, by the sample-median
    df_bin = df_bin.filter(regex=filter_genes).div(df_bin.Median, axis=0)
    # set values smaller than the median to 0 and 1 otherwise
    df_bin[df_bin.isna()] = 0
    df_bin[df_bin <= 1] = 0
    df_bin[df_bin > 1] = 1

    return df_bin





def preproc_new_data(cpm_gene_file, predictor_gene_file, sep='\t'):
    '''
    Read a CPM containing CSV-file, binarize it and return only the subset of genes relevant for the prediction

    :param cpm_gene_file: A file with non-duplicate WormbaseID gene identifiers in the first column and 1 or more columns with
    gene expression samples.
    For example:

    Gene_ID	SRR6207028	SRR6207027	SRR6207030	SRR6207029
    WBGene00007063	2.9784207155	1.4274445525	0.9728337164	3.2820503369
    WBGene00200402	0	0	0	0
    WBGene00007064	21.8278717092	16.4387110019	32.8859628004	30.5373379176
    WBGene00044951	0	0	0	0
    WBGene00195865	0	0	0	0
    WBGene00197051	0	0	0	0
    WBGene00199694	0	0	0	0
    WBGene00199940	0	0	0	0
    WBGene00201913	0	0	0	0
    WBGene00007065	4.6394630376	2.9846567916	3.0157845209	5.1371222665
    ....

    :param predictor_gene_file:  The file with the elastic net coefficients for the 576 genes (Predictor_Genes.csv):
    WormBaseID	ElasticNet_Coef
    WBGene00012747	-3.9183849126
    WBGene00011554	-3.6615313969
    WBGene00002259	-3.4419756727
    WBGene00018196	-3.3572860482
    ....

    :param sep: Delimiter to use, default '\t'
    :return: A DataFrame with binarized CPMs for the 576 genes relevant for the prediction and an additional last column with
    the elastic net coefficients
    '''

    # Read the predictor_gene_file
    predictor_gene_df = pd.read_csv(predictor_gene_file, index_col=0)

    # Read the new data and transpose it to make it fit for binarization
    cpm_df = pd.read_csv(cpm_gene_file, index_col=0, sep=sep)
    cpm_df = cpm_df.T
    binarized_df = make_binary(cpm_df)

    # take the subset of genes that is relevant for the prediction and transpose the DataFrame to join the elastic
    # net coefficients
    df_for_prediction = binarized_df[list(set(predictor_gene_df.index))]
    df_for_prediction = df_for_prediction.T
    # attach the elastic net coefficients to the DataFrame
    df_for_prediction = df_for_prediction.join(predictor_gene_df)

    return df_for_prediction


def predict(cpm_gene_file, predictor_gene_file, cpm_gene_file_sep='\t'):
    '''

    :param cpm_gene_file: A file with non-duplicate WormbaseID gene identifiers in the first column and 1 or more columns with
    gene expression samples.
    For example:
    Gene_ID	SRR6207028	SRR6207027	SRR6207030	SRR6207029
    WBGene00007063	2.9784207155	1.4274445525	0.9728337164	3.2820503369
    WBGene00200402	0	0	0	0
    WBGene00007064	21.8278717092	16.4387110019	32.8859628004	30.5373379176
    WBGene00044951	0	0	0	0
    WBGene00195865	0	0	0	0
    WBGene00197051	0	0	0	0
    WBGene00199694	0	0	0	0
    WBGene00199940	0	0	0	0
    WBGene00201913	0	0	0	0
    WBGene00007065	4.6394630376	2.9846567916	3.0157845209	5.1371222665
    ....

    :param predictor_gene_file:  The file with the elastic net coefficients for the 576 genes (Predictor_Genes.csv):
    WormBaseID	ElasticNet_Coef
    WBGene00012747	-3.9183849126
    WBGene00011554	-3.6615313969
    WBGene00002259	-3.4419756727
    WBGene00018196	-3.3572860482
    ....

    :param cpm_gene_file_sep: Delimiter to use, default '\t'

    :return: Pandas DataFrame with the sample names (column names) of cpm_gene_file as the index and
    the predicted biological age in hours in the second column
    '''

    # The intercept has to be added to get the final prediction
    intercept = 103.54631743289005

    # Read a CPM containing CSV-file (cpm_gene_file), binarize the gene expression and
    # get only the subset of genes relevant for the prediction (predictor_gene_file)
    # the last column of df_for_prediction contains the elastic net coefficients for the prediction
    df_for_prediction = preproc_new_data(cpm_gene_file, predictor_gene_file, sep=cpm_gene_file_sep)

    prediction_result_dict = {}
    prediction_result_dict['Sample'] = []
    prediction_result_dict['Predicted_Biological_Age'] = []
    # loop over the samples (columns) in the DataFrame, predict the biological age in hours and save it to the
    # dictionary
    for i in range(len(df_for_prediction.columns) - 1):
        prediction_result_dict['Sample'].append(df_for_prediction.columns[i])
        # Sum up all coefficients for genes that are 1 and thereby contribute to the biological age
        prediction_result_dict['Predicted_Biological_Age'].append(
            sum(df_for_prediction.iloc[:, i] * df_for_prediction.ElasticNet_Coef) + intercept)

    prediction_result_df = pd.DataFrame(prediction_result_dict)
    prediction_result_df = prediction_result_df.set_index('Sample')

    return prediction_result_df

def test_prediction():
    '''
    Usage Example.
    The CPM gene file and the predictor gene file can be found in data/
    :return: 
    '''
    cpm_gene_file = 'GSE65765_CPM.csv'
    predictor_gene_file = 'Predictor_Genes.csv'
    res = predict(cpm_gene_file, predictor_gene_file, cpm_gene_file_sep='\t')
    print(res)

