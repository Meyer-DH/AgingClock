import pandas as pd
import scipy.special
import numpy as np



def calculate_Bio_Age_correction(df, age_column='Bio_Age'):
    '''

    Calculate the second correction factor for every biological age of every sample (row) in df.


    :param df: Pandas DataFrame with a row for each Sample.
    Columns contain Genes, one column with the biological age, as well as possible further meta-data,i.e.
    the Strain, Treatment, RNAi, GEO accession number
    :return: A copy of df with the updated/corrected biological age
    '''
    df_corrected = df.copy()
    bio = df_corrected[age_column].values
    median = 372 # 15.5 days
    sd = (192 / 3) # 8/3 days

    # scipy.special.erf(((median - y) / sd) / np.sqrt(2)) calculates the number of SDs that the biological age y is
    # away from the median.
    # This value divided by 2 is the area that remains on the right side of the median.
    # We add 0.5 for the left side (50%) of the median to get the (approximated) total remaining area under the curve
    # To approximate at which percentage the new median of this truncated bell curve would be,
    # we divide the number by 2, in total:
    # (((scipy.special.erf(((median - y) / sd) / np.sqrt(2)) / 2) + 0.5) / 2)

    # If the biological age would be the median, i.e. y=372 hours, the right half of the bell curve would be trimmed.
    # This would leave just the left side of the curve, for which the new median can be approximated by the value
    # at which the area under the curve is 50%/2=25%.
    #
    # To calculate the biological age at which we get this percentage, we use the inverse error function
    # scipy.special.erfinv() to calculate how many SDs the new percentage value is apart from the original median at
    # 50%, i.e. scipy.special.erfinv(0.5 - X), where X is the approximated value explained above. This value needs to
    # be scaled to result in SDs by multiplying it with 2-times the square root of 2.
    # To calculate the new median, i.e. the rescaled biological age of the population,
    # we have to subtract the old median y by the number SDs we just calculated times 192/3,
    # which is the SD of the curve, i.e. y- numberOfSDs*sd

    corrected_bio_age = [(y - (scipy.special.erfinv(0.5 - (((scipy.special.erf(((median - y) / sd) /
                                                                  np.sqrt(2)) / 2) + 0.5) / 2)) * np.sqrt(2) * 2 * sd))
            for y in bio]
    df_corrected['Corrected_Biological_Age'] = corrected_bio_age

    return df_corrected


def test_correction():
    df = pd.read_csv('bio_age_example.csv')
    corrected_df = calculate_Bio_Age_correction(df)
    print(corrected_df)
