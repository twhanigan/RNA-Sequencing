# Protein level abundance associated with B819 Sensitivity
import numpy as np
import pandas as pd
import math
import numpy as d
import pandas as pd
import math
from pandas import ExcelWriter
from pandas import ExcelFile
import scipy
from scipy import stats
from scipy.stats import pearsonr
import warnings
from scipy.stats import ttest_ind
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics import jaccard_similarity_score
from more_itertools import unique_everseen
import numbers
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pyplot as pyplot
import matplotlib.animation as anim
from matplotlib.animation import FFMpegWriter
from matplotlib import style
import statistics

# Import Spreadsheets
Raw_Frame = pd.read_excel('joined.xlsx')
Raw_Frame.columns = ['gene', 'H460r1', 'H460r2', 'H460r3', 'H1975r1', 'H1975r2', 'H1975r3', 'MDAMB231r1', 'MDAMB231r2', 'MDAMB231r3', 'PC3r1', 'PC3r2', 'PC3r3', 'KP4r1', 'KP4r2',
                     'KP4r3', 'HCT116r1', 'HCT116r2', 'HCT116r3', 'A549r1', 'A549r2', 'A549r3', 'H2122r1', 'H2122r2', 'H2122r3', 'H2009r1', 'H2009r2', 'H2009r3', 'SW620r1', 'SW620r2', 'SW620r3']
# Drop rows with zero value
Raw_Frame = Raw_Frame.replace(0, np.nan).copy()
Raw_Frame = Raw_Frame.dropna(how='any', axis=0).copy()

# log2 transform
Log2_Frame = np.log2(Raw_Frame.loc[:, 'H460r1':'SW620r3']).copy()
# Get Protein Names from Description
Log2_Frame['Name'] = Raw_Frame['gene']
Log2_Frame.to_excel('Log2_Frame+B.xlsx')


Median_Norm = Log2_Frame.loc[:, 'H460r1':'SW620r3']
Median_Norm['Name'] = Log2_Frame['Name']

# Average Replicate Values and make same order as activity matrix
Median_Norm['H460'] = Median_Norm[['H460r1', 'H460r2', 'H460r3']].mean(axis=1)
Median_Norm['H1975'] = Median_Norm[[
    'H1975r1', 'H1975r2', 'H1975r3']].mean(axis=1)
Median_Norm['MDAMB231'] = Median_Norm[[
    'MDAMB231r1', 'MDAMB231r2', 'MDAMB231r3']].mean(axis=1)
Median_Norm['PC3'] = Median_Norm[['PC3r1', 'PC3r2', 'PC3r3']].mean(axis=1)
Median_Norm['KP4'] = Median_Norm[['KP4r1', 'KP4r2', 'KP4r3']].mean(axis=1)
Median_Norm['HCT116'] = Median_Norm[[
    'HCT116r1', 'HCT116r2', 'HCT116r3']].mean(axis=1)
Median_Norm['H2122'] = Median_Norm[[
    'H2122r1', 'H2122r2', 'H2122r3']].mean(axis=1)
Median_Norm['A549'] = Median_Norm[['A549r1', 'A549r2', 'A549r3']].mean(axis=1)
Median_Norm['H2009'] = Median_Norm[[
    'H2009r1', 'H2009r2', 'H2009r3']].mean(axis=1)
Median_Norm['SW620'] = Median_Norm[[
    'SW620r1', 'SW620r2', 'SW620r3']].mean(axis=1)
Median_Norm_Final = Median_Norm[['Name', 'KP4', 'MDAMB231', 'H2009',
                                 'H2122', 'PC3', 'H460', 'H1975', 'SW620', 'HCT116', 'A549']].reset_index()
Median_Norm_Final.to_excel('Median_Norm_B.xlsx')


Activity_Matrix = pd.read_csv('Activity_Matrix_B508.csv')
column = list(Activity_Matrix.columns)
column.append('Name')


Expression_Frame = Median_Norm_Final.copy()

Expression_Frame_Subset = Median_Norm.loc[:, 'KP4':'A549'].copy()
keys = Expression_Frame_Subset.index

# Make new activity matrix with same number of rows as your expression matrix
len_expression = len(Expression_Frame_Subset)
Activity_Matrix_New = []
Activity_Matrix_New = pd.DataFrame(
    pd.concat([Activity_Matrix]*(len_expression)))
Activity_Matrix_New.index = keys

# Replace null values with zero
Expression_Frame_Subset = Expression_Frame_Subset.dropna()
# print(Expression_Frame_Subset)

# Calculate correlation between Activity and Expresion
Activity_Expression_Correlation = Expression_Frame.loc[:, 'KP4':'A549'].corrwith(
    Activity_Matrix_New, axis=1)

# Calculate P-Value between groups
correlation = []
count = -1
for row in Expression_Frame.iterrows():
    count = count + 1
    correlation.append(ttest_ind(
        Expression_Frame.loc[count, 'KP4':'PC3'], Expression_Frame.loc[count, 'H460':'A549'], axis=0))
corr_df = pd.DataFrame(correlation)

Final_Table = pd.concat([Expression_Frame['Name'],
                         Activity_Expression_Correlation, corr_df], axis=1).dropna()
Final_Table.to_excel('Final_Table_B.xlsx')

#concatenated = pd.concat(Expression_Frame['KP4':'A549'],'')
ylabel = Expression_Frame.loc[0:, 'Name'].copy()
xlabel = ['KP4',	'MDAMB231',	'H2009',	'H2122',	'PC3',
          'H460',	'H1975',	'SW620',	'HCT116',	'A549']
g = sns.clustermap(Expression_Frame_Subset, cmap="RdYlBu", col_cluster=False, xticklabels=xlabel, figsize=(
    7, 7), method='average', robust=True, z_score=0)  # col_cluster=False #standard_scale=1, metric = 'correlation'
#f=sns.clustermap(Final_table, cmap='RdYlBu',col_cluster=False, yticklabels=Final_table_names,xticklabels=xlabel, method = 'centroid',robust = True,z_score = 0)
plt.show()

# plt.show()
