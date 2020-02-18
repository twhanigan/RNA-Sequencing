#Parse and compare RNASeq data
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
from sklearn.cluster import KMeans
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
import itertools
from pathlib import Path
import os
import glob
import functools
#Import data

# Step 1: get a list of all csv files in target directory
my_dir = "E:\\NGS\RNASeq\\"
filelist = []
filesList = []
os.chdir( my_dir )

files = glob.glob(os.path.join('E:\\NGS\RNASeq\*.txt'))
list_of_dfs = [pd.read_csv(f, sep = '\t', engine="python", skipfooter=1) for f in files]

col_names = ['gene', 'read']
for dataframe in list_of_dfs:
	dataframe.columns = col_names
	dataframe.set_index('gene')

#combined_df=pd.concat(list_of_dfs, axis = 1, sort = False)
joined_df = functools.reduce(lambda x,y: pd.merge(x,y, on = 'gene'),list_of_dfs)

joined_df.to_excel('joined.xlsx')

style.use('seaborn-muted')
joined_map = joined_df.iloc[0:(len(joined_df)-4),1:].copy()
joined_map = joined_map.replace(0,np.nan).copy()
joined_map = joined_map.dropna(how= 'any', axis = 0).copy()
joined_map = joined_map.astype(float).reset_index(drop=True).copy()
print(joined_map)
table_names = joined_df['gene']

xlabel = ['R1','R2','R3','R4','R5','R6']

#g=sns.clustermap(joined_map,cmap = "RdYlBu", col_cluster= False,yticklabels=table_names, xticklabels = xlabel,figsize = (7,7), method = 'average', robust = True,z_score = 0) #col_cluster=False #standard_scale=1, metric = 'correlation' 
#f=sns.clustermap(Final_table, cmap='RdYlBu',col_cluster=False, yticklabels=Final_table_names,xticklabels=xlabel, method = 'centroid',robust = True,z_score = 0)
#plt.show()