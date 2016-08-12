import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import scipy
from scipy import stats
%matplotlib inline
#read in data
motif_position_df = pd.read_csv('/home/jtao/for_shengnan/motif_start_frame_C57BL6J.tsv', sep='\t')
motif_position_df.index = motif_position_df['ID'].values
del motif_position_df['ID']
#seperate data into veh only, kla only, and unchanged
veh_indices=motif_position_df[motif_position_df['Factors'].str.contains('atac_veh')].index.values
kla_incices=motif_position_df[motif_position_df['Factors'].str.contains('atac_kla')].index.values
veh_indices=set(veh_indices)
kla_incices=set(kla_incices)
veh_only=veh_indices-kla_incices
kla_only=kla_incices-veh_indices
unchanged=veh_indices.intersection(kla_incices)
veh_only=np.array(list(veh_only))
kla_only=np.array(list(kla_only))
unchanged=np.array(list(unchanged))
veh_only_df=motif_position_df.loc[motif_position_df.index.isin(veh_only)]
del veh_only_df['Factors']
del veh_only_df['chr']
kla_only_df=motif_position_df.loc[motif_position_df.index.isin(kla_only)]
del kla_only_df['Factors']
del kla_only_df['chr']
unchanged_df=motif_position_df.loc[motif_position_df.index.isin(unchanged)]
del unchanged_df['Factors']
del unchanged_df['chr']
overall_motif_df=motif_position_df
del overall_motif_df['Factors']
del overall_motif_df['chr']
def co_occur(df):
    '''
    input is a dataframe contains all data of the genomic position of motif in chr1.
    output is a dataframe contains the frequency of co-occurance of every two motifs. 
    '''
    motifs = df.columns.values
    count_frame=np.zeros((df.shape[1],df.shape[1]),dtype=np.int)#a dataframe to store future data
    count_frame=pd.DataFrame(count_frame, columns=motifs)
    count_frame.index=motifs
    col_vector=np.zeros((df.shape[0],1),dtype=np.int) #a column vector use in loop
    for i in range (df.shape[1]):#use the column vector to store each column one by one
        col_vector=df.ix[:,i]!= -1# Find if the motifs exist
        for j in range (df.shape[1]): 
            log_v=df.ix[:,j]!= -1#Find if the motifs exist
            vector_sum = 1*col_vector + 1*log_v
            log_input = vector_sum == 2
            count_input=np.sum(1*log_input)#convert logical ou
            count_frame.ix[i,j]=count_frame.ix[i,j]+count_input
            np.fill_diagonal(count_frame.values, 0)# change diagonal back to zero
    return count_frame
veh_cooccurence_df=co_occur(veh_only_df)
kla_cooccurence_df=co_occur(kla_only_df)
unchanged_cooccurence_df=co_occur(unchanged_df)
overall_cooccurence_df=co_occur(overall_motif_df)
#plot veh co-occurence
veh_cooccurence_plot=[]
for i in range (veh_cooccurence_df.shape[1]-1):
    for j in range (i+1,veh_cooccurence_df.shape[1]):
        veh_cooccurence_plot.append(veh_cooccurence_df.ix[i,j])
sns.distplot(veh_cooccurence_plot)
plt.ylabel('Frequency')
plt.xlabel('co_occurence')
plt.title('veh') 
#plot kla co-occurence
kla_cooccurence_plot=[]
for i in range (kla_cooccurence_df.shape[1]-1):
    for j in range (i+1,kla_cooccurence_df.shape[1]):
        kla_cooccurence_plot.append(kla_cooccurence_df.ix[i,j])
sns.distplot(kla_cooccurence_plot)
plt.ylabel('Frequency')
plt.xlabel('co_occurence')
plt.title('kla') 
#plot unchanged co-occurence
unchanged_cooccurence_plot=[]
for i in range (unchanged_cooccurence_df.shape[1]-1):
    for j in range (i+1,unchanged_cooccurence_df.shape[1]):
        unchanged_cooccurence_plot.append(unchanged_cooccurence_df.ix[i,j])
sns.distplot(unchanged_cooccurence_plot)
plt.ylabel('Frequency')
plt.xlabel('co_occurence')
plt.title('unchanged') 
#motifpair occurence
def motifpair_occurence(df):
    '''
    input:motfis position dataframe
    outpit:a dataframe countain occurence count of each motif pair.
    '''
    motifs = df.columns.values
    motifs_occurence=sum(df.values!=-1)#a list contains occurence count of each motifs
    count_occurence=np.zeros((df.shape[1],df.shape[1]),dtype=np.int)#a dataframe to store future data
    count_occurence=pd.DataFrame(count_occurence, columns=motifs)
    count_occurence.index=motifs
    for i in range (len(motifs_occurence)):#use the column vector to store each column one by one
        for j in range (len(motifs_occurence)): 
            count_input=motifs_occurence[i]*motifs_occurence[j]
            count_occurence.ix[i,j]=count_input
            np.fill_diagonal(count_occurence.values, 1)# change diagonal, one with oneslef to one
    return count_occurence
#normazlie occurence by motifpair occurence
veh_motif_occurence=motifpair_occurence(veh_only_df)
veh_normalized_cooccurence=veh_cooccurence_df/veh_motif_occurence
kla_motif_occurence=motifpair_occurence(kla_only_df)
kla_normalized_cooccurence=kla_cooccurence_df/kla_motif_occurence
unchanged_motif_occurence=motifpair_occurence(unchanged_df)
unchanged_normalized_cooccurence=unchanged_cooccurence_df/unchanged_motif_occurence
#plot veh co-occurence
veh_normalized_cooccurence_plot=[]
for i in range (veh_normalized_cooccurence.shape[1]-1):
    for j in range (i+1,veh_normalized_cooccurence.shape[1]):
        veh_normalized_cooccurence_plot.append(veh_normalized_cooccurence.ix[i,j])
sns.distplot(veh_normalized_cooccurence_plot)
plt.ylabel('Frequency')
plt.xlabel('normalized_occurence')
plt.title('veh')
plt.xlim(0.00050, 0.00060)
#plot kla co-occurence
kla_normalized_cooccurence_plot=[]
for i in range (kla_normalized_cooccurence.shape[1]-1):
    for j in range (i+1,kla_normalized_cooccurence.shape[1]):
        kla_normalized_cooccurence_plot.append(kla_normalized_cooccurence.ix[i,j])
sns.distplot(kla_normalized_cooccurence_plot, bins=50)
plt.ylabel('Frequency')
plt.xlabel('normalized_occurence')
plt.title('kla')
plt.xlim(0.000034, 0.000036)
#plot unchanged co-occurence
unchanged_normalized_cooccurence_plot=[]
for i in range (unchanged_normalized_cooccurence.shape[1]-1):
    for j in range (i+1,unchanged_normalized_cooccurence.shape[1]):
        unchanged_normalized_cooccurence_plot.append(unchanged_normalized_cooccurence.ix[i,j])
sns.distplot(unchanged_normalized_cooccurence_plot, bins=50)
plt.ylabel('Frequency')
plt.xlabel('normalized_occurence')
plt.title('unchanged')
def Find_Z(n):
    '''
    For all pairs of motifs - is there a pair that co-occurs more or less often than you would expect.
    input: a dataframe contains frequency of co-occurence of every pair of motifs.
    output:a dataframe contains z score of each pair of motifs.
    '''
    motifs = n.columns.values
    #convert dataframe to matirx
    n=n.as_matrix(columns=None)
    z_matrix = np.zeros((n.shape[0],n.shape[1]-1),dtype=np.float)
    for i in range (n.shape[0]):
        co_motif = n[i,:]
        co_motif = np.delete(co_motif,i)#remove data of the motif co-occur with itself
        z_score=stats.zscore(co_motif)# find z socre 
        z_matrix[i,:]=z_score
    #convert z score matirx to dataframe
    zscore_all=np.zeros((z_matrix.shape[0],z_matrix.shape[0]),dtype=np.float)
    for i in range (z_matrix.shape[0]):
        z_motif_self=z_matrix[i,:]
        z_motif_self=np.insert(z_motif_self,i,100)
        zscore_all[i,:]=z_motif_self
    zscore_frame = pd.DataFrame(zscore_all, columns=motifs)
    zscore_frame.index = motifs
    return zscore_frame
# z socre of normalized cooccurence
veh_z_normalized=Find_Z(veh_normalized_cooccurence)
kla_z_normalized=Find_Z(kla_normalized_cooccurence)
unchanged_z_normalized=Find_Z(unchanged_normalized_cooccurence)
veh_nz_plot=veh_z_normalized.values
kla_nz_plot=kla_z_normalized.values
unchanged_nz_plot=unchanged_z_normalized.values
veh_nz_plot=veh_nz_plot.flatten()
kla_nz_plot=kla_nz_plot.flatten()
unchanged_nz_plot=unchanged_nz_plot.flatten()
#plot overall z score of normalized cooccurence
sns.distplot(veh_nz_plot[veh_nz_plot!=100])
plt.ylabel('Frequency')
plt.xlabel('z score')
plt.title('veh')
plt.show()
sns.distplot(kla_nz_plot[kla_nz_plot!=100])
plt.ylabel('Frequency')
plt.xlabel('z score')
plt.title('kla')
plt.show()
sns.distplot(unchanged_nz_plot[unchanged_nz_plot!=100])
plt.ylabel('Frequency')
plt.xlabel('z score')
plt.title('unchanged')
#make truth table
def make_truth_table(df):
    '''
    input: a 196*196 z score dataframe
    '''
    motifs = df.columns.values
    Truth_frame=np.zeros((df.shape[1],df.shape[1]),dtype=np.int)#a dataframe to store future data
    Truth_frame=pd.DataFrame(Truth_frame, columns=motifs)
    Truth_frame.index=motifs
    for i in range (df.shape[1]):
        for j in range (df.shape[1]):
            if motifs[i]==motifs[j]:
                Truth_frame.ix[i,j]=100
            elif abs(df.ix[i,j])>=2.5 and abs(df.ix[j,i])>=2.5:
                Truth_frame.ix[i,j]=1
            else:
                Truth_frame.ix[i,j]=0
    #reshape dataframe
    Pairs=[]
    Truth=[]
    #loop in part of count data that contain meaning counting
    for i in range (Truth_frame.shape[1]-1):
        for j in range (i+1,Truth_frame.shape[1]):
            #put motif pair and correlation into the empty list
            motif_pairs=(motifs[i],motifs[j])
            Pairs.append(motif_pairs)
            Truth.append(Truth_frame.ix[i,j])
    #reshape the dataframe
    reshaped_frame = pd.DataFrame({'Truth': Truth}, index=Pairs)
    return reshaped_frame
#Truth table of veh,remove the entry of one motif with itself.
veh_truth=make_truth_table(veh_z_normalized)
veh_truth=veh_truth[veh_truth.Truth !=100]
#Truth table of kla,remove the entry of one motif with itself.
kla_truth=make_truth_table(kla_z_normalized)
kla_truth=kla_truth[kla_truth.Truth !=100]
#Truth table of unchanged,remove the entry of one motif with itself.
unchanged_truth=make_truth_table(unchanged_z_normalized)
unchanged_truth=unchanged_truth[unchanged_truth.Truth !=100]
Truth_table_cooccurence=pd.concat([veh_truth, kla_truth,unchanged_truth], axis=1)
Truth_table_cooccurence.columns = ['veh', 'kla','unchanged']
Truth_table_cooccurence.to_csv('/home/shs038/veh_kla/Truth_table_cooccurence.tsv', sep='\t')
