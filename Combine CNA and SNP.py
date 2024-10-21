import pandas as pd

# Reading the file into a DataFrame
cna_path = 'data/raw data/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes'
cna_frame = pd.read_csv(cna_path,sep='\t', index_col=0, header=0)

svn_path = 'data/raw data/mc3.v0.2.8.PUBLIC.nonsilentGene.xena'
svn_frame = pd.read_csv(svn_path,sep='\t', index_col=0, header=0)

# Selecting common columns and rows
common_columns = cna_frame.columns.intersection(svn_frame.columns)
common_rows = cna_frame.index.intersection(svn_frame.index)

selected_cna_frame = cna_frame.loc[common_rows, common_columns]
selected_svn_frame = svn_frame.loc[common_rows, common_columns]
selected_svn_frame.to_csv("../data/source data/svn_somatic_mutation_matrix.txt", header=True)

# print(selected_cna_frame.head(5))
# print(selected_svn_frame.head(5))

selected_cna_frame.replace(2.0,1.0,inplace = True)
selected_cna_frame.replace(-2.0,1.0,inplace = True)
selected_cna_frame.replace(-1.0,1.0,inplace = True)
selected_cna_frame.to_csv("../data/source data/cna_somatic_mutation_matrix.txt", header=True)

# print(selected_cna_frame)

# add the tow dataframe
combined_frame = selected_cna_frame + selected_svn_frame
# replace 2 to 1
combined_frame.replace(2.0,1.0,inplace = True)
print(combined_frame)

def count_occurrences(df,num):
    count = (df == num).sum().sum()
    print(f"Count of occurrences of {num} in the dataframe:", count)
# # Counting occurrences of 2.0 in the dataframe
# count_2 = (combined_frame == 2.0).sum().sum()
#
# # Counting occurrences of 1.0 in the dataframe
# count_1 = (combined_frame == 1.0).sum().sum()
#
# # Displaying the count
# print("Count of occurrences of 2.0 in the dataframe:", count_2)
# print("Count of occurrences of 1.0 in the dataframe:", count_1)
def remove_null_rows_and_columns(df):
    # Removing rows where all values are 0.0
    df_rm_rows = df.loc[(df != 0).any(axis=1)]
    # Removing columns where all values are 0.0
    df_rm_rows_cols = df_rm_rows.loc[:,(df_rm_rows != 0).any(axis=0)]
    return df_rm_rows_cols
df_filtered = remove_null_rows_and_columns(combined_frame)

# write the result to file
somatic_mutation_file = '..\data\source data\somatic mutation data.txt'
df_filtered.to_csv(somatic_mutation_file, sep='\t', index=True)

