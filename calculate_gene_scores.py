import glob
import os
import pandas as pd
import numpy as np


def calculate_gene_scores(file_name, output_path,cancer_type):

    # Read the DataFrame
    # mutation_data = pd.read_csv(file_name, index_col=0) # 12 cancers
    mutation_data = pd.read_csv(file_name, sep='\t', index_col=0)
    # 计算每个样本中的突变基因总数m
    mutation_counts = mutation_data.sum(axis=1)

    # 计算每个sample的突变率（1/m），只对该基因在某样本中突变的情况计算
    sample_scores = 1/mutation_counts

    # Get the number of genes (columns in the mutation matrix)
    num_genes = mutation_data.shape[1]

    # Initialize an array to store gene mutation scores
    gene_mutation_scores = np.zeros(num_genes)

    # Multiply the mutation matrix by sample mutation scores
    weighted_matrix = mutation_data.multiply(sample_scores, axis=0)

    # Sum the weighted matrix across patients (axis=0) to get gene mutation scores
    gene_score = weighted_matrix.sum(axis=0)

    # 最大的m值
    max_m = mutation_counts.max()

    # 对于没有突变的基因，分数设置为1/(最大的m)
    gene_score[gene_score == 0] = 1 / max_m

    # 对基因得分进行排序
    sorted_gene_scores = gene_score.sort_values(ascending=False)

    output_filename = os.path.join(output_path, f'{cancer_type}_gene_score.csv')
    # 如果提供了输出文件名，将结果写入CSV文件
    if output_filename:
        sorted_gene_scores.to_csv(output_filename, header=False)

    return sorted_gene_scores

def pan_cancer_gene_score(cancers_path,output_path):
    # Define the pattern to match the cancer gene-patient matrix files
    # pattern = os.path.join(cancers_path, '*_gene_patient_matrix.csv')
    pattern = os.path.join(cancers_path, '*_somatic_mutation.csv')


    # Use glob to get a list of all matching file names
    file_names = glob.glob(pattern)

    # Loop through each file name and read the corresponding DataFrame
    for file_name in file_names:
        # Extract the cancer type from the file name (assuming the format '<cancer_type>_gene_patient_matrix.csv')
        cancer_type = os.path.basename(file_name).replace('_somatic_mutation.csv', '')
        # replace the backslashes (\\) with forward slashes (/)
        file_name = file_name.replace('\\', '/')
        calculate_gene_scores(file_name,output_path,cancer_type)


# 调用函数示例
cancers_path = "data/source data/pancancers"
output_path = "results/pancancers"
pan_cancer_gene_score(cancers_path,output_path)
