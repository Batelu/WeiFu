# coding=utf-8
import numpy as np
import pandas as pd
from numpy import *
from scipy.io import savemat
import scipy.io as scio
import copy
import itertools
import sys
import csv
#from function import *
from docplex.mp.model import Model
from docplex import *
from pyecharts import options as opts
from pyecharts.charts import Graph, Page
import time
from scipy.io import savemat

def txt2mat(txt_file_address, mat_file_address,matName):
    txt_list = list()
    with open(txt_file_address,'r') as f:
        for line in f:
            line = line.strip()
            txt_list.append([line.replace("'","")])
    txt_list = np.array(txt_list, dtype=object)
    savemat(mat_file_address, {matName: txt_list})

def txt2mat_matrix(txt_file_address, mat_file_address,matName):
    matrix = np.loadtxt(txt_file_address)
    savemat(mat_file_address, {matName: matrix.T})
# generate pan-cancer genes
#txt2mat('source/Data/largesize/analyzed_genes_large.txt','source/Data/largesize/GENE_PAN.mat','GENE_PAN')


#txt2mat_matrix('source/Data/largesize/patient_data_gen_matrix.txt','source/Data/largesize/MU_PAN.mat','MU_PAN')

def generate_gene2pathways_matrix(matGeneAddress,matGeneName,matPathAddress,matPathName,g2pAddress):
    top2000gene_path = list()
    gene_mat_file = scio.loadmat(matGeneAddress)
    temp = np.array(gene_mat_file[matGeneName])
    top2000genes = list()
    for i in temp[0]:
        top2000genes.append(i-1)

    gene_path_mat_file = scio.loadmat(matPathAddress)
    gene_path = np.array(gene_path_mat_file[matPathName])
    #print(gene_path[6432])
    top2000gene_path = gene_path[top2000genes]

    np.savetxt(g2pAddress,np.c_[top2000gene_path],fmt='%d', delimiter='\t')

def generate_patient_index(patient_address):
    patient2index = dict()
    index2patient = dict()
    i = 0
    with open(patient_address,'r') as f:
        for line in f:
            line = line.split()
            patient2index[line[0]] = i
            index2patient[i] = line[0]
            i += 1

    with open("source/Data/smallsize/BRCA_patient_index.txt",'w') as f:
        for item in patient2index.items():
            f.write(str(item[0]) + "\t" + str(item[1]) + '\n')

    return patient2index,index2patient

def get_genes(gene_file_address):
    genes = list()
    gene2index = dict()
    index2gene = dict()
    with open(gene_file_address, 'r') as f:
        print("Load genes...")
        index = 0
        for line in f:
            v = line.strip().split()
            genes.append(v[0].replace("'", ""))
            gene2index[v[0].replace("'", "")] = index
            index2gene[index] = v[0].replace("'", "")
            index += 1

    return genes,gene2index

def generate_gene_patients_index(mutation_file):
    all_samples = set()
    sample_mutatedGenes = dict()

    gene_mutatedSamples = dict()
    unMutGene = list()  # 非突变基因集
    GeneMutName = list()  # 用来存放发生了突变的基因名
    GeneMutScore = dict()  # 放基因名+突变分数
    gene2index = dict()
    patient2index, index2patient = generate_patient_index(mutation_file)


    # BRCA
    genes_original,g2i_original = get_genes('source/Data/smallsize/BRCA_genes.txt')

    genes_later, g2i_later = get_genes('source/Data/largesize/analyzed_genes_large.txt')
    #BRCA与泛癌基因的并集
    genes_union = list(set(genes_original).intersection(set(genes_later)))

    for gene in genes_union:
        gene_mutatedSamples[gene.replace("'", "")] = set()

    sample_mut_f = open(mutation_file, 'r')
    for line in sample_mut_f:
        v = line.split()
        sampleID = v[0]
        all_samples.add(sampleID)
        sample_mutatedGenes[sampleID] = set()
        for i in range(len(v) - 1):
            gene = v[i + 1].strip()
            if gene in genes_union:
                GeneMutName.append(gene)
                gene_mutatedSamples[gene].add(patient2index[sampleID])
                sample_mutatedGenes[sampleID].add(gene)
    sample_mut_f.close()

    with open("source/Data/smallsize/BRCA_genes_patientindex.txt",'w') as f:
        for item in gene_mutatedSamples:
            f.write(item)
            for i in gene_mutatedSamples[item]:
                f.write("\t"+str(i))
            f.write('\n')
    # generate patient index file
    # with open("source/Data/smallsize/BRCA_patient_index.txt",'w') as f:
    #     i = 0
    #     for item in sample_mutatedGenes:
    #         if len(sample_mutatedGenes[item]) != 0:
    #             f.write(str(item[0]) + "\t" + str(i) + '\n')
    #             i += 1

    SampMutNum = dict()  # 样本的突变分数=1除以sMutNum
    MaxSampGeneNum = 0  # 所有样本中，突变基因的最大数
    for sample in all_samples:  # 样本数
        sMutNum = 0  # 每个样本的突变基因数
        for geneName in GeneMutName:
            if patient2index[sample] in gene_mutatedSamples[geneName]:
                sMutNum = sMutNum + 1  # 计算在该样本中发生突变的基因数
        if sMutNum != 0:
            SampMutNum[sample] = 1.0 / sMutNum
            if sMutNum > MaxSampGeneNum:
                MaxSampGeneNum = sMutNum
        else:
            SampMutNum[sample] = 0  # 没有突变基因的样本突变分数为0
    # ##计算基因突变分数
    for geneName in GeneMutName:  # 遍历至少在一个样本中突变的基因
        GMS = 0  # 基因突变分数
        samples = gene_mutatedSamples[geneName]  # 发生这个基因的样本集合
        for sample in samples:  # 找每一个样本的突变分数
            GMS = GMS + SampMutNum[index2patient[sample]]
        GeneMutScore[geneName] = GMS

    #gene_sums = sorted(GeneMutScore.items(), key=lambda d: d[1], reverse=True)  # 排序,会以tuple的形式返回

    with open('source/Data/smallsize/BRCA_gene2freq.txt', 'w') as f:
        for item in GeneMutScore:
            f.write(item + "\t"+ str(GeneMutScore[item]) + '\n')

def add_outcome_column(score_file, validated_genes_file, output_file):
    # 读取基因得分文件
    scores_df = pd.read_csv(score_file)

    # 读取已验证基因文件
    with open(validated_genes_file, 'r') as file:
        validated_genes = set(line.strip() for line in file)

    scores_df['TotalScore'] = scores_df['TotalScore']/30
    # 添加Outcome列
    scores_df['Outcome'] = scores_df['Gene'].apply(lambda gene: 'Yes' if gene in validated_genes else 'No')
    new_scores_df = scores_df[['Gene','Outcome','TotalScore']]

    # 保存新的文件
    new_scores_df.to_csv(output_file, index=False)

# score_file = 'results/pancancer_gene_scores_alias.csv'
score_file = 'results/other methods/MaxMIF_results.csv'
validated_genes_file = 'data/integrated_true_genes.txt'
output_file = 'results/other methods/MaxMIF_results_for_ROC.csv'
add_outcome_column(score_file, validated_genes_file, output_file)

#patient_inde = generate_patient_index("source/Data/smallsize/BRCA_patient2gene.txt")
#generate_gene_patients_index("source/Data/smallsize/BRCA_patient2gene.txt")
#generate_gene2pathways_matrix('source/Data/largesize/gene_top_2000U.mat','gene_top_2000','source/CDPMiner/GENE_PATH.mat','GENE_PATH','source/Data/largesize/top2000gene_U_pathway.txt')

#generate_gene2pathways_matrix('source/Data/largesize/gene_top.mat','gene_topU','source/CDPMiner/GENE_PATH.mat','GENE_PATH','source/Data/largesize/gene_pathwayU.txt')

def get_GeneID_Name():
    info_path = '..\data\\raw data\gene_info'
    df_gene_info = pd.read_csv(info_path, sep='\t', header=0, low_memory=False)
    column_names = df_gene_info.columns
    print(column_names)

    GeneID_Name = df_gene_info[df_gene_info['#tax_id'] == 9606][['GeneID', 'Symbol']]

    print(GeneID_Name)
    # write the result to file
    GeneID_Name_file = '..\data\source data\GeneID_Name.txt'
    GeneID_Name.to_csv(GeneID_Name_file, sep='\t',index=False)

#get_GeneID_Name()

"""
Select rows that are not all 0.0
"""
def select_not_all_zero(path,cancer,filename, index_col):
    # 读取文件
    gene_samples = pd.read_csv(path+f'{cancer}'+filename, sep='\t', header=0, index_col=index_col)

    # 记录所有的行号
    all_rows = gene_samples.index.tolist()

    # 删除全0.0的行
    gene_samples = gene_samples.loc[(gene_samples != 0.0).any(axis=1)]

   # 记录保留下来的行号
    kept_rows = gene_samples.index.tolist()

    # 计算删除的行号
    deleted_rows = [row for row in all_rows if row not in kept_rows]

    # 将结果保存到新的文件
    gene_samples.to_csv(path+'filtered_'+f'{cancer}'+filename, sep='\t', index=False)

    # 将保留的行号写入文件
    with open(path+f'{cancer}_kept_genes'+filename, 'w') as file:
        file.write('\n'.join(map(str, kept_rows)))
        print(f"{cancer} Kept row numbers saved to 'kept_rows.txt'")

    # 将删除的行号写入文件
    with open(path+f'{cancer}_deleted_genes'+filename, 'w') as file:
        file.write('\n'.join(map(str, deleted_rows)))
        print(f"{cancer} Deleted row numbers saved to 'deleted_rows.txt'")

"""
行索引和对应的行号保存成一个文件
"""
def save_gene_index(input_file, output_file):
    """
    将输入文件中的 gene 名称和对应的行号保存到输出文件中。

    参数:
    input_file (str): 输入文件路径,每行包含一个 gene 名称。
    output_file (str): 输出文件路径,每行格式为 "gene_name,line_number"。

    返回:
    None
    """
    # 读取输入文件,提取 gene 名称
    with open(input_file, 'r') as f:
        next(f)
        gene_names = [line.strip().split('\t')[0] for line in f]

    # 创建 gene 名称和行号的字典
    gene_index = {gene_name: i for i, gene_name in enumerate(gene_names)}

    # 将结果保存到输出文件
    with open(output_file, 'w') as f:
        for gene_name, index in gene_index.items():
            f.write(f"{gene_name},{index}\n")

# save_gene_index('../data/source data/common_gene_mutation.txt', '../data/source data/gene_index.txt')

