import pandas as pd
import matplotlib.pyplot as plt
import os
import utilities

def get_common_gene_samples():
    # 读取基因突变文件,包含gene和samples的信息
    gene_mut = pd.read_csv('data/source data/somatic mutation data.txt', sep='\t', index_col=0)

    # load cancer-sample file
    cancer_samples = pd.read_csv('results/files/cancer_sample_type.csv')

    # 公共的基因和样本
    common_genes = gene_mut.index
    common_samples = list(set(gene_mut.columns) & set(cancer_samples['sample']))

    common_gene_mut_data = gene_mut.loc[common_genes,common_samples]

    # 根据公共样本修改病理数据
    common_cancer_samples = cancer_samples[cancer_samples['sample'].isin(common_samples)]

    # 保持文件
    common_gene_mut_data.to_csv('data/source data/common_gene_mutation.csv', sep='\t')
    common_cancer_samples.to_csv('data/source data/common_cancer_samples.csv', sep='\t')

def get_samples(sample_file):
    # 读取TSV文件
    samples = pd.read_csv(sample_file, sep='\t', header=0)

    # 选择需要的两列数据
    cancer_samples = samples[['cancer type abbreviation','sample']]

    # 将新表写入新的文件
    cancer_samples.to_csv('results/files/cancer_sample_type.csv', index=False)

    # 统计Cancer的重复值
    cancers_counts = samples['cancer type abbreviation'].value_counts(sort=False)

    # 将结果写入新的文件
    cancers_counts.to_csv('data/source data/cancers_counts.csv', index=True, index_label='cancer type',
                          header='count')

    return samples,cancer_samples,cancers_counts


# 比较两列的数值是否完全相同
# if (df['sample'] == df['samples']).all():
#     print("The two columns have the same values.")
# else:
#     print("The two columns have different values.")


def draw_bar():
    # 读取TSV文件
    samples = pd.read_csv('data/source data/common_cancer_samples.csv', sep='\t', header=0, index_col=None)

    # 统计Cancer的重复值
    cancers_counts = samples['cancer type abbreviation'].value_counts(sort=False)

    # 定义颜色列表
    #colors = ["#BCBD46", "#E8E8B9", "#86D3DE", "#A5DFE7", "#959595", "#D9D9D9"]
    colors = ["#E99C93", "#9FACD3", "#F1DBB9", "#D9D1E3", "#CAD4E7", "#D0E2E8"]

    # 绘制柱状图
    plt.figure(figsize=(10, 6))
    cancers_counts.plot(x='cancer type', y='count', kind='bar')
    ax = cancers_counts.plot(x='cancer type', y='count', kind='bar', color=colors)
    plt.title('Number of Every Cancer')
    plt.xlabel('cancer type')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig('results/figures/cancers_counts.png')
    plt.savefig('results/figures/cancers_counts.pdf', format='pdf')
    plt.show()

def get_cancer_sample_dict():

    cancer_samples = pd.read_csv('data/source data/common_cancer_samples.csv', sep='\t', header=0, index_col=None)

    # 读取mutation文件,包含gene和samples的信息
    gene_samples = pd.read_csv('data/source data/common_gene_mutation.csv', index_col=0,sep='\t')

    # 创建一个字典,用于存储每种cancer对应的samples
    cancer_sample_dict = {}
    for cancer, samples in zip(cancer_samples['cancer type abbreviation'], cancer_samples['sample']):
        if cancer not in cancer_sample_dict:
            cancer_sample_dict[cancer] = []
        cancer_sample_dict[cancer].append(samples)

    # 根据cancer将gene数据分开,并保存到不同的文件
    for cancer, samples in cancer_sample_dict.items():
        # 选择与当前cancer相关的samples列
        cancer_genes = gene_samples[gene_samples.columns.intersection(samples)].T
        cancer_genes.to_csv(base_dir + f'/{cancer}_somatic_mutation.csv', sep='\t')
        # generate_cancer_file(cancer, cancer_genes.T, base_dir)

    # Write the cancer_sample_dict to a file
    with open('data/source data/cancer_sample_dict.txt', 'w') as file:
        for cancer, samples in cancer_sample_dict.items():
            file.write(f"{cancer}: {', '.join(samples)}\n")

def check_if_duplicated(data):
    if len(data) == len(set(data)):
        print('列表里的元素互不重复！')
    else:
        print('列表里有重复的元素！')

def generate_cancer_file(cancer,cancer_genes,base_dir):


    # 创建文件夹
    # folder_path = os.path.join(base_dir, cancer)
    # if not os.path.exists(folder_path):
    #     os.makedirs(folder_path)

    # 保存到文件
    cancer_genes.to_csv(base_dir+f'/{cancer}_somatic_mutation.csv', index=True, sep='\t')

    #处理全0行
    # utilities.select_not_all_zero(base_dir+f'/{cancer}/', cancer)

def generate_gene_express_file(base_dir):

    # 读取公共基因表达文件
    gene_expr = pd.read_table('../data/source data/common_gene_expression.txt', sep='\t', index_col=0)
    # print(gene_expr.index)

    # 读取基因表达文件
    gene_expr_all = pd.read_table('../data/raw data/GeneExpression.xena', sep='\t', index_col=0)
    # print(gene_expr_all.index)
    gene_expr_all = gene_expr_all.loc[gene_expr.index]

    # print(set(gene_expr_all.index)-set(gene_expr.index))

    # if len(gene_expr_all.index) == len(set(gene_expr_all.index)):
    #     print('列表里的元素互不重复！')
    # else:
    #     print('列表里有重复的元素！')

    # 找出重复的索引值
    duplicates = gene_expr_all.index.duplicated(keep=False)  # keep=False标记所有重复项
    # 显示所有重复索引的行
    duplicate_rows = gene_expr_all[duplicates]
    print(duplicate_rows)

    # 删除重复的索引，保留第一个
    gene_expr_all = gene_expr_all[~gene_expr_all.index.duplicated(keep='first')]

    # # 找出df1中有但df2中没有的索引
    # unique_to_df1 = gene_expr.index.difference(gene_expr_all.index)
    # print("Indexes unique to df1:", unique_to_df1)

    # 读取cancer-samples字典文件
    with open('../data/source data/cancer_sample_dict.txt', 'r') as f:
        cancer_sample_dict = {}
        for line in f:
            cancer, samples = line.strip().split(': ')
            cancer_sample_dict[cancer] = samples.split(', ')

    # # 根据cancer和对应的samples,将基因表达数据分成不同的文件
    # for cancer, samples in cancer_sample_dict.items():
    #     cancer_gene_expr_data = gene_expr[samples]
    #     cancer_gene_expr_data.to_csv(base_dir+f'/{cancer}/{cancer}_gene_expression.txt', sep='\t')
    #     print(f"Saved {cancer} gene expression data to '{cancer}_gene_expression.txt'")

    # 读取差异表达基因数据
    with open('../data/source data/gene_express_cancer_normal.txt', 'r') as f:
        lines = [line.strip().split(',') for line in f.readlines()]
        normals = [x[1].strip() for x in lines]  # 提取第二列数据

    gene_diff_samples = list()
    for cancer, samples in cancer_sample_dict.items():
        for sample in samples:
            for normal in normals:
                if sample[:-3] == normal[:-3]:
                    if sample in gene_expr.columns and normal in gene_expr_all.columns:
                        gene_diff_samples.append(gene_expr[sample])
                        gene_diff_samples.append(gene_expr_all[normal])
        if len(gene_diff_samples) > 0:
            gene_diff_expr = pd.DataFrame(gene_diff_samples)
            gene_diff_expr.T.to_csv(base_dir+f'/{cancer}/{cancer}_gene_diff_expression.txt', sep='\t')
            print(f"Saved {cancer} gene expression data to '{cancer}_gene_diff_expression.txt'")

def del_gene_diff_expr_zero_rows():

    # 读取cancer-samples字典文件
    with open('../data/source data/cancer_sample_dict.txt', 'r') as f:
        cancer_sample_dict = {}
        for line in f:
            cancer, samples = line.strip().split(': ')
            cancer_sample_dict[cancer] = samples.split(', ')

    for cancer, samples in cancer_sample_dict.items():
        if cancer == 'TGCT':
            continue
        utilities.select_not_all_zero(base_dir+f'/{cancer}/', cancer, '_gene_diff_expression.txt',0)

def generate_gene_diff_expr_samples_file():
    # 读取cancer-samples字典文件
    with open('../data/source data/cancer_sample_dict.txt', 'r') as f:
        cancer_sample_dict = {}
        for line in f:
            cancer, samples = line.strip().split(': ')
            cancer_sample_dict[cancer] = samples.split(', ')

    for cancer, samples in cancer_sample_dict.items():
        if cancer == 'TGCT':
            continue
        with open(base_dir+f'/{cancer}/filtered_{cancer}_gene_diff_expression.txt', 'r') as f:
            first_line = f.readline().strip()

            # 将样本信息保存到另一个文件（例如 sample_names.txt）
            sample_names = first_line.split('\t')  # 样本信息以制表符分隔
            output_filename = base_dir+f'/{cancer}/filtered_{cancer}_kept_genes_diff_expression_samples.txt'

            with open(output_filename, 'w') as output_file:
                for sample_name in sample_names:
                    output_file.write(sample_name + '\n')



sample_file = 'data/raw/samples.tsv'
# get_samples(sample_file)
draw_bar()
# get_common_gene_samples()
# 指定要创建文件夹的目录
base_dir = 'data/source data/pancancers'
# get_cancer_sample_dict()
# geneExp_file = '../data/raw data/GeneExpression.xena'
#generate_gene_express_file(base_dir)

#del_gene_diff_expr_zero_rows()

# generate_gene_diff_expr_samples_file()

