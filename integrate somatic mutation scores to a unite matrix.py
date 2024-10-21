import pandas as pd

# Cancer weight dict
cancer_weights = {
    'ACC':0.005,
    'BLCA':0.031,
    'BRCA':0.116,
    'CESC':0.033,
    'CHOL':0.022,
    'COAD':0.048,
    'ESCA':0.026,
    'GBM':0.004,
    'HNSC':0.004,
    'KICH':0.006,
    'KIRC':0.006,
    'KIRP':0.005,
    'LGG':0.004,
    'LIHC':0.021,
    'LUAD':0.062,
    'LUSC':0.062,
    'MESO':0.002,
    'OV':0.016,
    'PAAD':0.026,
    'PCPG':0.004,
    'PRAD': 0.073,
    'READ':0.048,
    'SARC':0.002,
    'STAD':0.049,
    'TGCT':0.004,
    'THCA':0.041,
    'UCEC':0.021,
    'UCS':0.033,
    'UVM':0.017

}

# 初始化一个空的DataFrame来汇总数据
combined_scores = pd.DataFrame()

# 遍历每种癌症及其权重
for cancer_type, weight in cancer_weights.items():
    # 读取对应的癌症基因突变分数文件，指定分隔符为逗号
    # file_path = f'results/{cancer_type}_gene_score.csv'  # 12 cancers
    file_path = f'results/pancancers/{cancer_type}_gene_score.csv'
    mutation_scores = pd.read_csv(file_path, sep=',',header=None,names=['Gene', 'Score'])

    # 应用权重
    mutation_scores['WeightedScore'] = mutation_scores['Score'] * weight

    # 如果combined_scores为空，直接赋值
    if combined_scores.empty:
        combined_scores = mutation_scores[['Gene', 'WeightedScore']].rename(columns={'WeightedScore': cancer_type})
    else:
        # 否则，合并到总表中
        mutation_scores = mutation_scores[['Gene', 'WeightedScore']].rename(columns={'WeightedScore': cancer_type})
        combined_scores = pd.merge(combined_scores, mutation_scores, on='Gene', how='outer')

# 替换所有NaN为0，因为有些基因可能在某些癌症中没有突变分数
combined_scores.fillna(0, inplace=True)

# 汇总所有癌症的权重分数
combined_scores['TotalScore'] = combined_scores.drop('Gene', axis=1).sum(axis=1)
# 对combined_scores按照TotalScore进行降序排序
sorted_combined_scores = combined_scores.sort_values(by='TotalScore', ascending=False)

# print(sorted_combined_scores.head())  # 打印前几行以检查结果
sorted_combined_scores.to_csv('results/Pancancer_Combined_Mutation_Scores.csv', index=False)  # 保存到CSV文件