#!/usr/bin/env python3              # 指定使用Python 3解释器
import pandas as pd                 # 数据处理和分析库
import re                           # 正则表达式库
import concurrent.futures           # 并发执行库（多线程/多进程）
import os                           # 操作系统接口库
import logging                      # 日志记录库
from tqdm import tqdm               # 进度条显示库

# 配置日志系统
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# 常量设置
PEPTIDE_LENGTH_MIN = 3              # 生成肽段的最小长度
PEPTIDE_LENGTH_MAX = 7              # 生成肽段的最大长度
MAX_WORKERS = 25                    # 并行处理的最大线程数
OUTPUT_DIR = "Results"              # 输出结果目录

def load_data():
    """从CSV文件加载三组数据"""
    try:
        # 加载蛋白质结构域信息（gzip压缩格式）
        protein_domains = pd.read_csv('D:/pyitem/step2/CS_Matrikine_Discovery-main/Recourses/Domain_Info.csv.gz', index_col=False, compression='gzip')
        # 加载蛋白酶切割位点信息
        protein_protease_sites = pd.read_csv('D:/pyitem/step2/CS_Matrikine_Discovery-main/Recourses/Prosper.csv.gz', index_col=False, compression='gzip')
        # 加载蛋白质基本信息（含序列）
        protein_info = pd.read_csv('D:/pyitem/step2/CS_Matrikine_Discovery-main/Recourses/Summary_MSP_E.csv', index_col=False)
        return protein_domains, protein_protease_sites, protein_info
    except FileNotFoundError as e:
        logging.error(f"Error loading data: {e}")
        raise

def find_indices(data_list, target_value):
    """在列表中查找目标值出现的所有位置索引"""
    return [i for i, value in enumerate(data_list) if value == target_value]

def get_protein_protease_cleavage_sites(protein_protease_sites, gene_name):
    """获取指定基因的蛋白酶切割位点"""
    positions = find_indices(protein_protease_sites.iloc[:, 0].tolist(), gene_name)
    if not positions:
        return pd.DataFrame()  # 返回空DataFrame如果没有找到
    return protein_protease_sites.iloc[positions, :]

def get_protein_sequence(protein_info, gene_name):
    """获取指定基因的蛋白质序列"""
    positions = find_indices(protein_info.iloc[:, 0].tolist(), gene_name)
    if not positions:
        logging.warning(f"⚠️ Gene Not Found: {gene_name}")
        return ''  # 返回空字符串如果没有找到
    return protein_info.iloc[positions[0], 24]  # 第25列是序列列

def get_protein_domains(protein_domains, gene_name):
    """获取指定基因的结构域信息"""
    return protein_domains[protein_domains.iloc[:, 0] == gene_name]

def find_substring_index(sequence, substring):
    """在序列中查找子串的位置"""
    try:
        return sequence.index(substring)  # 返回起始索引
    except ValueError:
        return -1  # 找不到返回-1

def process_peptides_matrix(input_df, gene_name, protein_seq):
    """从切割位点生成肽段并筛选长度(3-7个氨基酸)"""
    peptide_list = []  # 存储生成的肽段
    seq_len = len(protein_seq)  # 蛋白质序列长度

    for i in range(len(input_df)):  # 遍历所有切割位点
        try:
            # 解析第i个切割位点信息
            row_i = input_df.iloc[i]
            pos_str = row_i['Position'].split(':')[1].strip()  # 提取位置字符串
            pos = int(round(float(pos_str.replace('"', ''))))  # 转换为整数位置
            score_i = float(row_i['Cleavage score'].split(':')[1].replace('"', ''))  # 提取切割得分
            enzyme_i = re.sub(r"\s+", " ", row_i['Enzyme Name'])  # 清理酶名称格式

            # 与后续切割位点(j)配对
            for j in range(i + 1, len(input_df)):
                # 解析第j个切割位点
                row_j = input_df.iloc[j]
                pos_str_j = row_j['Position'].split(':')[1].strip()
                pos_j = int(round(float(pos_str_j.replace('"', ''))))
                score_j = float(row_j['Cleavage score'].split(':')[1].replace('"', ''))
                enzyme_j = re.sub(r"\s+", " ", row_j['Enzyme Name'])

                # 确定两个切割位点的顺序
                if pos > pos_j:
                    start, end = pos_j, pos
                    enzyme1, enzyme2 = enzyme_j, enzyme_i
                else:
                    start, end = pos, pos_j
                    enzyme1, enzyme2 = enzyme_i, enzyme_j

                # 跳过无效位置
                if start < 0 or end > seq_len:
                    continue

                # 提取肽段序列
                peptide = protein_seq[start:end]
                length = len(peptide)

                # 长度筛选(3-7个氨基酸)
                if PEPTIDE_LENGTH_MIN <= length <= PEPTIDE_LENGTH_MAX:
                    peptide_list.append({
                        'GN': gene_name,
                        'Enzyme': f"{enzyme1} and {enzyme2}",  # 酶组合名称
                        'start_cleavage_site': start,  # 起始位置
                        'end_cleavage_site': end,      # 结束位置
                        'average_cleavage_score': (score_i + score_j) / 2,  # 平均得分
                        'length_of_peptide': length,   # 肽段长度
                        'peptide': peptide             # 肽序列
                    })

        except Exception as e:
            logging.error(f"Error processing row {i}: {e}")
    
    return pd.DataFrame(peptide_list)  # 返回肽段DataFrame

def find_matching_proteins(peptide, protein_info, protein_domains):
    """检查肽段在其他蛋白质中的存在情况并匹配结构域"""
    matches = []  # 存储匹配结果
    
    # 遍历所有蛋白质
    for idx, line in enumerate(protein_info[' The_original_sequence']):
        try:
            # 跳过非字符串或未找到的情况
            if not isinstance(line, str) or peptide not in line:
                continue
                
            gene = protein_info['Protein Gene Name'][idx]  # 当前蛋白质基因名
            seq = get_protein_sequence(protein_info, gene)  # 获取序列
            if not seq:  # 序列不存在则跳过
                continue
                
            # 确定肽段在蛋白质中的位置
            start = find_substring_index(seq, peptide)
            end = start + len(peptide)
            
            # 获取该蛋白质的结构域信息
            domains = get_protein_domains(protein_domains, gene)
            matched_domain = ""
            
            # 检查肽段是否位于结构域内
            for _, domain_row in domains.iterrows():
                d_start, d_end = domain_row.iloc[4], domain_row.iloc[5]  # 结构域起止位置
                # 检查重叠情况 (肽段部分进入结构域即视为匹配)
                if (d_start <= start < d_end) or (d_start < end <= d_end) or (start <= d_start and end >= d_end):
                    matched_domain = f"{gene}({domain_row.iloc[3]})"  # 记录基因(结构域名)
                    break
            
            # 存储匹配结果
            matches.append(matched_domain if matched_domain else f"{gene}( )")
                
        except Exception as e:
            logging.error(f"Error checking match at index {idx}: {e}")
    
    return ';'.join(matches)  # 用分号连接所有匹配结果

if __name__ == '__main__':
    # 加载三大数据集
    try:
        protein_domains, protein_protease_sites, protein_info = load_data()
    except Exception as e:
        logging.error(f"Failed to load data: {e}")
        exit(1)
    
    # 指定关注的蛋白酶ID列表
    protease_ids = ['C01.036', 'M10.003', 'M10.004', 'M10.005', 'M10.008',
                    'S01.131', 'S01.133', 'S01.010']
    
    # 待分析的基因列表(300+个基因)
    genes_to_analyze = [ 
        "HIST1H3A", "HGFAC", "HGF", "HADHA", "H3F3A", "GSN", "GRP", "GPX3", "GPI", "GPC2", "GPC1",
        "GNLY", "GLIPR2", "GLG1", "GLA", "GHR", "GGH", "GDNF", "GBP1", "GBA", "GAPDH", "GANAB",
        "FST", "FRAS1", "FOLR2", "FN1", "FMOD", "FLT1", "FLNB", "FLNA", "FKBP1A", "FGG", "FGFR4",
        "FGFR3", "FGFR2", "FGFR1", "FGF7", "FGF3", "FGF2", "FGB", "FGA", "FBN2", "FBN1", "FBLN5",
        "FBLN2", "FBLN1", "FAM3C", "F2R", "F2", "ESR2", "ERAP1", "EPO", "EPHB4", "ENDOD1", "EMILIN1",
        "EMCN", "ELN", "EGFR", "EFTUD2", "EFNA3", "EFNA1", "EFEMP2", "EFEMP1", "EEF2", "EDN1",
        "ECM1", "DYNC1H1", "DSG1", "DSC2", "DPT", "DPP7", "DMKN", "DMBT1", "DEFB4A", "DEFB1",
        "DEFA3", "DEFA1", "DCN", "DCD", "DAG1", "CXCL9", "CXCL8", "CXCL12", "CXCL1", "CXADR",
        "CTSV", "CTSK", "CTSG", "CTSD", "CTHRC1", "CTGF", "CSTA", "CST6", "CST4", "CST3", "CSPG4",
        "CRTAP", "CRIP2", "CRH", "CREG1", "CPQ", "CPE", "CPA3", "CP", "COPA", "COMP", "COLQ",
        "COL7A1", "COL6A6", "COL6A5", "COL6A3", "COL6A2", "COL6A1", "COL5A2", "COL5A1", "COL4A6",
        "COL4A5", "COL4A4", "COL4A3", "COL4A2", "COL4A1", "COL3A1", "COL2A1", "COL1A2", "COL1A1",
        "COL18A1", "COL17A1", "COL16A1", "COL15A1", "COL14A1", "COL12A1", "CMA1", "CLU", "CLTC",
        "CLCA4", "CKLF", "CKAP4", "CHGB", "CHGA", "CFI", "CFHR3", "CFHR1", "CFH", "CFD", "CFB",
        "CELA1", "CDSN", "CDH1", "CD8B", "CD8A", "CD5L", "CD55", "CD40", "CD34", "CD209", "CD163",
        "CD14", "CD109", "CCT6A", "CCT2", "CCL7", "CCL5", "CCL27", "CCL26", "CCL2", "CCL11",
        "CASP4", "CASP1", "CAPZA2", "CAPZA1", "CANX", "CAMP", "CALU", "CALR", "CALM1", "C9",
        "C8B", "C7", "C5", "C4BPA", "C4B", "C4A", "C3", "C1S", "C1R", "C1QC", "C1QBP", "C1QB",
        "BGN", "BDNF", "B2M", "AZGP1", "ATP5O", "ATP5B", "ATP5A1", "ASPN", "ARF4", "APOH",
        "APOE", "APOD", "APOB", "APOA4", "APOA2", "APOA1", "APCS", "ANXA2P2", "ANXA2", "ANXA1",
        "ANOS1", "ANGPT2", "ANGPT1", "AMBP", "ALDOA", "ALCAM", "ALB", "AIMP1", "AHSG", "AGTR2",
        "AGT", "AFM", "AEBP1", "ADM2", "ADM", "ADIPOQ", "ADAMTSL5", "ADAMTS17", "ACTN4", "ACTN2",
        "ACTN1", "ACPP", "ACE2", "ACE", "ABI3BP", "A2ML1", "A2M", "A1BG"
    ]

    def process_gene(gene):
        # 构造输出路径
        output_path = os.path.join(OUTPUT_DIR, f"{gene}.csv")
        
        # 检查结果是否已存在(避免重复计算)
        if os.path.exists(output_path):
            logging.info(f"Skipping {gene}: Result file already exists.")
            return
        
        logging.info(f"Analysing: {gene}")
        
        # 获取该基因的蛋白酶切割位点
        cleavage_sites = get_protein_protease_cleavage_sites(protein_protease_sites, gene)
        # 调试：输出列名
        logging.info(f"Columns in cleavage_sites: {cleavage_sites.columns.tolist()}")
        
        # 查找ID列(兼容不同大小写格式)
        id_column = [col for col in cleavage_sites.columns 
                     if 'ID' in col or 'id' in col or 'Id' in col]
        if not id_column:
            logging.error(f"No ID column found in cleavage_sites for {gene}")
            return
            
        # 筛选指定蛋白酶ID的切割位点
        filtered_sites = cleavage_sites[cleavage_sites[id_column[0]].isin(protease_ids)]
        
        # 没有匹配位点则跳过
        if filtered_sites.empty:
            logging.info(f"No matching protease sites for {gene}.")
            return
        
        # 获取蛋白质序列
        protein_seq = get_protein_sequence(protein_info, gene)
        if not protein_seq:
            logging.info(f"No sequence found for {gene}. Skipping.")
            return
        
        # 生成符合条件的肽段(3-7aa)
        peptides = process_peptides_matrix(filtered_sites, gene, protein_seq)
        if peptides.empty:
            logging.info(f"No peptides generated for {gene}.")
            return
        
        
        # 为每个肽段寻找在其他蛋白质中的匹配情况
        domain_annotations = []
        # 使用进度条显示匹配进度
        for peptide in tqdm(peptides['peptide'],
                            total=len(peptides),
                            desc=f"Matching peptides for {gene}",
                            leave=False):
            result = find_matching_proteins(peptide, protein_info, protein_domains)
            domain_annotations.append({'peptide_proteins_domains': result})
        
        # 将匹配结果合并到肽段DataFrame
        peptides = pd.concat([peptides.reset_index(drop=True),
                             pd.DataFrame(domain_annotations)], axis=1)
        
        # 保存结果到CSV
        try:
            peptides.to_csv(output_path, index=False)
            logging.info(f"Saved results to: {output_path}")
        except Exception as e:
            logging.error(f"Error saving results for {gene}: {e}")

    # 确保输出目录存在
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    
    # 使用线程池并行处理所有基因
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # 用tqdm显示总进度条
        list(tqdm(executor.map(process_gene, genes_to_analyze),
                  total=len(genes_to_analyze),
                  desc="Processing genes"))