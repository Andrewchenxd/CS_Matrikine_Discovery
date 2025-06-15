#!/usr/bin/env python3
import pandas as pd
import re
import os
import logging
from collections import OrderedDict
from Bio.SeqIO import parse

# 配置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# 全局阈值
CUT_THRESHOLD = 0.9
min_length = 4  # 最小片段长度
max_length = 50  # 最大片段长度

def read_fasta(fastafile):
    """读取FASTA文件到有序字典"""
    fasta_dict = OrderedDict()
    for record in parse(fastafile, "fasta"):
        seq = str(record.seq)
        seq = re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '-', seq.upper())
        fasta_dict[record.id] = seq
    return fasta_dict

def process_protein(seq_id, protein_seq, protein_cuts):
    """处理单个蛋白质，生成切割片段"""
    fragments = []
    length = len(protein_seq)
    
    # 获取所有切割位置并排序
    cut_positions = sorted({int(pos) for pos, _, _ in protein_cuts})  # 确保位置是整数
    all_cuts = [0] + cut_positions + [length]
    
    # 构建位置信息字典 {位置: (所有酶名, 最大概率)}
    pos_info = {}
    for pos, proteases, pro in protein_cuts:
        pos_int = int(pos)  # 转换为整数
        if pos_int not in pos_info:
            pos_info[pos_int] = (proteases, pro)
        else:
            # 合并相同位置的酶和概率
            existing_proteases, existing_pro = pos_info[pos_int]
            new_proteases = existing_proteases + "," + proteases
            new_pro = max(existing_pro, pro)
            pos_info[pos_int] = (new_proteases, new_pro)
    
    # 生成切割片段
    for i in range(1, len(all_cuts)):
        start = all_cuts[i-1]  # 上一个切割点位置
        end = all_cuts[i]      # 当前切割点位置
        
        # 片段序列
        frag_seq = protein_seq[start:end]
        frag_length = end - start
        
        # 处理起始点信息
        if start == 0:
            front_protease = "None"
            front_pro = 1.0
        else:
            front_protease, front_pro = pos_info[start]
        
        # 处理终止点信息
        if end == length:
            back_protease = "None"
            back_pro = 1.0
        else:
            back_protease, back_pro = pos_info[end]
        
        # 计算切割概率
        frag_prob = front_pro * back_pro
        if frag_length < min_length or frag_length > max_length:
            continue
        else:
            fragments.append({
                'sequence_id': seq_id,
                'fragment_seq': frag_seq,
                'front_protease': front_protease,
                'back_protease': back_protease,
                'start': start + 1,  # 转换为1-based起始位置
                'end': end,          # 1-based终止位置
                'length': frag_length,
                'probability': frag_prob
            })
    
    return fragments

def main():
    # 读取数据
    fasta_dict = read_fasta('./data/15个可能具有治疗作用.fasta')
    csv_path = './data/20250528235831_jiBTUqyb_results.csv'
    
    # 读取CSV文件并转换列类型
    df = pd.read_csv(csv_path, header=None, names=[
        'protease', 'sequence_id', 'position', 'seqs', 'prediction', 'pro'
    ])
    
    # 转换列类型，先处理缺失值
    df['pro'] = pd.to_numeric(df['pro'], errors='coerce')
    df['position'] = pd.to_numeric(df['position'], errors='coerce')
    
    # 删除包含缺失值的行
    df = df.dropna(subset=['position', 'pro'])
    
    # 现在可以安全地转换为整数
    df['position'] = df['position'].astype(int)
    
    # 过滤阈值以上的切割点
    df = df[df['pro'] > CUT_THRESHOLD]
    
    # 按sequence_id分组收集切割信息
    protein_data = {}
    for seq_id, group in df.groupby('sequence_id'):
        # 只处理存在于FASTA中的蛋白质
        if seq_id in fasta_dict:
            cuts = []
            for _, row in group.iterrows():
                cuts.append((
                    row['position'], 
                    row['protease'], 
                    row['pro']
                ))
            protein_data[seq_id] = cuts
    
    # 处理所有蛋白质
    all_fragments = []
    for seq_id, protein_seq in fasta_dict.items():
        if seq_id in protein_data:
            fragments = process_protein(seq_id, protein_seq, protein_data[seq_id])
            all_fragments.extend(fragments)
        else:
            if  len(protein_seq) < min_length or len(protein_seq) > max_length:
                logger.warning(f"蛋白质 {seq_id} 的长度不符合要求，跳过处理")
                continue
            else:
                # 无切割点的蛋白质作为完整片段
                all_fragments.append({
                    'sequence_id': seq_id,
                    'fragment_seq': protein_seq,
                    'front_protease': "None",
                    'back_protease': "None",
                    'start': 1,
                    'end': len(protein_seq),
                    'length': len(protein_seq),
                    'probability': 1.0
                })
    
    # 转换为DataFrame并保存
    result_df = pd.DataFrame(all_fragments)
    output_path = './output/protein_fragments_results.csv'
    result_df.to_csv(output_path, index=False)
    logger.info(f"结果已保存至: {output_path}")
    logger.info(f"共处理 {len(fasta_dict)} 个蛋白质，生成 {len(result_df)} 个片段")

if __name__ == "__main__":
    main()