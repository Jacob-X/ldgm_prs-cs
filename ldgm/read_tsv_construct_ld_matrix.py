# --*-- conding:utf-8 --*--
# @Time : 2025/4/16 11:01
# @Author : Jacob Xu
# @Email : xu.tianxin@foxmail.com
# @File : read_tsv_construct_ld_matrix.py
# @Software : PyCharm
# @Description : This script reads a TSV file containing genetic data and constructs a linkage disequilibrium (LD) matrix.

import pandas as pd
import numpy as np
from scipy.sparse import lil_matrix

def read_ldgm_tsv(file_path):

    ldgm_blocks = {}
    # "/Volumes/data_files/LDGM/example/height/ldgms.GRCh38/ldgm_EUR_rsid.tsv"
    df = pd.read_csv(file_path, sep="\t")

    # 这里[-20:]是用chr22来测试的，去掉就是全部chr的数据用来计算
    first_blocks = sorted(df['LD_block'].unique(), key=int)[-20:]

    for block_id in first_blocks:
        blk_name = f"blk_{block_id}"
        blk_df = df[df['LD_block'] == block_id].copy()

        # 1. 获取所有出现过的 node，并重编号（如原编号是 [3, 7, 9]，我们变成 [0, 1, 2]）
        node_ids = sorted(blk_df['LD_node'].unique())
        node_index = {node_id: i for i, node_id in enumerate(node_ids)}  # 映射 old_id -> new_id
        num_nodes = len(node_ids)
        # print(blk_name, num_nodes)
        # print(node_index)

        node_to_snp = dict(zip(blk_df['LD_node'], blk_df['SNP']))
        snplist = [node_to_snp[n] for n in node_ids]

        # 2. 初始化 LD 矩阵
        ld_matrix = lil_matrix((num_nodes, num_nodes), dtype=np.float32)

        for _, row in blk_df.iterrows():
            if row['LD_neighbors'] == '.' or pd.isna(row['LD_neighbors']):
                node = node_index[int(row['LD_node'])]
                ld_matrix[node, node] = float(row['LD_diagonal'])
                continue

            node = node_index[int(row['LD_node'])]
            neighbors = list(map(int, row['LD_neighbors'].split(',')))
            weights = list(map(float, row['LD_weights'].split(',')))

            for n, w in zip(neighbors, weights):
                if n in node_index:  # 只处理在当前 block 中的节点
                    ni = node_index[n]
                    ld_matrix[node, ni] = w
                    ld_matrix[ni, node] = w

            ld_matrix[node, node] = float(row['LD_diagonal'])

        chrom = blk_df["CHROM"].iloc[0]

        ldgm_blocks[blk_name] = {
            "ldblk": ld_matrix.tocsr(),
            "snplist": snplist,
            "node_map": node_index,
            "CHROM": chrom
        }

    return ldgm_blocks
