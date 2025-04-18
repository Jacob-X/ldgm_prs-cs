# --*-- conding:utf-8 --*--
# @Time : 2025/4/14 15:52
# @Author : Jacob Xu
# @Email : xu.tianxin@foxmail.com
# @File : 2.bcf_to_hdf5.py
# @Software : PyCharm
import pysam
import numpy as np
import h5py
from collections import defaultdict

def parse_ldgm_bcf_to_hdf5(bcf_path, output_h5_path):
    bcf = pysam.VariantFile(bcf_path)

    # Block-wise storage
    block_snps = defaultdict(list)
    block_nodes = defaultdict(list)
    block_diag = defaultdict(dict)
    block_neighbors = defaultdict(dict)
    block_weights = defaultdict(dict)

    for rec in bcf.fetch():
        chrom = rec.chrom
        pos = rec.pos
        site_id = f"{chrom}:{pos}"

        info = rec.info
        block = info.get("LD_block")
        node = info.get("LD_node")
        diag = info.get("LD_diagonal", 1.0)
        neighbors = info.get("LD_neighbors", [])
        weights = info.get("LD_weights", [])

        block_snps[block].append(site_id)
        block_nodes[block].append(node)
        block_diag[block][node] = diag
        block_neighbors[block][node] = neighbors
        block_weights[block][node] = weights

    # 写入 HDF5 文件
    with h5py.File(output_h5_path, "w") as h5f:
        for block_id in block_snps:
            snps = block_snps[block_id]
            nodes = block_nodes[block_id]
            node_to_idx = {node: i for i, node in enumerate(nodes)}
            N = len(nodes)

            # 初始化 LD matrix
            ld_matrix = np.zeros((N, N))
            for node in nodes:
                i = node_to_idx[node]
                ld_matrix[i, i] = block_diag[block_id][node]
                neighbors = block_neighbors[block_id].get(node, [])
                weights = block_weights[block_id].get(node, [])
                for neighbor, weight in zip(neighbors, weights):
                    if neighbor in node_to_idx:
                        j = node_to_idx[neighbor]
                        ld_matrix[i, j] = weight

            # 对称化（PRS-CS 要求 LD 是对称矩阵）
            ld_matrix = (ld_matrix + ld_matrix.T) - np.diag(np.diag(ld_matrix))

            grp = h5f.create_group(f"blk_{block_id}")
            grp.create_dataset("ldblk", data=ld_matrix)
            grp.create_dataset("snplist", data=np.array(snps, dtype='S'))  # 字符串保存为字节串

    print(f"转换完成: {output_h5_path}")


if __name__ == '__main__':

    bcf_path="/Volumes/data_files/LDGM/example/height/split/LD_blocks/1kg_ldgm.EUR.chr10:100331627:104378781.bcf"
    output_h5_path="test.hdf5"

    parse_ldgm_bcf_to_hdf5(bcf_path, output_h5_path)