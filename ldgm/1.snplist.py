# --*-- conding:utf-8 --*--
# @Time : 2025/4/14 15:02
# @Author : Jacob Xu
# @Email : xu.tianxin@foxmail.com
# @File : 1.snplist.py
# @Software : PyCharm

import pandas as pd
import glob
import os

if __name__ == '__main__':

    input_path = "/Volumes/data_files/LDGM/example/height/snplist/"
    output_file = "snpinfo_mult_1kg_hm3.txt"

    snplist_files = glob.glob(os.path.join(input_path, "*.snplist"))

    df_list = []

    for filepath in snplist_files:
        filename = os.path.basename(filepath)

        chr_str = filename.split('_')[1]  # 'chr2'
        chr_num = chr_str.replace('chr', '')  # '2'

        df = pd.read_csv(filepath)

        new_df = pd.DataFrame({
            'CHR': chr_num,
            'SNP': df['site_ids'],
            'A1': df['deriv_alleles'],
            'A2': df['anc_alleles'],
            'FRQ_AFR': df['AFR'],
            'FRQ_AMR': df['AMR'],
            'FRQ_EAS': df['EAS'],
            'FRQ_EUR': df['EUR'],
            'FRQ_SAS': df['SAS']
        })

        df_list.append(new_df)

    merged_df = pd.concat(df_list, ignore_index=True)

    merged_df['CHR'] = merged_df['CHR'].astype(int)
    merged_df = merged_df.sort_values(by='CHR')

    merged_df.to_csv(output_file, sep="\t", index=False)

    print(f"合并完成，输出文件：{output_file}")