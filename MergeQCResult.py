#!/usr/bin/python3
import re
import os
import sys
import pandas as pd
import numpy as np

class MergeQCResult:
    '''
    所有样本质控结果合并
    名称：MergeQCResult
    输入文件：样本质控信息
    输出：所有样本合并结果
    '''

    def __init__(self, sample_list, gender_dict, gender_calculate, out_path):
        self.sample_list = sample_list
        self.gender_dict = gender_dict
        self.gender_calculate = gender_calculate
        self.out_path = out_path

    def merge_qc(self):
        '''
        合并不同样本质控结果
        Args:
            sample_list: 样本列表
            gender_dict：输入性别信息
            gender_calculate：性别计算结果
            out_path: 输出结果路径
        Returns: 所有样本质控结果合并文件
        '''
        n = 1
        merge_raw_df = pd.DataFrame()
        merge_reali_df = pd.DataFrame()
        raw_qc_file = os.path.join(self.out_path, 'QC_raw.xlsx')
        reali_qc_file = os.path.join(self.out_path, '质控结果.xlsx')
        for sample in self.sample_list:
            raw_file = os.path.join(self.out_path, 'QualityControl', '{}'.format(sample), 'raw_qc', 'qc_statistic.txt')
            reali_file = os.path.join(self.out_path, 'QualityControl', '{}'.format(sample), 're_ali_qc', 'qc_statistic.txt')
            if os.path.exists(raw_file) and os.path.exists(reali_file):
                raw_df = pd.read_table(raw_file, sep = '\t', header = 0)
                reali_df = pd.read_table(reali_file, sep = '\t', header = 0)
                if n == 1:
                    merge_raw_df = raw_df
                    merge_reali_df = reali_df
                else:
                    merge_raw_df = pd.concat([merge_raw_df, raw_df], axis = 0)
                    merge_reali_df = pd.concat([merge_reali_df, reali_df], axis = 0)
                n += 1
            else:
                print ('Warning: Could not find {} sample qc file, skipping this QC'.format(sample))

        merge_raw_df.to_excel(raw_qc_file, index = False, header = True)

        merge_reali_df['输入性别'] = merge_reali_df.样本名称.map(self.gender_dict)
        merge_reali_df['判断性别'] = merge_reali_df.样本名称.map(self.gender_calculate)

        merge_reali_df['性别核对结果'] = np.where(merge_reali_df['输入性别'].eq(merge_reali_df['判断性别']), '正确', '错误')

        merge_reali_df.to_excel(reali_qc_file, index=False, header=True)

