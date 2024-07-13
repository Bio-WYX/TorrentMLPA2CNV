#!/usr/bin/python3
import re
import os
import sys
import numpy as np
import pandas as pd

class FilterGene:
    '''
    过滤掉阴性样本中覆盖度低的基因
    名称：FilterGene
    输入文件：样本合并结果
    输出：过滤后的结果
    '''

    def __init__(self, merge_df, contron_list, out_path):
        self.merge_df = merge_df
        self.contron_list = contron_list
        self.out_path = out_path

    def filters(self):
        '''
        过滤阴性样本中覆盖度太低的检测位点
        Args:
            merge_df: 检测样本合并结果
            contron_list: 阴性样本列表
            out_path: 输出结果路径
        Returns: 过滤位点后的检测结果
        '''

        autosome_df = self.merge_df[(self.merge_df['CHR'] != 'chrX') & (self.merge_df['CHR'] != 'chrY')]
        genderX_df = self.merge_df[self.merge_df['CHR'] == 'chrX']
        genderY_df = self.merge_df[self.merge_df['CHR'] == 'chrY']

        n = 1
        filtered_df = pd.DataFrame()
        for contron_sample in self.contron_list:
            filtered_out = os.path.join(self.out_path, '{}.filtered.txt'.format(contron_sample))
            autosome_filtered_df = autosome_df[autosome_df['CORRECT_READS.{}'.format(contron_sample)] <= \
                                               autosome_df['CORRECT_READS.{}'.format(contron_sample)].median() * 0.05]
            genderX_filtered_df = genderX_df[genderX_df['CORRECT_READS.{}'.format(contron_sample)] <= \
                                               genderX_df['CORRECT_READS.{}'.format(contron_sample)].median() * 0.05]
            genderY_filtered_df = genderY_df[genderY_df['CORRECT_READS.{}'.format(contron_sample)] <= \
                                               genderY_df['CORRECT_READS.{}'.format(contron_sample)].median() * 0.05]
            sample_filtered_df = pd.concat([autosome_filtered_df, genderX_filtered_df, genderY_filtered_df], axis=0)
            sample_filtered_df.to_csv(filtered_out, sep='\t', index=False, header=True)

            if n == 1:
                filtered_df = sample_filtered_df
            else:
                filtered_df = pd.merge(filtered_df, sample_filtered_df, \
                                       on=['CHR', 'START', 'END', 'GENE'], how = 'outer', suffixes=('', '_drop_right'))
            n += 1

        diff_df = self.merge_df.merge(filtered_df, how='outer', indicator=True, \
                                 suffixes=('', '_drop_right')).loc[lambda x: x['_merge'] == 'left_only']
        diff_df = diff_df.drop(diff_df.filter(regex='_drop_right').columns, axis=1)
        diff_df = diff_df.drop('_merge', axis=1)

        diff_df.reset_index(drop=True, inplace=True)
        return diff_df