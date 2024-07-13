#!/usr/bin/python3
import re
import os
import sys
import pandas as pd

class MergeStatistic:
    '''
    所有样本的比对信息结果合并
    名称：MergeStatistic
    输入文件：样本比对信息
    输出：所有样本合并结果
    '''

    def __init__(self, contron_list, case_list, out_path):
        self.contron_list = contron_list
        self.case_list = case_list
        self.out_path = out_path

    def merge_file(self):
        '''
        合并不同样本分析结果
        Args:
            contron_list: 阴性样本列表
            case_list: 检测样本列表
            out_path: 输出结果路径
        Returns: 所有样本检测结果合并文件
        '''
        n = 1
        merge_df = pd.DataFrame()
        for contron in self.contron_list:
            stat_file = os.path.join(self.out_path, '{}.stat.txt'.format(contron))
            contron_df = pd.read_table(stat_file, sep = '\t', header = 0)
            if n == 1:
                merge_df = contron_df
            else:
                merge_df = pd.merge(merge_df, contron_df, on = ['CHR', 'START', 'END', 'GENE'], \
                                    suffixes = ('', '.{}'.format(contron)))
            n += 1
        for case in self.case_list:
            stat_file = os.path.join(self.out_path, '{}.stat.txt'.format(case))
            case_df = pd.read_table(stat_file, sep='\t', header=0)
            merge_df = pd.merge(merge_df, case_df, on=['CHR', 'START', 'END', 'GENE'], \
                                suffixes=('', '.{}'.format(case)))
            n += 1

        merge_df = merge_df.rename(columns={'ALL_READS': 'ALL_READS.{}'.format(self.contron_list[0]), \
                                            'CORRECT_READS': 'CORRECT_READS.{}'.format(self.contron_list[0]), \
                                            'FAILED_READS': 'FAILED_READS.{}'.format(self.contron_list[0]), \
                                            'CORRECT_RATIO(%)': 'CORRECT_RATIO(%).{}'.format(self.contron_list[0])})

        out_file = os.path.join(self.out_path, 'merge_out.txt')
        merge_df.to_csv(out_file, sep='\t', index=False, header = True)
        return merge_df