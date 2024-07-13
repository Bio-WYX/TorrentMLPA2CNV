#!/usr/bin/python3
import re
import os
import sys
import numpy as np
import pandas as pd

class FilterPlotGene:
    '''
    过滤掉作图中的内参基因
    名称：FilterPlotGene
    输入文件：样本计算结果
    输出：删除内参基因的结果
    '''

    def __init__(self, result_df, refgene_list):
        self.result_df = result_df
        self.refgene_list = refgene_list

    def ref_filters(self):
        '''
        过滤掉内参基因作为作图使用结果
        Args:
            result_df: 样本计算结果
            refgene_list: 内参基因列表
        Returns: 过滤内参基因后的检测结果
        '''

        ref_gene_df = pd.read_table(self.refgene_list, header=None, names=['GENE'])
        refgene_gender_df = pd.merge(left=ref_gene_df[['GENE']], right=self.result_df[['GENE', 'CHR']], on='GENE')
        autosome_refgene_df = refgene_gender_df.loc[(refgene_gender_df['CHR'] != 'chrX') & \
                                                      (refgene_gender_df['CHR'] != 'chrY'), 'GENE']
        filter_gene_df = self.result_df.merge(autosome_refgene_df, on = ['GENE'], how = 'outer', indicator = True, \
                                 suffixes=('', '_drop_right')).loc[lambda x: x['_merge'] == 'left_only']
        filter_gene_df = filter_gene_df.drop(filter_gene_df.filter(regex='_drop_right').columns, axis=1)
        filter_gene_df = filter_gene_df.drop('_merge', axis=1)

        filter_gene_df.reset_index(drop=True, inplace=True)
        return filter_gene_df