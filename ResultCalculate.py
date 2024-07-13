#!/usr/bin/python3
import re
import os
import sys
import numpy as np
import pandas as pd

class ResultCalculate:
    '''
    计算拷贝数结果
    名称：ResultCalculate
    输入文件：样本信息统计结果
    输出：计算结果
    '''

    def __init__(self, merge_df, contron_list, case_list, out_path, refgene_list):
        self.merge_df = merge_df
        self.contron_list = contron_list
        self.case_list = case_list
        self.out_path = out_path
        self.refgene_list = refgene_list

    def arrange(self):
        '''
        提取基本信息与样本计算reads的列
        Args:
            merge_df: 样本合并或过滤结果
        Returns: 提取的信息列
        '''

        common_df = self.merge_df[['CHR', 'START', 'END', 'GENE']]
        filtered_df = self.merge_df.filter(regex='CORRECT_READS')
        arrange_df = pd.concat([common_df, filtered_df], axis=1)
        arrange_df = arrange_df.rename(columns=lambda x: re.sub(r'CORRECT_READS\.', '', x))
        return arrange_df

    def calculate_pro(self, analyses_df, contron_relative_list, case_relative_list, refgene_list):
        '''
        计算检测样本各个位点的相对比值
        Args:
            analyses_df: 相关信息提取结果列
            contron_relative_list: 相对阴性样本列表
            case_relative_list: 相对检测样本列表
        Returns: 归一化结果及其参考均值
        '''

        analyses_df = analyses_df.copy()
        refgene_df = analyses_df.loc[analyses_df['GENE'].isin(refgene_list)]
        relative_list = []
        for contron_sample in contron_relative_list:
            analyses_df['{}.relative'.format(contron_sample)] = analyses_df[contron_sample] / refgene_df[contron_sample].median()
            relative_list.append('{}.relative'.format(contron_sample))

        if len(relative_list) > 0:
            analyses_df['relative_means'] = analyses_df[relative_list].mean(axis=1)
        else:
            analyses_df['relative_means'] = 1

        for case_sample in case_relative_list:
            analyses_df['{}.relative'.format(case_sample)] = analyses_df[case_sample] / refgene_df[case_sample].median()
            analyses_df['{}.pro'.format(case_sample)] = analyses_df['{}.relative'.format(case_sample)] / analyses_df['relative_means']
        return analyses_df

    def calculate_CN(self, calculate_df):
        '''
        计算检测样本各个位点的CN值
        Args:
            calculate_df: 相对比值计算结果
        Returns: 各位点CN值计算结果
        '''

        for sample in self.case_list:
            calculate_df['{}.pro'.format(sample)] = \
                calculate_df['{}.pro'.format(sample)].apply(lambda x: x if x <= 2 else 2)
            calculate_df['{}_CN'.format(sample)] = 'CN2'
            calculate_df.loc[(calculate_df['{}.pro'.format(sample)] >= 0) & (calculate_df['{}.pro'.format(sample)] < 0.4), \
                             '{}_CN'.format(sample)] = 'CN0'
            calculate_df.loc[(calculate_df['{}.pro'.format(sample)] >= 0.4) & (calculate_df['{}.pro'.format(sample)] <= 0.65), \
                             '{}_CN'.format(sample)] = 'CN1'
            calculate_df.loc[(calculate_df['{}.pro'.format(sample)] >= 1.3) & (calculate_df['{}.pro'.format(sample)] < 1.65), \
                             '{}_CN'.format(sample)] = 'CN3'
            calculate_df.loc[(calculate_df['{}.pro'.format(sample)] >= 1.65) & (calculate_df['{}.pro'.format(sample)] <= 2), \
                             '{}_CN'.format(sample)] = 'CN4'
        return calculate_df

    def result(self):
        '''
        性别识别及检测样本CN值结果
        Args:
            merge_df: 检测样本合并结果
            contron_list: 阴性样本列表
            case_list: 检测样本列表
            out_path: 输出结果路径
            refgene_list：内参基因列表
        Returns: 性别结果及样本检测结果
        '''

        arrange_df = self.arrange()
        ref_gene_df = pd.read_table(self.refgene_list, header=None, names=['GENE'])
        autosome_df = arrange_df[(arrange_df['CHR'] != 'chrX') & (arrange_df['CHR'] != 'chrY')]
        genderX_df = arrange_df[arrange_df['CHR'] == 'chrX']
        genderY_df = arrange_df[arrange_df['CHR'] == 'chrY']
        sample_list = self.contron_list + self.case_list
        refgene_gender_df = pd.merge(left=ref_gene_df[['GENE']], right=arrange_df[['GENE', 'CHR']], on='GENE')
        autosome_refgene_list = refgene_gender_df.loc[(refgene_gender_df['CHR'] != 'chrX') & \
                                                      (refgene_gender_df['CHR'] != 'chrY'), 'GENE'].to_list()
        chrX_refgene_list = refgene_gender_df.loc[refgene_gender_df['CHR'] == 'chrX', 'GENE'].to_list()
        chrY_refgene_list = refgene_gender_df.loc[refgene_gender_df['CHR'] == 'chrY', 'GENE'].to_list()
        gender_calculate = {}
        male_list = []
        female_list = []
        for sample in sample_list:
            if genderY_df[sample].median() >= autosome_df[sample].median() * 0.05:
                gender_calculate[sample] = '男'
                male_list.append(sample)
            else:
                gender_calculate[sample] = '女'
                female_list.append(sample)
        gender_df = pd.DataFrame.from_dict(gender_calculate, orient='index', columns=['性别'])
        gender_df = gender_df.reset_index().rename(columns={'index': '样本名称'})
        gender_file = os.path.join(self.out_path, '样本性别.xlsx')
        gender_df.to_excel(gender_file, index = False, header = True)

        male_row = ['CHR', 'START', 'END', 'GENE'] + male_list
        female_row = ['CHR', 'START', 'END', 'GENE'] + female_list
        maleX_df = genderX_df[male_row]
        femaleX_df = genderX_df[female_row]
        maleY_df = genderY_df[male_row]
        femaleY_df = genderY_df[female_row]

        autosome_df = self.calculate_pro(autosome_df, self.contron_list, self.case_list, autosome_refgene_list)
        femaleX_df = self.calculate_pro(femaleX_df, list(set(self.contron_list) - set(male_list)), \
                                        list(set(self.case_list) - set(male_list)), chrX_refgene_list)
        maleX_df = self.calculate_pro(maleX_df, list(set(self.contron_list) - set(female_list)), \
                                      list(set(self.case_list) - set(female_list)), chrX_refgene_list)
        maleY_df = self.calculate_pro(maleY_df, list(set(self.contron_list) - set(female_list)), \
                                      list(set(self.case_list) - set(female_list)), chrY_refgene_list)
        femaleY_df = femaleY_df.copy()
        for contron_female in list(set(self.contron_list) - set(male_list)):
            femaleY_df['{}.relative'.format(contron_female)] = 0
        femaleY_df['relative_means'] = 0
        for case_female in list(set(self.case_list) - set(male_list)):
            femaleY_df['{}.relative'.format(case_female)] = 0
            femaleY_df['{}.pro'.format(case_female)] = 0

        female_df = pd.concat([femaleX_df, femaleY_df], axis=0)
        male_df = pd.concat([maleX_df, maleY_df], axis=0)
        female_df.drop('relative_means', axis=1, inplace=True)
        male_df.drop('relative_means', axis=1, inplace=True)


        gender_df = pd.merge(female_df, male_df, on=['CHR', 'START', 'END', 'GENE'], how = 'left')
        calculate_df = pd.concat([autosome_df, gender_df], axis=0)

        calculate_df = self.calculate_CN(calculate_df)
        calculate_df.reset_index(drop=True, inplace=True)

        out_file = os.path.join(self.out_path, 'result.txt')
        calculate_df.to_csv(out_file, sep='\t', index=False, header=True)
        return calculate_df, gender_calculate
