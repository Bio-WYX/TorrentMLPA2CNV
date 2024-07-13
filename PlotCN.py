#!/usr/bin/python3
import re
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class PlotCN:
    '''
    分析结果作图
    名称：PlotCN
    输入文件：分析结果
    输出：检测样本CN比例图
    '''

    def __init__(self, result_df, out_path):
        self.result_df = result_df
        self.out_path = out_path

    def plot(self):
        '''
        分析结果作图
        Args:
            result_df: 分析结果
            out_path: 输出路径
        Returns: 检测样本CN比例图
        '''

        #result_df = pd.read_table(result_file, sep = '\t', header = 0)
        case_df = self.result_df.filter(like='.pro')
        case_list = case_df.columns.tolist()
        width = self.result_df.shape[0] / 100 * 20

        colors = {'CN0':'darkorchid', 'CN1':'orangered', 'CN2':'limegreen', 'CN3':'cornflowerblue', 'CN4':'mediumblue'}

        for sample in case_list:
            name = re.sub("\.pro", "", sample)
            self.result_df['colors'] = 'limegreen'
            self.result_df.loc[:, 'colors'] = self.result_df['{}_CN'.format(name)].apply(lambda x: colors[x])
            color_dict = self.result_df[['colors']].to_dict()

            plt.figure(figsize=(width, 8))
            ax = plt.axes()
            plt.scatter(self.result_df['GENE'], self.result_df[sample], c=self.result_df['colors'])
            ax.set_ylim(bottom = -0.05, top = 2.05)
            ax.set_xlim(left = -1.5, right = self.result_df.shape[0] + 1.5)
            plt.xticks(rotation=90, fontsize=8)
            plt.yticks(list(np.arange(0, 2.2, 0.2)))
            #plt.legend(loc='upper right')  #图例位置
            plt.xlabel('GENE')
            plt.ylabel('ratio')
            plt.title('{} reads ratio'.format(name))
            #plt.grid(True)  #背景网格线
            for n in [0.4, 0.65, 1.3, 1.65]:
                plt.axhline(n, color="gray", linestyle='--', label='CN')
            plt.axhspan(0.65, 0.8, alpha=0.25, color='gray')  # 灰区1
            plt.axhspan(1.2, 1.3, alpha=0.25, color='gray')  # 灰区2
            #plt.axhline([n for n in [0.2, 0.65, 1.3, 1.65]])
            # plt.xticks(CN1_df['GENE'], color='orangered')
            # plt.xticks(CN2_df['GENE'], color='orangered')

            for key in color_dict['colors'].keys():
                color = color_dict['colors'][key]
                ax.get_xticklabels()[key].set_color(color)

            ax.spines['top'].set_visible(False)  #去掉上框线
            ax.spines['right'].set_visible(False)  #去掉右框线
            #ax.spines['bottom'].set_visible(False)
            #ax.spines['left'].set_visible(False)

            plot_file = os.path.join(self.out_path, '{}_result.pdf'.format(name))
            plt.savefig(plot_file, format='pdf', bbox_inches='tight', pad_inches=0.2)
#
# if __name__ == '__main__':
#     result_file = 'E:\\ansaisi\\MLPA_NGS\\test (2).txt'
#     plot(result_file)