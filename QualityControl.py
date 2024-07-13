#!/usr/bin/python3
import re
import os
import sys
import pandas as pd

class QualityControl:
    '''
    Bam文件质量值计算与统计
    名称：QualityControl
    输入文件：bam_list
    输出：QC_statistics
    '''
    # def __init__(self, bamdst, bed_file, in_bam, out_path):
    #     self.bamdst = bamdst
    #     self.bed_file = bed_file
    #     self.in_bam = in_bam
    #     self.out_path = out_path

    @staticmethod
    def qc_calculate(bamdst, bed_file, out_path, in_bam):
        '''
        bam文件质量值计算
        Args:
            bamdst: bamdst软件路径
            bed_file: 目标区域文件
            in_bam: 输入bam文件
            out_path: 输出结果路径
        Returns: 样本质量值计算结果
        '''
        bam2qc_cmd = '{} -p {} -o {} {}'.format(bamdst, bed_file, out_path, in_bam)
        os.system(bam2qc_cmd)

    @staticmethod
    def qc_statistic(sample, out_path):
        qc_result = os.path.join(out_path, 'coverage.report')
        dp_result = os.path.join(out_path, 'depth.tsv.gz')
        dp_file = pd.read_csv(dp_result, compression='gzip', sep='\t')
        if os.path.exists(qc_result) and os.path.exists(dp_result):
            stat_file = os.path.join(out_path, 'qc_statistic.txt')
            stat_dict = {}
            stat_dict['样本名称'] = sample

            with open(qc_result, 'r') as qc:
                for q in qc:
                    info = q.strip().split('\t')
                    if '##' in info[0]:
                        continue
                    elif info[0] == '[Total] Raw Data(Mb)':
                        stat_dict['测序数据量(Mb)'] = info[1]
                    elif info[0] == '[Total] Fraction of Mapped Data(Mb)':
                        stat_dict['比对率'] = info[1]
                    elif info[0] == '[Total] Fraction of MapQ reads in all reads':
                        stat_dict['MapQReadsRatio'] = info[1]
                    elif info[0] == '[Total] Mapped Reads':
                        stat_dict['比对上的reads'] = info[1]
                    elif info[0] == '[Target] Target Reads':
                        stat_dict['目标区域reads'] = info[1]
                    elif info[0] == '[Target] Len of region':
                        stat_dict['目标区域长度'] = info[1]
                    elif info[0] == '[Target] Target Data(Mb)':
                        stat_dict['目标区域数据量(Mb)'] = info[1]
                    elif info[0] == '[Target] Target Data Rmdup(Mb)':
                        stat_dict['TargetDataRmdup(Mb)'] = info[1]
                    elif info[0] == '[Target] Average depth':
                        stat_dict['平均深度'] = info[1]
                    elif info[0] == '[Target] Average depth(rmdup)':
                        stat_dict['AverageDepth(rmdup)'] = info[1]
                    elif info[0] == '[Target] Fraction of Target Data in all data':
                        stat_dict['目标区域数据比率'] = info[1]
                    elif info[0] == '[Target] Coverage (>0x)':
                        stat_dict['覆盖度(>0x)'] = info[1]
                    elif info[0] == '[Target] Coverage (>=10x)':
                        stat_dict['覆盖度(>=10x)'] = info[1]
                    elif info[0] == '[Target] Coverage (>=100x)':
                        stat_dict['覆盖度(>=100x)'] = info[1]
                    else:
                        continue

            stat_dict['均一性'] = (dp_file["Raw Depth"] >= int(float(stat_dict['平均深度'])) * 0.2).mean().round(4)
            stat_dict['捕获效率'] = round(float(stat_dict['目标区域reads']) / float(stat_dict['比对上的reads']), 4)
            stat_dict["覆盖度(>4x)"] = (dp_file["Raw Depth"] >= 4).mean().round(4)
            stat_dict["覆盖度(>=20x)"] = (dp_file["Raw Depth"] >= 20).mean().round(4)
            stat_dict["覆盖度(>=500x)"] = (dp_file["Raw Depth"] >= 500).mean().round(4)

            stat_df = pd.DataFrame(data=[stat_dict], columns=['样本名称', '测序数据量(Mb)', '比对率', \
                                                              '目标区域长度', '目标区域数据量(Mb)', \
                                                              '捕获效率', '平均深度', '均一性', \
                                                              '覆盖度(>4x)', '覆盖度(>=20x)', '覆盖度(>=500x)'])

            stat_df.to_csv(stat_file, sep='\t', index=False, header=True)
        else:
            print ('Warning: Could not find {} file, skipping this QC'.format(qc_result))

    def run_qc_calculate(self, sample, bamdst, bed_file, out_path, in_bam):
        '''
        运行程序
        Args:
            bamdst: bamdst软件路径
            bed_file: 目标区域文件
            in_bam: 输入bam文件
            out_path: 输出结果路径
        Returns: 样本质量值计算结果
        '''
        self.qc_calculate(bamdst, bed_file, out_path, in_bam)
        self.qc_statistic(sample, out_path)