#!/usr/bin/python3
import re
import os
import sys
import math
import logging
import argparse
import json
import configparser
import multiprocessing
import pandas as pd
from InputResolve import InputResolve
from BamReAlignment import BamReAlignment
from QualityControl import QualityControl
from ReadsStatisticReAlignment import ReadsStatisticReAlignment
from MergeStatistic import MergeStatistic
from FilterGene import FilterGene
from ResultCalculate import ResultCalculate
from MergeQCResult import MergeQCResult
from FilterPlotGene import FilterPlotGene
from PlotCN import PlotCN
from PackResult import PackResult

def read_config(config_file):
    config_info = configparser.ConfigParser()
    config_info.read(config_file)
    return config_info
    #config_info.get('pathkey', 'obs_path')

def main():
    parser = argparse.ArgumentParser(description='Torrent_MLPA使用手册', prog='Torrent_MLPA', \
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-c', '--config', metavar = 'FILE', dest='config', required=True, help='配置文件')
    parser.add_argument('-i', '--input', metavar='FILE', dest='input', required=True, help='输入文件')
    parser.add_argument('-b', '--bed', metavar='FILE', dest='bed', required=True, help='目标区域文件')
    parser.add_argument('-g', '--genelist', metavar='FILE', dest='genelist', required=True, help='内参基因列表文件')
    parser.add_argument('-o', '--outdir', dest='outdir', type=str, required=True, help='输出文件夹')
    parser.add_argument('-f', '--filter', dest='filter', type=str, required=False, default='False', help='过滤参数')
    parser.add_argument('-r', '--refgene', dest='refgene', type=str, required=False, default='True', help='图显基因过滤参数')

    args = parser.parse_args()
    config_file = args.config
    input_file = args.input
    target_file = args.bed
    refgene_file = args.genelist
    out_path = args.outdir
    filter = args.filter
    refgene = args.refgene

    config_info = read_config(config_file)
    fasta_file = config_info.get('refrence', 'fasta_file')
    adapter1 = config_info.get('adapter', 'adapter1')
    adapter2 = config_info.get('adapter', 'adapter2')
    CPU = int(config_info.get('resource', 'threads'))
    samtools = config_info.get('program', 'samtools')
    bgzip = config_info.get('program', 'bgzip')
    bwa = config_info.get('program', 'bwa')
    cutadapt = config_info.get('program', 'cutadapt')
    bamdst = config_info.get('program', 'bamdst')

    log_file = os.path.join(out_path, 'log.log')
    logging.basicConfig(filename=log_file, \
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s-%(funcName)s', \
                        level=logging.INFO)

    if os.system("which {} >/dev/null 2>&1".format(samtools)) != 0:
        logging.error('samtools命令未找到，请核查系统是否已安装samtools以及配置文件指定路径是否正确！')
        sys.exit("Command {} does not exist, stop running!".format(samtools))
    if os.system("which {} >/dev/null 2>&1".format(bgzip)) != 0:
        logging.error('bgzip命令未找到，请核查系统是否已安装bgzip以及配置文件指定路径是否正确！')
        sys.exit("Command {} does not exist, stop running!".format(bgzip))
    if os.system("which {} >/dev/null 2>&1".format(bwa)) != 0:
        logging.error('bwa命令未找到，请核查系统是否已安装bwa以及配置文件指定路径是否正确！')
        sys.exit("Command {} does not exist, stop running!".format(bwa))
    if os.system("which {} >/dev/null 2>&1".format(cutadapt)) != 0:
        logging.error('cutadapt命令未找到，请核查系统是否已安装cutadapt以及配置文件指定路径是否正确！')
        sys.exit("Command {} does not exist, stop running!".format(cutadapt))

    logging.info('分析任务开始')
    logging.info('开始解析输入文件')
    resolve_input = InputResolve(input_file)
    control_list, case_list, gender_dict, bam_dict = resolve_input.resolve_json()
    sample_list = control_list + case_list
    sample_num = len(sample_list)

    cpu_count = int(multiprocessing.cpu_count())
    if sample_num <= math.ceil(cpu_count * 0.8):
        threads = sample_num
    else:
        threads = math.ceil(cpu_count * 0.8)

    if CPU <= math.ceil(cpu_count * 0.8 / threads):
        CPU = CPU
    else:
        CPU = math.ceil(cpu_count * 0.8 / threads)

    try:
        logging.info('开始进行样本重比对')
        p = multiprocessing.Pool(threads)
        for s in sample_list:
            in_bam = bam_dict[s]
            bam_re_alignment = BamReAlignment()
            p.apply_async(func=bam_re_alignment.run_re_alingment, args=(samtools, bgzip, cutadapt, bwa, s, in_bam, \
                                                                        adapter1, adapter2, fasta_file, CPU, out_path))
        p.close()
        p.join()

        logging.info('开始进行样本质量值计算')
        p = multiprocessing.Pool(threads)
        for s in sample_list:
            in_bam1 = bam_dict[s]
            in_bam2 = os.path.join(out_path, '{}.cutadapt.bam'.format(s))
            out_path1 = os.path.join(out_path, 'QualityControl', '{}'.format(s), 'raw_qc')
            out_path2 = os.path.join(out_path, 'QualityControl', '{}'.format(s), 're_ali_qc')
            if not os.path.exists(out_path1):
                os.makedirs(out_path1)
            if not os.path.exists(out_path2):
                os.makedirs(out_path2)
            run_qc = QualityControl()
            p.apply_async(func=run_qc.run_qc_calculate, args=(s, bamdst, target_file, out_path1, in_bam1))
            p.apply_async(func=run_qc.run_qc_calculate, args=(s, bamdst, target_file, out_path2, in_bam2))
        p.close()
        p.join()

        logging.info('开始进行样本基本信息计算')
        p = multiprocessing.Pool(threads)
        for s in sample_list:
            bam_file = os.path.join(out_path, '{}.cutadapt.bam'.format(s))
            s_out = os.path.join(out_path, '{}.stat.txt'.format(s))
            reads_statistic = ReadsStatisticReAlignment()
            p.apply_async(func=reads_statistic.reads_stat, args=(target_file, bam_file, fasta_file, s_out))
        p.close()
        p.join()

        logging.info('开始进行样本计算结果合并')
        merge_stat = MergeStatistic(control_list, case_list, out_path)
        merge_df = merge_stat.merge_file()

        if filter == 'True':
            logging.info('开始过滤低覆盖位点')
            filter_gene = FilterGene(merge_df, control_list, out_path)
            merge_df = filter_gene.filters()

        logging.info('开始进行结果分析计算')
        result_calculate = ResultCalculate(merge_df, control_list, case_list, out_path, refgene_file)
        result_df, gender_calculate = result_calculate.result()

        logging.info('开始进行质控结果合并')
        merge_qc_stat = MergeQCResult(sample_list, gender_dict, gender_calculate, out_path)
        merge_qc_stat.merge_qc()

        logging.info('开始进行分析结果作图')
        if refgene == 'False':
            result_df = result_df
        else:
            refgene = refgene_file
            refgene_filter = FilterPlotGene(result_df, refgene)
            result_df = refgene_filter.ref_filters()
        plot_CN = PlotCN(result_df, out_path)
        plot_CN.plot()

        logging.info('开始打包分析结果')
        zipresult = PackResult(case_list, out_path)
        zipresult.zip_file()

        logging.info('分析任务结束')

    except Exception as e:
        logging.error(e)

if __name__ == '__main__':
    main()
