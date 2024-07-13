#!/usr/bin/python3
import re
import os
import sys
import pysam
import difflib
import multiprocessing

class BamReAlignment:
    '''
    测序结果Bam文件去掉adapter重新比对
    名称：BamReAlignment
    输入文件：bam
    输出：rebam
    '''
    # def __init__(self, bed_file, bam_file, out_file):
    #     self.bed_file = bed_file
    #     self.bam_file = bam_file
    #     self.out_file = out_file

    @staticmethod
    def rev_com_seq(seq):
        '''
        获得反向、互补、反向互补序列
        Args:
            seq: sequence
            Returns: reverse complement sequence
        '''

        transtable = str.maketrans('ATGCatcgNn', 'TACGtagcNn')
        seqreverse = seq[::-1]
        seqcomplement = seq.translate(transtable)
        seqrevcom = seqreverse.translate(transtable)
        return seqreverse, seqcomplement, seqrevcom

    @staticmethod
    def bam_to_fq(samtools, bgzip, sample_name, in_bam, CPU, out_path):
        '''
        bam转化为fastq文件
        Args:
            samtools: samtools命令路径
            bgzip: bgzip命令路径
            sample_name: 样本名称
            in_bam: 输入bam文件
            CPU: 调用线程数
            out_path: 输出路径
        Returns: fastq文件
        '''
        out_fq = os.path.join(out_path, '{}.fq'.format(sample_name))
        bam2fq_cmd = '{} fastq -@ {} {} -c 9 >{} 2>/dev/null'.format(samtools, CPU, in_bam, out_fq)
        bgzip_cmd = '{} -f -@ {} {} >/dev/null 2>&1'.format(bgzip, CPU, out_fq)
        os.system(bam2fq_cmd)
        os.system(bgzip_cmd)

    @staticmethod
    def cut_adapter(cutadapt, sample_name, adapter1, adapter2, \
                    seqreverse1, seqreverse2, seqcomplement1, seqcomplement2, seqrevcom1, seqrevcom2, CPU, out_path):
        '''
        去除接头序列
        Args:
            cutadapt: cutadapt命令路径
            sample_name: 样本名称
            adapter1: 接头1
            adapter2: 接头2
            seqreverse1: 接头1反向序列
            seqreverse2: 接头2反向序列
            seqcomplement1: 接头1互补序列
            seqcomplement2: 接头2互补序列
            seqrevcom1: 接头1反向互补序列
            seqrevcom2: 接头2反向互补序列
            CPU: 调用线程数
            out_path: 输出路径
        Returns: 去除接头序列后的fastq文件
        '''
        out_fq = os.path.join(out_path, '{}.fq.gz'.format(sample_name))
        out_adapt = os.path.join(out_path, '{}.cutadapt.fq.gz'.format(sample_name))
        cutadapt_cmd = '{} -b {} -b {} -b {} -b {} -b {} -b {} -b {} -b {} -e 0.1 -n 2 -o {} {} -j {} >/dev/null 2>&1'\
            .format(cutadapt, adapter1, adapter2, seqreverse1, seqreverse2, seqcomplement1, seqcomplement2, \
                    seqrevcom1, seqrevcom2, out_adapt, out_fq, CPU)
        os.system(cutadapt_cmd)

    @staticmethod
    def re_alignments(bwa, samtools, sample_name, CPU, fasta_file, out_path):
        '''
        重新比对
        Args:
            bwa: bwa命令路径
            samtools: samtools命令路径
            sample_name: 样本名称
            CPU: 调用线程数
            fasta_file: fasta文件
            out_path: 输出路径
        Returns: 重新比对后的bam文件
        '''
        out_adapt = os.path.join(out_path, '{}.cutadapt.fq.gz'.format(sample_name))
        adapt_bam = os.path.join(out_path, '{}.cutadapt.bam'.format(sample_name))
        bwa_cmd = '{} mem -R "@RG\\tID:{}\\tSM:{}\\tPL:Iontorrent" -t {} {} {} 2>/dev/null | \
        samtools view -bS -@ {} - 2>/dev/null | samtools sort -@ {} -o {} >/dev/null 2>&1'.format\
            (bwa, sample_name, sample_name, CPU, fasta_file, out_adapt, CPU, CPU, adapt_bam)
        samtoos_cmd = '{} index {} >/dev/null 2>&1'.format(samtools, adapt_bam)
        os.system(bwa_cmd)
        os.system(samtoos_cmd)

    def run_re_alingment(self, samtools, bgzip, cutadapt, bwa, sample_name, in_bam, \
                         adapter1, adapter2, fasta_file, CPU, out_path):
        '''
        运行程序
        Args:
            samtools: samtools命令路径
            bgzip: bgzip命令路径
            cutadapt: cutadapt命令路径
            bwa: bwa命令路径
            sample_name: 样本名称
            adapter1: 接头1
            adapter2: 接头2
            fasta_file: fasta文件
            CPU: 调用线程数
            out_path: 输出路径
        Returns: 运行结果的bam文件
        '''
        seqreverse1, seqcomplement1, seqrevcom1 = self.rev_com_seq(adapter1)
        seqreverse2, seqcomplement2, seqrevcom2 = self.rev_com_seq(adapter2)
        self.bam_to_fq(samtools, bgzip, sample_name, in_bam, CPU, out_path)
        self.cut_adapter(cutadapt, sample_name, adapter1, adapter2, \
                         seqreverse1, seqreverse2, seqcomplement1, seqcomplement2, seqrevcom1, seqrevcom2, CPU, out_path)
        self.re_alignments(bwa, samtools, sample_name, CPU, fasta_file, out_path)
