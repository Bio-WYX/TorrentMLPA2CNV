#!/usr/bin/python3
import re
import os
import sys
import pysam
import difflib
import multiprocessing

class ReadsStatisticReAlignment:
    '''
    统计样本测序结果中reads的比对信息
    名称：ReadsStatistic
    输入文件：bam
    输出：txt
    '''
    # def __init__(self, bed_file, bam_file, out_file):
    #     self.bed_file = bed_file
    #     self.bam_file = bam_file
    #     self.out_file = out_file

    @staticmethod
    def rev_com_seq(seq):
        '''
        获得反向互补序列
        Args:
            seq: sequence
            Returns: reverse complement sequence
        '''
        seqreverse = seq[::-1]
        transtable = str.maketrans('ATGCatcgNn', 'TACGtagcNn')
        finalseq = seqreverse.translate(transtable)
        return finalseq

    @staticmethod
    def mapping_info(reads, sequence):
        '''
        比对序列分析分类及统计
        Args:
            reads: sequence
            tag: alignment info
            length: target length
        Returns: correct_reads, failed_reads
        '''

        reads_match = difflib.get_close_matches(reads, sequence, 1, cutoff=0.9)
        correct_reads = 0
        failed_reads = 0
        if reads_match:
            correct_reads = 1
        else:
            failed_reads = 1
        return correct_reads, failed_reads

    #@staticmethod
    def reads_stat(self, bed_file, bam_file, fasta_file, out_file):
        '''
        统计样本目标区域的比对结果
        Returns: out file
        '''
        OUT = open(out_file, 'w')
        OUT.write('CHR\tSTART\tEND\tGENE\tALL_READS\tCORRECT_READS\tFAILED_READS\tCORRECT_RATIO(%)\n')
        samfile = pysam.AlignmentFile(bam_file, "rb")
        fastafile = pysam.Fastafile(fasta_file)
        with open(bed_file, 'r') as bed:
            for i in bed:
                info  = i.strip().split('\t')
                chr = info[0]
                start = int(info[1])
                end = int(info[2])
                gene = info[3]
                length = end - start + 1

                sequence = fastafile.fetch(chr, start, end + 1)
                sequence = [sequence.upper()]
                reads_num = samfile.count(contig = chr, start = start, stop = end)
                allreads = samfile.fetch(contig = chr, start = start, stop = end)
                correct_reads = 0
                failed_reads = 0
                for read in allreads:
                    tag = read.cigarstring
                    reads1 = read.seq.upper()
                    tag_match1 = re.match('(\d+)S\S+', tag)
                    tag_match2 = re.match('\S+?(\d+)S', tag)
                    if tag_match1:
                        left_num = int(tag_match1.group(1))
                        reads1 = reads1[left_num:]
                    if tag_match2:
                        right_num = int(tag_match2.group(1))
                        reads1 = reads1[:-right_num]

                    reads2 = self.rev_com_seq(reads1).upper()
                    correct_reads1, failed_reads1 = self.mapping_info(reads1, sequence)
                    correct_reads2, failed_reads2 = self.mapping_info(reads2, sequence)
                    if correct_reads1 >= correct_reads2:
                        correct_reads += correct_reads1
                        failed_reads += failed_reads1
                    else:
                        correct_reads += correct_reads2
                        failed_reads += failed_reads2

                if reads_num > 0:
                    correct_ratio = round(correct_reads/reads_num * 100, 2)
                else:
                    correct_ratio = 0
                OUT.write(f'{chr}\t{start}\t{end}\t{gene}\t{reads_num}\t{correct_reads}\t{failed_reads}\t{correct_ratio}\n')
        OUT.close()