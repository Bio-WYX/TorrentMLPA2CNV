#!/usr/bin/python3
import re
import os
import sys
import zipfile

class PackResult:
    '''
    分析结果文件打包
    名称：PackResult
    输入文件：样本列表与结果路径
    输出：所有结果打包文件
    '''

    def __init__(self, sample_list, out_path):
        self.sample_list = sample_list
        self.out_path = out_path

    def zip_file(self):
        '''
        分析结果文件打包
        Args:
            sample_list: 检测样本列表
            out_path: 输出结果路径
        Returns: 所有结果打包文件
        '''

        # 打包后的文件名
        zip_file_name = os.path.join(self.out_path, '分析结果.zip')

        # 创建一个新的ZIP文件
        with zipfile.ZipFile(zip_file_name, 'w') as my_zip:
            for sample in self.sample_list:
                file_path = os.path.join(self.out_path, '{}_result.pdf'.format(sample))
                my_zip.write(file_path, os.path.basename(file_path))

            file_qc = os.path.join(self.out_path, '质控结果.xlsx')
            my_zip.write(file_qc, os.path.basename(file_qc))

# s = ['sample10', 'sample12', 'sample13', 'sample15', 'sample16', 'sample1', 'sample3', 'sample4', 'sample6', 'sample7']
# o = sys.argv[1]
# zipresult = PackResult(s, o)
# zipresult.zip_file()