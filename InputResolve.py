#!/usr/bin/python3
import re
import os
import sys
import json
import pandas as pd

class InputResolve:
    '''
    解析输入样本文件
    名称：InputResolve
    输入文件：输入信息json文件
    输出：样本文件信息
    '''

    def __init__(self, json_file):
        self.json_file = json_file

    def resolve_json(self):
        with open(self.json_file, 'r') as j:
            data_dict = json.load(j)

        control_list = []
        case_list = []
        gender_dict = {}
        bam_dict = {}
        for info in data_dict['samples']:
            name = info['name']
            bam = info['bam']
            gender = info['gender']
            label = info['label']
            bam_dict[name] = bam
            gender_dict[name] = gender
            if label == 'control':
                control_list.append(name)
            elif label == 'case':
                case_list.append(name)
            else:
                print ('样本分类选择错误，请重新选择')

        return control_list, case_list, gender_dict, bam_dict
