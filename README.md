# TorrentMLPA2CNV

## 说明文档：

- 参数：
  
  - -c 输入的配置文件（config.init）
  
  - -i 输入文件（json格式）
  
  - -b 目标区域文件（bed格式）
  
  - -g 内参基因列表文件（基因list）
  
  - -o 输出结果文件夹（文件夹路径）
  
  - -f 是否过滤覆盖度较低的基因区域，默认为False（不过滤）
  
  - -r 图显去掉内参基因，默认为True（去掉）

- 文件内容格式见示例文件

- 程序说明
  
  - Torrent_MLPA_ReAlignment.py #软件主程序
  
  - InputResolve.py #输入文件解析程序
  
  - BamReAlignment.py #bam文件重比对程序
  
  - QualityControl.py #质控程序
  
  - ReadsStatisticReAlignment.py #样本基本信息计算程序
  
  - MergeStatistic.py #合并每个样本计算结果程序
  
  - FilterGene.py #低覆盖度位点过滤程序
  
  - ResultCalculate.py #最终结果计算程序
  
  - MergeQCResult.py #质控结果合并程序
  
  - FilterPlotGene.py #内参基因图显不显示过滤程序
  
  - PlotCN.py #分析结果作图程序
  
  - PackResult.py #分析结果打包程序

- 运行命令示例：`python3 Torrent_MLPA_ReAlignment.py -c config.init -i sample.json -b torrentMLPA_hg19-v2.0.bed -g refgene-v2.0.txt -o .`
