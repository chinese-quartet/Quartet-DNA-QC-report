# Dnaseqc
Guideline for Quartet DNA and Methylation QC Pipeline

# Quartet DNA Quality Control Evaluation

This repository contains the quality control (QC) evaluation pipelines for DNA sequencing data and DNA methylation sequencing data generated by the Quartet project. Our aim is to assess and ensure the quality of the sequencing data through comprehensive QC processes. The project is divided into two main parts: the DNA Data QC Pipeline and the DNA Methylation Data QC Pipeline.

## Table of Contents

- [Introduction](#introduction)
- [DNA Data QC Pipeline](#dna-data-qc-pipeline)
  - [Requirements](#requirements)
  - [Usage](#usage)
- [DNA Methylation Data QC Pipeline](#dna-methylation-data-qc-pipeline)
  - [Requirements](#requirements-1)
  - [Usage](#usage-1)
- [Contributing](#contributing)

## Introduction

This project focuses on the quality control of DNA and DNA methylation sequencing data. Quality control is essential to validate the sequencing data for further analysis and research. Our pipelines leverage various software tools to analyze and generate QC reports, ensuring the reliability and accuracy of the data.

## DNA Data QC Pipeline

The DNA Data QC Pipeline starts with VCF files, using hap.py and VBT software to analyze and transform the data format, followed by dnaseqc to generate quality control reports.

### Requirements
- hap.py
- VBT
- dnaseqc
### installation
- hap.py
(https://github.com/Illumina/hap.py)
- VBT
(https://github.com/sbg/VBT-TrioAnalysis)
- rtg-tools 下载地址
（https://github.com/RealTimeGenomics/rtg-tools）
- Chinese Quartet 标准数据集下载地址
（https://chinese-quartet.org/#/reference-datasets/download）

### Usage

```bash
# Step 1: Analyze and transform data format with hap.py and VBT
## 利用标准数据集计算F1 score
### truth.vcf为Quartet标准数据集，confident.bed为高置信区间，reference.fa为参考基因组文件
hap.py truth.vcf query.vcf -f confident.bed -o output_prefix -r reference.fa

### wes数据需要先与探针取交集
### 推荐使用samtools 计算bed文件交集
samtools intersect -a target.bed -b reference_dataset.bed > intersect.bed
hap.py truth.vcf query.vcf -f intersect.bed -o output_prefix -r reference.fa

### 得到的输出文件.summary.csv文件

## 使用multiqc整合同一批次D5、D6、F7和M8，hap计算结果
multiqc ./dir_to_four_summary_csv_file
```
## 文件格式修改参考 extract_hap_result.ipynb文件，输出variants.calling.qc.txt文件

### variants.calling.qc.txt文件示例
| Sample  | SNV number | INDEL number | SNV precision | INDEL precision | SNV recall | INDEL recall |
| :---: | :--: | :------: | :------:|  :------:|  :------:|  :------:|
| LCL5_UU_Illumina_D5_20230629_20171028_EATRISPLUS_UU_LCL5_hc  |  3855821  | 980430  | 99.73| 98.44 | 99.25 | 98.59 |
| LCL6_UU_Illumina_D6_20230629_20171028_EATRISPLUS_UU_LCL6_hc  |  3861023  | 976804  | 99.74| 98.51 | 99.38 | 98.68 |

```bash
## 计算孟德尔遗传率，需要先安装VBT

## 利用rtgtools 合并同一批次的D5、D6、F7、M8文件，用于后续计算孟德尔遗传率
rtg vcfmerge --force-merge-all -o ${project}.family.vcf.gz ${D5_vcf} ${D6_vcf} ${F7_vcf} ${M8_vcf}
### 示例
## /Volumes/移动硬盘/FD/software/rtg-tools/rtg-tools-3.12.1/rtg vcfmerge --force-merge-all -o Quartet_DNA_ILM_Nova_WUX_1.family.vcf.gz Quartet_DNA_ILM_Nova_WUX_LCL5_1_20171024_RAW.vcf.gz Quartet_DNA_ILM_Nova_WUX_LCL6_1_20171024_RAW.vcf.gz Quartet_DNA_ILM_Nova_WUX_LCL7_1_20171024_RAW.vcf.gz Quartet_DNA_ILM_Nova_WUX_LCL8_1_20171024_RAW.vcf.gz
## 解压缩，用于vbt输入
gunzip ${project}.family.vcf.gz

## 运行vbt.sh脚本（测试使用，可根据具体需求进行修改）
bash vbt.sh
## 输出 ${family_name}.D5.txt、${family_name}.D6.txt和${family_name}.consensus.txt 三个文件

## 运行merge_two_family_with_genotype.py脚本
python merge_two_family_with_genotype.py -LCL5 ${family_name}.D5.txt -LCL6 ${family_name}.D6.txt -genotype ${family_name}.consensus.txt -family {family_name}
## 输出${family_name}.summary.txt文件

```
### 输出${family_name}.summary.txt文件示例
| Family  | Total_Variants | Mendelian_Concordant_Variants | Mendelian_Concordance_Rate |
| :---: | :--: | :------: | :------:|
| EATRISPLUS_UU.INDEL  |  1294054  | 1178598  | 0.910779611979|
| EATRISPLUS_UU.SNV  |  5034285  | 4868250   | 0.967019149691|


# Step 2: Generate QC report with dnaseqc
```R
##下载并安装R包dnaseqc和相关依赖
library(devtools)
devtools::install_github("markx945/Dnaseqc/dnaseqc")
library(dnaseqc)

## 读取F1_score计算和孟德尔遗传率计算结果，以历史数据为例
variant_qc <- system.file("example","variants.calling.qc.txt",package = "dnaseqc")
mendelian_qc <- system.file("example","EATRISPLUS_UU.summary.txt",package = "dnaseqc")

## 输入测序类型，“WGS”或“WES”，计算得到DNAseq QC指标
result = Dnaseqc(variant_qc_file = variant_qc, mendelian_qc_file = mendelian_qc, data_type = "WGS")

## 生成报告
### 从R包中读取报告模板路径
doc_path <- system.file("extdata","Quartet_temp.docx",package = "dnaseqc")
### 指定路径生成DNAseq
GenerateDNAReport(DNA_result = result,doc_file_path = doc_path,output_path = './DNAseq/' )

```


## DNA Methylation Data QC Pipeline
The DNA Methylation Data QC Pipeline begins with processed methylation sequencing data to generate quality control reports.

### Usage

```R
<command to generate QC report>
```

## Contributing
We welcome contributions to improve the QC pipelines. Please feel free to fork the repository, make changes, and submit pull requests. For major changes, please open an issue first to discuss what you would like to change.





