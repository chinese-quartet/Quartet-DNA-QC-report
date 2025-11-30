#! /bin/bash

export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

## 指定测试数据批次
family_name="Quartet_DNA_ILM_Nova_WUX_1"
family_vcf="Quartet_DNA_ILM_Nova_WUX_1.family.vcf"


nt=16

## 构建D5、F7、M8 ped文件
echo -e "${family_name}\tLCL8\t0\t0\t2\t-9\n${family_name}\tLCL7\t0\t0\t1\t-9\n${family_name}\tLCL5\tLCL7\tLCL8\t2\t-9" > ${family_name}.D5.ped

mkdir VBT_D5
## 计算D5家系孟德尔遗传率
/home/cfff_r2636/data/software/VBT/VBT-TrioAnalysis/vbt mendelian -ref /home/cfff_r2636/data/reference/hg38/genome/hg38.fa -mother ${family_vcf} -father ${family_vcf} -child ${family_vcf} -pedigree ${family_name}.D5.ped -outDir VBT_D5 -out-prefix ${family_name}.D5 --output-violation-regions -thread-count $nt

## 文件重命名
cat VBT_D5/${family_name}.D5_trio.vcf > ${family_name}.D5.vcf

## 构建D6、F7、M8 ped文件
echo -e "${family_name}\tLCL8\t0\t0\t2\t-9\n${family_name}\tLCL7\t0\t0\t1\t-9\n${family_name}\tLCL6\tLCL7\tLCL8\t2\t-9" > ${family_name}.D6.ped

mkdir VBT_D6
## 计算D6家系的孟德尔遗传率
/home/cfff_r2636/data/software/VBT/VBT-TrioAnalysis/vbt mendelian -ref /home/cfff_r2636/data/reference/hg38/genome/hg38.fa -mother ${family_vcf} -father ${family_vcf} -child ${family_vcf} -pedigree ${family_name}.D6.ped -outDir VBT_D6 -out-prefix ${family_name}.D6 --output-violation-regions -thread-count $nt

cat VBT_D6/${family_name}.D6_trio.vcf > ${family_name}.D6.vcf

## 文件重命名
cat ${family_name}.D5.vcf | grep -v '##' > ${family_name}.D5.txt
cat ${family_name}.D6.vcf | grep -v '##' > ${family_name}.D6.txt
cat ${family_vcf} | grep -v '##' | awk '
		BEGIN { OFS = "\t" }
		NF > 2 && FNR > 1 { 
		    for ( i=9; i<=NF; i++ ) { 
		        split($i,a,":") ;$i = a[1];
		    } 
		} 
		{ print }
' > ${family_name}.consensus.txt
