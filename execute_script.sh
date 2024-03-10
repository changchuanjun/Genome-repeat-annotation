
#step1 单独运行
#=================================
#查看数据库中Malus分类的family数目
#singularity exec /data/Erick_Tong/02Data_analysis_project/09Genome_assembly/03geneke/TETools202309.sif famdb.py lineage -ad Malus

#统计Malus families数目
#singularity exec /data/Erick_Tong/02Data_analysis_project/09Genome_assembly/03geneke/TETools202309.sif famdb.py lineage -ad --format totals Malus
#================================


#================================global variable
threads=8
pa=2 #1pa=4threads
TEsingularity=/data/Erick_Tong/02Data_analysis_project/09Genome_assembly/03geneke/TETools202309.sif
genome_file=Fuji_Ral.fasta #这里结尾必须是.fasta
genome_prefix=Fuji_Ral
timesj=$(date +"%Y_%m_%d")
#Malus(基于已有数据库RepBase+Dfam  运行repeatmasker)
#===============================global variable


#step2 RepeatModeler
#=====================================================================

# 建库
singularity exec -B $PWD:$PWD -B /data:/data -B /biodata:/biodata ${TEsingularity} \
	BuildDatabase -name ${genome_prefix}_BuildDatabase_${timesj} ${genome_file}

# 运行 RepeatModeler
singularity exec -B $PWD:$PWD -B /data:/data -B /biodata:/biodata ${TEsingularity} \
	RepeatModeler \
		-database ${genome_prefix}_BuildDatabase_${timesj} \
  		-threads ${threads} \
  		-LTRStruct  1>repeatmodeler.log 2>&1

mkdir 00_${genome_prefix}_BuildDatabase_${timesj}
mv ${genome_prefix}_BuildDatabase_${timesj}.* 00_${genome_prefix}_BuildDatabase_${timesj}

mkdir 00_${genome_prefix}_denovo_repeatlibrary_${timesj}
mv ${genome_prefix}_BuildDatabase_${timesj}-* 00_${genome_prefix}_denovo_repeatlibrary_${timesj}

# 如果运行终止可以设置recoverDir自动检查运行
# singularity exec ../container/TETools202309.sif RepeatModeler -database Bra \
#    -threads 20  -LTRStruct -recoverDir RM_xxxx  1>repeatmodeler.log 2>&1
#=====================================================================


#step3 RepeatMasker
#=====================================================================


## 基于已有数据库RepBase+Dfam  运行repeatmasker
## -e 比对软件；-pa 并行任务数，实际占用线程为pa*4；-gff,输出gff；-s,缓慢模式；-a,输出align文件；-species,物种范围；-dir,结果目录；
singularity exec -B $PWD:$PWD -B /data:/data -B /biodata:/biodata ${TEsingularity} \
	RepeatMasker \
		-e rmblast \
		-pa ${pa} \
		-gff \
		-s \
		-a \
		-species Malus \
		-dir ./01_${genome_prefix}_repeatmasker_result1_RepBaseDfam_${timesj} \
		${genome_file} \
		> ${genome_prefix}_repeatmasker_result1_RepBaseDfam.log 2>&1


# 基于denovo lib 运行repeatmasker
#-pa 5 5个并行（32线程）-s 慢速* -q 快速 -qq 高速 -gff/-a 生成.gff/.align文件 .tbl结果汇总 -nolow 不包括低复杂重复序列（后面基因结构注释）
singularity exec -B $PWD:$PWD -B /data:/data -B /biodata:/biodata ${TEsingularity} \
	RepeatMasker \
		-e rmblast \
		-pa ${pa} \
		-gff \
		-s \
		-a \
		-lib 00_${genome_prefix}_denovo_repeatlibrary_${timesj}/${genome_prefix}_BuildDatabase_${timesj}-families.fa \
		-dir 02_${genome_prefix}_repeatmasker_result2_denovolib_${timesj} \
		01_${genome_prefix}_repeatmasker_result1_RepBaseDfam_${timesj}/${genome_prefix}.fasta.masked \
		> ${genome_prefix}_repeatmasker_result2_denovolib.log 2>&1


# 合并结果repeatmasker_result1、repeatmasker_result2----------------------------(全部重复)
mkdir 03_${genome_prefix}_combine_repeatmasker_result12_${timesj}
#合并
singularity exec -B $PWD:$PWD -B /data:/data -B /biodata:/biodata ${TEsingularity} \
	perl /opt/RepeatMasker/util/combineRMFiles.pl \
		01_${genome_prefix}_repeatmasker_result1_RepBaseDfam_${timesj}/${genome_prefix}.fasta \
		02_${genome_prefix}_repeatmasker_result2_denovolib_${timesj}/${genome_prefix}.fasta.masked \
		03_${genome_prefix}_combine_repeatmasker_result12_${timesj}/${genome_prefix}


# 根据.out文件生成生成.gff文件
singularity exec -B $PWD:$PWD -B /data:/data -B /biodata:/biodata ${TEsingularity} \
	rmOutToGFF3.pl \
		03_${genome_prefix}_combine_repeatmasker_result12_${timesj}/${genome_prefix}.out \
		> 03_${genome_prefix}_combine_repeatmasker_result12_${timesj}/${genome_prefix}.out.gff

#软屏蔽
singularity exec -B $PWD:$PWD -B /data:/data -B /biodata:/biodata ${TEsingularity} \
	maskFile.pl \
		-fasta ${genome_file} \
		-annotations 03_${genome_prefix}_combine_repeatmasker_result12_${timesj}/${genome_prefix}.out \
		-softmask

mv ${genome_file}.masked 03_${genome_prefix}_combine_repeatmasker_result12_${timesj}/${genome_prefix}.soft_masked.fasta


#硬屏蔽
singularity exec -B $PWD:$PWD -B /data:/data -B /biodata:/biodata ${TEsingularity} \
	maskFile.pl \
		-fasta ${genome_file} \
		-annotations 03_${genome_prefix}_combine_repeatmasker_result12_${timesj}/${genome_prefix}.out


mv ${genome_file}.masked 03_${genome_prefix}_combine_repeatmasker_result12_${timesj}/${genome_prefix}.hard_masked.fasta

# 绘图
#calcDivergenceFromAlign.pl -s  hap02-repeat.div  rmout_combine02/hap02.align
#createRepeatLandscape.pl  -div hap02-repeat.div -g 660000000 > hap02-repeat.div.html
#-g 基因组大小

# 合并结果repeatmasker_result1、repeatmasker_result2----------------------------(全部重复)
####################################################

####################################################
# 合并结果repeatmasker_result1、repeatmasker_result2----------------------------(去除低复杂重复重复)
mkdir 04_delete03_Simple_repeat.Low_complexity_masked_${timesj}

##01去除低复杂重复

grep -v "Simple_repeat" 03_${genome_prefix}_combine_repeatmasker_result12_${timesj}/${genome_prefix}.out | grep -v "Low_complexity" > 04_delete03_Simple_repeat.Low_complexity_masked_${timesj}/${genome_prefix}_delete_sample_repeat_masked.out

# 根据.out文件生成生成.gff文件
singularity exec -B $PWD:$PWD -B /data:/data -B /biodata:/biodata ${TEsingularity} \
	rmOutToGFF3.pl \
		04_delete03_Simple_repeat.Low_complexity_masked_${timesj}/${genome_prefix}_delete_sample_repeat_masked.out \
		> 04_delete03_Simple_repeat.Low_complexity_masked_${timesj}/${genome_prefix}_delete_sample_repeat_masked.out.gff


#软屏蔽
singularity exec -B $PWD:$PWD -B /data:/data -B /biodata:/biodata ${TEsingularity} \
        maskFile.pl \
                -fasta ${genome_file} \
                -annotations 04_delete03_Simple_repeat.Low_complexity_masked_${timesj}/${genome_prefix}_delete_sample_repeat_masked.out \
                -softmask

mv ${genome_file}.masked 04_delete03_Simple_repeat.Low_complexity_masked_${timesj}/${genome_prefix}.soft_masked.fasta


#硬屏蔽
singularity exec -B $PWD:$PWD -B /data:/data -B /biodata:/biodata ${TEsingularity} \
        maskFile.pl \
                -fasta ${genome_file} \
                -annotations 04_delete03_Simple_repeat.Low_complexity_masked_${timesj}/${genome_prefix}_delete_sample_repeat_masked.out


mv ${genome_file}.masked 04_delete03_Simple_repeat.Low_complexity_masked_${timesj}/${genome_prefix}.hard_masked.fasta

# 合并结果repeatmasker_result1、repeatmasker_result2----------------------------(去除低复杂重复重复)

# 绘图
#calcDivergenceFromAlign.pl -s  hap02-repeat.div  rmout_combine02/hap02.align
#createRepeatLandscape.pl  -div hap02-repeat.div -g 660000000 > hap02-repeat.div.html
#-g 基因组大小
###################################################




