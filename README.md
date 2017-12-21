# CWAP(Cancer WES Analysis Pipeline)
# A WES pipline for BGISEQ-500

## Bioinformatics pipelines
生物信息流程（Bioinformatics pipelines）指二代测序数据分析过程中的一连串处理步骤，为了将原始数据转化为解读结果。例如，基因组测序后，将下机数据FASTQ文件转化为突变数据VCF文件的分析过程（Reads-to-variants workflows）。在肿瘤基因组研究中，为了准确、快速得到突变结果，一个科学、有效的生物信息流程非常重要。目前，已经有不少公开发表或华大内部开发的肿瘤基因组生物信息流程，例如Broad研究所开发的GATK最佳实践（GATK Best Practices）流程，但是还鲜有基于BGISEQ-500测序平台数据特点的流程，为了适应新形势的需要，我们开发了一个新的肿瘤基因组信息分析流程。

## Cancer WES Analysis Pipeline，CWAP
本流程全称为“肿瘤全外显子测序信息分析流程”（Cancer WES Analysis Pipeline，CWAP），版本为2017a。本流程的分析目标为将肿瘤和对照样品的全外显子测序下机数据分析得到突变结果，主要包括样本测序数据处理、突变检测及注释两个分析步骤。具体而言，分析过程为：（1）利用SOAPnuke对原始FASTQ文件过滤和质控，利用Edico或者BWA进行比对，利用Picard工具包对比对结果进行标记重复（BWA比对时），然后用GATK进行局部重比对和碱基质量值分数矫正。（2）突变检测包括遗传突变检测（包括SNP（GATK）和Indel（GATK））和体细胞突变检测（包括SNV（Mutect，Varscan，Mutect2），Sindel（GATK，Platypus，Mutect2）和SCNV（ADTEx）），突变注释支持AnnoDB、ANNOVAR和VEP三个软件，突变过滤是注释完成后的可选步骤，据项目实际情况而定。本流程基于GATK最佳实践框架，使用DRAGEN处理器（Edico Genomics公司）加速比对、标记重复读段和突变检测等过程，本流程生成的任务脚本适用于华大基因的高性能计算集群，可以通过Pymonitor软件投递有依赖关系的任务列表并自动监控任务状态。
