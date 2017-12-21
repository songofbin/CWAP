#!/usr/bin/perl -w
#A WES pipline for BGISEQ-500
#author: Song Bin, songbin@genomics.cn
use strict;
use warnings;
use FindBin '$Bin';
use FileHandle;
use File::Basename qw/dirname basename/;

die "Usage:\tperl $0 <Fastq list> <Config file> <Outdir> \n"if(@ARGV!=3);

my($fqlist,$config,$outdir)=@ARGV;
`mkdir -p $outdir`;

################################### main parameters ###################################
my(%configs,%Splitbeds,%SplitbedCount);
my(%rawfqs,%cleanfqs,%cleanList);
my(%bwaBams,%rmdupBams,%realnBams,%brecalBams,%sampleBams,%mergeBams);
my(%patientIds,%samplePairs,%patients2samples,%edge_sampair); 
my(%maxmems); ## max memory use for co-realnment
my @chrs=();
foreach my $i("M",1..22,"X","Y"){push @chrs,"chr$i";}

################################ Process sequence #######################################
&readConfig($config);

&readFqlist($fqlist);

if ($configs{FastqClean} =~ /true/i){
    &FastqClean(\%rawfqs,\%cleanfqs);
}

if ($configs{Aligment} ne "false"){
    if($configs{Aligment} eq "edico"){
        &Aligment_Edico(\%cleanList);
    }elsif($configs{Aligment} =~ /bwa/){
        &Aligment_BWA(\%cleanfqs);
    }
}

if ($configs{MarkDuplicates} eq "true"){
    &MarkDuplicates(\%bwaBams);
}

if(defined $configs{Corealn} && $configs{RealnType} eq "co-realn"){
    &PatientRealn($configs{Corealn});
    &BaseRecalibrator(\%realnBams);
}elsif($configs{RealnType} eq "realn"){
    &GATKrealnRecal(\%rmdupBams);
}

if ($configs{mergeBrecalBAM} =~ /true/i){
    &mergeBrecalBAM(\%brecalBams,\%mergeBams);
}

if ($configs{bamdstQC} =~ /true/i){
    &bamdstQC(\%sampleBams);
}

if ($configs{TargetQC} =~ /true/i){
    &TargetQC(\%sampleBams);
}

#if ($configs{DNA_damage} =~ /true/i){
#    &DNA_damage(\%sampleBams);
#}

######################### Germline variant calling: snps,indels ##########################
if($configs{gatkDISCOVERY}=~/true/i){
    foreach my $samp (sort keys %brecalBams){
        &makedir("$outdir/$samp/gatkDISCOVERY");
        foreach my $chr (sort keys %{$brecalBams{$samp}}){
            &gatkDISCOVERY($samp,$brecalBams{$samp}{$chr},$chr);
        }
        &mergeGATKVars($samp);
    }
}

if($configs{gatkGVCF}=~/true/i){
    foreach my $samp (sort keys %brecalBams){
	&makedir("$outdir/$samp/gatkGVCF");
        if($configs{Aligment} eq "edico"){
            &gatkGVCF_GenotypeGVCFs($samp);
        }else{
            foreach my $chr (sort keys %{$brecalBams{$samp}}){
                &gatkGVCF_sample($samp,$brecalBams{$samp}{$chr},$chr);
            }
            &gatkGVCF_merge($samp);
        }
        &gatkGVCF_process($samp);        
    }

    &makedir("$outdir/combine/gatkGVCF") if(!-e "$outdir/combine/gatkGVCF");
    if($configs{Aligment} eq "edico"){
        &gatkGVCF_combine_bysample("$outdir/combine/gatkGVCF");
    }else{
        &gatkGVCF_combine_bychr("$outdir/combine/gatkGVCF");
        &gatkGVCF_merge("combine");
    }
    &gatkGVCF_VQSR("combine");
    &gatkGVCF_process("combine");

    if($configs{patients}){
        open IN,"$configs{patients}" or die $!;
        while(<IN>){
                chomp;
                my @F = split /\t/;
                my $p = shift @F;
                @{$patients2samples{$p}} = @F;
        }
        close IN;

        gatkGVCF_combineGTs(\%patients2samples,\%configs,$outdir); 
    }
}

######################## Somatic variant calling: Ssnv, Sindel #########################
if($configs{somatic}){
    open SMTF,"$configs{somatic}" or die $!;
    while(<SMTF>){
        chomp;
        my ($patient,$normal,$tumor)=(split /\s+/)[0,1,2];
        die "bam files of $normal don't exists!\n" if(!exists $brecalBams{$normal});
        die "bam files of $tumor don't exists!\n" if(!exists $brecalBams{$tumor});	
        my $samplePair="$normal-VS-$tumor";      
        `mkdir -p $outdir/somatic/$samplePair` if(!-e "$outdir/somatic/$samplePair");
        $samplePairs{$samplePair}{normal}=$normal;
        $samplePairs{$samplePair}{tumor} =$tumor;
    }
    close SMTF;

    if($configs{mutectSNV}=~/true/i){
	&somatic_batch_byChr(\%brecalBams,\%samplePairs,\&mutectSNV,\&mutectSNV_process,"mutectSNV");
    }

    if($configs{mutect2}=~/true/i){
    	if($configs{seqType} =~ /WES/i){
                &somatic_batch_byRegion(\%brecalBams,\%samplePairs,\&RunMutect2,\&mutect2_process,"mutect2");
        }elsif($configs{seqType} =~ /WGS/i){
                &somatic_batch_byChr(\%brecalBams,\%samplePairs,\&RunMutect2,\&mutect2_process,"mutect2");
        }
    }

    if($configs{varscanSNV}=~/true/i){
	&somatic_batch_byChr(\%brecalBams,\%samplePairs,\&varscanSNV,\&varscanSNV_process,"varscanSNV");
    }

    if($configs{gatkSindel}=~/true/i){
	&somatic_batch_byChr(\%brecalBams,\%samplePairs,\&gatkSindel,\&gatkSindel_process,"gatkSindel");
    }

    if($configs{platypusSindel}=~/true/i){
	&somatic_batch_byChr(\%brecalBams,\%samplePairs,\&platypusSindel,\&platypusSindel_process,"platypusSindel");
    }

    if($configs{ADTExSCNV}=~/true/i){
	&somatic_batch_bySamplePair(\%brecalBams,\%samplePairs,\&ADTExSCNV);
    }
    
    if($configs{ABSOLUTE}=~/true/i){
	&somatic_batch_bySamplePair(\%brecalBams,\%samplePairs,\&ABSOLUTE);
    }
}

############################### monitor to control the pipeline #######################
if($configs{edgeList}=~/true/i){
	&edgeList();  ##creat scrips to monitor
	&edgeList2();
	&edgeList_Mutect2();
}
if($configs{shlist}=~/true/i){
        &shlist();
}
if($configs{rmlist}=~/true/i){
	&rmlist();
}

################################          sub function          #########################
## This Function is using to generate script. 
sub generateShell{
    my ($output_shell, $content, $finish_string) = @_;
    unlink glob "$output_shell.*" if($configs{rmlog} eq "true");
    $finish_string ||= "Still_waters_run_deep";
    open OUT,">$output_shell" or die "Cannot open file $output_shell:$!";
    print OUT "#!/bin/bash\n";
    print OUT "echo ==========start at : `date` ==========\n";
    print OUT "$content && \\\n";
    print OUT "echo ==========end at : `date` ========== && \\\n";
    print OUT "echo $finish_string 1>&2 && \\\n";
    print OUT "echo $finish_string > $output_shell.sign\n";
    close OUT;
}

sub makedir{
    my ($dir)=@_;
    if(!-e $dir){
    `mkdir -p $dir`;
    `mkdir $dir/shell`;
    `mkdir $dir/tmp`;
    `mkdir $dir/raw`;
    `mkdir $dir/result`;
    }
}

sub readConfig{
    my ($config) = @_;
    open CFG,"$config" or die $!;
    while(<CFG>){
        chomp;
        if($_=~/^\#/){next;}
        if($_=~/^(.*)=(.+);.?/){
		$configs{$1}=$2;
        }
    }
    close CFG;

    $configs{home}         ||=$Bin;
    $configs{bin}          ||="$configs{home}/bin";
    $configs{database}     ||="$configs{home}/database";
    $configs{lib}          ||="$configs{home}/lib"; 
    $configs{reference}    ||="$configs{database}/reference/hg19/hg19.fa";
    $configs{targetRegion} ||="$configs{database}/targetRegion/BGI_exome_V4_region.bed";
    $configs{Splitbed}     ||="$configs{database}/targetRegion/BGIv4_TR/Splitbed/";

    ## GATK_bundle
    $configs{bundle}         ||="$configs{database}/GATK_bundle";
    $configs{GATKdbsnp}      ||="$configs{bundle}/dbsnp_138.hg19.vcf.gz";
    $configs{GATKdbsnp150}   ||="$configs{bundle}/dbsnp_150.hg19.vcf.gz";
    $configs{GATKomni}       ||="$configs{bundle}/1000G_omni2.5.hg19.vcf.gz";
    $configs{GATKkgSNP}      ||="$configs{bundle}/1000G_phase1.snps.high_confidence.hg19.vcf.gz";
    $configs{GATKkgIndel}    ||="$configs{bundle}/1000G_phase1.indels.hg19.vcf.gz";
    $configs{GATKmillsIndel} ||="$configs{bundle}/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz";
    $configs{GATKhapmap}     ||="$configs{bundle}/hapmap_3.3.hg19.vcf.gz";

    ## mutect dbsnp and cosmic
    $configs{mutectDbsnp}    ||="$configs{bundle}/dbsnp_138.hg19.vcf.gz";
    #$configs{mutectCosmic}   ||="$configs{bundle}/cosmic_v79.hg19.vcf.gz";
    $configs{mutectCosmic} ||="$configs{bundle}/cosmic_v82.hg19.vcf.gz";

    ##annovar
    $configs{annovar}     ||="$configs{bin}/ANNOVAR/v20170716";
    $configs{annovardb}   ||="$configs{bin}/ANNOVAR/humandb20170901";
    $configs{annovarPara} ||="-remove -protocol refGene,ensGene,knownGene,cytoBand,genomicSuperDups,gwasCatalog,phastConsElements100way,phastConsElements46way,targetScanS,tfbsConsSites,wgRna,avsnp147,clinvar_20170130,cosmic82_coding,cosmic82_noncoding,dbnsfp33a,exac03nonpsych,exac03nontcga,gnomad_exome,gnomad_genome,popfreq_all_20150413,intervar_20170202 -operation g,g,g,r,r,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput";

    ## tools
    $configs{java8}    ||="$configs{bin}/java";
    $configs{java}     ||="/usr/java/latest/bin/java";
    $configs{SOAPnuke} ||="$configs{bin}/SOAPnuke";
    $configs{fqcheck}  ||="$configs{bin}/fqcheck33";
    $configs{bwa}      ||="$configs{bin}/bwa";
    $configs{samtools} ||="$configs{bin}/samtools";
    $configs{bcftools} ||="$configs{bin}/bcftools";
    $configs{picard}   ||="$configs{bin}/picard.jar";
    $configs{GATK}     ||="$configs{bin}/GenomeAnalysisTK.jar";
    $configs{GATK23}   ||="$configs{bin}/GenomeAnalysisTK-2.3.jar";
    $configs{mutect}   ||="$configs{bin}/mutect-1.1.7.jar";
    $configs{Platypus} ||="$configs{bin}/Platypus_0.8.1";
    $configs{varscan}  ||="$configs{bin}/VarScan.v2.4.3.jar";
    $configs{monitor}  ||="$configs{bin}/monitor";
    $configs{PanCanQC} ||="$configs{bin}/PanCanQC-1.2.2";

    ## python path
    $configs{pythonBin}="$configs{bin}/Anaconda2-4.3.0/bin";
    $configs{pythonLib}="$configs{bin}/Anaconda2-4.3.0/lib";
    $configs{htslib}="$configs{lib}/htslib";
    $configs{pythonPath}="$configs{bin}/Anaconda2-4.3.0/lib/python2.7/site-packages";

    ## clean ,bwa-mem and picard parameters
    $configs{platformPara} ||="COMPLETE";
    #$configs{SOAPnukePara} ||= "-n 0.1 -q 0.5 -l 5 -Q 2 -E 35 -G -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r CAACTCCTTGGCTCACAGAACGACATGGCTACGATCCGACTT -M 2";
    $configs{SOAPnukePara} ||= "-n 0.1 -q 0.5 -l 5 -Q 2 -E 35 -G -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG -M 2"; #modify at 20171017
    $configs{bwaMemPara} ||= "-t 6 -k 32";
    $configs{bwaAlnPara} ||= "-o 1 -e 50 -m 100000 -l 32 -k 2 -t 4 -L -i 15";
    $configs{bwaSamPara} ||= "";
    $configs{picardRmdupPara} ||= "REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true PROGRAM_RECORD_ID=null";
    $configs{covdepPara} ||= "--maxdepth 1000 -q 1";
    ## Varscan somatic parameters
    $configs{samtoolsMpileupPara}||="-q 20 -Q 20 -B";
    $configs{varscanPara}||="--output-vcf 1 --strand-filter 1 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.001 --min-freq-for-hom 0.75 --somatic-p-value 0.05 --normal-purity 1.0 --tumor-purity 1.0 --p-value 0.99 --mpileup 1";
    ## Platypus sindel parameters
    $configs{PlatypusPara} ||= "--filterDuplicates=1 --nCPU=1 --logFileName=/dev/null";
    $configs{PlatypussmtStrictPara} ||= "--minPosterior 5";
    ## gatk sindel parameters
    $configs{GATKsomaticIndelPara} ||= "--window_size 300";
    $configs{GATKsmtFlitPara} ||= "--RStartTh 5 --RendTh 5";

    $configs{seqType}        ||= "WES";   
    $configs{readLen}        ||= "PE50"; 
    $configs{rmlog}           ||= "true";
   
    $configs{FastqClean}     ||= "true"; 
    $configs{Aligment}        ||= "edico";
    $configs{rmCleanReads}   ||= "true";
    $configs{MarkDuplicates} ||= "true";
    $configs{RealnType}      ||= "realn";
    $configs{mergeBrecalBAM} ||= "false"; 
    $configs{bamdstQC}       ||= "true";
    $configs{TargetQC}       ||= "true";
    $configs{DNA_damage}     ||= "true";
    $configs{gatkDISCOVERY}  ||= "false";
    $configs{gatkGVCF}       ||= "false";
    $configs{mutectSNV}      ||= "true";
    $configs{mutect2}        ||= "true";
    $configs{varscanSNV}     ||= "true";
    $configs{gatkSindel}     ||= "true";
    $configs{platypusSindel} ||= "true";
    $configs{ADTExSCNV}      ||= "true";
    $configs{ABSOLUTE}       ||= "true";

    $configs{edgeList}       ||= "true";
    $configs{shlist}         ||= "true";
    $configs{rmlist}         ||= "true";

    # file exist check
    foreach my $key (sort keys %configs){
        next if($key=~/Para/);
        if($configs{$key}=~/\// && !-e $configs{$key}){
             print "Warning: Value of $key don't exist! $configs{$key}\n";
        }
    }

    # value check    
    if($configs{targetRegion} =~ /BGI_exome_V4/){
        `if [ -e $outdir/TR ];then rm -rf $outdir/TR;fi`;
        `ln -s $configs{database}/targetRegion/BGIv4_TR $outdir/TR`;
    }
    elsif($configs{USER_TR}){
	`if [ -e $outdir/TR ];then rm -rf $outdir/TR;fi`;
	`ln -s $configs{USER_TR} $outdir/TR`;
    }	
    else{
        die "targetRegion only support BGI_exome_V4 now\n";
    }
    if ($configs{somatic}){
	foreach my $chr(@chrs){
		my @files = glob "$configs{Splitbed}/$chr.list.*";
		@{$Splitbeds{$chr}}=@files;
		$SplitbedCount{$chr} = scalar @files;
	}
    }	
}

## parse targetRegion for gatk
sub parseTR{
    my $tr=shift;
    `mkdir -p $outdir/TR`;
    `mkdir $outdir/TR/interval`;
    `mkdir $outdir/TR/bed`;
    open TR,"$tr" or die $!;
    my (%outInt,%outBed);
    foreach my $i (@chrs){
        open $outInt{$i},">$outdir/TR/interval/$i.intervals" or die $!;
        open $outBed{$i},">$outdir/TR/bed/$i.bed" or die $!; 
    }
    while(<TR>){
        chomp;
        my ($c,$start,$end)=split /\s+/,$_;
        $outInt{$c}->print("$c:$start-$end\n");
        $outBed{$c}->print("$c\t$start\t$end\n");
    }
    close TR;
    foreach my $chr(@chrs){
        undef $outInt{$chr};
        undef $outBed{$chr};
    }
}

sub readFqlist{
    my ($fqlist) = @_;
    open LL,$fqlist or die $!;
    while(<LL>){
        chomp;
        my @F=split /\t/;
        die "Wrong format of $fqlist!" if(scalar(@F)!=7);
        my ($samp,$lib,$laneIndex,$fqfiles,$readLen,$insert,$baseNum)=@F[0,1,2,3,4,5,6];
        my $dir="$outdir/$samp/$laneIndex";
	my ($fq1,$fq2)=(split /,/,$fqfiles)[0,1];
	$rawfqs{$samp}{$lib}{$laneIndex}{1}=$fq1;
	$rawfqs{$samp}{$lib}{$laneIndex}{2}=$fq2;
        my $cleanfq1="$dir/raw/$laneIndex\_1.clean.fq.gz";
        my $cleanfq2="$dir/raw/$laneIndex\_2.clean.fq.gz";
        $cleanfqs{$samp}{$lib}{$laneIndex}{1}=$cleanfq1;
        $cleanfqs{$samp}{$lib}{$laneIndex}{2}=$cleanfq2;
	my ($chip,$lane,$index)=(split /_/,$laneIndex)[0,1,2];
	my $rgid="$chip\_$lane";
	push @{$cleanList{$samp}}, "$rgid,$samp,$lib,$lane,$configs{platformPara},BGI,$cleanfq1,$cleanfq2";

        $bwaBams{$samp}{$lib}{$laneIndex}="$dir/result/$laneIndex.sort.bam";
        if(!exists $rmdupBams{$samp}){
            $rmdupBams{$samp}  = "$outdir/$samp/rmdup/result/$samp.rmdup.bam";
            foreach my $chr (@chrs){
                $realnBams{$samp}{$chr}  = "$outdir/$samp/realn/result/$samp.$chr.realn.bam";
                $brecalBams{$samp}{$chr} = "$outdir/$samp/brecal/result/$samp.$chr.brecal.bam"
            }
            if($configs{Aligment} eq "edico"){
                $sampleBams{$samp} = "$outdir/$samp/edico/result/$samp.mkdup.bam";
            }else{
                $sampleBams{$samp} = "$outdir/$samp/brecal/result/$samp.bam";
            }
	    if($configs{mergeBrecalBAM} eq "true"){
		$mergeBams{$samp} = "$outdir/$samp/brecal/result/$samp.bam";
	    }
        }
    }
    close LL;
}

sub FastqClean{
    my ($hashref1,$hashref2) = @_;
    my %rawfqs   = %$hashref1;
    my %cleanfqs = %$hashref2;
    foreach my $samp (sort keys %rawfqs){
        foreach my $lib (sort keys %{$rawfqs{$samp}}){
	   foreach my $laneIndex (sort keys %{$rawfqs{$samp}{$lib}}){
               my $dir="$outdir/$samp/$laneIndex";
               makedir($dir) if(!-e $dir);;
               my $fq1=$rawfqs{$samp}{$lib}{$laneIndex}{1};
               my $fq2=$rawfqs{$samp}{$lib}{$laneIndex}{2};
               my $cleanfq1=$cleanfqs{$samp}{$lib}{$laneIndex}{1};
               my $cleanfq2=$cleanfqs{$samp}{$lib}{$laneIndex}{2};
               &peLaneClean($samp,$laneIndex,$dir,$fq1,$fq2,$cleanfq1,$cleanfq2);
           } 
       }
   }
}

sub Aligment_BWA{
    my ($hashref) = @_;
    my %cleanfqs = %$hashref;
    foreach my $samp (sort keys %cleanfqs){
        foreach my $lib (sort keys %{$cleanfqs{$samp}}){
	   foreach my $laneIndex (sort keys %{$cleanfqs{$samp}{$lib}}){
                my $dir="$outdir/$samp/$laneIndex";
		my $cleanfq1=$cleanfqs{$samp}{$lib}{$laneIndex}{1};
		my $cleanfq2=$cleanfqs{$samp}{$lib}{$laneIndex}{2};
                my ($chip,$lane,$index)=(split /_/,$laneIndex)[0,1,2];
                my $rgid="$chip\_$lane";
                if($configs{Aligment} eq "bwa_aln"){
                    &BWA_aln($samp,$cleanfq1,$cleanfq2,$laneIndex,$lib,$rgid,$dir);
                }
                elsif($configs{Aligment} eq "bwa_mem"){
                    &BWA_mem($samp,$cleanfq1,$cleanfq2,$laneIndex,$lib,$rgid,$dir);
                }
                if($configs{Aligment} =~ /bwa/){
	            &laneQC($samp,$laneIndex,$lib,$dir);
                }
            }
        }
    }
}

# from /hwfssz1/ST_CANCER/POL/SHARE/CancerPipeline/Edico_BGISEQ_v1.0/CancerPipeline_Edico_BGISEQ_v1.0.pl
sub Aligment_Edico{
    my ($hashref) = @_;
    my %cleanList = %$hashref;
 
    # default /staging/human/reference/hg19/hg19.fa
    $configs{vc_reference}  = $configs{reference}; 
    if ($configs{readLen} =~ /PE50/){
        $configs{hash_table} = "$configs{database}/reference/hg19.fa.k_19";
    }elsif($configs{readLen} =~ /PE100/){
        #$hash_table = "/staging/examples/reference/hg19/hg19.fa.k_21.f_16.m_149";
        $configs{hash_table} = "$configs{database}/reference/hg19.fa.k_21";
    }

    foreach my $samp (sort keys %cleanList){
        my $sampdir="$outdir/$samp/edico";
        makedir($sampdir);
	my $infqcsv = "$sampdir/result/01fastqList4edico.csv";
	open FQCSV, ">$infqcsv" || die $!;
	print FQCSV "RGID,RGSM,RGLB,Lane,RGPL,RGCN,Read1File,Read2File\n";
	print FQCSV join ("\n", @{$cleanList{$samp}}), "\n";
	close FQCSV;
	
	&RunEdico($samp, $infqcsv);
    }
}

## fastq Clean by SOAPnuke
sub peLaneClean{
    my ($samp,$laneIndex,$dir,$fq1,$fq2,$cleanfq1,$cleanfq2)=@_;
    my $cleandir=dirname $cleanfq1;
    my $cleanfq1name = basename $cleanfq1;
    my $cleanfq2name = basename $cleanfq2;

    my $shell1="$dir/shell/clean.sh";
    my $cleanCmd="if [ -e \"$cleanfq1\" ];then rm -rf $cleanfq1;fi && \\\n";
    $cleanCmd.="if [ -e \"$cleanfq2\" ];then rm -rf $cleanfq2;fi && \\\n";
    $cleanCmd.="$configs{SOAPnuke} filter -1 $fq1 -2 $fq2 $configs{SOAPnukePara} -o $cleandir -C $cleanfq1name -D $cleanfq2name && \\\n"; 
    $cleanCmd.="gzip -t $cleanfq1 $cleanfq2";
    generateShell($shell1, $cleanCmd);

    my $shell2="$dir/shell/cleanQC.sh";    
    my $cleanQCCmd.="perl $configs{bin}/soapnuke_stat.pl $cleandir/Basic_Statistics_of_Sequencing_Quality.txt $cleandir/Statistics_of_Filtered_Reads.txt > $cleandir/$laneIndex.clean.stat && \\\n";
    #$cleanQCCmd.="rm -rf $cleandir/*.txt && \\\n";
    $cleanQCCmd.="$configs{fqcheck} -r $cleanfq1 -c $cleandir/$laneIndex\_1.clean.fqcheck && \\\n";
    $cleanQCCmd.="$configs{fqcheck} -r $cleanfq2 -c $cleandir/$laneIndex\_2.clean.fqcheck && \\\n";
    $cleanQCCmd.="export GNUPLOT_PS_DIR=/share/app/gnuplot-4.6.7/share/gnuplot/4.6/PostScript && \\\n";
    $cleanQCCmd.="perl $configs{bin}/fqcheck_distribute.pl $cleandir/$laneIndex\_1.clean.fqcheck $cleandir/$laneIndex\_2.clean.fqcheck -o $cleandir/$laneIndex.clean.";
    generateShell($shell2, $cleanQCCmd);
}

## Alignment by bwa mem
sub BWA_mem{
    my ($samp,$cleanfq1,$cleanfq2,$laneIndex,$lib,$rgid,$dir)=@_;
    my $shell="$dir/shell/bwa.sh";
    my $memCmd="export LD_LIBRARY_PATH=/hwfssz1/ST_PRECISION/PUB/softwares/packages/lib:\$LD_LIBRARY_PATH && \\\n";
    $memCmd.= "$configs{bwa} mem $configs{bwaMemPara} -M -R \'\@RG\\tID:$rgid\\tSM:$samp\\tLB:$lib\\tPU:$laneIndex\\tPL:$configs{platformPara}\\tCN:BGI\' $configs{reference} $cleanfq1 $cleanfq2 | $configs{samtools} view -b -S -F 256 -t $configs{reference}.fai -o $dir/result/$laneIndex.bam - && \\\n";
    $memCmd.= "$configs{java8} -Xmx6g -XX:-UseGCOverheadLimit -jar $configs{picard} SortSam VALIDATION_STRINGENCY=SILENT TMP_DIR=$dir/result/tmp CREATE_INDEX=true SORT_ORDER=coordinate I=$dir/result/$laneIndex.bam O=$dir/result/$laneIndex.sort.bam";
    generateShell($shell, $memCmd);
}

sub BWA_aln{
    my ($samp,$cleanfq1,$cleanfq2,$laneIndex,$lib,$rgid,$dir)=@_;
    
    my $shell1="$dir/shell/bwa_aln1.sh";
    my $aln1Cmd="$configs{bwa} aln $configs{bwaAlnPara} -f $dir/result/$laneIndex\_1.sai $configs{reference} $cleanfq1";
    generateShell($shell1,$aln1Cmd);

    my $shell2="$dir/shell/bwa_aln2.sh";
    my $aln2Cmd="$configs{bwa} aln $configs{bwaAlnPara} -f $dir/result/$laneIndex\_2.sai $configs{reference} $cleanfq2";
    generateShell($shell2,$aln2Cmd);

    my $shell3="$dir/shell/bwa.sh";
    my $sampeCmd="$configs{bwa} sampe $configs{bwaSamPara} -r \"\@RG\\tID:$rgid\\tSM:$samp\\tLB:$lib\\tPU:$laneIndex\\tPL:$configs{platformPara}\\tCN:BGI\" $configs{reference} $dir/result/$laneIndex\_1.sai $dir/result/$laneIndex\_2.sai $cleanfq1 $cleanfq2 | $configs{samtools} view -b -S -t $configs{reference}.fai - >$dir/result/$laneIndex.bam && \\\n";
    $sampeCmd.="$configs{java8} -Xmx6g -XX:-UseGCOverheadLimit -jar $configs{picard} SortSam VALIDATION_STRINGENCY=SILENT TMP_DIR=$dir/result/tmp CREATE_INDEX=true SORT_ORDER=coordinate I=$dir/result/$laneIndex.bam O=$dir/result/$laneIndex.sort.bam";
    generateShell($shell3,$sampeCmd);
}

sub laneQC{
    my ($samp,$laneIndex,$lib,$dir)=@_;
    my $bwabam="$dir/result/$laneIndex.bam";
    my $sortbam="$dir/result/$laneIndex.sort.bam";

    my $shell="$dir/shell/laneQC.sh";
    my $cmd="export LD_LIBRARY_PATH=/hwfssz1/ST_PRECISION/PUB/softwares/packages/lib:\$LD_LIBRARY_PATH && \\\n";
    $cmd.="$configs{samtools} flagstat $bwabam > $bwabam.stat && \\\n";
    $cmd.="$configs{samtools} flagstat $sortbam > $sortbam.stat && \\\n";
    $cmd.="diff1=\$(diff $bwabam.stat $sortbam.stat | wc -l) && \\\n";
    $cmd.="[ \$diff1 -eq 0 ] && [ -e $dir/raw/$laneIndex\_1.clean.fq.gz ] && rm $dir/raw/$laneIndex\_1.clean.fq.gz && \\\n";
    $cmd.="[ \$diff1 -eq 0 ] && [ -e $dir/raw/$laneIndex\_2.clean.fq.gz ] && rm $dir/raw/$laneIndex\_2.clean.fq.gz && \\\n";
    $cmd.="[ \$diff1 -eq 0 ] && [ -e $dir/result/$laneIndex\_1.sai ] && rm $dir/result/$laneIndex\_1.sai && \\\n";
    $cmd.="[ \$diff1 -eq 0 ] && [ -e $dir/result/$laneIndex\_2.sai ] && rm $dir/result/$laneIndex\_2.sai && \\\n";
    $cmd.="[ \$diff1 -eq 0 ] && [ -e $bwabam ] && rm -r $bwabam $dir/result/tmp";
    generateShell($shell,$cmd);   
}

# from /hwfssz1/ST_CANCER/POL/SHARE/CancerPipeline/Edico_BGISEQ_v1.0/CancerPipeline_Edico_BGISEQ_v1.0.pl
sub RunEdico{
	my ($samp, $infqcsv) = @_;
	my $sampdir="$outdir/$samp/edico";
	my $tmpdir = "$sampdir/tmp";
	my $addpara = "--dbsnp=$configs{GATKdbsnp150}";
	if ($configs{seqType} =~ /WES/i){ # WES & Panel
		$addpara= " --vc-target-bed=$outdir/TR/CallVariantRegion/ex_region.sort.bed";
	}
	
	open ECSH, ">$sampdir/shell/edico.$samp.sh" || die $!;
	print ECSH "#!/bin/bash\necho ==========start at : `date` ==========\n";
	print ECSH "/opt/edico/bin/dragen_reset && \\\n";
	#print ECSH "if [ ! -d \"$tmpdir\" ];then mkdir $tmpdir;fi && \\\n";
	print ECSH "/opt/edico/bin/dragen --force --ref-dir=$configs{hash_table} --fastq-list=$infqcsv --fastq-list-sample-id=$samp --output-file-prefix=$samp.mkdup --output-directory=$sampdir/result --intermediate-results-dir=$tmpdir --output-format=bam --enable-map-align-output=true --enable-sort=true --enable-duplicate-marking=true --enable-bam-indexing=true --enable-variant-caller=true --vc-emit-ref-confidence=GVCF --vc-reference=$configs{reference} --vc-sample-name=$samp $addpara && \\\n";
	print ECSH "/opt/edico/bin/dragen --ref-dir=$configs{hash_table} --variant=$sampdir/result/$samp.mkdup.hard-filtered.gvcf.gz --output-file-prefix=$samp.mkdup.joint --output-directory=$sampdir/result --intermediate-results-dir=$tmpdir --enable-joint-genotyping=true $addpara && \\\n";
	#print ECSH "rm -rf $sampdir/result/$samp.mkdup.gvcf.gz* $sampdir/result/$samp.mkdup.joint.vcf.gz* && \\\n";
	print ECSH "if [ -d \"$tmpdir\" ];then rm -rf $tmpdir;fi && \\\n";
	print ECSH echostring("$sampdir/shell/edico.$samp.sh");
	close ECSH;


       my $shell1="$sampdir/shell/edico_process.sh";
       my $cmd1.= "/usr/bin/gunzip $sampdir/result/$samp.mkdup.gvcf.gz && \\\n";
       $cmd1.= << "HERE";
sed -i '2s#this location">#this location">\\n\\#\\#FILTER=<ID=LowQual,Description="Low quality">#' $sampdir/result/$samp.mkdup.gvcf && \\
HERE
        $cmd1.="$configs{bin}/bgzip $sampdir/result/$samp.mkdup.gvcf && \\\n";
        $cmd1.="$configs{bin}/tabix -p vcf $sampdir/result/$samp.mkdup.gvcf.gz && \\\n";
        $cmd1.="ln -s $sampdir/result/$samp.mkdup.gvcf.gz $outdir/$samp/gatkGVCF/raw/$samp.gvcf.gz && \\\n";
	  $cmd1.="ln -s $sampdir/result/$samp.mkdup.gvcf.gz.tbi $outdir/$samp/gatkGVCF/raw/$samp.gvcf.gz.tbi && \\\n";
	  $cmd1.="mkdir -p $outdir/$samp/rmdup/result && \\\n";
	  $cmd1.="ln -s $sampdir/result/$samp.mkdup.bam $outdir/$samp/rmdup/result/$samp.rmdup.bam && \\\n";
   	  $cmd1.="ln -s $sampdir/result/$samp.mkdup.bam.bai $outdir/$samp/rmdup/result/$samp.rmdup.bam.bai";
	  generateShell($shell1,$cmd1);
	
	open ANNSH, ">$sampdir/shell/variant_Ann.$samp.sh" || die $!;
	print ANNSH "#!/bin/bash\necho ==========start at : `date` ==========\n";
	print ANNSH "export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
	print ANNSH "zcat $sampdir/result/$samp.mkdup.joint.hard-filtered.vcf.gz | awk '\$1~/^#/ || \$7==\"PASS\" {print}' | $configs{bcftools} norm -m-both | $configs{bcftools} norm -f $configs{reference} -o $sampdir/result/$samp.mkdup.joint.hard-filtered.PASS.vcf.gz -O z && \\\n";
	print ANNSH "perl $configs{annovar}/table_annovar.pl $sampdir/result/$samp.mkdup.joint.hard-filtered.PASS.vcf.gz $configs{annovardb} -buildver hg19 -out $sampdir/result/$samp.mkdup.joint.hard-filtered.PASS -remove -protocol refGene,ensGene,knownGene,cytoBand,genomicSuperDups,gwasCatalog,phastConsElements100way,phastConsElements46way,targetScanS,tfbsConsSites,wgRna,avsnp147,clinvar_20170130,cosmic82_coding,cosmic82_noncoding,dbnsfp33a,exac03nonpsych,exac03nontcga,gnomad_exome,gnomad_genome,popfreq_all_20150413,intervar_20170202 -operation g,g,g,r,r,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput && \\\n";
	print ANNSH "gzip -f $sampdir/result/$samp.mkdup.joint.hard-filtered.PASS.hg19_multianno.txt && \\\n";
	print ANNSH "rm -rf $sampdir/result/$samp.mkdup.joint.hard-filtered.PASS.vcf.gz $sampdir/result/$samp.mkdup.joint.hard-filtered.PASS.hg19_multianno.vcf $sampdir/result/$samp.mkdup.joint.hard-filtered.PASS.avinput && \\\n";
	print ANNSH echostring("$sampdir/shell/variant_Ann.$samp.sh");
	close ANNSH;
	
	# read count per 1000bp
	open RCSH, ">$sampdir/shell/bam_readcount.$samp.sh" || die $!;
	print RCSH "#!/bin/bash\necho ==========start at : `date` ==========\n";
	if (`grep chr $configs{reference}\.fai`){
		print RCSH "$configs{bin}/readCounter -w 1000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY $sampleBams{$samp} | gzip > $sampdir/result/$samp.mkdup.win1k.wig.gz && \\\n";
	}
	else{
		print RCSH "$configs{bin}/readCounter -w 1000 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $sampleBams{$samp} | gzip > $sampdir/result/$samp.mkdup.win1k.wig.gz && \\\n";
	}
	print RCSH echostring("$sampdir/shell/bam_readcount.$samp.sh");
	close RCSH;
	
	# calculate the quality control values used in the ICGC PanCan project. Generally for WGS data.
	# modify form https://github.com/eilslabs/PanCanQC
	my $qcdir = "$outdir/$samp/QC";
	my $ACEseq = "$outdir/$samp/QC/ACEseq";
	my $qcshell = "$outdir/$samp/QC/shell";
	unless(-e $qcdir){system("mkdir -p $qcdir")}
	unless(-e $ACEseq){system("mkdir -p $ACEseq")}
	unless(-e $qcshell){system("mkdir -p $qcshell")}
	open STEP1, ">$qcshell/flagstat.$samp.sh" || die $!;
	print STEP1 "#!/bin/bash\necho ==========start at : `date` ==========\n";
	print STEP1 "export LD_LIBRARY_PATH=/hwfssz1/ST_CANCER/POL/SHARE/tools/samtools/xz/v5.2.3/lib:\$LD_LIBRARY_PATH && \\\n";
	print STEP1 "$configs{samtools} flagstat $sampleBams{$samp} > $qcdir/$samp.flagstat && \\\n";
	print STEP1 echostring("$qcshell/flagstat.$samp.sh");
	close STEP1;
	
	open STEP2, ">$qcshell/coverageQc.$samp.sh" || die $!;
	print STEP2 "#!/bin/bash\necho ==========start at : `date` ==========\n";
	print STEP2 "$configs{PanCanQC}/coverageQc --alignmentFile=$sampleBams{$samp} --outputFile=$qcdir/$samp.genome_coverage --processors=1 --basequalCutoff=0 --ungappedSizes=$configs{PanCanQC}/hg19.chrLenOnlyACGT_realChromosomes.tab && \\\n";
	print STEP2 echostring("$qcshell/coverageQc.$samp.sh");
	close STEP2;
	
	open STEP3, ">$qcshell/genomeCoverage.$samp.sh" || die $!;
	print STEP3 "#!/bin/bash\necho ==========start at : `date` ==========\n";
	print STEP3 "$configs{PanCanQC}/genomeCoverage --alignmentFile=$sampleBams{$samp} --outputFile=/dev/stdout --processors=4 --mode=countReads --windowSize=1 | perl $configs{PanCanQC}/filter_readbins.pl - $configs{PanCanQC}/hg19.chrLenOnlyACGT_realChromosomes.tab | gzip > $qcdir/$samp.readbin_coverage.gz && \\\n";
	print STEP3 "gzip -dc $qcdir/$samp.readbin_coverage.gz | awk '{print \$1,\$2,\$2+999,\$3}' | sed 's/ /\\t/g' |  sed '1i\\#chr\\tpos\\tend\\tcoverage' | perl $configs{PanCanQC}/annotate_vcf.pl -a - --aFileType=custom --aChromColumn chr --aPosColumn pos --aEndColumn end -b $configs{PanCanQC}/wgEncodeCrgMapabilityAlign100mer_chr.bedGraph.gz --tabix_bin /$configs{bin}/tabix --bFileType=bed --reportBFeatCoord --columnName map | $configs{pythonBin}/python $configs{PanCanQC}/addMappability.py -o $ACEseq/$samp.readbin_coverage.Mappability.gz && \\\n";
	print STEP3 "$configs{pythonBin}/python $configs{PanCanQC}/merge_and_filter_cnv.py --inputfile $ACEseq/$samp.readbin_coverage.Mappability.gz --output $ACEseq/$samp.readbin_coverage.Mappability.filtered.gz --coverage 0 --mappability 1000 --NoOfWindows 5 && \\\n";
	print STEP3 "/share/app/R-3.3.2/bin/Rscript $configs{PanCanQC}/correctGCBias.R --windowFile $ACEseq/$samp.readbin_coverage.Mappability.filtered.gz --timefile $configs{PanCanQC}/ReplicationTime_10cellines_mean_10KB.Rda --chrLengthFile $configs{PanCanQC}/chrlengths.txt --pid $samp --outfile $samp.corrected.txt --corPlot $samp.gc_corrected.png --corTab $samp.qc_gc_corrected.tsv --qcTab $samp.qc_gc_corrected.slim.txt --gcFile $configs{PanCanQC}/hg19_GRch37_100genomes_gc_content_10kb.txt --outDir $ACEseq --lowess_f 0.1 --scaleFactor 0.9 --coverageYlims 4 && \\\n";
	print STEP3 echostring("$qcshell/genomeCoverage.$samp.sh");
	close STEP3;
	
	open STEP4, ">$qcshell/bam_stats.$samp.sh" || die $!;
	print STEP4 "#!/bin/bash\necho ==========start at : `date` ==========\n";
	print STEP4 "/share/app/glibc-2.17/lib/ld-2.17.so --library-path /share/app/glibc-2.17/lib:/share/app/libz/zlib-1.2.11 $configs{PanCanQC}/bam_stats -i $sampleBams{$samp} -o $qcdir/$samp.read_edits && \\\n";
	print STEP4 echostring("$qcshell/bam_stats.$samp.sh");
	close STEP4;
	
	open STEP5, ">$qcshell/summary.$samp.sh" || die $!;
	print STEP5 "#!/bin/bash\necho ==========start at : `date` ==========\n";
	print STEP5 "perl $configs{PanCanQC}/summaryQC.pl $qcdir/$samp.flagstat $qcdir/$samp.genome_coverage $ACEseq/$samp.qc_gc_corrected.slim.txt $qcdir/$samp.read_edits $qcdir/$samp.readbin_coverage.gz $configs{PanCanQC}/hg19_reducedGenome.n300l5M.sorted.bed $qcdir/$samp.PanCanQC.summary.txt && \\\n";
	print STEP5 "tar -czf $qcdir/ACEseq.tar.gz $ACEseq && rm -rf $ACEseq   && \\\n";
	print STEP5 echostring("$qcshell/summary.$samp.sh");
	close STEP5;
	
	# remove clean reads
	if ($configs{rmCleanReads}=~/true/i){
		open RMSH, ">$sampdir/shell/rmCleanReads.$samp.sh" || die $!;
		print RMSH "#!/bin/bash\necho ==========start at : `date` ==========\n";
		print RMSH "perl $configs{bin}/rm_cleanReads_edico_pipeline.pl $sampleBams{$samp} $infqcsv && \\\n";
		print RMSH echostring("$sampdir/shell/rmCleanReads.$samp.sh");
		close RMSH;
	}
}


## mark duplication and merge 
sub MarkDuplicates{
    my ($hashref) = @_;
    my %bwaBams = %$hashref;
    foreach my $samp (sort keys %bwaBams){
        my $dir="$outdir/$samp/rmdup";
        makedir($dir);
        my $shell1="$dir/shell/rmdup.$samp.sh";
        my $cmd1="$configs{java8} -Djava.io.tmpdir=$dir/tmp -Xmx6g -XX:-UseGCOverheadLimit -jar $configs{picard} MarkDuplicates CREATE_INDEX=true";
        my $bwaBam='';
        foreach my $lib (sort keys %{$bwaBams{$samp}}){
            foreach my $laneIndex (sort keys %{$bwaBams{$samp}{$lib}}) {
                $cmd1.=" I=$bwaBams{$samp}{$lib}{$laneIndex}";
                $bwaBam.="$bwaBams{$samp}{$lib}{$laneIndex} ";
            }
        }
        $cmd1.=" O=$dir/result/$samp.rmdup.bam METRICS_FILE=$dir/result/$samp.rmdup.bam.met TMP_DIR=$dir/tmp $configs{picardRmdupPara}";
        generateShell($shell1,$cmd1);

        my $shell2="$dir/shell/rmdupQC.sh";
        my $cmd2.="$configs{samtools} flagstat $dir/result/$samp.rmdup.bam > $dir/result/$samp.rmdup.bam.stat && \\\n";
        $cmd2.="$configs{samtools} idxstats $dir/result/$samp.rmdup.bam > $dir/result/$samp.rmdup.bam.idxstats && \\\n";
        chop($bwaBam);
        $cmd2.="#rm -rf $bwaBam";
        generateShell($shell2,$cmd2);
    }
}

## Perform local realnment of reads around indels on patient level
sub PatientRealn{
	my ($corealnFile)=@_;
        makedir("$outdir/Co-realn");
        `mkdir -p $outdir/Co-realn/intervals`;
        open CORE,"$corealnFile" or die $!;
        while(<CORE>){
            chomp;
            my @RealnSamples=split /\s+/,$_;
            foreach my $chr (@chrs){
                corealnShell(\@RealnSamples,$chr);
            }
        }
        close CORE;
}

sub GATKrealnRecal{
    my ($ref)=@_;
    my %rmdupBams = %$ref;
    foreach my $samp(sort keys %rmdupBams){
        my $rmdupBam=$rmdupBams{$samp};
        my $realndir="$outdir/$samp/realn";
        makedir($realndir);
        my $brecaldir="$outdir/$samp/brecal";
        makedir($brecaldir);
        foreach my $chr (@chrs){    
            my $realnBam  = $realnBams{$samp}{$chr};
            my $brecalBam = $brecalBams{$samp}{$chr};
            my $shell="$realndir/shell/GATKrealnRecal.$samp.$chr.sh";
            my $cmd= "$configs{java8} -Djava.io.tmpdir=$realndir/tmp -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T RealignerTargetCreator -l INFO -I $rmdupBam -R $configs{reference} -L $outdir/TR/interval/$chr.intervals -o $realndir/raw/$samp.$chr.realn.intervals -known $configs{GATKkgIndel} -known $configs{GATKmillsIndel} -nt 4 && \\\n";
            $cmd.= "$configs{java8} -Djava.io.tmpdir=$realndir/tmp -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T IndelRealigner -l INFO -maxReads 60000 -maxInMemory 400000 -I $rmdupBam -R $configs{reference} ";
            $cmd.="-L $chr " if($chr ne "allchr");
            $cmd.="-o $realnBam -targetIntervals $realndir/raw/$samp.$chr.realn.intervals -known $configs{GATKkgIndel} -known $configs{GATKmillsIndel} && \\\n";
            $cmd.="$configs{java8} -Djava.io.tmpdir=$brecaldir/tmp -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T BaseRecalibrator -l INFO -I $realnBam -R $configs{reference} -L $outdir/TR/interval/$chr.intervals -knownSites $configs{GATKdbsnp} -knownSites $configs{GATKkgIndel} -knownSites $configs{GATKmillsIndel} -o $brecaldir/raw/$samp.$chr.brecal.grp -nct 4 && \\\n";
            $cmd.="$configs{java8} -Djava.io.tmpdir=$brecaldir/tmp -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T PrintReads -l INFO -I $realnBam -R $configs{reference} -BQSR $brecaldir/raw/$samp.$chr.brecal.grp -o $brecalBam -nct 4 && \\\n";
            $cmd.="rm $realnBam";
            generateShell($shell,$cmd);
        }
    }
}

sub corealnShell{
    my ($ref,$chr)=@_;
    my @RealnSamples=@$ref;
    my $patientId= shift @RealnSamples;
    my $realIn='';
    my @outbamList=();
    for my $samp(@RealnSamples){
        makedir("$outdir/$samp/realn") if(!-e "$outdir/$samp/realn");
        my $rmdupBam=$rmdupBams{$samp};
        my $realnBam=$realnBams{$samp}{$chr};
        my $rmdupBamName=(split /\//,$rmdupBam)[-1];
        push @outbamList,"$rmdupBamName\t$realnBam";
        $realIn.="-I $rmdupBam ";
        $patientIds{$samp}=$patientId;
    }
    chop($realIn);
    
    my $sampleNum=scalar(@RealnSamples);
    if(!exists $maxmems{$patientId}){
        $maxmems{$patientId} = 4 * ($sampleNum/4 + 1) > 10 ? 14 : 10 * ($sampleNum/10 + 1);
        $maxmems{$patientId} .= "g";
    }
    my $maxR = 20000 * $sampleNum > 150000 ? 150000 : 20000 * $sampleNum;
    my $maxRiM = 150000 * $sampleNum > 1000000 ? 1000000 : 150000 * $sampleNum;

    my $outbam="$outdir/Co-realn/result/$patientId.$chr.co-realn.bam.map"; 
    open COL,">$outbam";print COL join("\n",@outbamList);close COL;
    
    my $shell="$outdir/Co-realn/shell/realn.$patientId.$chr.sh";
    my $realnCmd="$configs{java8} -Djava.io.tmpdir=$outdir/Co-realn/tmp -Xmx$maxmems{$patientId} -XX:-UseGCOverheadLimit -jar $configs{GATK} -T RealignerTargetCreator -l INFO $realIn -R $configs{reference} -L $outdir/TR/interval/$chr.intervals -o $outdir/Co-realn/intervals/$patientId.$chr.realn.intervals -known $configs{GATKkgIndel} -known $configs{GATKmillsIndel} -nt 4 && \\\n";
    $realnCmd.="$configs{java8} -Djava.io.tmpdir=$outdir/Co-realn/tmp -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T IndelRealigner -l INFO -maxReads $maxR -maxInMemory $maxRiM $realIn -R $configs{reference} ";
    $realnCmd.="-L $chr " if($chr ne "allchr");
    $realnCmd.="-nWayOut $outbam -targetIntervals $outdir/Co-realn/intervals/$patientId.$chr.realn.intervals -known $configs{GATKkgIndel} -known $configs{GATKmillsIndel}";
    foreach(@outbamList){
        my $realnBam=(split /\t/,$_)[-1];
        $realnCmd.=" && \\\n$configs{samtools} index $realnBam";
    }
    generateShell($shell,$realnCmd);
}

sub BaseRecalibrator{
    my ($hashref) = @_;
    my %realnBams = %$hashref;
    foreach my $samp(sort keys %realnBams){
        my $dir="$outdir/$samp/brecal";
        makedir($dir);
        foreach my $chr( sort keys %{$realnBams{$samp}}){
            my $shell="$dir/shell/brecal.$samp.$chr.sh";
            my $brecalCmd="$configs{java8} -Djava.io.tmpdir=$dir/tmp -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T BaseRecalibrator -l INFO -I $realnBams{$samp}{$chr} -R $configs{reference} -L $outdir/TR/interval/$chr.intervals -knownSites $configs{GATKdbsnp} -knownSites $configs{GATKkgIndel} -knownSites $configs{GATKmillsIndel} -o $dir/raw/$samp.$chr.brecal.grp -nct 4 && \\\n";
            $brecalCmd.="$configs{java8} -Djava.io.tmpdir=$dir/tmp -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T PrintReads -l INFO -I $realnBams{$samp}{$chr} -R $configs{reference} -BQSR $dir/raw/$samp.$chr.brecal.grp -o $dir/result/$samp.$chr.brecal.bam -nct 4";
            generateShell($shell,$brecalCmd);
        }
    }
}

sub mergeBrecalBAM{
    my ($hashref,$hashref_merge) = @_;
    my %brecalBams = %$hashref;
    my %mergeBams = %$hashref_merge;
    foreach my $samp(sort keys %brecalBams){
	my $sampleBAM=$sampleBams{$samp};
        my $shell="$outdir/$samp/brecal/shell/mergeBam.$samp.sh";
        my $mergeCmd = "export LD_LIBRARY_PATH=/hwfssz1/ST_PRECISION/PUB/softwares/packages/lib:\$LD_LIBRARY_PATH && \\\n";
	$mergeCmd.="$configs{java8} -Djava.io.tmpdir=$outdir/$samp/brecal/tmp -Xmx2g -XX:-UseGCOverheadLimit -jar $configs{picard} MergeSamFiles CREATE_INDEX=true";
        foreach my $chr( sort keys %{$brecalBams{$samp}}){
            $mergeCmd.=" I=$brecalBams{$samp}{$chr}";   
        }
        $mergeCmd.=" O=$mergeBams{$samp} SO=coordinate AS=true VALIDATION_STRINGENCY=SILENT && \\\n";
        $mergeCmd.="$configs{samtools} flagstat $mergeBams{$samp} > $outdir/$samp/brecal/result/$samp.bam.stat";
        generateShell($shell,$mergeCmd);
    }
}

## calculate the coverage and depth information for each sampleM
sub bamdstQC{
    my ($hashref) = @_;
    my %sampleBams = %$hashref;
    foreach my $samp(sort keys %sampleBams){ 
        my $dir="$outdir/$samp/coverage";
        makedir($dir);
        my $shell="$dir/shell/covdep.$samp.sh";
        my $cmd="export LD_LIBRARY_PATH=$configs{bin}/bamtools-2.4.1/lib:\$LD_LIBRARY_PATH && \\\n";
        $cmd.="$configs{bin}/bamdst $configs{covdepPara} -p $outdir/TR/bed/allchr.bed -o $dir/result $sampleBams{$samp} && \\\n";
        $cmd.="perl $configs{bin}/bamdstPlot.pl -i $dir/result/depth_distribution.plot -c $dir/result/coverage.report -o $dir/result && \\\n";
  	$cmd.="$configs{samtools} view -bu -F 4 $sampleBams{$samp} -L $outdir/TR/bed/allchr.bed | $configs{bin}/bamtools stats -in /dev/stdin -insert > $dir/result/Target.report";
        generateShell($shell,$cmd);
    }
}

## Analyze coverage distribution and validate read mates per interval and per sample
## Collect quality metrics for a set of intervals
sub TargetQC{
    my ($hashref) = @_;
    my %sampleBams = %$hashref;
    foreach my $samp(sort keys %sampleBams){
        my $dir="$outdir/$samp/coverage";
        my $bam=$sampleBams{$samp};
        my $shell="$dir/shell/targetqc.$samp.sh";
        my $targCmd="$configs{java8} -Djava.io.tmpdir=$dir/tmp -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T DiagnoseTargets -l INFO -I $bam -R $configs{reference} -L $outdir/TR/interval/allchr.intervals -o $dir/result/$samp.TqrgetQC.vcf.gz -missing $dir/result/$samp.TqrgetQC.missing.intervals && \\\n";
        $targCmd.="$configs{java8} -Djava.io.tmpdir=$dir/tmp -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T QualifyMissingIntervals -l INFO -I $bam -R $configs{reference} -targets $outdir/TR/interval/allchr.intervals -L $dir/result/$samp.TqrgetQC.missing.intervals -o $dir/result/$samp.TqrgetQC.missing.grp";
        generateShell($shell,$targCmd);
    }
}

sub DNA_damage{ 
    my ($hashref) = @_;
    my %sampleBams = %$hashref;
    foreach my $samp(sort keys %sampleBams){
		my $dir = "$outdir/$samp/DNA_damage";
		makedir($dir);
		my $bam=$sampleBams{$samp};
		my $shell="$dir/shell/DNA_damage.$samp.sh";
		my $Cmd = "perl $configs{bin}/Damage-estimator.pl $bam $configs{reference} $dir/result $samp";
		generateShell($shell,$Cmd);
    }
}

## calling SNP and INDEL by GATK34
sub gatkDISCOVERY{
    my ($samp,$bam,$chr)=@_;
    my $dir="$outdir/$samp/gatkDISCOVERY"; 
    my $shell1 = "$dir/shell/variant.$samp.$chr.sh";
    my $gatkCmd="$configs{java8} -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T HaplotypeCaller -l INFO -R $configs{reference} -I $bam -L $outdir/TR/CallVariantRegion/ex_region.sort.$chr.bed --genotyping_mode DISCOVERY -stand_call_conf 30 -o $dir/raw/$samp.$chr.HC.variant.vcf.gz && \\\n";
    $gatkCmd.="$configs{java8} -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T UnifiedGenotyper -l INFO -R $configs{reference} -I $bam -L $outdir/TR/CallVariantRegion/ex_region.sort.$chr.bed --genotyping_mode DISCOVERY -stand_call_conf 30 -o $dir/raw/$samp.$chr.UG.variant.vcf.gz";
    generateShell($shell1,$gatkCmd);
 
    #hard-filtering for SNP (single WES unable to use VQSR)
    my $shell2="$dir/shell/variant.$samp.$chr.process.sh";
    my $gatkFiltcmd="$configs{java8} -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T SelectVariants -l INFO -V $dir/raw/$samp.$chr.UG.variant.vcf.gz -R $configs{reference} -selectType SNP -o $dir/raw/$samp.$chr.snp_raw.vcf.gz && \\\n";
    $gatkFiltcmd.="$configs{java8} -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T VariantFiltration -l INFO -V $dir/raw/$samp.$chr.snp_raw.vcf.gz -R $configs{reference} -filter \"QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 4.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" -filterName \"SNP_filter\" -o $dir/raw/$samp.$chr.snp_filter.vcf.gz && \\\n";
    #Vhard-filtering for INDEL
    $gatkFiltcmd.="$configs{java8} -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T SelectVariants -l INFO -V $dir/raw/$samp.$chr.HC.variant.vcf.gz -R $configs{reference} -selectType INDEL -o $dir/raw/$samp.$chr.indel_raw.vcf.gz && \\\n";
    $gatkFiltcmd.="$configs{java8} -Xmx4g -XX:-UseGCOverheadLimit -jar $configs{GATK} -T VariantFiltration -l INFO -V $dir/raw/$samp.$chr.indel_raw.vcf.gz -R $configs{reference} -filter \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0\" -filterName \"INDEL_filter\" -o $dir/raw/$samp.$chr.indel_filter.vcf.gz && \\\n";
    #annovar
    $gatkFiltcmd.="zcat $dir/raw/$samp.$chr.snp_filter.vcf.gz | awk '\$1~/^#/ || \$7==\"PASS\" {print}' | $configs{bcftools} norm -m-both -f $configs{reference} - | $configs{bcftools} norm -f $configs{reference} -o $dir/result/$samp.$chr.snp.vcf.gz -O z - && \\\n";
    $gatkFiltcmd.="zcat $dir/raw/$samp.$chr.indel_filter.vcf.gz | awk '\$1~/^#/ || \$7==\"PASS\" {print}' | $configs{bcftools} norm -m-both -f $configs{reference} - | $configs{bcftools} norm -f $configs{reference} -o $dir/result/$samp.$chr.indel.vcf.gz -O z - && \\\n";
#    $gatkFiltcmd.="perl $configs{annovar}/table_annovar.pl $dir/result/$samp.$chr.snp.vcf.gz $configs{annovardb} -buildver hg19 -out $dir/result/$samp.$chr.snp $configs{annovarPara} && \\\n";
#    $gatkFiltcmd.="perl $configs{annovar}/table_annovar.pl $dir/result/$samp.$chr.indel.vcf.gz $configs{annovardb} -buildver hg19 -out $dir/result/$samp.$chr.indel $configs{annovarPara} && \\\n";
#    $gatkFiltcmd.="gzip -f $dir/result/$samp.$chr.snp.hg19_multianno.txt $dir/result/$samp.$chr.indel.hg19_multianno.txt && \\\n";
#    $gatkFiltcmd.="rm -rf $dir/raw/$samp.$chr.snp_raw.vcf.gz* $dir/raw/$samp.$chr.indel_raw.vcf.gz* && \\\n";
#    $gatkFiltcmd.="rm -rf $dir/result/$samp.$chr.snp.avinput $dir/result/$samp.$chr.snp.hg19_multianno.vcf $dir/result/$samp.$chr.indel.avinput $dir/result/$samp.$chr.indel.hg19_multianno.vcf";
    generateShell($shell2,$gatkFiltcmd);
}

sub gatkGVCF_sample{
    my ($samp,$bam,$chr)=@_;
    my $dir = "$outdir/$samp/gatkGVCF";
    my $shell = "$dir/shell/gatkGVCF.$samp.$chr.sh";
    my $cmd;
    if ($configs{seqType} =~ /WES/i){
    $cmd="$configs{java8} -Xmx5G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T HaplotypeCaller -R $configs{reference} -I $bam -L $outdir/TR/CallVariantRegion/ex_region.sort.$chr.bed --emitRefConfidence GVCF -o $dir/raw/$samp.$chr.g.vcf.gz && \\\n";
    }elsif ($configs{seqType} =~ /WES/i){
    $cmd="$configs{java8} -Xmx5G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T HaplotypeCaller -R $configs{reference} -I $bam -L $chr --emitRefConfidence GVCF -o $dir/raw/$samp.$chr.g.vcf.gz && \\\n";
    }
    $cmd.="$configs{java8} -Xmx2G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T GenotypeGVCFs -R $configs{reference} --variant $dir/raw/$samp.$chr.g.vcf.gz -o $dir/result/$samp.$chr.vcf.gz -stand_call_conf 30 -allSites";
    generateShell($shell,$cmd);
}

sub gatkGVCF_GenotypeGVCFs{
    my ($samp)=@_;
    my $dir = "$outdir/$samp/gatkGVCF";
    my $shell = "$dir/shell/gatkGVCF.$samp.sh";
    my $cmd="$configs{java8} -Xmx2G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T GenotypeGVCFs -R $configs{reference} --variant $dir/raw/$samp.gvcf.gz -o $dir/result/$samp.vcf.gz -stand_call_conf 30 -allSites";
    generateShell($shell,$cmd);
}

sub gatkGVCF_merge{
    my ($samp) = @_;
    my $dir = "$outdir/$samp/gatkGVCF";
    my $shell1 = "$dir/shell/gatkGVCF.$samp.merge.sh";
    my $cmd1   = "$configs{java8} -cp $configs{GATK} org.broadinstitute.gatk.tools.CatVariants -R $configs{reference}";
    for my $chr (@chrs){
        $cmd1 .= " -V $dir/result/$samp.$chr.vcf.gz";
    }
    $cmd1 .= " -out $dir/result/$samp.vcf.gz --assumeSorted";
    generateShell($shell1,$cmd1);
}

sub gatkGVCF_process{
    my ($samp) = @_;
    my $dir = "$outdir/$samp/gatkGVCF";
    my $shell2 = "$dir/shell/gatkGVCF.$samp.SNP.sh";
    my $cmd2   = "export PERL5LIB=\"$configs{lib}:\$PERL5LIB\"\n";
    $cmd2 .= "export LIBRARY_PATH=\"/hwfssz1/ST_PRECISION/PUB/softwares/packages/lib:\$LIBRARY_PATH\"\n";
    $cmd2 .= "$configs{java8} -Xmx3G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T SelectVariants -R $configs{reference} -V $dir/result/$samp.vcf.gz -selectType SNP --excludeNonVariants -o $dir/raw/$samp.raw.snp.vcf.gz && \\\n";
    $cmd2 .= "$configs{java8} -Xmx3G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T VariantFiltration -R $configs{reference} \\\n";
    $cmd2 .= "-V $dir/raw/$samp.raw.snp.vcf.gz \\\n";
    $cmd2 .= "--filterExpression \"QD < 2.0 || FS > 60.0 || MQ <40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filterName \"filter\"  -o $dir/raw/$samp.filtered_snp.vcf.gz && \\\n";
    $cmd2 .= "$configs{bcftools} view -e 'ALT==\"*\"' -f \"PASS\" -o $dir/result/$samp.filtered_snp.vcf.gz -O z $dir/raw/$samp.filtered_snp.vcf.gz && \\\n";
    $cmd2 .= "$configs{bin}/annodb/annodb.pl --mode snp --if vcf --optdb \'Conservative,Cancer,Disease,Functional,ENCODE,Population\' --genedb \'proteinatlas,omim,pharmgkb,dbnsfp,cgc\' --remove --verdbsnp v149_hg19 --buildver hg19 --outfile $dir/result/$samp.filtered_snp.vcf $dir/result/$samp.filtered_snp.vcf.gz && \\\n";
    $cmd2 .= "perl $configs{bin}/annodb/Annodb_stat_for_all.pl -i $dir/result/$samp.filtered_snp.vcf.genome_summary.csv -o $dir/result/$samp.snp.stat -s $samp -v snp && \\\n";
    $cmd2 .= "gzip $dir/result/$samp.filtered_snp.vcf.genome_summary.csv $dir/result/$samp.filtered_snp.vcf.exome_summary.csv $dir/result/$samp.filtered_snp.vcf.gene_summary.csv";
    generateShell($shell2,$cmd2);

    my $shell3 = "$dir/shell/gatkGVCF.$samp.INDEL.sh";
    my $cmd3   = "export PERL5LIB=\"$configs{lib}:\$PERL5LIB\"\n";
    $cmd3 .= "export LIBRARY_PATH=\"/hwfssz1/ST_PRECISION/PUB/softwares/packages/lib:\$LIBRARY_PATH\"\n";
    $cmd3 .= "$configs{java8} -Xmx3G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T SelectVariants -R $configs{reference} -V $dir/result/$samp.vcf.gz -selectType INDEL --excludeNonVariants -o $dir/raw/$samp.raw.indel.vcf.gz && \\\n";
    $cmd3 .= "$configs{java8} -Xmx3G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T VariantFiltration -R $configs{reference} \\\n";
    $cmd3 .= "-V $dir/raw/$samp.raw.indel.vcf.gz \\\n";
    $cmd3 .= "--filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filterName \"filter\" -o $dir/raw/$samp.filtered_indel.vcf.gz && \\\n";
    $cmd3 .= "$configs{bcftools} view -f \"PASS\" -o $dir/result/$samp.filtered_indel.vcf.gz -O z $dir/raw/$samp.filtered_indel.vcf.gz && \\\n";
    $cmd3 .= "$configs{bin}/annodb/annodb.pl --mode indel --if vcf --optdb \'Conservative,Cancer,Disease,Functional,ENCODE,Population\' --genedb \'proteinatlas,omim,pharmgkb,dbnsfp,cgc\' --remove --verdbsnp v149_hg19 --buildver hg19 --outfile $dir/result/$samp.filtered_indel.vcf $dir/result/$samp.filtered_indel.vcf.gz && \\\n";
    $cmd3 .= "perl $configs{bin}/annodb/Annodb_stat_for_all.pl -i $dir/result/$samp.filtered_indel.vcf.genome_summary.csv -o $dir/result/$samp.indel.stat -s $samp -v indel && \\\n";
    $cmd3 .= "perl $configs{bin}/indel_stat.pl $dir/result/$samp.filtered_indel.vcf.genome_summary.csv $dir/result/$samp.indel_len && \\\n";
    $cmd3 .= "perl $configs{bin}/indel_stat.pl $dir/result/$samp.filtered_indel.vcf.exome_summary.csv $dir/result/$samp.indel_cds_len && \\\n";
    $cmd3 .= "perl $configs{bin}/indel_lenght_R.pl $samp $dir/result/$samp.indel_len.xls $dir/result/$samp.indel_cds_len.xls $dir/result/ && \\\n";
    $cmd3 .= "gzip $dir/result/$samp.filtered_indel.vcf.genome_summary.csv $dir/result/$samp.filtered_indel.vcf.exome_summary.csv $dir/result/$samp.filtered_indel.vcf.gene_summary.csv";
    generateShell($shell3,$cmd3);
}

sub gatkGVCF_VQSR{
    my ($sample) = @_;
    my $dir = "$outdir/$sample/gatkGVCF";
    my $shell1 = "$dir/shell/gatkGVCF.SNP.VQSR.sh";
    my $cmd1   = "export PERL5LIB=\"$configs{lib}:\$PERL5LIB\"\n";
    $cmd1 .= "$configs{java8} -Xmx3G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T SelectVariants -R $configs{reference} -V $dir/result/$sample.vcf.gz -selectType SNP --excludeNonVariants -o $dir/raw/$sample.raw.snp.vcf.gz && \\\n";
    $cmd1 .= "$configs{java8} -Xmx5G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T VariantRecalibrator -R $configs{reference} -input $dir/raw/$sample.raw.snp.vcf.gz \\\n";
    $cmd1 .= "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $configs{GATKhapmap} \\\n";
    $cmd1 .= "-resource:omni,known=false,training=true,truth=true,prior=12.0 $configs{GATKomni} \\\n";
    $cmd1 .= "-resource:1000G,known=false,training=true,truth=false,prior=10.0 $configs{GATKkgSNP} \\\n";
    $cmd1 .= "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $configs{GATKdbsnp} \\\n";
    # SNP specific recommendations for Exome data (not use DP)
    $cmd1 .= "-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \\\n";
    $cmd1 .= "-mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\\n"; 
    $cmd1 .= "-recalFile $dir/raw/$sample.recalibrate_SNP.recal -tranchesFile $dir/raw/$sample.recalibrate_SNP.tranches -rscriptFile $dir/raw/VQSR_SNP_plots.R && \\\n";
    $cmd1 .= "$configs{java8} -Xmx5G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T ApplyRecalibration -R $configs{reference} -input $dir/raw/$sample.raw.snp.vcf.gz -mode SNP --ts_filter_level 99.0 -recalFile $dir/raw/$sample.recalibrate_SNP.recal -tranchesFile $dir/raw/$sample.recalibrate_SNP.tranches -o $dir/raw/$sample.filtered_snp.vcf.gz && \\\n";
    $cmd1 .= "$configs{java8} -Xmx10G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T SelectVariants -R $configs{reference} -V $dir/raw/$sample.filtered_snp.vcf.gz --excludeFiltered -o $dir/result/$sample.filtered_snp.vcf.gz";
    generateShell($shell1,$cmd1);
    
    my $shell2 = "$dir/shell/gatkGVCF.Indel.VQSR.sh";
    my $cmd2 = "export PERL5LIB=\"$configs{lib}:\$PERL5LIB\"\n";
    $cmd2 .= "$configs{java8} -Xmx5G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T SelectVariants -R $configs{reference} -V $dir/result/$sample.vcf.gz -selectType INDEL --excludeNonVariants -o $dir/raw/$sample.raw.indel.vcf.gz && \\\n";
    $cmd2 .= "$configs{java8} -Xmx5G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T VariantRecalibrator -R $configs{reference} -input $dir/raw/$sample.raw.indel.vcf.gz \\\n";
    $cmd2 .= "-resource:mills,known=true,training=true,truth=true,prior=12.0 $configs{GATKmillsIndel} \\\n";
    $cmd2 .= "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $configs{GATKdbsnp} \\\n";
    $cmd2 .= "-an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff -mode INDEL \\\n";
    $cmd2 .= "-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 --maxGaussians 4 \\\n";
    $cmd2 .= "-recalFile $dir/raw/$sample.recalibrate_INDEL.recal -tranchesFile $dir/raw/$sample.recalibrate_INDEL.tranches -rscriptFile $dir/raw/$sample.recalibrate_INDEL_plots.R && \\\n";
    $cmd2 .= "$configs{java8} -Xmx5G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T ApplyRecalibration -R $configs{reference} -input $dir/raw/$sample.raw.indel.vcf.gz -mode INDEL --ts_filter_level 99.0 -recalFile $dir/raw/$sample.recalibrate_INDEL.recal -tranchesFile $dir/raw/$sample.recalibrate_INDEL.tranches -o $dir/raw/$sample.filtered_indel.vcf.gz && \\\n";
    $cmd2 .= "$configs{java8} -Xmx5G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T SelectVariants -R $configs{reference} -V $dir/raw/$sample.filtered_indel.vcf.gz --excludeFiltered -o $dir/result/$sample.filtered_indel.vcf.gz";
    generateShell($shell2,$cmd2);
}

sub gatkGVCF_combine_bychr{
    my ($dir) = @_;
    for my $chr (@chrs){    
        my $cmd="$configs{java8} -Xmx2G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T GenotypeGVCFs -R $configs{reference}";
        foreach my $samp(sort keys %brecalBams){
            $cmd.=" --variant $outdir/$samp/gatkGVCF/raw/$samp.$chr.g.vcf.gz";
        }
        $cmd.=" -o $dir/result/combine.$chr.vcf.gz";
        my $shell = "$dir/shell/gatkGVCF.combine.$chr.sh";
        generateShell($shell,$cmd);
    }
}

sub gatkGVCF_combine_bysample{
    my ($dir) = @_;
    my $cmd="$configs{java8} -Xmx2G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T GenotypeGVCFs -R $configs{reference} \\\n";
    foreach my $samp(sort keys %brecalBams){
        $cmd.="--variant $outdir/$samp/gatkGVCF/result/$samp.gvcf.gz \\\n";
    }
    $cmd.=" -o $dir/result/combine.vcf.gz";
    my $shell = "$dir/shell/gatkGVCF.combine.sh";
    generateShell($shell,$cmd);
}

sub gatkGVCF_combineGTs{
	my ($ref1,$ref2,$outdir) = @_;
	my %patients2samples = %$ref1;
	my %configs          = %$ref2;

	my $dir = "$outdir/combine/gatkGVCF";
	foreach my $p (sort keys %patients2samples){
		my @samples = @{$patients2samples{$p}};
		my $SNcmd;
		for my $s (@samples){
			$SNcmd.="-sn $s ";
		}
		my $shell = "$dir/shell/gatkGVCF.$p.combineGTs.sh";
		my $cmd = "$configs{java8} -Xmx3G -Djava.io.tmpdir=$dir/tmp -jar $configs{GATK} -T SelectVariants -R $configs{reference} -V $dir/result/combine.filtered.vep.vcf.gz -o $dir/raw/$p.VQSR.vep.vcf $SNcmd --restrictAllelesTo BIALLELIC && \\\n";
		$cmd .= "$configs{bin}/02SNP/Germline_combineGTs.pl $dir/raw/$p.VQSR.vep.vcf $dir/raw/$p.VQSR.vep.combineGTs.vcf && \\\n";
		$cmd .= "ln -s $dir/raw/$p.VQSR.vep.combineGTs.vcf $dir/raw/$p.VQSR.vep.combineGTs.vep.vcf && \\\n";
		$cmd .= "sh $configs{bin}/vcf2maf.sh $dir/raw/$p.VQSR.vep.combineGTs.vcf $p-A $p-T";
		generateShell($shell,$cmd);
	}
}


sub somatic_batch_byChr{
	my ($hashref1,$hashref2,$func,$func_pro,$type)=@_;
	my %brecalBams = %$hashref1;
	my %samplePairs = %$hashref2;
	foreach my $samplePair (sort keys %samplePairs){
		my $normal=$samplePairs{$samplePair}{normal};
		my $tumor=$samplePairs{$samplePair}{tumor};
		foreach my $chr(sort keys %{$brecalBams{$normal}}){
			my $nbam=$brecalBams{$normal}{$chr};
			my $tbam=$brecalBams{$tumor}{$chr};
			$func->($normal,$tumor,$nbam,$tbam,$chr);
		}
		&mergeSoVars($normal,$tumor,$type);
		$func_pro->($normal,$tumor);
	}
}

sub somatic_batch_byRegion{
	my ($hashref1,$hashref2,$func,$func_pro,$type)=@_;
	my %brecalBams = %$hashref1;
	my %samplePairs = %$hashref2;
	foreach my $samplePair (sort keys %samplePairs){
		my $normal=$samplePairs{$samplePair}{normal};
		my $tumor=$samplePairs{$samplePair}{tumor};
		foreach my $chr(sort keys %{$brecalBams{$normal}}){
			my $bedNum = $SplitbedCount{$chr};
			for my $bedID (1..$bedNum){
				my $nbam=$brecalBams{$normal}{$chr};
				my $tbam=$brecalBams{$tumor}{$chr};
				$func->($normal,$tumor,$nbam,$tbam,$chr,$bedID);
			}
		}
		&mergeSoVars($normal,$tumor,$type);
		$func_pro->($normal,$tumor);
	}
}

sub somatic_batch_bySamplePair{
        my ($hashref1,$hashref2,$func)=@_;
        my %brecalBams = %$hashref1;
        my %samplePairs = %$hashref2;
        foreach my $samplePair (sort keys %samplePairs){
                my $normal=$samplePairs{$samplePair}{normal};
                my $tumor=$samplePairs{$samplePair}{tumor};
		$func->($normal,$tumor);
        }
}


## calling somatic SNV by mutect
sub mutectSNV{
    my ($normal,$tumor,$nbam,$tbam,$chr)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $dir="$outdir/somatic/$samplePair/mutectSNV";
    makedir($dir) if(!-e $dir);
    
    my $shell1="$dir/shell/mutect.$samplePair.$chr.sh";
    my $cmd="$configs{java} -Xmx2g -Djava.io.tmpdir=$dir/tmp -XX:-UseGCOverheadLimit -jar $configs{mutect} -T MuTect -R $configs{reference} --input_file:normal $nbam --input_file:tumor $tbam --normal_sample_name $normal --tumor_sample_name $tumor --vcf $dir/raw/$samplePair.$chr.mutect.vcf.gz --out $dir/raw/$samplePair.$chr.mutect.txt";
    $cmd.=" --dbsnp $configs{mutectDbsnp}" if(-f $configs{mutectDbsnp});
    $cmd.=" --cosmic $configs{mutectCosmic}" if(-f $configs{mutectCosmic});
    if ($configs{seqType} =~ /WES/i){
        $cmd.=" -L $outdir/TR/CallVariantRegion/ex_region.sort.$chr.bed";
    }elsif($configs{seqType} =~ /WGS/i){
        $cmd.=" -L $chr";
    }
    $cmd.=" --enable_extended_output && \\\n";
    $cmd.="gzip -f $dir/raw/$samplePair.$chr.mutect.txt";
    generateShell($shell1,$cmd);
}

sub mutectSNV_process{
    my ($normal,$tumor)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $dir="$outdir/somatic/$samplePair/mutectSNV";
    my $shell1="$dir/shell/mutect.$samplePair.process.sh";
    my $cmd1="perl $configs{bin}/Mutect/mutectFilter.pl $dir/raw/$samplePair.mutect.txt.gz $dir/raw/$samplePair.mutect.vcf.gz 0.02 0.05 $dir/result/$samplePair.snv.txt.gz $dir/result/$samplePair.snv.vcf && \\\n";
    $cmd1.="perl $configs{annovar}/table_annovar.pl $dir/result/$samplePair.snv.vcf $configs{annovardb} -buildver hg19 -out $dir/result/$samplePair.snv $configs{annovarPara} && \\\n";
    $cmd1.="gzip -f $dir/result/$samplePair.snv.hg19_multianno.txt && \\\n";
    $cmd1.="rm -rf $dir/result/$samplePair.snv.avinput $dir/result/$samplePair.snv.hg19_multianno.vcf";
    generateShell($shell1,$cmd1);

    my $shell2="$dir/shell/mutect.$samplePair.vcf2maf.sh";
    my $cmd2 = "sh $configs{bin}/vcf2maf.sh $dir/result/$samplePair.snv.vcf $dir/result/$samplePair.snv.maf $normal $tumor && \\\n";
    $cmd2 .= "gzip -f $dir/result/$samplePair.snv.vep.vcf $dir/result/$samplePair.snv.maf";
    generateShell($shell2,$cmd2);
}

sub RunMutect2{
    my ($normal,$tumor,$nbam,$tbam,$chr,$bed)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $dir="$outdir/somatic/$samplePair/mutect2";
    makedir($dir) if(!-e $dir);
    my $shell="$dir/shell/mutect.$samplePair.$chr.$bed.sh";
    my $cmd="$configs{java8} -Xmx4g -Djava.io.tmpdir=$dir/tmp -XX:-UseGCOverheadLimit -jar $configs{GATK} -T MuTect2 -R $configs{reference} -I:normal $nbam -I:tumor $tbam -o $dir/raw/$samplePair.$chr.$bed.mutect2.vcf.gz";
    $cmd.=" --dbsnp $configs{mutectDbsnp}";
    $cmd.=" --cosmic $configs{mutectCosmic}";
    $cmd.=" -nct 4";
    if ($configs{seqType} =~ /WES/i){
    	$cmd.=" -L $configs{Splitbed}/$chr.list.$bed.bed";
    }elsif($configs{seqType} =~ /WGS/i){
        $cmd.=" -L $chr";
    }
    generateShell($shell,$cmd);
}

sub mutect2_process{
    my ($normal,$tumor)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $dir="$outdir/somatic/$samplePair/mutect2";
    my $shell1="$dir/shell/mutect2.$samplePair.process.sh";
    my $cmd1 = "sh $configs{bin}/vcf2maf_vcfID.sh $dir/result/$samplePair.mutect2.vcf $dir/result/$samplePair.mutect2.maf $normal $tumor NORMAL TUMOR && \\\n";
    $cmd1 .= "gzip -f $dir/result/$samplePair.mutect2.vcf $dir/result/$samplePair.mutect2.vep.vcf $dir/result/$samplePair.mutect2.maf";
    generateShell($shell1,$cmd1);
}

## calling somatic SNV by Varscan
sub varscanSNV{
    my ($normal,$tumor,$nbam,$tbam,$chr)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $dir="$outdir/somatic/$samplePair/varscanSNV";
    makedir($dir) if(!-e $dir);

    my $shell1 = "$dir/shell/varscan.$samplePair.$chr.sh";
    my $cmd="mkfifo $dir/raw/$samplePair.$chr.fifo ;";
    $cmd.=" $configs{samtools} mpileup $configs{samtoolsMpileupPara}";
    if ($configs{seqType} =~ /WES/i){
        $cmd.=" -l $outdir/TR/CallVariantRegion/ex_region.sort.$chr.bed";
    }elsif($configs{seqType} =~ /WGS/i){
        $cmd.="";#do nothing in samtools mpileup -l 
    }
    $cmd.=" -f $configs{reference} $nbam $tbam >$dir/raw/$samplePair.$chr.fifo & ";
    $cmd.="$configs{java} -Xmx3g -jar $configs{varscan} somatic $dir/raw/$samplePair.$chr.fifo --output-snp $dir/raw/$samplePair.$chr.snp --output-indel $dir/raw/$samplePair.$chr.indel $configs{varscanPara} && \\\n";
#    $cmd.="gzip -f $dir/raw/$samplePair.$chr.snp $dir/raw/$samplePair.$chr.indel && \\\n";
    $cmd.="rm -f $dir/raw/$samplePair.$chr.fifo";
    generateShell($shell1,$cmd);    
}

sub varscanSNV_process{
    my ($normal,$tumor)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $dir="$outdir/somatic/$samplePair/varscanSNV";
    my $shell1 = "$dir/shell/varscan.$samplePair.process.sh";
    my $cmd1="#perl $configs{annovar}/table_annovar.pl $dir/raw/$samplePair.varscan.snp.vcf $configs{annovardb} -buildver hg19 -out $dir/result/$samplePair.varscan.snp $configs{annovarPara} && \\\n";
    $cmd1.="#perl $configs{annovar}/table_annovar.pl $dir/raw/$samplePair.varscan.indel.vcf $configs{annovardb} -buildver hg19 -out $dir/result/$samplePair.varscan.indel $configs{annovarPara} && \\\n";
    $cmd1.="#gzip -f $dir/result/$samplePair.varscan.snp.hg19_multianno.txt && \\\n";
    $cmd1.="#rm -rf $dir/result/$samplePair.varscan.snp.avinput $dir/result/$samplePair.varscan.snp.hg19_multianno.vcf && \\\n";
    $cmd1.="#gzip -f $dir/result/$samplePair.varscan.indel.hg19_multianno.txt && \\\n";
    $cmd1.="#rm -rf $dir/result/$samplePair.varscan.indel.avinput $dir/result/$samplePair.varscan.indel.hg19_multianno.vcf";
    generateShell($shell1,$cmd1);

    my $shell2 = "$dir/shell/varscan.$samplePair.vcf2maf.snp.sh";
    my $cmd2 = "sh $configs{bin}/vcf2maf_vcfID.sh $dir/raw/$samplePair.varscan.snp.vcf $dir/result/$samplePair.varscan.snp.maf $normal $tumor NORMAL TUMOR && \\\n";
    $cmd2 .= "gzip -f $dir/raw/$samplePair.varscan.snp.vep.vcf $dir/result/$samplePair.varscan.snp.maf";
    generateShell($shell2,$cmd2);

    my $shell3 = "$dir/shell/varscan.$samplePair.vcf2maf.indel.sh";
    my $cmd3 = "sh $configs{bin}/vcf2maf_vcfID.sh $dir/raw/$samplePair.varscan.indel.vcf $dir/result/$samplePair.varscan.indel.maf $normal $tumor NORMAL TUMOR && \\\n";
    $cmd3 .= "gzip -f $dir/raw/$samplePair.varscan.indel.vep.vcf $dir/result/$samplePair.varscan.indel.maf"; 
    generateShell($shell3,$cmd3);
}

sub varscanSNV_native_process{
    my ($normal,$tumor)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $dir="$outdir/somatic/$samplePair/varscanSNV";
    my $shell2 = "$dir/shell/varscan.$samplePair.process.sh";
    my $cmd = << "HERE";
zcat $dir/raw/$samplePair.snp.gz | perl -ane '\$F[1]=\$F[1]."\\t".\$F[1];print join("\\t",\@F),"\\n";' | sed '1d' > $dir/result/$samplePair.snp.annovar && \\
perl $configs{annovar}/table_annovar.pl $dir/result/$samplePair.snp.annovar $configs{annovardb} -buildver hg19 -out $dir/result/$samplePair.snp -remove -protocol refGene,ensGene,knownGene,cytoBand,genomicSuperDups,gwasCatalog,phastConsElements100way,phastConsElements46way,targetScanS,tfbsConsSites,wgRna,avsnp147,clinvar_20170130,cosmic82_coding,cosmic82_noncoding,dbnsfp33a,exac03nonpsych,exac03nontcga,gnomad_exome,gnomad_genome,popfreq_all_20150413,intervar_20170202 -operation g,g,g,r,r,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f -nastring . -otherinfo && \\
gzip -f $dir/result/$samplePair.snp.hg19_multianno.txt && \\
rm $dir/result/$samplePair.snp.annovar
HERE
    chomp($cmd);
    generateShell($shell2,$cmd);
}

# merge GATK and Varscan's SNP results
sub varMerge{
    my ($samp,$normal,$tumor)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $dir="$outdir/somatic/$samplePair/varMerge";
    my $varscanDir="$outdir/somatic/$samplePair/varscanSNV";
    makedir($dir);
    
    my $shell1 = "$dir/shell/SNPmerge.$samplePair.sh";
    my $cmd1="perl $configs{bin}/02SNP/SNP_GATK_union.pl $samp $outdir/$normal/gatkDISCOVERY/result/$normal.snp.hg19_multianno.txt.gz $outdir/$tumor/gatkDISCOVERY/result/$tumor.snp.hg19_multianno.txt.gz $dir/result/$samp.GATK_union.stat > $dir/result/$samp.GATK_union.snp && \\\n";
    $cmd1.="perl $configs{bin}/02SNP/SNP_GATK_Varscan_union.pl $samp $varscanDir/result/$samplePair.snp.hg19_multianno.txt.gz $dir/result/$samp.GATK_union.snp $dir/result/$samp.GATK_Varscan_union.stat > $dir/result/$samp.GATK_Varscan_union.snp && \\\n";
    $cmd1.="gzip $dir/result/$samp.GATK_union.snp $dir/result/$samp.GATK_Varscan_union.snp";
    generateShell($shell1,$cmd1);    
}

## calling somatic indel by Platypus
sub platypusSindel{
    my ($normal,$tumor,$nbam,$tbam,$chr)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $dir="$outdir/somatic/$samplePair/platypusSindel";
    makedir($dir) if(!-e $dir);

    my $shell1="$dir/shell/sindel.$samplePair.$chr.sh";
    my $platCmd= "export PATH=$configs{pythonBin}:\$PATH\n";
    $platCmd.="export LD_LIBRARY_PATH=$configs{pythonLib}:$configs{htslib}:\$LD_LIBRARY_PATH\n";
    $platCmd.= "export PYTHONPATH=$configs{pythonPath}:\$PYTHONPATH\n";
    if($configs{seqType} =~ /WES/i){
    $platCmd.= "$configs{pythonBin}/python $configs{Platypus}/Platypus.py callVariants $configs{PlatypusPara} --refFile=$configs{reference} --bamFiles=$nbam,$tbam --regions=$outdir/TR/CallVariantRegion/ex_region.sort.$chr.bed --output=$dir/raw/$samplePair.$chr.platypus.vcf";
    }elsif($configs{seqType} =~ /WGS/i){
    $platCmd.= "$configs{pythonBin}/python $configs{Platypus}/Platypus.py callVariants $configs{PlatypusPara} --refFile=$configs{reference} --bamFiles=$nbam,$tbam  --output=$dir/raw/$samplePair.$chr.platypus.vcf";
    }
    generateShell($shell1,$platCmd);
}

sub platypusSindel_process{
    my ($normal,$tumor)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $dir="$outdir/somatic/$samplePair/platypusSindel";
    my $shell1="$dir/shell/sindel.$samplePair.process.sh";

    my $smtIndelCmd = "$configs{pythonBin}/python $configs{Platypus}/somaticMutationDetector.py $configs{PlatypussmtStrictPara} --inputVCF $dir/raw/$samplePair.platypus.vcf --outputVCF $dir/result/$samplePair.platypus.somatic.vcf --tumourSample $tumor --normalSample $normal && \\\n";
    $smtIndelCmd.="sh $configs{bin}/vcf2maf_vcfID.sh $dir/result/$samplePair.platypus.somatic.vcf $dir/result/$samplePair.platypus.somatic.maf $normal $tumor NORMAL TUMOR && \\\n";
    $smtIndelCmd .= "gzip -f $dir/result/$samplePair.platypus.somatic.vcf $dir/result/$samplePair.platypus.somatic.maf";
    generateShell($shell1,$smtIndelCmd);
}

## calling somatic indel by gatk
sub gatkSindel{
    my ($normal,$tumor,$nbam,$tbam,$chr)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $gatkSindelDir="$outdir/somatic/$samplePair/gatkSindel";
    makedir($gatkSindelDir) if(!-e $gatkSindelDir);
    my $gatkCmd;
    my $shell1 ="$gatkSindelDir/shell/sindel.$samplePair.$chr.sh";
    if($configs{seqType} =~ /WES/i){
    $gatkCmd="$configs{java} -Xmx4g -jar $configs{GATK23} -T SomaticIndelDetector -l INFO -R $configs{reference} -I:normal $nbam -I:tumor $tbam -o $gatkSindelDir/raw/$chr.indel.vcf -L $outdir/TR/CallVariantRegion/ex_region.sort.$chr.bed $configs{GATKsomaticIndelPara}";
    }elsif($configs{seqType} =~ /WGS/i){
    $gatkCmd="$configs{java} -Xmx4g -jar $configs{GATK23} -T SomaticIndelDetector -l INFO -R $configs{reference} -I:normal $nbam -I:tumor $tbam -o $gatkSindelDir/raw/$chr.indel.vcf $configs{GATKsomaticIndelPara}";
    }
    generateShell($shell1,$gatkCmd);
}

sub gatkSindel_process{
    my ($normal,$tumor)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $gatkSindelDir="$outdir/somatic/$samplePair/gatkSindel";
    my $shell2="$gatkSindelDir/shell/sindel.$samplePair.process.sh";
    my $filtCmd= "perl $configs{bin}/gatk/GATKSomaticIndel_filt_fisher.pl $configs{GATKsmtFlitPara} -i $gatkSindelDir/raw/indel.vcf -normal $normal -o $gatkSindelDir/result/$samplePair.sindel && \\\n";
    $filtCmd.= "perl $configs{bin}/gatk/GATKSomaticIndel_filtanno.pl -vcf $gatkSindelDir/result/$samplePair.sindel.filt-repeat_overlap.xls $gatkSindelDir/result/$samplePair.sindel.filt.vcf $gatkSindelDir/result/$samplePair.sindel.filt.fisher.vcf && \\\n";
    $filtCmd.= "perl $configs{annovar}/table_annovar.pl $gatkSindelDir/result/$samplePair.sindel.filt.fisher.vcf $configs{annovardb} -buildver hg19 -out $gatkSindelDir/result/$samplePair.sindel.filt.fisher $configs{annovarPara} && \\\n";
    $filtCmd.= "gzip -f $gatkSindelDir/result/$samplePair.sindel.filt.fisher.hg19_multianno.txt $gatkSindelDir/raw/indel.vcf && \\\n";
    $filtCmd.= "rm -rf $gatkSindelDir/result/$samplePair.sindel.filt.fisher.avinput";
    generateShell($shell2,$filtCmd);
}

sub mergeGATKVars{
    my ($samp)=@_;
    my $mergeVarsDir="$outdir/$samp/gatkDISCOVERY"; 

    my $shell="$mergeVarsDir/shell/mergeVars.$samp.sh";
    my $mergeVarsCmd="export PERL5LIB=$configs{bin}/vcftools/lib && \\\n";
    $mergeVarsCmd.="perl $configs{bin}/mergeVar.pl $mergeVarsDir/result snp.hg19_multianno.txt.gz 1 $samp.snp.hg19_multianno.txt && \\\n";
    $mergeVarsCmd.="perl $configs{bin}/mergeVar.pl $mergeVarsDir/result indel.hg19_multianno.txt.gz 1 $samp.indel.hg19_multianno.txt && \\\n";
    my ($files1,$files2);
    for my $chr (@chrs){
        $files1.="$mergeVarsDir/result/$samp.$chr.snp.vcf.gz ";
        $files2.="$mergeVarsDir/result/$samp.$chr.indel.vcf.gz ";
    }
    $mergeVarsCmd.="perl $configs{bin}/vcf-concat $files1 > $mergeVarsDir/result/$samp.snp.vcf && \\\n";
    $mergeVarsCmd.="perl $configs{bin}/vcf-concat $files2 > $mergeVarsDir/result/$samp.indel.vcf && \\\n";
    $mergeVarsCmd.="gzip -f $mergeVarsDir/result/$samp.snp.vcf && \\\n";
    $mergeVarsCmd.="gzip -f $mergeVarsDir/result/$samp.indel.vcf && \\\n";
    $mergeVarsCmd.="rm $files1 $files2";
    generateShell($shell,$mergeVarsCmd);
}

sub mergeSoVars{
    my ($normal,$tumor,$type)=@_;
    my $samplePair="$normal-VS-$tumor";
    my $mergeVarsDir="$outdir/somatic/$samplePair/$type";

    my $shell="$mergeVarsDir/shell/mergeVars.$samplePair.sh";
    my $mergeVarsCmd="export PERL5LIB=$configs{bin}/vcftools/lib && \\\n"; 
    if($type eq "mutectSNV"){
        my ($files1);
        for my $chr (@chrs){
            $files1.="$mergeVarsDir/raw/$samplePair.$chr.mutect.vcf.gz ";
        }
	if($type eq "mutectSNV"){
	        $mergeVarsCmd.="perl $configs{bin}/mergeVar.pl $mergeVarsDir/raw mutect.txt.gz 2 $samplePair.mutect.txt && \\\n";
	}
        $mergeVarsCmd.="perl $configs{bin}/vcf-concat $files1 > $mergeVarsDir/raw/$samplePair.mutect.vcf && \\\n";
        $mergeVarsCmd.="gzip -f $mergeVarsDir/raw/$samplePair.mutect.vcf";
        #$mergeVarsCmd.="rm $files1 $mergeVarsDir/raw/*mutect.vcf.gz.idx";
    }
    if($type eq "mutect2"){
	my ($files1);
        for my $chr (@chrs){
            for my $n(1..$SplitbedCount{$chr}){
                 $files1.="$mergeVarsDir/raw/$samplePair.$chr.$n.mutect2.vcf.gz \\\n";
            }
        }
        $mergeVarsCmd.="perl $configs{bin}/vcf-concat $files1 > $mergeVarsDir/result/$samplePair.mutect2.vcf";
	#$mergeVarsCmd.="gzip -f $mergeVarsDir/result/$samplePair.mutect2.vcf";
    }
    elsif($type eq "varscan_native"){
        $mergeVarsCmd.="perl $configs{bin}/mergeVar.pl $mergeVarsDir/raw indel 1 $samplePair.indel && \\\n";
        $mergeVarsCmd.="perl $configs{bin}/mergeVar.pl $mergeVarsDir/raw snp 1 $samplePair.snp";
    }
    elsif($type eq "varscanSNV"){
        my ($files1,$files2);
        for my $chr (@chrs){
            $files1.="$mergeVarsDir/raw/$samplePair.$chr.snp.vcf ";
            $files2.="$mergeVarsDir/raw/$samplePair.$chr.indel.vcf ";
        }
        $mergeVarsCmd.="perl $configs{bin}/vcf-concat $files1 > $mergeVarsDir/raw/$samplePair.varscan.snp.vcf && \\\n";
        #$mergeVarsCmd.="gzip -f $mergeVarsDir/raw/$samplePair.varscan.snp.vcf && \\\n";
        $mergeVarsCmd.="perl $configs{bin}/vcf-concat $files2 > $mergeVarsDir/raw/$samplePair.varscan.indel.vcf";
        #$mergeVarsCmd.="gzip -f $mergeVarsDir/raw/$samplePair.varscan.indel.vcf";
    }
    elsif($type eq "gatkSindel"){
        my ($files1);
        for my $chr (@chrs){
            $files1.="$mergeVarsDir/raw/$chr.indel.vcf ";
        }
        $mergeVarsCmd.="perl $configs{bin}/vcf-concat $files1 > $mergeVarsDir/raw/indel.vcf";
        #$mergeVarsCmd.="rm $files1 $mergeVarsDir/raw/*indel.vcf.idx";
    }
    elsif($type eq "platypusSindel"){
        my ($files1,$files2);
        for my $chr (@chrs){
            $files1.="$mergeVarsDir/raw/$samplePair.$chr.platypus.vcf.gz ";
            #$files2.="$mergeVarsDir/raw/$samplePair.$chr.platypus.somatic.vcf.gz ";
        }
        $mergeVarsCmd.="perl $configs{bin}/vcf-concat $files1 > $mergeVarsDir/raw/$samplePair.platypus.vcf";
        #$mergeVarsCmd.="perl $configs{bin}/vcf-concat $files2 > $mergeVarsDir/raw/$samplePair.platypus.somatic.vcf";
        #$mergeVarsCmd.="gzip -f $mergeVarsDir/raw/$samplePair.platypus.vcf $mergeVarsDir/raw/$samplePair.platypus.somatic.vcf";
        #$mergeVarsCmd.="rm $files1 $files2";
    } 
    generateShell($shell,$mergeVarsCmd);
}

sub ADTExSCNV{
	my ($normal,$tumor) = @_;
	my $samplePair="$normal-VS-$tumor";
	my $sampledir="$outdir/somatic/$samplePair/ADTEx";
	makedir($sampledir) if(!-e $sampledir);
        
        my $ADTExbin = "$configs{bin}/05CNV/ADTEx/ADTEx";
	my $shell    = "$sampledir/shell/ADTEx.$samplePair.sh";
	my $ADTExCmd = "export PATH=$configs{bin}:\$PATH\n";
	$ADTExCmd.= "perl $ADTExbin/ADTEx_vcf2baf.pl $outdir/$normal/gatkGVCF/result/$normal.filtered_snp.vcf.gz $outdir/$tumor/gatkGVCF/result/$tumor.filtered_snp.vcf.gz > $sampledir/raw/$samplePair.baf && \\\n";
	$ADTExCmd.= "perl $ADTExbin/depth4ADTEx.pl $outdir/$normal/coverage/result/depth.tsv.gz $configs{targetRegion} > $sampledir/raw/$normal.cov && \\\n";
	$ADTExCmd.= "perl $ADTExbin/depth4ADTEx.pl $outdir/$tumor/coverage/result/depth.tsv.gz $configs{targetRegion} > $sampledir/raw/$tumor.cov && \\\n";
	$ADTExCmd.="rm -r $sampledir/result && \\\n";
	$ADTExCmd.="$configs{pythonBin}/python $ADTExbin/ADTEx.py --normal $sampledir/raw/$normal.cov --tumor $sampledir/raw/$tumor.cov --bed $configs{targetRegion} --out $sampledir/result --DOC --baf $sampledir/raw/$samplePair.baf --estimatePloidy && \\\n";
	$ADTExCmd.="$configs{bin}/Rscript $ADTExbin/plot_results.R $sampledir/result && \\\n";
        $ADTExCmd.="$configs{bin}/Rscript $ADTExbin/summary.R $sampledir/result/cnv.result $sampledir/result/extracted.baf $sampledir/result/$samplePair\_summary.pdf && \\\n";
        $ADTExCmd.="$configs{bin}/Rscript $ADTExbin/ADTEx2AllelicCapsegPNG.r $sampledir/result/cnv.result $sampledir/result/zygosity/zygosity.res $sampledir/result/ $samplePair && \\\n";
        $ADTExCmd.="rm $sampledir/raw/$normal.cov $sampledir/raw/$tumor.cov";
	generateShell($shell,$ADTExCmd)
}

sub ABSOLUTE{
        my ($normal,$tumor) = @_;
        my $samplePair="$normal-VS-$tumor";
        my $sampledir="$outdir/somatic/$samplePair/ABSOLUTE";
        makedir($sampledir) if(!-e $sampledir);
        my $shell = "$sampledir/shell/ABSOLUTE.$samplePair.sh";
        my $ABSOLUTECmd = "perl $configs{bin}/ADTEx2seg.pl $outdir/somatic/$samplePair/ADTEx/result/cnv.result $sampledir/raw/$samplePair.seg && \\\n";
        #$ABSOLUTECmd.= "perl $configs{bin}/Mutect/mutect2maf.pl $samplePair $outdir/somatic/$samplePair/mutectSNV/result/$samplePair.snv.hg19_multianno.txt.gz $sampledir/raw && \\\n";
        $ABSOLUTECmd.= "perl $configs{bin}/ABSOLUTE/maf_transher.v2.pl /hwfssz1/ST_PRECISION/PMO/F14ZF1QSSY1656/Breast/01WES/Zebra/result/somatic/$samplePair/mutectSNV/result/$samplePair.snv.maf.gz $sampledir/raw/$samplePair.maf $samplePair && \\\n";
	$ABSOLUTECmd.= "$configs{bin}/Rscript $configs{bin}/ABSOLUTE/run_ABSOLUTE.R $sampledir/raw/$samplePair.seg $sampledir/raw/$samplePair.maf $sampledir/result $samplePair Illumina_WES";
        generateShell($shell,$ABSOLUTECmd)

}


sub edgeList_Mutect2{
	`mkdir -p $outdir/shell_run/`;
        chomp(my $day=`date +%Y%m%d`);
        open EDGE,">$outdir/shell_run/edge_mutect2.$day.list" or die $!;
	if(exists $configs{somatic}){
		foreach my $samplePair (sort keys %samplePairs){
                	my $normal=$samplePairs{$samplePair}{normal};
         		my $tumor =$samplePairs{$samplePair}{tumor};
			foreach my $chr(@chrs){
				for my $bed (1..$SplitbedCount{$chr}){
				print EDGE "$outdir/somatic/$samplePair/mutect2/shell/mutect.$samplePair.$chr.$bed.sh:12G:4cpu $outdir/somatic/$samplePair/mutect2/shell/mergeVars.$samplePair.sh:1G:4cpu \n";
				}
			}
		print EDGE "$outdir/somatic/$samplePair/mutect2/shell/mergeVars.$samplePair.sh:1G:1cpu $outdir/somatic/$samplePair/mutect2/shell/mutect2.$samplePair.process.sh:1G:1cpu \n";
		}
	}
	close EDGE;
}

sub edgeList2{
	`mkdir -p $outdir/shell_run/`;
	chomp(my $day=`date +%Y%m%d`);
	open EDGE1,">$outdir/shell_run/edge1.$day.list" or die $!;
    	foreach my $samp(sort keys %bwaBams){
        	foreach my $lib(sort keys %{$bwaBams{$samp}}){
            		foreach my $laneIndex(sort keys %{$bwaBams{$samp}{$lib}}){
            		print EDGE1 "$outdir/$samp/$laneIndex/shell/clean.sh:1G:4cpu $outdir/$samp/$laneIndex/shell/cleanQC.sh:1G:1cpu\n";
            		}
        	}
   	}
	close EDGE1;
	#######edico-edge2
	open EDGE2,">$outdir/shell_run/edge2_edico.$day.list" or die $!;
	foreach my $samp(sort keys %bwaBams){
		print EDGE2 "$outdir/$samp/edico/shell/edico.$samp.sh\n";
	}
	close EDGE2;
	open EDGE2,">$outdir/shell_run/edge2_edico_process.$day.list" or die $!;
	foreach my $samp(sort keys %bwaBams){
                print EDGE2 "$outdir/$samp/edico/shell/edico_process.sh\n";
        }
	close EDGE2;
	###########edge3 no co-realn and M2 need to be added
	open EDGE3,">$outdir/shell_run/edge3.$day.list" or die $!;
	foreach my $samp(sort keys %bwaBams){
		print EDGE3 "$outdir/$samp/rmdup/shell/rmdupQC.sh:1G:1cpu\n";
       		print EDGE3 "$outdir/$samp/coverage/shell/covdep.$samp.sh:1G:1cpu\n";
        	print EDGE3 "$outdir/$samp/coverage/shell/targetqc.$samp.sh:7G:1cpu\n";
        	if($configs{gatkGVCF}=~/true/i){
        		print EDGE3 "$outdir/$samp/gatkGVCF/shell/gatkGVCF.$samp.sh:10G:1cpu $outdir/$samp/gatkGVCF/shell/gatkGVCF.$samp.INDEL.sh:10G:1cpu\n";
        		print EDGE3 "$outdir/$samp/gatkGVCF/shell/gatkGVCF.$samp.sh:10G:1cpu $outdir/$samp/gatkGVCF/shell/gatkGVCF.$samp.SNP.sh:10G:1cpu\n";
        	}
		foreach my $chr(@chrs){
			print EDGE3 "$outdir/$samp/realn/shell/GATKrealnRecal.$samp.$chr.sh:15G:4cpu\n";
		}
	}
	######before somatic more model ???
	if(exists $configs{somatic}){
		foreach my $samplePair (sort keys %samplePairs){
        		my $normal=$samplePairs{$samplePair}{normal};
        		my $tumor =$samplePairs{$samplePair}{tumor};
        		foreach my $chr(@chrs){
				my $NormalBrecalEdge = "$outdir/$normal/realn/shell/GATKrealnRecal.$normal.$chr.sh:15G:4cpu";
            			my $TumorBrecalEdge  = "$outdir/$tumor/realn/shell/GATKrealnRecal.$tumor.$chr.sh:15G:4cpu";
        			if($configs{mutectSNV}=~/true/i){
        				my $shellDir="$outdir/somatic/$samplePair/mutectSNV/shell";
        				print EDGE3 "$NormalBrecalEdge $shellDir/mutect.$samplePair.$chr.sh:8G:1cpu\n";
                			print EDGE3 "$TumorBrecalEdge $shellDir/mutect.$samplePair.$chr.sh:8G:1cpu\n";
                			print EDGE3 "$shellDir/mutect.$samplePair.$chr.sh:8G:1cpu $shellDir/mergeVars.$samplePair.sh:1G:1cpu\n";
                			if(!exists $samplePairs{$samplePair}{mutectSNV}){
                    				print EDGE3 "$shellDir/mergeVars.$samplePair.sh:1G:1cpu $shellDir/mutect.$samplePair.process.sh:1G:1cpu\n";
                    				print EDGE3 "$shellDir/mutect.$samplePair.process.sh:1G:1cpu $shellDir/mutect.$samplePair.vcf2maf.sh:2G:1cpu\n";
                    				$samplePairs{$samplePair}{mutectSNV}="done";
                			}
                		}
			#
				if($configs{gatkSindel}=~/true/i){
                			my $shellDir="$outdir/somatic/$samplePair/gatkSindel/shell";
                			print EDGE3 "$NormalBrecalEdge $shellDir/sindel.$samplePair.$chr.sh:6G:1cpu\n";
                			print EDGE3 "$TumorBrecalEdge $shellDir/sindel.$samplePair.$chr.sh:6G:1cpu\n";
                			print EDGE3 "$shellDir/sindel.$samplePair.$chr.sh:6G:1cpu $shellDir/mergeVars.$samplePair.sh:1G:1cpu\n";
                			if(!exists $samplePairs{$samplePair}{gatkSindel}){
                    				print EDGE3 "$shellDir/mergeVars.$samplePair.sh:1G:1cpu $shellDir/sindel.$samplePair.process.sh:1G:1cpu\n";
                    				$samplePairs{$samplePair}{gatkSindel}="done";
                			}
            			}
            			if($configs{platypusSindel}=~/true/i){
                			my $shellDir="$outdir/somatic/$samplePair/platypusSindel/shell";
                			print EDGE3 "$NormalBrecalEdge $shellDir/sindel.$samplePair.$chr.sh:5G:1cpu\n";
                			print EDGE3 "$TumorBrecalEdge $shellDir/sindel.$samplePair.$chr.sh:5G:1cpu\n";
                			print EDGE3 "$shellDir/sindel.$samplePair.$chr.sh:5G:1cpu $shellDir/mergeVars.$samplePair.sh:1G:1cpu\n";
                			if(!exists $samplePairs{$samplePair}{platypusSindel}){
                    				print EDGE3 "$shellDir/mergeVars.$samplePair.sh:1G:1cpu $shellDir/sindel.$samplePair.process.sh:1G:1cpu\n";
                    				$samplePairs{$samplePair}{platypusSindel}="done";
                			}
            			}
            			if($configs{varscanSNV}=~/true/i){
                			my $shellDir="$outdir/somatic/$samplePair/varscanSNV/shell";
               				print EDGE3 "$NormalBrecalEdge $shellDir/varscan.$samplePair.$chr.sh:10G:1cpu\n";
                			print EDGE3 "$TumorBrecalEdge $shellDir/varscan.$samplePair.$chr.sh:10G:1cpu\n";
                			print EDGE3 "$shellDir/varscan.$samplePair.$chr.sh:10G:1cpu $shellDir/mergeVars.$samplePair.sh:1G:1cpu\n";
                			if(!exists $samplePairs{$samplePair}{varscanSNV}){
                    				print EDGE3 "$shellDir/mergeVars.$samplePair.sh:1G:1cpu $shellDir/varscan.$samplePair.process.sh:2G:1cpu\n";
                    				print EDGE3 "$shellDir/mergeVars.$samplePair.sh:1G:1cpu $shellDir/varscan.$samplePair.vcf2maf.snp.sh:10G:1cpu\n";
                    				print EDGE3 "$shellDir/mergeVars.$samplePair.sh:1G:1cpu $shellDir/varscan.$samplePair.vcf2maf.indel.sh:10G:1cpu\n";
                    				$samplePairs{$samplePair}{varscanSNV}="done";
                			}
            			}
			}
		}
	}
	close EDGE3;

}

## generate the execution script 
sub edgeList{
    `mkdir -p $outdir/shell_run/`;
    chomp(my $day=`date +%Y%m%d`);
    open EDGE,">$outdir/shell_run/edge.$day.list" or die $!;
    foreach my $samp(sort keys %bwaBams){
        foreach my $lib(sort keys %{$bwaBams{$samp}}){
            foreach my $laneIndex(sort keys %{$bwaBams{$samp}{$lib}}){
                if($configs{Aligment} eq "bwa_aln"){
                    print EDGE "$outdir/$samp/$laneIndex/shell/clean.sh:1G:4cpu $outdir/$samp/$laneIndex/shell/bwa_aln1.sh:4G:1cpu\n";
                    print EDGE "$outdir/$samp/$laneIndex/shell/clean.sh:1G:4cpu $outdir/$samp/$laneIndex/shell/bwa_aln2.sh:4G:1cpu\n";
                    print EDGE "$outdir/$samp/$laneIndex/shell/bwa_aln1.sh:4G:1cpu $outdir/$samp/$laneIndex/shell/bwa.sh:6G:1cpu\n";
                    print EDGE "$outdir/$samp/$laneIndex/shell/bwa_aln2.sh:4G:1cpu $outdir/$samp/$laneIndex/shell/bwa.sh:6G:1cpu\n";
                }
                elsif($configs{Aligment} eq "bwa_mem"){
                    print EDGE "$outdir/$samp/$laneIndex/shell/clean.sh:1G:4cpu $outdir/$samp/$laneIndex/shell/bwa.sh:6G:1cpu\n";
                }
                print EDGE "$outdir/$samp/$laneIndex/shell/clean.sh:1G:4cpu $outdir/$samp/$laneIndex/shell/cleanQC.sh:1G:1cpu\n";
                print EDGE "$outdir/$samp/$laneIndex/shell/bwa.sh:6G:1cpu $outdir/$samp/$laneIndex/shell/laneQC.sh:1G:1cpu\n";
                print EDGE "$outdir/$samp/$laneIndex/shell/bwa.sh:6G:1cpu $outdir/$samp/rmdup/shell/rmdup.$samp.sh:8G:1cpu\n"; 
            }
        }
        print EDGE "$outdir/$samp/rmdup/shell/rmdup.$samp.sh:8G:1cpu $outdir/$samp/rmdup/shell/rmdupQC.sh:1G:1cpu\n";
        print EDGE "$outdir/$samp/rmdup/shell/rmdup.$samp.sh:8G:1cpu $outdir/$samp/coverage/shell/covdep.$samp.sh:1G:1cpu\n";
        print EDGE "$outdir/$samp/rmdup/shell/rmdup.$samp.sh:8G:1cpu $outdir/$samp/coverage/shell/targetqc.$samp.sh:7G:1cpu\n";

        my ($realnEdge,$brecalEdge,$patientId);
        if($configs{Corealn} && $configs{RealnType} eq "co-realn"){
            $patientId=$patientIds{$samp};
            $maxmems{$patientId} =~ s/g/G/;
        }
        foreach my $chr(@chrs){
            if($configs{Corealn} && $configs{RealnType} eq "co-realn"){
                $realnEdge  = "$outdir/Co-realn/shell/realn.$patientId.$chr.sh:$maxmems{$patientId}:1cpu";
                $brecalEdge = "$outdir/$samp/brecal/shell/brecal.$samp.$chr.sh:8G:1cpu";
                print EDGE "$realnEdge $brecalEdge\n";
            }else{
                $realnEdge  = "$outdir/$samp/realn/shell/GATKrealnRecal.$samp.$chr.sh:15G:4cpu";
                $brecalEdge = "$outdir/$samp/realn/shell/GATKrealnRecal.$samp.$chr.sh:15G:4cpu";
            	#print EDGE "$realnEdge $outdir/$samp/brecal/shell/mergeBam.$samp.sh:1G:1cpu\n";
	    }

            print EDGE "$outdir/$samp/rmdup/shell/rmdup.$samp.sh:8G:1cpu $realnEdge\n";


            if($configs{gatkDISCOVERY}=~/true/i){
                my $shellDir="$outdir/$samp/gatkDISCOVERY/shell"; 
                print EDGE "$brecalEdge $shellDir/variant.$samp.$chr.sh:6G:1cpu\n";
                print EDGE "$shellDir/variant.$samp.$chr.sh:6G:1cpu $shellDir/variant.$samp.$chr.process.sh:1G:1cpu\n";
                print EDGE "$shellDir/variant.$samp.$chr.process.sh:1G:1cpu $shellDir/mergeVars.$samp.sh:1G:1cpu\n";
            }
            if($configs{gatkGVCF}=~/true/i){
                my $shellDir="$outdir/$samp/gatkGVCF/shell";
                print EDGE "$brecalEdge $shellDir/gatkGVCF.$samp.$chr.sh:10G:1cpu\n";
                print EDGE "$shellDir/gatkGVCF.$samp.$chr.sh:10G:1cpu $shellDir/gatkGVCF.$samp.merge.sh:10G:1cpu\n";
                print EDGE "$shellDir/gatkGVCF.$samp.$chr.sh:10G:1cpu $outdir/combine/gatkGVCF/shell/gatkGVCF.combine.$chr.sh:10G:1cpu\n";
            }
        }
        if($configs{gatkGVCF}=~/true/i){
            my $shellDir="$outdir/$samp/gatkGVCF/shell";
            print EDGE "$shellDir/gatkGVCF.$samp.merge.sh:10G:1cpu $shellDir/gatkGVCF.$samp.SNP.sh:11G:1cpu\n";
            print EDGE "$shellDir/gatkGVCF.$samp.merge.sh:10G:1cpu $shellDir/gatkGVCF.$samp.INDEL.sh:11G:1cpu\n";
        }
    }

    foreach my $chr(@chrs){
        print EDGE "$outdir/combine/gatkGVCF/shell/gatkGVCF.combine.$chr.sh:10G:1cpu $outdir/combine/gatkGVCF/shell/gatkGVCF.combine.merge.sh:10G:1cpu\n";    
    }
    print EDGE "$outdir/combine/gatkGVCF/shell/gatkGVCF.combine.merge.sh:10G:1cpu $outdir/combine/gatkGVCF/shell/gatkGVCF.combine.SNP.sh:10G:1cpu\n";
    print EDGE "$outdir/combine/gatkGVCF/shell/gatkGVCF.combine.merge.sh:10G:1cpu $outdir/combine/gatkGVCF/shell/gatkGVCF.combine.INDEL.sh:10G:1cpu\n";

    if(exists $configs{somatic}){
    foreach my $samplePair (sort keys %samplePairs){
        my $normal=$samplePairs{$samplePair}{normal};
        my $tumor =$samplePairs{$samplePair}{tumor};
        foreach my $chr(@chrs){
            my ($NormalBrecalEdge,$TumorBrecalEdge);
            if($configs{Corealn} && $configs{RealnType} eq "co-realn"){
                $NormalBrecalEdge = "$outdir/$normal/brecal/shell/brecal.$normal.$chr.sh:6G:1cpu";
                $TumorBrecalEdge  = "$outdir/$tumor/brecal/shell/brecal.$tumor.$chr.sh:6G:1cpu";
            }else{
                $NormalBrecalEdge = "$outdir/$normal/realn/shell/GATKrealnRecal.$normal.$chr.sh:15G:4cpu";
                $TumorBrecalEdge  = "$outdir/$tumor/realn/shell/GATKrealnRecal.$tumor.$chr.sh:15G:4cpu";
            }
           
            if($configs{mutectSNV}=~/true/i){
                my $shellDir="$outdir/somatic/$samplePair/mutectSNV/shell";
                print EDGE "$NormalBrecalEdge $shellDir/mutect.$samplePair.$chr.sh:8G:1cpu\n";
                print EDGE "$TumorBrecalEdge $shellDir/mutect.$samplePair.$chr.sh:8G:1cpu\n";
                print EDGE "$shellDir/mutect.$samplePair.$chr.sh:8G:1cpu $shellDir/mergeVars.$samplePair.sh:1G:1cpu\n";
                if(!exists $samplePairs{$samplePair}{mutectSNV}){
                    print EDGE "$shellDir/mergeVars.$samplePair.sh:1G:1cpu $shellDir/mutect.$samplePair.process.sh:1G:1cpu\n";
                    print EDGE "$shellDir/mutect.$samplePair.process.sh:1G:1cpu $shellDir/mutect.$samplePair.vcf2maf.sh:2G:1cpu\n";
                    $samplePairs{$samplePair}{mutectSNV}="done";
                }
            }
            if($configs{gatkSindel}=~/true/i){
                my $shellDir="$outdir/somatic/$samplePair/gatkSindel/shell";
                print EDGE "$NormalBrecalEdge $shellDir/sindel.$samplePair.$chr.sh:6G:1cpu\n";
                print EDGE "$TumorBrecalEdge $shellDir/sindel.$samplePair.$chr.sh:6G:1cpu\n";
                print EDGE "$shellDir/sindel.$samplePair.$chr.sh:6G:1cpu $shellDir/mergeVars.$samplePair.sh:1G:1cpu\n";
                if(!exists $samplePairs{$samplePair}{gatkSindel}){
                    print EDGE "$shellDir/mergeVars.$samplePair.sh:1G:1cpu $shellDir/sindel.$samplePair.process.sh:1G:1cpu\n";
                    $samplePairs{$samplePair}{gatkSindel}="done";
                }
            }
            if($configs{platypusSindel}=~/true/i){
		my $shellDir="$outdir/somatic/$samplePair/platypusSindel/shell";
                print EDGE "$NormalBrecalEdge $shellDir/sindel.$samplePair.$chr.sh:5G:1cpu\n";
                print EDGE "$TumorBrecalEdge $shellDir/sindel.$samplePair.$chr.sh:5G:1cpu\n";
                print EDGE "$shellDir/sindel.$samplePair.$chr.sh:5G:1cpu $shellDir/mergeVars.$samplePair.sh:1G:1cpu\n";
                if(!exists $samplePairs{$samplePair}{platypusSindel}){
                    print EDGE "$shellDir/mergeVars.$samplePair.sh:1G:1cpu $shellDir/sindel.$samplePair.process.sh:1G:1cpu\n";
                    $samplePairs{$samplePair}{platypusSindel}="done";
                }
            }
            if($configs{varscanSNV}=~/true/i){
                my $shellDir="$outdir/somatic/$samplePair/varscanSNV/shell";
                print EDGE "$NormalBrecalEdge $shellDir/varscan.$samplePair.$chr.sh:10G:1cpu\n";
                print EDGE "$TumorBrecalEdge $shellDir/varscan.$samplePair.$chr.sh:10G:1cpu\n";
                print EDGE "$shellDir/varscan.$samplePair.$chr.sh:10G:1cpu $shellDir/mergeVars.$samplePair.sh:1G:1cpu\n";
                if(!exists $samplePairs{$samplePair}{varscanSNV}){
                    print EDGE "$shellDir/mergeVars.$samplePair.sh:1G:1cpu $shellDir/varscan.$samplePair.process.sh:2G:1cpu\n";
                    print EDGE "$shellDir/mergeVars.$samplePair.sh:1G:1cpu $shellDir/varscan.$samplePair.vcf2maf.snp.sh:10G:1cpu\n";
                    print EDGE "$shellDir/mergeVars.$samplePair.sh:1G:1cpu $shellDir/varscan.$samplePair.vcf2maf.indel.sh:10G:1cpu\n";
                    $samplePairs{$samplePair}{varscanSNV}="done";
                }
            }
        }
    }
    }
    close EDGE;
    
    my @r = ('a'..'z','A'..'Z');
    $configs{subProjectName} ||= join '', map { $r[int rand @r] } 0..6;

    open RUNSH,">$outdir/shell_run/run.$day.sh" or die $!;
    print RUNSH "$configs{monitor} taskmonitor -q $configs{queue} -P $configs{priority} -p $configs{subProjectName} -i $outdir/shell_run/edge.$day.list -f 1\n";
    close RUNSH;
}

sub shlist{
        `mkdir -p $outdir/shell_run/shlist`;
        open Cleanshlt,">$outdir/shell_run/shlist/clean.list" or die $!;
	open Cleanqcshlt,">$outdir/shell_run/shlist/cleanQC.list" or die $!;
        open Edicoshlt,">$outdir/shell_run/shlist/edico.list" or die $!;
        open Edico_processshlt,">$outdir/shell_run/shlist/edico_process.list" or die $!;
        open Covdepshlt,">$outdir/shell_run/shlist/covdep.list" or die $!;
        open Targetqcshlt,">$outdir/shell_run/shlist/targetqc.list" or die $!;
        open GatkGVCFshlt,">$outdir/shell_run/shlist/gvcf2vcf.list" or die $!;
        open VCF2SNPshlt,">$outdir/shell_run/shlist/vcf2snp.list" or die $!;
        open VCF2INDELshlt,">$outdir/shell_run/shlist/vcf2indel.list" or die $!;
        open Brecal_Realnshlt, ">$outdir/shell_run/shlist/brecal_realn.list" or die $!;
        foreach my $samp(sort keys %bwaBams){
                foreach my $lib(sort keys %{$bwaBams{$samp}}){
                        foreach my $laneIndex(sort keys %{$bwaBams{$samp}{$lib}}){
                        print Cleanshlt "$outdir/$samp/$laneIndex/shell/clean.sh\n";
			print Cleanqcshlt "$outdir/$samp/$laneIndex/shell/cleanQC.sh\n";
			}
                }
                print Edicoshlt "$outdir/$samp/edico/shell/edico.sh\n";
                print Edico_processshlt "$outdir/$samp/edico/shell/edico_process.sh\n";
                print Covdepshlt "$outdir/$samp/coverage/shell/covdep.$samp.sh\n";
                print Targetqcshlt "$outdir/$samp/coverage/shell/targetqc.$samp.sh\n";
                print GatkGVCFshlt "$outdir/$samp/gatkGVCF/shell/gatkGVCF.$samp.sh\n";
                print VCF2SNPshlt "$outdir/$samp/gatkGVCF/shell/gatkGVCF.$samp.SNP.sh\n";
                print VCF2INDELshlt "$outdir/$samp/gatkGVCF/shell/gatkGVCF.$samp.INDEL.sh\n";
                foreach my $chr(@chrs){
                        print Brecal_Realnshlt "$outdir/$samp/realn/shell/GATKrealnRecal.$samp.$chr.sh\n";
                }
        }
	close Cleanshlt;close Edicoshlt;close Edico_processshlt;close Covdepshlt;close Targetqcshlt;close GatkGVCFshlt;close VCF2SNPshlt;close VCF2INDELshlt;close Brecal_Realnshlt;
        if(exists $configs{somatic}){
        open Mutect,">$outdir/shell_run/shlist/mutectSNV.list" or die $!;
        open GatkSindel,">$outdir/shell_run/shlist/gatkSindel.list" or die $!;
        open PlatypusSindel,">$outdir/shell_run/shlist/platypusSindel.list" or die $!;
        open VarscanSNV,">$outdir/shell_run/shlist/varscanSNV.list" or die $!;
        open Mutect_process,">$outdir/shell_run/shlist/mutectSNV_process.list" or die $!;
        open GatkSindel_process,">$outdir/shell_run/shlist/gatkSindel_process.list" or die $!;
        open VarscanSNV_process,">$outdir/shell_run/shlist/varscanSNV_process.list" or die $!;
        open PlatypusSindel_process,">$outdir/shell_run/shlist/platypusSindel_process.list" or die $!;
        open VarscanSNV_process,">$outdir/shell_run/shlist/varscanSNV_process.list" or die $!;
        open Mutect2,">$outdir/shell_run/shlist/mutect2.list" or die $!;
	open Mutect2_process,">$outdir/shell_run/shlist/mutect2_process.list" or die $!;
	foreach my $samplePair (sort keys %samplePairs){
                my $normal=$samplePairs{$samplePair}{normal};
                my $tumor =$samplePairs{$samplePair}{tumor};
		foreach my $chr(@chrs){
                print Mutect "$outdir/somatic/$samplePair/mutectSNV/shell/mutect.$samplePair.$chr.sh\n";
                print GatkSindel "$outdir/somatic/$samplePair/gatkSindel/shell/sindel.$samplePair.$chr.sh\n";
                print PlatypusSindel "$outdir/somatic/$samplePair/platypusSindel/shell/sindel.$samplePair.$chr.sh\n";
                print VarscanSNV "$outdir/somatic/$samplePair/varscanSNV/shell/varscan.$samplePair.$chr.sh\n";
                for my $bed (1..$SplitbedCount{$chr}){
                	print Mutect2 "$outdir/somatic/$samplePair/mutect2/shell/mutect.$samplePair.$chr.$bed.sh\n";
                	}
		}
                print Mutect_process "$outdir/somatic/$samplePair/mutectSNV/shell/mutect.$samplePair.process.sh\n";
                print GatkSindel_process "$outdir/somatic/$samplePair/gatkSindel/shell/sindel.$samplePair.process.sh\n";
                print PlatypusSindel_process "$outdir/somatic/$samplePair/platypusSindel/shell/sindel.$samplePair.process.sh\n";
                print VarscanSNV_process "$outdir/somatic/$samplePair/varscanSNV/shell/varscan.$samplePair.process.sh\n";
        	}
        close Mutect;close GatkSindel;close PlatypusSindel;close VarscanSNV;close Mutect_process;close GatkSindel_process;close VarscanSNV_process;close PlatypusSindel_process;close VarscanSNV_process;close Mutect2;close Mutect2_process
	}
}

sub rmlist{
	`mkdir -p $outdir/shell_run/rmlist`;
	open RMcleanfq,">$outdir/shell_run/rmlist/rmcleanfq.list";
	open RMedicobam,">$outdir/shell_run/rmlist/rmedicobam.list";
	foreach my $samp(sort keys %bwaBams){
                foreach my $lib(sort keys %{$bwaBams{$samp}}){
                        foreach my $laneIndex(sort keys %{$bwaBams{$samp}{$lib}}){
                        print RMcleanfq "$outdir/$samp/$laneIndex/raw/$laneIndex\_1.clean.fq.gz\n";
			print RMcleanfq "$outdir/$samp/$laneIndex/raw/$laneIndex\_2.clean.fq.gz\n";
                        }
                }
		print RMedicobam "$outdir/$samp/edico/result/$samp.bam\n";
	}
	close RMcleanfq;
	close RMedicobam;

}

## echo sign for monitor  
sub echostring{
    my $sh=shift;
    my $ostr="echo ==========end at : `date` ========== && \\\n";
    $ostr.="echo Still_waters_run_deep 1>&2 && \\\n";
    $ostr.="echo Still_waters_run_deep > $sh.sign\n";
    return $ostr;
}
