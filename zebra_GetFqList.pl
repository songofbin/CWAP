#!/usr/bin/perl

=head1 Name

	zebra_GetFqList.pl -- 	Get fastq list of zebra for analysis pipeline.

=head1 Description

	This script aims at finding relationship of these three files and 
	printing input file for WES or RNAseq pipeline.
	sample - library     -                   index  (pool.info)
	            ||                            ||
            library(+/-\d)   -  chip-lane               (seq.info)
	                           ||             ||  
	                        chip-lane   -    index  (fqfile.list)

=head1 Version

	Author: Song Bin, songbin@genomics.cn
	Version: 1.0	Date: 2016-9-20
	Version: 1.1	Date: 2017-2-7
	Version: 1.2	Date: 2017-2-20
	Version: 1.3	Date: 2017-3-5 
	changing key of hash %lib2lane from $rawlib-$library to $rawlib-$id.
	Version: 1.4    Date: 2017-5-4
	Make sure that fqfile.list contains all samples/librarys' chip in pool.info and seq.info file.
	Version: 1.5	Date: 2017-8-8
	Fix this line: if($library =~ /(.*-N)[+-]\d/){	

=head1 Usage
	
	perl $0 <pool.info with header> <seq.info with header> <fqfile.list> <SeqType:new/old> <output>

	Result examples:
	for Zebra RNA-seq SE50+10 pipeline input  (new)
	BC008-A1	RHUMfddNAAAARAAZebra-N	CL100006572_L01_6	/lfs1/zebra03/F14ZF1QSSY1656_Temp/CL100006572_L01/CL100006572_L01_6.fq.gz	50	200	1580787850

	for Zebra WES PE50+10 pipeline input      (new)
	BC008-N	RHUMbmaXDAAZRAAZebra-N	CL100006996_L01_17	/szhwfs1/zebra_bgiseq500/zebra01/F14ZF1QSSY1656_Temp/CL100006996_L01/clean_data/CL100006996_L01_17_1.fq.gz.clean.gz,/szhwfs1/zebra_bgiseq500/zebra01/F14ZF1QSSY1656_Temp/CL100006996_L01/clean_data/CL100006996_L01_17_2.fq.gz.clean.gz	50;50	150	4890829790

	for Zebra WES PE50 pipeline input          (new)
	BC044-A	RHUMbmaXJALARAAZebra-N	CL100010118_L02	/szhwfs1/zebra_bgiseq500/zebra01/F14ZF1QSSY1656_Temp/CL100010118_L02/CL100010118_L02_read_1.fq.gz,/szhwfs1/zebra_bgiseq500/zebra01/F14ZF1QSSY1656_Temp/CL100010118_L02/CL100010118_L02_read_2.fq.gz	50;50	150	68002894800

	for Zebra WES PE50+10 pipeline input      (old)
	BC018-N	33	zebra	/lfs1/zebra03/F14ZF1QSSY1656_Temp/CL300002657_L02/CL300002657_L02_45	library
=cut

use strict;
use warnings;
use File::Basename qw/dirname basename/;
use Data::Dumper;
die "perl $0 <pool.info> <seq.info> <fqfile.list> <SeqType:new/old> <output>\n" if(@ARGV != 5);

my ($poolfile,$seqfile,$fqlist,$seqtype,$output)=@ARGV;
die "SeqType should be 'new' or 'old' !" if($seqtype ne "new" and $seqtype ne "old");

my $report_backup_dir="/ifs3/yjy_data/Zebra_HTML_report_backup";
my (%pool,%libCount,%lib2lane,%lane2lib,%fqfiles,%fqStats);

open IN1,"$seqfile" || die $!;
<IN1>;
while(<IN1>){
	chomp;
	my ($id,$chip,$lane,$library,$seqType,$startTime,$endTime,$info)=(split /\t/)[0,1,2,3,4,6,8,9];
	if(defined $info && $info =~ /作废/){
		next;
	}
	my $rawlib=$library;
	
	# Library sequenced for several times, eg. RHUMbmaXFACTRAAZebra-N+3
	if($library =~ /(.*-N)[+-]\d/){
		$rawlib=$1;
	}

	if(exists $lib2lane{$rawlib}){
		$libCount{$rawlib}++;
	}else{
		$libCount{$rawlib}=1;
	}
	
	$lib2lane{$rawlib}{$libCount{$rawlib}}="$chip\_$lane";
	$lane2lib{"$chip\_$lane"}=$library;
}
close IN1;

#print Dumper(\%lib2lane);

open IN2,"$poolfile" || die $!;
<IN2>;
while(<IN2>){
	chomp;
	my @F = split /\t/;
	my ($taskID,$sample,$sampleID,$poolNum,$indexRange,$rawlib)=@F[0,1,2,3,4,5];
	#Fix sample name if needed.
	#$sample=~s/(.*)-(\d+)-(.*)A/$1-$3/ if $sample=~/(.*)-(\d+)-(.*)/;
	
	if(!exists $lib2lane{$rawlib}){
		die "Can't find the chip_lane of library $rawlib !\n";
	}

	foreach my $id (sort keys %{$lib2lane{$rawlib}}){
		my $chipLaneID=$lib2lane{$rawlib}{$id};	
		#Index could be "41-44" or "41"
		my @indexs;
		if($indexRange=~/-/){
			my ($from,$to)=(split /-/,$indexRange)[0,1];
			for my $index ($from..$to){
				push @indexs,$index;
			}
		}
		else{push @indexs,$indexRange;}
		for my $index (@indexs){
			my $indexID="$chipLaneID\_$index";
			$pool{$indexID}=$sample;
		}
		
		my $indexID1="$chipLaneID\_read";
		$pool{$indexID1}=$sample;
	}
}
close IN2;

open IN3,"$fqlist" || die $!;
while(<IN3>){
	chomp;
	my ($id,$file)=(split /\t/)[0,1];
	
	# id could be "CL100004346_L01_41_1","CL100004346_L01_41_2"(PE with Index) 
	# or "CL100011896_L01_read_1","CL100011896_L01_read_2"(PE without Index)
	# or "CL100006572_L01_10"(SE with Index)
	my @ids=split /_/,$id;
	my ($chip,$lane,$index)=@ids[0,1,2];
	my $fqid=$#ids==3?$ids[3]:1;
	
	my $indexID    = "$chip\_$lane\_$index";
	my $sample;
	if(exists $pool{$indexID}){
		$sample=$pool{$indexID};
		$fqfiles{$sample}{$indexID}{$fqid}=$file;
	}
	else{
		warn "Can't find sample name of the Chip_lane_index $indexID !\n";
		next;
	}
	
	next if ($seqtype eq "old");
	my $dir  = dirname($file);
	my $name = basename($file);
	$name =~ s/(.*).fq.gz$/$1/;

	#get fqStat file
	my $fqStatFile = "$dir/$name.fq.fqStat.txt";
	if(!-e $fqStatFile){
		$fqStatFile = "$dir/$chip\_$lane\_Report/FqReport/$name.fq.fqStat.txt";
		if(!-e $fqStatFile){
			$fqStatFile = "$report_backup_dir/$chip\_$lane/FqReport/$name.fq.fqStat.txt";
		}
	}

	if(-e $fqStatFile){
		my $info = `cat $fqStatFile | head -13`;
		map {	my ($key,$value) = (split /\t/,$_)[0,1];
			$key =~ s/^#(.*)/$1/;
			$fqStats{$sample}{$indexID}{$fqid}{$key} = $value
		} split (/\n/,$info);
	}
	else{
		warn "Can't find fqStat file of $name !\n";
	}
}
close IN3;

open OUT,">$output" || die $!;
foreach my $sample (sort keys %fqfiles){
	foreach my $indexID (sort keys %{$fqfiles{$sample}}){
		my ($chip,$lane,$index)=(split /_/,$indexID)[0,1,2]; 
		my $library    = $lane2lib{"$chip\_$lane"};
		my $fqfile     = $fqfiles{$sample}{$indexID}{1};
		my $lanePath   = $fqfile;
		$lanePath=~s/(.*).fq.gz$/$1/;
		$lanePath=~s/(.*)_read_1$/$1/;

		if($seqtype eq "old"){
			print OUT "$sample\t33\tzebra\t$lanePath\t$library\n";
		}
		elsif($seqtype eq "new"){
			my $fqPrefix=(split /\//,$lanePath)[-1];
			if(exists $fqfiles{$sample}{$indexID}{2}){
				$fqfile.=",$fqfiles{$sample}{$indexID}{2}";
			}
			my $insertSize = "150";
			my $readLength = "NA";
			my $BaseNum    = "NA";
			if(exists $fqStats{$sample}{$indexID}{1}{"BaseNum"}){
				$BaseNum=$fqStats{$sample}{$indexID}{1}{"BaseNum"};
				$readLength = $fqStats{$sample}{$indexID}{1}{"row_readLen"};
				if(exists $fqStats{$sample}{$indexID}{2}{"BaseNum"}){
					$BaseNum+=$fqStats{$sample}{$indexID}{2}{"BaseNum"};
					$readLength.=";$fqStats{$sample}{$indexID}{2}{row_readLen}";
				}
			}
			print OUT "$sample\t$library\t$fqPrefix\t$fqfile\t$readLength\t$insertSize\t$BaseNum\n";
		}
	}
}

