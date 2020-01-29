#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use File::Slurp;
use Statistics::R;
use List::MoreUtils qw(uniq);
use Cwd;

my ($dataFile) = ("");

my $var_pat = getcwd();
GetOptions(	"in=s"	=> \$dataFile);

my $inPath = <$var_pat/$dataFile>;
open my $inFile, '<', "$inPath" or die "Can't open $inPath: $! ";
my @fileArray = <$inFile>;


my @allcdr3aa;
foreach my $line (@fileArray){
		chomp $line;
		my @lineArray = split /\t/,$line;
		my $cdr3aa = $lineArray[0];
		if($cdr3aa =~ /\*/){}
		elsif($cdr3aa =~ /\_/){}
		elsif($cdr3aa =~ /^CA/ ){push @allcdr3aa, $cdr3aa."\n";}
}

my @unique_cdr3aa = uniq @allcdr3aa;

my $R = Statistics::R -> new();
$R -> run("library(ggplot2)");
$R -> run("myData <- read.table('//data/evilla/cdr3_ggplot/vaccines/vaccines_table.txt', header=TRUE)");
my $plotNum=1;
my ($a,$b,$c,$d)=("") x 4;
for (my $i=0;$i<24;$i++){
		chomp $unique_cdr3aa[$i];
		my $a = $unique_cdr3aa[$i];
		if($unique_cdr3aa[$i+1]){chomp $unique_cdr3aa[$i+1];	$b = $unique_cdr3aa[$i+1];}
		if($unique_cdr3aa[$i+2]){chomp $unique_cdr3aa[$i+2];	$c = $unique_cdr3aa[$i+2];}
		if($unique_cdr3aa[$i+3]){chomp $unique_cdr3aa[$i+3];	$d = $unique_cdr3aa[$i+3];}

		$R -> run("subData <- subset(myData,cdr3aa=='".$a."'|cdr3aa=='".$b."'|cdr3aa=='".$c."'|cdr3aa=='".$d."')");
		$R -> run("subDay0 <- data.frame(Day=rep('0',nrow(subData)))");
		$R -> run("subDay0 <- cbind(subData\$cdr3aa,subDay0,subData\$Day_0)");
		$R -> run("colnames(subDay0)[1]<- 'cdr3aa'");
		$R -> run("colnames(subDay0)[3]<- 'Value'");
		$R -> run("subDay7 <- data.frame(Day=rep('7',nrow(subData)))");
		$R -> run("subDay7 <- cbind(subData\$cdr3aa,subDay7,subData\$Day_7)");
		$R -> run("colnames(subDay7)[1]<- 'cdr3aa'");
		$R -> run("colnames(subDay7)[3]<- 'Value'");
		$R -> run("subDay28 <- data.frame(Day=rep('28',nrow(subData)))");
		$R -> run("subDay28 <- cbind(subData\$cdr3aa,subDay28,subData\$Day_28)");
		$R -> run("colnames(subDay28)[1]<- 'cdr3aa'");
		$R -> run("colnames(subDay28)[3]<- 'Value'");
		$R -> run("dataCombined <- rbind(subDay0,subDay7,subDay28)");
		$R -> run("ggplot(dataCombined,aes(x=Day,y=Value))+geom_jitter(width=0.3,aes(colour=Day))");
		$i=$i+4;
		$plotNum++;
}