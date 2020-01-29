#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use File::Slurp;
use List::MoreUtils qw(uniq);
use Cwd;

my $data_path = "";
my $var_path = getcwd();
require "/".$var_path."/_scripts/basic_scripts.pl";
GetOptions("in|i=s" => \$data_path);

my $outpath = 
my @filesArray = <$data_path/*txt>;
my @vdjFile;

foreach my $fileIn(@filesArray){
		my $fileName = removeFileType(returnFile($fileIn));
		my @fileArray = read_file($fileIn);
		shift @fileArray;
		my @finalArray;
		my $totalCount = 0;
		foreach my $line(@fileArray){
				chomp $line;
				my ($v,$nt,$aa,$d,$j,$length,$fromV,$fromJ,$dlength,$add,$count) = split /\t/,$line;
				$totalCount = ($totalCount+$count);
		}
		my $header = "count\tfrequency\tCDR3nt\tCDR3aa\tV\tD\tJ\n";
		push @finalArray, $header;
		foreach my $line(@fileArray){
				$line =~ s/\r//;
				chomp $line;
				my ($v,$nt,$aa,$d,$j,$length,$fromV,$fromJ,$dlength,$add,$count) = split /\t/,$line;
				my $freq = ($count/$totalCount);
				chomp $count;
				my $vdjToolsLine = $count."\t".$freq."\t".$nt."\t".$aa."\t".$v."\t".$d."\t".$j."\n";
				push @finalArray, $vdjToolsLine;
				
		}
write_file($var_pat."/".$fileName.".vdjtools", @finalArray);
}
