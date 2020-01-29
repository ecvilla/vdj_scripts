#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use File::Slurp;
use Statistics::R;
use List::MoreUtils qw(uniq);

use Cwd;
my $var_path = getcwd();
require"/".$var_path."/_scripts/basic_script.pl";
my @outputFiles = <$var_pat/_reports/*.all.aln>;
my @oligoFiles;
foreach my $file (@outputFiles){
		print $file."\n";
		push @oligoFiles, $file;
}
my @txt_file;
push @txt_file, "FILE_NUM\tOLIGO_NUM\tOLIGO\tPOINT\tVIRUS\tEXP_NUM\tEXPRESSION\tFILE_NAME\n";
my $fileNum = 0;
my @namesArray;

foreach my $oligoFile (@oligoFiles){
		my ($name) = removeFileType(removeFileType(bk_basic_returnFile($oligoFile)));
		push @namesArray, $name;
		$fileNum++;
		my @clusterOligos = read_file($oligoFile);
		@clusterOligos = grep {$_ !~ /00000/} @clusterOligos;
		@clusterOligos = grep {$_ !~ /OLIGO/} @clusterOligos;
		@clusterOligos = grep {$_ !~ /Consensus/} @clusterOligos;
		@clusterOligos = grep {$_ =~ /\-GP\-/ or /\-GP\_/} @clusterOligos;

		my @clusterDataOnly;
		foreach my $clus (@clusterOligos){
				my @tempClus = split "\t", $clus;
				my @oligoDetails = split "_", $tempClus[0];
				my $length = @oligoDetails;
				push @clusterDataOnly, [$tempClus[0]."\t".$tempClus[1], $oligoDetails[$length-5]];
		}

		my @uniqueClusOligos = uniqueArray(@clusterDataOnly);
		my @sortUniqueClusOligos = map{$_->[0]}sort{ $a->[1] <=> $b->[1]} @uniqueClusOligos;
		@sortUniqueClusOligos = reverse @sortUniqueClusOligos;

		my @array;
		my @EBOV;
		my @SUDV;
		my @other;
		foreach my $clus(@sortUniqueClusOligos){
				if($clus =~ /EBOV/){push @EBOV, $clus."\tEBOV";}
				elsif($clus =~ /SUDV/){push @SUDV, $clus."\tSUDV";}
				else{push @other, $clus."\tother";}
				@array = (@other,@SUDV,@EBOV);
		}


		my $oligoNum = 0;
		my @expressionArray;
		foreach my $cluster (@array){
				$oligoNum++;
				my ($oligo, $exp,$virus) = split "\t", $cluster;
				my @oligoDetails = split "_", $oligo;
				my $exp_class = "";
				my $length = @oligoDetails;
				if($exp>=1 and $exp<2){$exp_class = "1-2";push @expressionArray, "1-2";}
				elsif($exp>=2 and $exp<10){$exp_class = "2-9.9";push @expressionArray, "2-10";}
				elsif($exp>=10 and $exp<50){$exp_class = "9.9-50";push @expressionArray, "10-50";}
				elsif($exp>=50){$exp_class = ">50";push @expressionArray, ">50";}
				push @txt_file, $fileNum."\t".$oligoNum."\t".$oligo."\t".$oligoDetails[$length-5]."\t".$virus."\t".$exp."\t".$exp_class."\t".$name."\n";
				push @txt_file, $fileNum."\t".$oligoNum."\t".$oligo."\t".$oligoDetails[$length-4]."\t".$virus."\t".$exp."\t".$exp_class."\t".$name."\n";
		}
		push @expressionArray, ">50";


		my @virusArray;
		foreach my $clus(@array){
				my @checkAllVariables = split "\t", $clus;
				push @virusArray, $checkAllVariables[2];
		}

		my @uniquevirus = bk_basic_uniqueArray(@virusArray);
		my @uniqueExp = bk_basic_uniqueArray(@expressionArray);
		my ($ebovFound, $sudvFound, $otherFound, $GPFound) = (0) x 4;
		foreach my $vir(@uniquevirus){
				if($vir =~ /EBOV/){$ebovFound++;}
				elsif($vir =~ /SUDV/){$sudvFound++;}
				elsif($vir =~ /other/){$otherFound++;}
				elsif($vir =~ /GP/){$GPFound++;}
		}

		my ($oneFound, $twoFound, $threeFound, $fourFound) = (0) x 4;
		foreach my $exp(@uniqueExp){
				if($exp =~ /1-2/){$oneFound++;}
				elsif($exp =~ /2-9\.9/){$twoFound++;}
				elsif($exp =~ /9\.9-50/){$threeFound++;}
				elsif($exp =~ />50/){$fourFound++;}
		}

		if ($ebovFound ==0){$oligoNum++;push @txt_file, ($fileNum)."\t".($oligoNum)."\tmissing\t0\tEBOV\t1\t1-2\t".$name."\n";push @txt_file, ($fileNum)."\t".($oligoNum)."\tmissing\t0\tEBOV\t1\t1-2\t".$name."\n";}
		if ($sudvFound ==0){$oligoNum++;push @txt_file, ($fileNum)."\t".($oligoNum)."\tmissing\t0\tSUDV\t1\t1-2\t".$name."\n";push @txt_file, ($fileNum)."\t".($oligoNum)."\tmissing\t0\tSUDV\t1\t1-2\t".$name."\n";}
		if ($otherFound ==0){$oligoNum++;push @txt_file, ($fileNum)."\t".($oligoNum)."\tmissing\t0\tother\t1\t1-2\t".$name."\n";push @txt_file, ($fileNum)."\t".($oligoNum)."\tmissing\t0\tother\t1\t1-2\t".$name."\n";}

		if ($oneFound ==0){$oligoNum++;push @txt_file, ($fileNum)."\t".($oligoNum)."\tmissing\t0\tEBOV\t1\t1-2\t".$name."\n";push @txt_file, ($fileNum)."\t".($oligoNum)."\tmissing\t0\tEBOV\t1\t1-2\t".$name."\n";}
		if ($twoFound ==0){$oligoNum++;push @txt_file, ($fileNum)."\t".($oligoNum)."\tmissing\t0\tEBOV\t1\t2-9.9\t".$name."\n";push @txt_file, ($fileNum)."\t".($oligoNum)."\tmissing\t0\tEBOV\t1\t2-9.9\t".$name."\n";}
		if ($threeFound ==0){$oligoNum++;push @txt_file, ($fileNum)."\t".($oligoNum)."\tmissing\t0\tEBOV\t1\t9.9-50\t".$name."\n";push @txt_file, ($fileNum)."\t".($oligoNum)."\tmissing\t0\tEBOV\t1\t9.9-50\t".$name."\n";}

		$oligoNum++;
		push @txt_file, $fileNum."\t".($oligoNum)."\tGP\t0\tGP\t50\t>50\t".$name."\n";
		push @txt_file, $fileNum."\t".($oligoNum)."\tGP\t700\tGP\t50\t>50\t".$name."\n";

}
write_file($var_pat."/toPlot.txt", @txt_file);

my $R = Statistics::R -> new();
$R -> run("library(ggplot2)");
$R -> run("myData <- read.table('".$var_pat."/toPlot.txt', header=TRUE)");

for (my $i = 1; $i<=$fileNum;$i++){
		$R -> run("subData <- subset(myData,FILE_NUM=='".$i."')");
		$R -> run("ggplot(subData,aes(x=POINT,y=OLIGO_NUM,group=OLIGO_NUM))+xlab('Position')+ylab('')+ggtitle('".$namesArray[$i-1]."')+geom_line(aes(size = EXPRESSION, colour = VIRUS ))+scale_size_manual(values =c(0.1,0.5,1.1,1.9))+facet_grid(.~FILE_NUM)+ theme_bw() + theme(panel.border = element_blank(),strip.background = element_blank(),strip.text.x = element_blank())");
}
