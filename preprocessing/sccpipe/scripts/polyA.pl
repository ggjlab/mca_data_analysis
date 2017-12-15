#!/usr/bin/perl
#Yincong Zhou
#2017-03-12

use strict;
use warnings;
use Switch;
use vars qw/$opt_b $opt_i $opt_o $opt_h/;

#use re::engine::TRE  max_cost =>1;
my $Mode=$ARGV[4];
my $Alength=$ARGV[3];
my $fastout=$ARGV[1];
my $fastbad=$ARGV[2];
my $line1;
my $line2;
my $line3;
my $line4;
my $count=0;
my $polyA = "A"x$Alength;
open(my $inputfile,$ARGV[0]) || die;
open(OUTPUT,">".$fastout);
open(BADOUT,">".$fastbad);

switch($Mode){
	case 1 {Mode_1($inputfile);}
	case 2 {Mode_2($inputfile);}
	case 3 {Mode_3($inputfile);}
	case 4 {Mode_4($inputfile);}
	else { my $outhash=polyAhash();
		foreach my $key (keys %$outhash){
			print $key."\n";
		}
	}
}

close($inputfile);
close(OUTPUT);
close(BADOUT);


sub Mode_1{
	my $input = shift;
	while($line1=<$input>){
		$line2=<$input>;
  		$line3=<$input>;
  		$line4=<$input>;
  		if($line2 =~ /($polyA)/g){
  			if(pos($line2) == $Alength){
				print BADOUT  join("",$line1,$line2,$line3,$line4);
				next;
	  		}
	 		else{
 		 		$line2 = substr($line2,0,pos($line2)-$Alength);
 		 		$line4 = substr($line4,0,length($line2))."\n";	
				print OUTPUT  join("",$line1,$line2."\n",$line3,$line4);
		     	}	
		}
  		else {
			print  OUTPUT join("",$line1,$line2,$line3,$line4);
		}


	}
}

sub Mode_2{
#	use re::engine::TRE  max_cost =>1;
	my $polyhash=polyAhash();
	my $input = shift;
	while($line1=<$input>){
		$line2=<$input>;
		$line3=<$input>;
		$line4=<$input>;
		chomp($line2);
		#	if($line2 =~ /\($polyA\)/i){
			#if(pos($line2) == $Alength){
			#	next;					  		
			#}
			#else{
			#	$line2 = substr($line2,0,pos($line2)-$Alength)."\n";	
			#print OUTPUT  join("",$line1,$line2,$line3,$line4);			
			OUTERLOOP:for(my $i=0;$i<length($line2)-$Alength;$i++){
				my $kmer=substr($line2,$i,$Alength);
				{
					if(exists $polyhash->{$kmer}){
						if($i<=4){
							print BADOUT join("",$line1,$line2."\n",$line3,$line4);
							last OUTERLOOP;
							#print substr($line2,$i,$Alength)."\n";
						}
						else {
							my $outfind=substr($line2,0,$i)."\n";
							$line4 = substr($line4,0,$i)."\n";
							print OUTPUT join("",$line1,$outfind,$line3,$line4);
							last OUTERLOOP;
						}
						
					}
					else{
						#print OUTPUT join("",$line1,$line2."\n",$Line4,$line4);
						if($i==(length($line2)-$Alength-1)){
							print OUTPUT join("",$line1,$line2."\n",$line3,$line4);
						}

					}
				}
			#print $kmer."\n";
			}
			#   }			
			#}
		#	else {
			#	print  OUTPUT join("",$line1,$line2,$line3,$line4);												
			#	}										
	}
}

sub Mode_3 {
	my $input = shift;
	my $start;
	while($line1=<$input>){
		$line2=<$input>;
		$line3=<$input>;
		$line4=<$input>;
		chomp($line2);
		my $tailA=substr($line2,length($line2)-$Alength,$Alength);
		if($tailA eq $polyA){
			
			for(my $i=(length($line2)-$Alength-1);$i>=0;$i--){
				if(substr($line2,$i,1) eq "A"){
					if($i==0){
						print BADOUT join("",$line1,$line2."\n",$line3,$line4);
						last;
					}
					else {
						next;
					}
				
				}
				else {
					if($i<=4){
						print BADOUT join("",$line1,$line2."\n",$line3,$line4);
						last;
					}
					else{
					print OUTPUT join("",$line1,substr($line2,0,length($line2)-$i)."\n",$line3,substr($line4,0,length($line2)-$i)."\n");
					last;

				}

			}
		}
		}
		else {
			print OUTPUT join("",$line1,$line2."\n",$line3,$line4);
		}
}
}
sub Mode_4 {
	print "mode_4\n";
}


sub kmercomform {
	my $strings = shift;

}

sub polyAhash {
	my %hash=();
	for(my $i=0;$i<$Alength;$i++){
		foreach my $base (qw/A T C G N/){
			my $new_poly = $polyA;
			substr($new_poly,$i,1)=$base;
			#	print $new_poly."\n";
			$hash{$new_poly}="1";
		}
	}
	return(\%hash);
}
