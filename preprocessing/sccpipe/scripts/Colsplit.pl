#!/usr/bin/perl -w

#split the raw file to 20 cols
#2017-03-08
#modify from  mRNASeqHT_demultiplex
#Any infringement please contact yczhou@zju.edu.cn


use strict;

use vars qw/$opt_b $opt_i $opt_o $opt_h/;
use IO::File;
my $blengthin=$ARGV[4];
my $mismatch=$ARGV[3];
my $input = $ARGV[2];
my $indexcol = index_data($input);
my $indexlen;
my $indexhash=index_hash($indexcol);


my $fastq_data = get_fastq_data($ARGV[0]);
my $out_dir = $ARGV[1];
my %summary = ();
system("mkdir $out_dir") if not -e $out_dir;
foreach my $key (keys %$indexcol){
	print $key."\n";
}

foreach my $sample(keys %$fastq_data){
		print  "process fastq data for sample $sample\n";
		my $r1_fastq = $fastq_data->{$sample}->{R1};
		my $r2_fastq = $fastq_data->{$sample}->{R2};
		print $r1_fastq, "\n";
		print $r2_fastq, "\n";
		my $r1_in = IO::File->new("$r1_fastq") or die $!;
		my $r2_in = IO::File->new("$r2_fastq") or die $!;

		my %out_handles = ();
		my ($sample_custom_name) = split("_", $sample);	
		print $sample_custom_name."\n";			
		foreach my $index(keys %$indexcol){
			my $id = $indexcol->{$index};

			my $cell_sample_name = join("_", $id, $sample_custom_name);
			my $out_sample_dir=join("/", $out_dir, $id);
			system("mkdir $out_sample_dir") if not -e $out_sample_dir;
			my $r1_sample_file = join("/", $out_dir, $id,join("_", $cell_sample_name, "R1.fastq"));
			my $r2_sample_file = join("/", $out_dir, $id,join("_", $cell_sample_name, "R2.fastq"));
			# print $r1_sample_file."\n";

			# my $r1_sample_file = join("/", $out_dir, join("_", $cell_sample_name, "R1.fastq"));
			# my $r2_sample_file = join("/", $out_dir, ,join("_", $cell_sample_name, "R2.fastq"));
			my $r1_sample_out = IO::File->new(">$r1_sample_file") or die $!;
			my $r2_sample_out = IO::File->new(">$r2_sample_file") or die $!;
			
			$out_handles{$cell_sample_name}{R1} = $r1_sample_out;
			$out_handles{$cell_sample_name}{R2} = $r2_sample_out;
			
		}
	
		my $r1_index_out;
		my $r2_index_out;
		while (my $temp1_1 = <$r1_in>){
		
			my $temp1_2 = <$r1_in>;
    		my $temp1_3 = <$r1_in>;
    		my $temp1_4 = <$r1_in>;
    		
    		my $temp2_1 = <$r2_in>;
    		my $temp2_2 = <$r2_in>;
    		my $temp2_3 = <$r2_in>;
    		my $temp2_4 = <$r2_in>;
    	
    		my $bc = substr($temp1_1,-($indexlen+1),$indexlen);
    		
    		if (exists $indexhash->{$bc}){
    			my $id = $indexhash->{$bc};
    			my $cell_name = join("_", $id, $sample_custom_name);
    			
    			# print $cell_name."\n";
    			$r1_index_out  = $out_handles{$cell_name}{R1};
    			$r2_index_out  = $out_handles{$cell_name}{R2};

    			
    			print $r1_index_out  join("", $temp1_1, $temp1_2,  $temp1_3, $temp1_4);	
    			print $r2_index_out  join("", $temp2_1, $temp2_2, $temp2_3, $temp2_4);
    			$summary{$cell_name}++;
    		}
    		else{
    			$summary{Undetermined}++;
    		}
	
    		
    	}
}


sub index_hash{
	my ($indexs) = @_;
    my %revised_indexs = ();
	foreach my $index(keys %$indexs){
		
		my $id = $indexs->{$index};
		if ($mismatch==0){
		$revised_indexs{$index} = $id;
}
		else {
		for(my $i =0; $i < length($index); $i++){

			foreach my $base (qw/A T G C N/){
				my $new_index = $index;

				substr($new_index, $i, 1) = $base;
				$revised_indexs{$new_index} = $id;
		}
		
	}
}
	$indexlen=length($index);

	}

	print $indexlen;
	
	return(\%revised_indexs);
}

sub index_data {
	my ($input) = @_;
	my $data = get_tab_delimited_data_from_file($input);
	my %indexs = ();
	
	foreach my $row(@$data){
		my ($index, $id) = @$row;
		$index=substr($index,0,$blengthin);
		$indexs{$index} = $id;
	}
	
	return(\%indexs);
}

sub get_tab_delimited_data_from_file {
	my ($input) = @_;

	my $in = IO::File->new($input) or die $!;
	my @data = ();
	
	while(<$in>){
		$_ =~ s/\s+$//;
		next if not $_;
		my @line = split("\t", $_);

		push @data, \@line;
		
	}

	return \@data;
}	

sub get_fastq_data {
	my ($dir) = @_;
    my %fastq = ();
    opendir(my $dh, $dir) || die;
    while(my $file = readdir $dh) {
       next if $file !~ /fastq/;
       	my $command;
       
       if($file =~ /fastq\.gz/){
       	   my $file_str = join("/", $dir, $file);
       	  $command = "gunzip $file_str";
       	  `$command`;
       	  $file =~ s/\.gz//;
       }
       elsif($file =~ /fastq\.zip/ ){
       	 my $file_str = join("/", $dir, $file);
          $command = "unzip $file_str";
       	  `$command`;
       	  $file =~ s/\.zip//;
       }	
       
       my $sample_name;
       my $strand;
       
       if($file =~ /^(\S+)\_L001\_(R\d)/){
   	 		$sample_name = $1;
   	 		$strand = $2;
   		}
   		elsif($file =~ /^(\S+)\_(R\d)/){
                $sample_name = $1;
                $strand = $2;
   		}
   		elsif($file =~ /^(\S+)\_(\d)\.fastq/){
                $sample_name = $1;
                $strand = $2;
   		}
   		else{
   			die "fastq file name is not supported\n ";
   		}
   		
   		$fastq{$sample_name}{$strand} = join("/", $dir, $file);
   		
    }
    closedir $dh;
   	return(\%fastq);
}

