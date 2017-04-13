#!/usr/bin/perl

while(<>){

	$line = $_;

	if($line =~ /^#CHROM/){
		chomp($line);

		@samples = split(/\t/,$line);

		for($i= 9; $i < @samples; $i++){
			$samples[$i] =~ s/\s/_/g;
			print "$i\t$samples[$i]\t0\t0\t0\t1\n";

		}


		last;
	}



}
