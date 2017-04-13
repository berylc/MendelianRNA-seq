#!/usr/bin/perl


%purcellSite = ();

open(SITES,"<../data/purcell5k.intervals");
while(<SITES>){
	$line = $_;
	chomp $line;
	$purcellSite{$line} = 1;
}
close(SITES);



#$vcfFile = $ARGV[0];
#open(VCF,"<$vcfFile");

while(<>){
	
	$line = $_;
	chomp $line;

	if($line !~ /^#/){
		@data = split(/\t/,$line);		
		$ref = $data[3];
		$alt = $data[4];		

		@alt_alleles = split(/,/,$alt);


		$indel = 0;
		$snp = 0;
		$refLen = length($ref);

		foreach $a (@alt_alleles){
			if(length($a) == $refLen){
				$snp = 1;				
			}
			else{
				$indel = 1;
			}
		}

		$alt_n = 1;

		if($alt =~ /,/){
			$data[7] =~ /AC=(.*?);/;
			@acs = split(/,/,$1);
			$alt_ac = 0;
			$alt_n = 1;			

			for($i= 0; $i < @acs; $i++){
				if($acs[$i] > $alt_ac){
					$alt_ac = $acs[$i];
					$alt = $alt_alleles[$i];
					$alt_n = $i+1;
				}
			}
		
		}

		#SNPs ONLY. Ignores Indels that overlap these sites and SNPs that overlap with Indels.
		if($snp == 1 && $indel == 0 && defined($purcellSite{"$data[0]:$data[1]"})){
			
			if($data[2] eq "."){
				$varID = sprintf("var_%d_%d",$data[0],$data[1])
			}
			else{
				$varID = $data[2];
			}

			printf "%d\t%s\t%.6f\t%d", $data[0],$varID,$data[1]/1000000,$data[1];			
			

			for($i= 9; $i < @data; $i++){
				@gt = split(/:/,$data[$i]);
				$geno = $gt[0];

				if($geno eq "0/0"){
					print " $ref $ref";
				}
				elsif($geno eq "0/$alt_n"){
					print " $alt $ref";				
				}
				elsif($geno eq "$alt_n/$alt_n"){
					print " $alt $alt";				
				}
				#Very rare alt will print here! Need to go back and check the allele order.
				else{
					print " 0 0";				
				}
			}
			print "\n";

		}

	}

}

close(VCF);
