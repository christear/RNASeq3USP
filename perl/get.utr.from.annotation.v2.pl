#!/usr/bin/perl -w
use strict;
=use
get the utr from annotation ...
new version, based on txn start/end and start codon, end codon ...
=cut

sub getUTRCDS {
	my ($method,$gtf,$output) = @_;
	my %Txns = ();
	my %CDS = ();
	open IN,"$gtf" or die $!;
	while(<IN>){
		chomp;
		next if(/^\#/);
		my @Columns = split("\t");
		my @Arra = split(";",$Columns[8]);
		my %Inf = ();
		foreach my $ele (@Arra){
			$ele =~ s/^\s//g;
			my ($k,$v) = split(/[\s\=]/,$ele);
			$k =~ s/\"//g;
			$v =~ s/\"//g;
			$Inf{$k} = $v;
		}
		if(exists $Inf{transcript_id}){
			my $txnid = $Inf{transcript_id};
			if($Columns[2] eq "transcript"){
				$Txns{$txnid}{info} = join("\t",@Columns);
#				$Txns{$txnid}{start} = $Columns[3];
#				$Txns{$txnid}{end} = $Columns[4];
			}
			if(($Columns[2] eq "start_codon") or ($Columns[2] eq "stop_codon")){
				$CDS{$txnid}{start} = $Columns[3] if($Columns[2] eq "start_codon");
				$CDS{$txnid}{end} = $Columns[4] if($Columns[2] eq "stop_codon");
				$Txns{$txnid}{$Columns[2]} = $Columns[3]."-".$Columns[4];
			}
		}
	}
	close IN;
	open OUT,">$output" or die $!;
	foreach my $txnid (sort keys %Txns){
		if($method eq "UTR"){
			my @Columns = split("\t",${Txns{$txnid}{info}});
			my $eo = "";
			if(exists $Txns{$txnid}{stop_codon}){
				my ($sp1,$sp2) = split("-",$Txns{$txnid}{stop_codon}); 
				my @Utr3 = @Columns;
				$Utr3[2] = "UTR3";
				$Utr3[3] = $sp1 if($Utr3[6] eq "+");
				$Utr3[4] = $sp2 if($Utr3[6] eq "-");
				$eo = join("\t",@Utr3);
				print OUT "$eo\n";
			}else{
				print "transcript $txnid does not have stop codon information\n";
			}
			if(exists $Txns{$txnid}{start_codon}){
				my ($st1,$st2) = split("-",$Txns{$txnid}{start_codon}); 
				my @Utr5 = @Columns;
				$Utr5[2] = "UTR5";
				$Utr5[4] = $st2 if($Utr5[6] eq "+");
				$Utr5[3] = $st1 if($Utr5[6] eq "-");
				$eo = join("\t",@Utr5);
				print OUT "$eo\n";
			}else{
				print "transcript $txnid does not have start codon information\n";
			}
		}else{
			if((exists $CDS{$txnid}{start}) and (exists $CDS{$txnid}{end})){
				my @Columns = split("\t",${Txns{$txnid}{info}});
				$Columns[2] = "CDS";
				$Columns[3] = ($Columns[6] eq "+") ? $CDS{$txnid}{start} : $CDS{$txnid}{end};
				$Columns[4] = ($Columns[6] eq "+") ? $CDS{$txnid}{end} : $CDS{$txnid}{start};
				my $l = join("\t",@Columns);
				print OUT "$l\n";
			}
		}	
	}
	close OUT;
}

if(@ARGV < 3){
	print "Usage:method\tGTF_input\toutput\nmethod should be UTR/CDS\n";
	exit(0);
}else{
	if($ARGV[0] =~ /[Uu][Tt][Rr]/){
		&getUTRCDS("UTR",$ARGV[1],$ARGV[2]);
	}else{
		&getUTRCDS("CDS",$ARGV[1],$ARGV[2]);
	}
	
}