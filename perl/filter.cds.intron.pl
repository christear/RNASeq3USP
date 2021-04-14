#!/usr/bin/perl -w
use strict;
=use
filter those 3' UTR splicing events overlapped with intron linked two coding exons ...
...
=cut

sub filterEvents {
	my ($events,$gtf,$output) = @_;
	my %CDS = ();
	my %UTR = ();
	open IN,"$gtf" or die $!;
	my $txnid = ""; my $cdsid = "CDSNID";my $gn = "";my $gid = ""; 
	while(<IN>){
		chomp;
		next if(/^\#/);
		my @Columns = split("\t");
		my @Attr = split(";",$Columns[8]);
		my %Inf = ();
		foreach my $ele (@Attr){
			$ele =~ s/^\s+//g;
			my ($k,$v) = split(/\s+/,$ele);
			$k =~ s/\"//g;
			$v =~ s/\"//g;
			#print "$k\t$v";
			$Inf{$k} = $v;
		}
		$gid = $Inf{transcript_id} if(exists $Inf{gene_id});
		$txnid = $Inf{transcript_id} if(exists $Inf{transcript_id});
		$gn = $Inf{gene_name} if(exists $Inf{gene_name});
		if($Columns[2] eq "CDS"){
			$cdsid = $Inf{ccdsid} if(exists $Inf{ccdsid});
			my $left = $Columns[4] + 1;my $right = $Columns[3] - 1;
			$CDS{$Columns[0]}{$left} = $cdsid.":".$txnid.":".$gn;
			$CDS{$Columns[0]}{$right} = $cdsid.":".$txnid.":".$gn;
#			print "$txnid\t$gn\t$cdsid\n";
		}elsif($Columns[2] eq "UTR"){
			push @{$UTR{$txnid}{coord}},join(":",@Columns);
			$UTR{$txnid}{gene} = $gn;
		}
	}
	close IN;
	my %Introns = ();
	foreach my $txnid(keys %UTR){
		my @Coords = @{$UTR{$txnid}{coord}};
		my @C1 = split(":",$Coords[0]);
		my $gn = $UTR{$txnid}{gene};
		next if($#Coords < 1);
		if($C1[6] eq "+"){
			my $left = $C1[4] + 1;my $right = 0;
			for(1..$#Coords){
				my @Ec = split(":",$Coords[$_]);
				$right = $Ec[3] - 1;
	#			print "$txnid\t$_\t$Ids[$_]\t$gn\n";
				$Introns{$C1[0]}{$left}{$right} = "UTR:UTR:".$txnid.":".$gn;
				$left = $Ec[4] + 1;
			}
		}else{
			my $right = $C1[3] - 1;my $left = 0;
			for(1..$#Coords){
				my @Ec = split(":",$Coords[$_]);
				$left = $Ec[4] + 1;
				$Introns{$C1[0]}{$left}{$right} = "UTR:UTR:".$txnid.":".$gn;
				$right = $Ec[3] - 1;
			}
		}
	}
	open IN,"$events" or die $!;
	open OUT,">$output" or die $!;
	while(<IN>){
		chomp;
		my @Columns = split("\t");
		my ($chr,$l,$r) = split(/[:\-]/,$Columns[0]);
		$chr = "chr".$chr if ($chr !~ /chr/);
		$chr =~ s/MT/M/;
		my $id = ""; my $lid = "UNA1";my $rid = "UNA2";
		$lid = $CDS{$chr}{$l} if(exists $CDS{$chr}{$l});
		$rid = $CDS{$chr}{$r} if(exists $CDS{$chr}{$r});
		$id = $lid.":".$rid;
		$id = $Introns{$chr}{$l}{$r} if((exists $Introns{$chr}{$l}{$r}) and (($lid eq "UNA1") and ($rid eq "UNA2")));
		print OUT "$_\t$id\n";
	}
	close IN;
	close OUT;
}

if(@ARGV < 3){
	print "Usage:events_list\tgtf\toutput\n";
	exit(0);
}else{
	&filterEvents($ARGV[0],$ARGV[1],$ARGV[2]);
}