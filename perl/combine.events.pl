#!/usr/bin/perl -w
=use
combine the splicing events ... here are the intron left and right gtf ... together ...
... 
=cut

sub combineEv {
	my @Gtfs = @_;
	my %Sout;
	my %Results = ();
	my %SoutSums = ();
	my $SpNums = ();
	my $id = ""; my $c = 0; my $i = 0;
	foreach my $gtffile (@Gtfs){
		open IN,"$gtffile" or die $!;
		while(<IN>){
			chomp;
			my @Columns = split("\t");
			my $line = $_;
			my @Attr = split(";",$Columns[8]);
			foreach my $ele (@Attr){
				my ($k,$v) = split(" ",$ele);
				$k =~ s/ //g;
				$v =~ s/ //g;
				$Info{$k} = $v;
			}
			$id = $Info{intron_id};
			if(!(exists $SpNums{$id})){
				$SpNums{$id} = 1;
				$SoutSums{$id} = $Info{sout_count};
			}else{
				my $n = $SpNums{$id};
				my $s = $SoutSums{$id};
				$n++;
				$s+= $Info{sout_count};
				$SpNums{$id} = $n;
				$SoutSums{$id} = $s;
			}
			$Columns[8] = $Columns[8]." ;total_sout ".$SoutSums{$id}." ;total_spnum ".$SpNums{$id};
			$Results{$id} = join("\t",@Columns);
			$Sout{$id}{$i} = $Info{sout_count};
		}
		close IN;
		$i++;
	}
	foreach my $id (sort keys %Sout){
		print "$id";
		for (0..($i - 1)) {
			my $order = $_;
			if(exists $Sout{$id}{$order}){
				print "\t$Sout{$id}{$order}";
			}else{
				print "\t0";
			}
		}
		print "\n";
	}
	return(\%Results);
}

if(@ARGV < 2){
	print "Usage: output\tgtf1\tgtf2...\t gtf... \n";
	exit(0);
}else{
	my $outfile = shift @ARGV;
	my %Results = %{&combineEv(@ARGV)};
	open OUT,">$outfile" or die $!;
	foreach my $id (sort keys %Results){
		print OUT "$Results{$id}\n";
	}
	close OUT;
}


