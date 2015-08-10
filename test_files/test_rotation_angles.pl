#!/usr/local/bin/perl -w
use Statistics::Basic;
use Math::Matrix;
use Math::Trig;
use strict;

my $dir = '/home/lixr/Human_1k/temp/';
my $loading_file = $dir.'PC123_loading';
#my $hmp_file = $dir.'H1k_hq_f0.1.hmp';
my $pc_score_file = $dir.'PCA_by_EIG';
my $pc_score_hashref = PC_score();

#my $exome_snp_hashref = Exome_SNP_list();

my @angles_1 = (5); ##(-20, -10, 0, 10, 20); ## , 
my @angles_2 = (0);
my @angles_3 = (-30);
my $outfile = $dir.'pc123_base_content_correlation_10';
open (O, '>>'.$outfile) || die;

#print O "Angle_1\tAngle_2\tAngle_3\tr1\tr2\tr3\n";
foreach my $an_1 (@angles_1) {
#	next unless $an_1 == 10;
	foreach my $an_2 (@angles_2) {
#		next unless $an_2 == 0;
		foreach my $an_3 (@angles_3) {
#			next unless $an_3 == -20;
#			my $rB_matrix = Rotation_matrix($an_1, $an_2, $an_3);
#			my $pc_snp_hashref = Rotation_back_spherical($rB_matrix);

			my $pc_snp_hashref = Parse_rB_PC_tag_file();
			my $pc_snp_base_content_hashref = Calculate_PC_snps_base_content($pc_snp_hashref);
			my $rs = Correlation_PC_score_base_content($pc_snp_base_content_hashref, $pc_score_hashref);
			print O "$an_1\t$an_2\t$an_3".$rs."\n";
			}
		}
	}

##############################################
sub Parse_rB_PC_tag_file {
#	my ($exome_snp_hashref) = @_;
	my $f = $dir.'rB_10_0_-30';
	my %hash;
	open (F, $f) || die;
	my ($i, $j);
	while (<F>) {
		chomp;
		my @t = split /\t/;
		my $id = $t[1];
			$id =~ s/\"//g;
		my ($pc1, $pc2, $pc3) = ($t[-3], $t[-2], $t[-1]); 
		next if ($pc1 + $pc2 + $pc3) == 0;
			 $hash{$id}{1} = 1 if $pc1 == 1;
			 $hash{$id}{2} = 1 if $pc2 == 1;
			 $hash{$id}{3} = 1 if $pc3 == 1;
		
#		if (exists $$exome_snp_hashref{$id}) {
#			
#			$i ++;
#			print $id." 1\n" if $i <= 2;	
#			 $hash{$id}{1} = 1 if $pc1 == 1;
#			 $hash{$id}{2} = 1 if $pc2 == 1;
#			 $hash{$id}{3} = 1 if $pc3 == 1;
#			}
#			else { $j ++; print $id." 2\n" if $j <= 10 }
		}
	close F;
	return \%hash;	
	}
sub Exome_SNP_list {
	my %hash;
	my $f = '/home/lixr/Human_1k/PC_as_trait/bed_files/Exome_snp_list';
	open (F, $f) || die;
	while (<F>) {
		chomp;
		my @t = split /\s+/;
		$hash{$t[1]} = 1;
		last if $t[2] == 2;
		}
	close F;
	return \%hash;	
	}	

sub Correlation_PC_score_base_content {
	my ($base_content_hashref, $score_hashref) = @_;
	my $rs ;
	for (my $pc = 1; $pc <= 3; $pc ++) {
		my @base_content = @{ $$base_content_hashref{$pc} };
		my $a = Statistics::Basic::vector(@base_content);
		my @score = @{ $$score_hashref{$pc} };
		my $b = Statistics::Basic::vector(@score);
		my $correlation = Statistics::Basic::correlation( $a, $b);
		  $rs .= "\t".$correlation;
		}
	return $rs;	
	}
	
sub PC_score {
	my %hash;
	open (Eig, $pc_score_file) || die;
	while (<Eig>) {
		chomp;
		my @t = split /\t/;
		next if /ind/;
		for (my $i = 1; $i <= 3; $i ++) {
			push @{ $hash{$i} }, $t[$i];
			}
		}
	close Eig;
	return \%hash;
	}	
sub Calculate_PC_snps_base_content {
	my ($pc_snp_hashref) = @_;
	my (%strain_base_count, %strain_snp_count); ## , %density
	for (my $ch = 1; $ch <= 22; $ch ++) {
		my $hmp_file = '/home/lixr/Human_1k/data_4_GAPIT/hmp_files/chunk_0.05/hq_ch'.$ch.'.hmp';
		
#		my $hmp_file = $dir.'H1k_hq_f0.1.hmp';
		open (Hmp, $hmp_file) || die;
		while (<Hmp>) {
			chomp;
			my @t = split /\t/;
			my $id = $t[0];
			next unless exists $$pc_snp_hashref{$id};
			my @pc_tags = keys %{ $$pc_snp_hashref{$id} };
			for (my $i = 11; $i <= $#t; $i ++) { ## 
				my $geno = $t[$i] =~ /N/ ? 'NN' : $t[$i];
				next if $geno =~ /N/;
				my ($g1, $g2) = $geno =~ /(\S)(\S)/;
#				my $strain = $$strain_id_arrayref[$i];
				my $strain = $i ;
				foreach my $pc_tag (@pc_tags) {
					$strain_snp_count{$strain}{$pc_tag} ++;
					$strain_base_count{$strain}{$pc_tag}{$g1} ++ if $g1 eq 'A';
					$strain_base_count{$strain}{$pc_tag}{$g2} ++ if $g2 eq 'A';
					
					}
				}
#			$k ++; last if $k > 100;	
			}
		close Hmp;
		}
#	my $k;
	my $out = $dir.'base_5_0_-30';
	open (OUT, '>'.$out) || die;
	print OUT "Strain\tA_PC1\tTotal_PC1_2n\tA_PC2\tTotal_PC2_2n\tA_PC3\tTotal_PC3_2n\n";
	my %hash;
	for (my $i = 11; $i < 1092 + 11; $i ++)  {
#		my $strain = $$strain_id_arrayref[$i];
		my $strain = $i;
		print OUT $strain;
		for (my $pc = 1; $pc <= 3; $pc ++) {
			my $total_num = exists $strain_snp_count{$strain}{$pc} ? 2 * $strain_snp_count{$strain}{$pc} : 1;
			my $a = exists $strain_base_count{$strain}{$pc}{'A'} ? sprintf "%.3f",$strain_base_count{$strain}{$pc}{'A'} / $total_num : 0;
			push @{ $hash{$pc} }, $a ;
			print OUT "\t".$a."\t".$total_num;
			}
		print OUT "\n";	
		}	
	close OUT;	
	return \%hash;
	}

sub Rotation_back_spherical {
	my ($W) = @_;
	my %hash;
	open (OL, $loading_file) || die;
	my $out = $dir.'rB_5_0_-30_f0.1';
	open (O2, '>'.$out) || die;
#	my $k;
	while (<OL>) {
		chomp;
		my $line = $_;
		my @t = split /\t/;
		next if /Loading/;
		my ($pc1, $pc2, $pc3) = @t[2..4];
		my $xyz = new Math::Matrix([$pc1], [$pc2], [$pc3]);
		my $xyz_rB = $W->multiply($xyz)->transpose;
		my ($n, $x, $y, $z) = split /\s+/, $xyz_rB;
		  $x = ($x * 350) ; $y = $y * 350 ; $z = ($z * 350) ;
		my $r = sqrt($x**2 + $y**2 + $z**2) ;
		my $theta =  acos_angle($z, $r) ;
		my $phi = atan_angle($x,$y);
		my ($pc1_tag, $pc2_tag, $pc3_tag) = PC_SNPs_tag($x, $y, $z, $r, $phi, $theta);
#		next if ($pc1_tag + $pc2_tag + $pc3_tag) == 0;
		my $id = $t[1];
	   $id =~ s/\"//g;
		if ($pc1_tag == 1) { $hash{$id}{1} = 1};
		if ($pc2_tag == 1) { $hash{$id}{2} = 1};
		if ($pc3_tag == 1) { $hash{$id}{3} = 1};
		print O2 $line."\t"."$x\t$y\t$z\t".$pc1_tag."\t".$pc2_tag."\t".$pc3_tag."\n";
#		$k ++; last if $k > 1000;
		}
	close OL;
	close O2;
	return \%hash;	
	}

sub PC_SNPs_tag {
	my ($x, $y, $z, $r, $phi, $theta) = @_;
#	my $xy = sqrt($x **2 + $y **2);
	my ($pc1, $pc2, $pc3) = (0, 0, 0);
	if ($theta <= 30 && $r >= 0) {$pc3 = 1; };
	if ($r >= 0.0 && $theta >= 30 && $phi >= 90 && $phi <= 210 ) {
#		$pc2 = 1 unless ( $z < 0 && $xy <= 0.5/sqrt(3)) ## for b
#		$pc2 = 1 unless ( $z < 0 && $xy <= 0.1) ## for c
		$pc2 = 1 unless $theta > 150 ; ## for d;
#		$pc2 = 1;
		}
	if ($r >= 0.0 && $theta >= 30 && ($phi <= 90 || $phi >= 330))  {
#		$pc1 = 1 unless ($z < 0 && $xy <= 0.5/sqrt(3)); ## for b
#		$pc1 = 1 unless ( $z < 0 && $xy <= 0.1) ## for c
		$pc1 = 1 unless $theta > 150 ; ## for d;
#		$pc1 = 1;
		}	
	return ($pc1, $pc2, $pc3);	
	}	
	
sub Rotation_matrix {
	my ($B, $C, $D) = @_;
	my $d = $D / 180 * pi;
	my $b = $B / 180 * pi;
	my $c = $C / 180 * pi;
	my $B_matrix = new Math::Matrix ([1, 0,       0       ],
	                                 [0, cos($b), -sin($b)],
	                                 [0, sin($b), cos($b) ] );
	my $C_matrix = new Math::Matrix ([cos($c),  0, sin($c)],
	                                 [0,        1, 0      ],
	                                 [-sin($c), 0, cos($c)] );
	                                 
	my $D_matrix = new Math::Matrix ([cos($d), -sin($d), 0],
	                                 [sin($d), cos($d),  0],
	                                 [0,       0,        1] );
	                                
	my $DC_matrix = $D_matrix->multiply($C_matrix);  
	my $DCB_matrix = $DC_matrix -> multiply($B_matrix);
	my $DCB_inv = $DCB_matrix->invert;
	return $DCB_inv;
	}	
	
sub acos_angle {
	my ($z, $r) = @_;
	my $theta;
	if (abs($z) >= $r) {
		$theta = $z >= 0 ? 0 : 180;
		}
		else {
			$theta = sprintf "%.2f", acos($z/$r) * 180 / pi;
			}
	return $theta;	
	}

sub atan_angle {
	my ($x, $y) = @_;
	my $phi;
	if ($y == 0 && $x >= 0) {
		$phi = 0;
		}
		elsif ( $y == 0 && $x < 0) {
			$phi = 180;
			}
			elsif ($x == 0 && $y > 0) {
				$phi = 90
				}
				elsif ($x == 0 && $y < 0) {
					$phi = 270
					}
					else {
						my $atan = atan($y / $x) * 180 / pi;
						
						if ($x > 0 && $y > 0) {
							$phi = $atan;
							}
							elsif ( $x < 0 && $y < 0) {
								$phi = 180 + $atan;
								}
								elsif ($x > 0 && $y < 0) {
									$phi = 360 + $atan
									}
									elsif ($x < 0 && $y > 0) {
										$phi = 180 + $atan
										}
						}
		$phi = sprintf "%.2f", $phi;
		return $phi;	 
	}	

	