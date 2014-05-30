#! /usr/bin/env perl
#
# Andrew Janke - rotor@cmr.uq.edu.au
# Center for Magnetic Resonance
# The University of Queensland
# http://www.cmr.uq.edu.au/~rotor
# 
# Copyright Andrew Janke, Mark Griffin The University of Queensland.
# Permission to use, copy, modify, and distribute this software and its
# documentation for any purpose and without fee is hereby granted,
# provided that the above copyright notice appear in all copies.  The
# author and the University of Queensland make no representations about the
# suitability of this software for any purpose.  It is provided "as is"
# without express or implied warranty.
# 
# decompose an affine matrix into translation, rotation, scale and shear

use strict;
use warnings "all";

use Math::Trig;
use Math::MatrixReal;
use Getopt::Tabular;
use File::Basename;

my($me) = &basename($0);
my($clobber) = '';
my($verbose) = 0;

my(@opt_table) = (
   ["-verbose", "boolean",  0, \$verbose, 
      "be verbose"],
   );

my($Usage) = "Usage: $me [options] <in.xfm>\n".
             "       $me -help to list options\n\n";

my($Help) = <<HELP;
   decompose a MNI .xfm file
HELP

# Check arguments
&Getopt::Tabular::SetHelp ($Help, $Usage);
&GetOptions (\@opt_table, \@ARGV) || exit 1;
die $Usage if ($#ARGV != 0);
my($xfmfile) = $ARGV[0];

# check for input file
die "$me: Couldn't find $xfmfile\n" if (!-e "$xfmfile");


my($Mtext, $M, @Trans, @Scale, @Shear, @Rotn);

# perl line noise to read in a text file Matrix into MatrixReal
my @MNIdump = `cat $xfmfile`;

my $line;
do{ $line = shift(@MNIdump) or die "$me: Couldn't find a MNI Transform in: $xfmfile\n\n";} 
until $line =~ m/Linear_Transform/;

foreach (@MNIdump[0..2]){
   chomp;
   s/\;//;
   $Mtext .= "[ $_ ]\n";
   }
$Mtext .= "[ 0 0 0 1 ]\n";

$M = Math::MatrixReal->new_from_string($Mtext);
decompose_transformation_matrix($M, \@Trans, \@Scale, \@Shear, \@Rotn);

my $c = 0;
foreach (@Rotn){
   $Rotn[$c] = rad2deg($_);
   $c ++;
   }

printf("Translation: %10.5f %10.5f %10.5f\n", @Trans);
printf("Rotation:    %10.5f %10.5f %10.5f\n", @Rotn);
printf("Scale:       %10.5f %10.5f %10.5f\n", @Scale);
printf("Shear:       %10.5f %10.5f %10.5f\n", @Shear);


# Subroutines #######################################################

# Takes a 4x4 transformation matrix $M (using the Math:MatrixReal module)
# Returns the translation, scale, shear and rotations encoded in the input
# matrix.
#
# Andrew Janke - rotor@cmr.uq.edu.au
# Losely based on Louis Collins' make_rots.c from the mni_autoreg package
# With substantial help from Mark Griffin - mark.griffin@cmr.uq.edu.au
sub decompose_transformation_matrix{
   my($M, $Trans, $Scale, $Shear, $Rotn) = @_;

   my ($Sx, $Sy, $Sz, $SHa, $SHb, $SHc, $Rx, $Ry, $Rz);

   # TRANSLATIONS - [M] = [T][S][SH][R]
   # As of yet I am assuming the center of rotation is (0,0,0) as
   # I am not exactly sure as to what SPM does here.
   @$Trans[0] = $M->element(1, 4);
   @$Trans[1] = $M->element(2, 4);
   @$Trans[2] = $M->element(3, 4);
   
   my $T = $M->shadow(); $T->one();         # Create and zero the Translation Matrix
   $T->assign(1, 4, @$Trans[0]);
   $T->assign(2, 4, @$Trans[1]);
   $T->assign(3, 4, @$Trans[2]);
   
   # SCALES - [M] = inv[T][T][S][SH][R] = [S][SH][R]
   # Here we use an identical method to Louis Collins's in mni_autoreg
   # Namely multiply a unit vector in each direction and measure the length
   # after the transformation.
   $M = $T->decompose_LR()->invert_LR() * $M;
   my $SSHRinv = $M->decompose_LR()->invert_LR();
   
   my $Unit = Math::MatrixReal->new(4, 1);
   $Unit->zero(); $Unit->assign(1, 1, 1); $Sx = ($SSHRinv * $Unit)->length();
   $Unit->zero(); $Unit->assign(2, 1, 1); $Sy = ($SSHRinv * $Unit)->length();
   $Unit->zero(); $Unit->assign(3, 1, 1); $Sz = ($SSHRinv * $Unit)->length();
   
   my $Sinv = $M->shadow(); $Sinv->zero();   # Create and zero the inverse Scaling Matrix
   $Sinv->assign(1, 1, $Sx);     $Sx = 1/$Sx;
   $Sinv->assign(2, 2, $Sy);     $Sy = 1/$Sy;
   $Sinv->assign(3, 3, $Sz);     $Sz = 1/$Sz;
  
  
   # SHEARS - [M] = inv[T][T]inv[S][S][SH][R] = [SH][R]
   # We assume the shear matrix:      SH [ 1 0 0 0 ]
   # where x' = x                        [ a 1 0 0 ]
   #       y' = ax + y                   [ b c 1 0 ]
   #       z' = bx + cy + z              [ 0 0 0 1 ]
   # 
   # However as M at this point is in fact [SH][R]
   # we can extract a, b and c as such:
   # 
   # let [ M1 ]
   #     [ M2 ]  =  [SH][R]
   #     [ M3 ] 
   #
   # thus:
   # 
   # a = (M2 . M1) / |M1|^2
   # b = (M3 . M1) / |M1|^2
   # c = (M3 . T)  / |T|^2     where T = M2 - (a . M1)
   #
   # We could also use the determinant to determine whether we have an 
   # Orthogonal matrix and thus don;t have shears, but we don't do this yet....
   
   $M = $Sinv * $M;
   
   # check determinant for "sheariness" if det != 0 shears exist.
   # my $det = $M->decompose_LR()->det_LR();
   
   my $M1 = ~$M->row(1);
   my $M2 = ~$M->row(2);
   my $M3 = ~$M->row(3);
   
   $SHa = $M2->scalar_product($M1) / $M1->scalar_product($M1);
   $SHb = $M3->scalar_product($M1) / $M1->scalar_product($M1);
   my $TMP = $M2 - ($SHa * $M1);
   $SHc = $M3->scalar_product($TMP) / $TMP->scalar_product($TMP);
   
   my $SH = $M->shadow(); $SH->one();      # Create and zero the Shear Matrix
   $SH->assign(2, 1, $SHa);
   $SH->assign(3, 1, $SHb);
   $SH->assign(3, 2, $SHc);
   
   
   # ROTATIONS - [M] = inv[T][T]inv[S][S]inv[SH][SH][R] = [R]
   # We assume cy is positive to ensure we get one of the 2 possible solutions
   # where rotations are between -pi and pi.
   #
   # Here we deduce Rx, Ry and Rz by virtue or the fact that the rotation
   # matrix is as follows. 
   #
   # R = [ cos(Ry)*cos(Rz)  <stuff>          <stuff>          0 ]
   #     [ cos(Ry)*sin(Rz)  <stuff>          <stuff>          0 ]
   #     [ sin(Ry)          sin(Rx)*cos(Ry)  cos(Rx)*cos(Ry)  0 ]
   #     [ 0                0                0                0 ]
   #
   # Then the quadrant of the angle must be deduced by the sign of
   # cos and sin for the particular rotation.
   
   $M = $SH->decompose_LR()->invert_LR() * $M;
   
   # Get Y Rotation and check that we aren't up a creek without a paddle
   my $sy = $M->element(3, 1);
   if (abs($sy) == 1) { die "cos X = 0. I haven't solved this yet...\n"; }
   $Ry = asin($sy);
   
   # Get X Rotation
   my $cy = cos($Ry);
   my $sx = $M->element(3, 2) / $cy;
   my $cx = $M->element(3, 3) / $cy;
   $Rx = asin($sx);
   if ($cx < 0){
      if ($sx > 0){ $Rx =  pi() - $Rx; }   # quadrant 2
      else        { $Rx = -pi() - $Rx; }   # quadrant 3
      }
   
   # Get Z Rotation
   my $cz = $M->element(1, 1) / $cy;
   my $sz = $M->element(2, 1) / $cy;
   $Rz = asin($sz);
   if ($cz < 0){ 
      if ($sz > 0){ $Rz =  pi() - $Rz; }   # quadrant 2
      else        { $Rz = -pi() - $Rz; }   # quadrant 3
      }
   
   
   # If verbose do a bit of checking and output the remainder which should
   # be the identity matrix or close to it  
   if ($verbose){ 
      my $RX = $M->shadow(); $RX->one();     # Create and zero the X Rotation Matrix
      $RX->assign(2, 2,  cos($Rx));   $RX->assign(2, 3, -sin($Rx));
      $RX->assign(3, 2,  sin($Rx));   $RX->assign(3, 3,  cos($Rx));
   
      my $RY = $M->shadow(); $RY->one();     # Create and zero the Y Rotation Matrix
      $RY->assign(1, 1,  cos($Ry));   $RY->assign(1, 3, -sin($Ry));
      $RY->assign(3, 1,  sin($Ry));   $RY->assign(3, 3,  cos($Ry));
   
      my $RZ = $M->shadow(); $RZ->one();     # Create and zero the Z Rotation Matrix
      $RZ->assign(1, 1,  cos($Rz));   $RZ->assign(1, 2, -sin($Rz));
      $RZ->assign(2, 1,  sin($Rz));   $RZ->assign(2, 2,  cos($Rz));
   
      my $RZYXinv = ($RZ * $RY * $RX)->decompose_LR()->invert_LR();
      print "Remainder:\n" . ($RZYXinv * $M);
      }
   
   @$Scale[0] = $Sx;   @$Scale[1] = $Sy;   @$Scale[2] = $Sz;
   @$Shear[0] = $SHa;  @$Shear[1] = $SHb;  @$Shear[2] = $SHc;
   @$Rotn[0]  = $Rx;   @$Rotn[1]  = $Ry;   @$Rotn[2]  = $Rz;
   }

# Creates a 4x4 transformation matrix $M using the input
# Translations, scales, shears and rotations (or not)
sub create_transformation_matrix{
   my($Trans, $Scale, $Shear, $Rotn) = @_;

   # set a few defaults
   my $c;
   for ($c = 0; $c < 3; $c ++){
      if (!defined @$Trans[$c]){ @$Trans[$c] = 0; }
      if (!defined @$Scale[$c]){ @$Scale[$c] = 1; }
      if (!defined @$Shear[$c]){ @$Shear[$c] = 0; }
      if (!defined @$Rotn[$c] ){ @$Rotn[$c]  = 0; }
      }

   my $M = Math::MatrixReal->new(4, 4); $M->one();
   $M->assign(1, 4, @$Trans[0]);
   $M->assign(2, 4, @$Trans[1]); 
   $M->assign(3, 4, @$Trans[2]);
   
   $M->assign(1, 1, @$Scale[0]);
   $M->assign(2, 2, @$Scale[1]);
   $M->assign(3, 3, @$Scale[2]);
  
   my $SH = $M->shadow(); $SH->one();      # Create and zero the Shear Matrix
   $SH->assign(2, 1, @$Shear[0]);
   $SH->assign(3, 1, @$Shear[1]);
   $SH->assign(3, 2, @$Shear[2]);
   
   my $RX = $M->shadow(); $RX->one();     # Create and zero the X Rotation Matrix
   $RX->assign(2, 2,  cos(@$Rotn[0]));   $RX->assign(2, 3, -sin(@$Rotn[0]));
   $RX->assign(3, 2,  sin(@$Rotn[0]));   $RX->assign(3, 3,  cos(@$Rotn[0]));
   
   my $RY = $M->shadow(); $RY->one();     # Create and zero the Y Rotation Matrix
   $RY->assign(1, 1,  cos(@$Rotn[1]));   $RY->assign(1, 3, -sin(@$Rotn[1]));
   $RY->assign(3, 1,  sin(@$Rotn[1]));   $RY->assign(3, 3,  cos(@$Rotn[1]));
   
   my $RZ = $M->shadow(); $RZ->one();     # Create and zero the Z Rotation Matrix
   $RZ->assign(1, 1,  cos(@$Rotn[2]));   $RZ->assign(1, 2, -sin(@$Rotn[2]));
   $RZ->assign(2, 1,  sin(@$Rotn[2]));   $RZ->assign(2, 2,  cos(@$Rotn[2]));
   
   return $M * $SH * $RZ * $RY * $RX;
   }

sub print_transformation_matrix{
   my($M, $name) = @_;
   my (@Trans, @Scale, @Shear, @Rotn, @Rotnd);
   
   print "$name:\n$M";
   decompose_transformation_matrix($M, \@Trans, \@Scale, \@Shear, \@Rotn);
   foreach (@Rotn){ push(@Rotnd, rad2deg($_)); }
   printf("Translation: %10.5f %10.5f %10.5f\n", @Trans);
   printf("Rotation:    %10.5f %10.5f %10.5f\n", @Rotnd);
   printf("Scale:       %10.5f %10.5f %10.5f\n", @Scale);
   printf("Shear:       %10.5f %10.5f %10.5f\n", @Shear);
   }
