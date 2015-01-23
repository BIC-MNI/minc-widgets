#! /usr/bin/env perl


use warnings "all";

$opt{'verbose'} = 1;
$opt{'fake'} = 0;


$infile = $ARGV[0];
$outfile = $ARGV[1];

$hdrfile = "$infile.hdr";

open(FH, "<$hdrfile");
foreach (<FH>){
   chomp;
   if(m/x_dimension/){
		(undef, $xsize) = split(/\ /, $_);
	   }
   if(m/y_dimension/){
		(undef, $ysize) = split(/\ /, $_);
	   }
   if(m/z_dimension/){
		(undef, $zsize) = split(/\ /, $_);
	   }
   }
close(FH);

# create config file
$conffile = $outfile;
$conffile =~ s/\.h5$/\.conf/;
open(FH, ">$conffile");
$conftxt =
   "INPUT-CLASS IN\n" .
   "INPUT-SIZE 16\n" .
   "RANK 3\n" .
   "DIMENSION-SIZES $zsize $ysize $xsize\n";

print "++++CONFTXT++++\n" . $conftxt . "+++++++++\n" if $opt{'verbose'};

print FH $conftxt;
close(FH);


# convert the sucker
&do_cmd('h5import', $infile, '-c', $conffile, '-o', $outfile);


sub do_cmd { 
   print STDOUT "@_\n" if $opt{'verbose'};
   if(!$opt{'fake'}){
      system(@_) == 0 or die;
      }
   }
