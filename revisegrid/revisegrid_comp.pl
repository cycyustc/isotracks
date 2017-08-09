#!/usr/bin/perl
# takes a file with tracks and revise the format for trilegal
# WARNING CHANGED TO USE DBERT_COMP FILES
use File::Basename;
use sort 'stable';		# guarantee stability
#use CGI qw/:standard/;

$argc = @ARGV;          # get the number of arguments
if ($argc<1) {
  die "Usage: $0 <revisegrid_exe> <grids_dir> <tppattern> <newgridname.dat> \n";
} else {
  $exe = $ARGV[0];
  $inpdir = $ARGV[1];
  $set = $ARGV[2];
  $tppattern = $ARGV[3];
  $outfile = $ARGV[4];
  print "Will check files in $inpdir and list correct files into $outfile\n";
  print "Track files will be looked for in $inpdir/dbert_comp\n";
  print "1TP files should have names $inpdir/${tppattern}*_1TP.INP\n";
}

#####at begin
$pwd = `pwd`;

open(TMPFILE, ">tmp");

print "gunzipping files...\n";
@inpfiles = <$inpdir/dbert_comp$set/*.gz>;
foreach $file (@inpfiles) {
    `gunzip -f $file`;
}

@inpfiles = <$inpdir/dbert_comp$set/*.HB>;
#find z,y pairs
foreach $file (@inpfiles) {
  $pz = index($file,'Z')+1;
  $py = index($file,'Y')+1;
  $z = substr($file,$pz); $zz=1.0*$z;
  $y = substr($file,$py); $yy=1.0*$y;
  $mh = log(($zz/(1.0-$zz-$yy))/0.0207)/log(10.0);
  push(@zs, $zz);
  print "found $zz $yy in $file ([M/H]=$mh)\n";
  $Ys{$zz}=$yy;
  $fileHB{$zz}=$file;
  $_=$file; s/HB/INT/; $fileINT{$zz}= $_ ;
  $_=$file; s/HB/LOW/; $fileLOW{$zz}= $_ ;
}

@tpfiles = <$inpdir/${tppattern}*1TP.INP>;
#find z,y pairs
foreach $file (@tpfiles) {
  $tpz = index($file,'Z')+1;
  $tpy = index($file,'Y')+1;
  $tz = substr($file,$tpz); $tzz=1.0*$tz;
  $ty = substr($file,$tpy); $tyy=1.0*$ty;
  push(@tzs, $tzz);
  print "1TP: found with $tzz $tyy in $file\n";
}

#sort by value does not work!!!!
@zsorted = sort {$a <=> $b} @zs;
#sort lexically then, at least it works
@zsorted = sort @zs;
#print @zsorted, "\n";
foreach $z (@zsorted)  {
    print "$z $Ys{$z}\n"."$fileLOW{$z} \n"."$fileINT{$z} \n"."$fileHB{$z} \n";
}

$nmet = keys(%Ys);
print "directory has $nmet metallicities\n";
print TMPFILE "$nmet\n";
foreach $z (@zsorted) {
  $y = $Ys{$z};
  print "doing Z=$z Y=$y\n";
  print TMPFILE "$z $y\n";
  ##### each track file parsed and corrected
  $tracklowfile = $fileLOW{$z};
  $trackhbfile = $fileHB{$z};
  $trackintfile = $fileINT{$z};
  $cutfile="";
  foreach $i (0..$#tpfiles) { #check  if there's a TP file
      if ($tzs[$i]==$z) {
	  print "with 1TP at $tpfiles[$i]\n";
	  $cutfile=$tpfiles[$i];
      }
  }
  print TMPFILE "$tracklowfile\n".
      "$trackhbfile $cutfile\n".
      "$trackintfile $cutfile\n";

  foreach $file ($tracklowfile,$trackhbfile) {
    open(TRACK, "<$file");
    #first pass to capture masses:
    @masses=();
    while ($inputline=<TRACK>) {
      chomp ($inputline);
      $_= $inputline;
      #print "read $inputline\n";
      if (/M= /) {
	$inputline =~ s#.*M= ##;
	($mass) = split(' ',$inputline);
	#print "found mass M=$mass\n";
	push(@masses,$mass);
      }
    }
    close TRACK;
    print "$file: found $#masses masses, last one is $mass Msun\n";
  }
}
close TMPFILE;
#print "Now proceed with: revisegrid/main tmp $outfile 0.02 0.005\n";
print "Now proceed with: $exe tmp $outfile 0.04 0.01\n";
#  0.04 0.01 was used for public database!!!!!

# WARNING: REVISEGRID MAIN SHOULD BE REVISED TO INCLUDE NEW VARIABLES, BOTH FROM DBERT_COMP AND 1tp FILES
# DO IT BEFORE REACTIVATING THIS. TEST WITH 
# revisegrid_comp.pl CAF09_V1.2S_M36_S12D_NS_MAS CAF09_V1.2S_M36_S12D_NS_MAS_comp.dat
#$a=`revisegrid/main tmp $outfile 0.04 0.01\n`;
# END OF WARNING

# this one instead is for the overshooting-test tracks:
print "executing: $exe tmp $outfile 0.04 0.01\n";
$a=`$exe tmp $outfile 0.04 0.01\n`;
#$a=`$exe tmp $outfile 0.08 0.02\n`;
#print "$a\n";
#sleep(10);
#print "gzipping files...\n";
#@inpfiles = <$inpdir/*.HB $inpdir/*.LOW $inpdir/*.INT $inpdir/*.INP>;
#foreach $file (@inpfiles) {
#    `gzip -f $file`; #discouraged, has caused errors in the past
#}
unlink tmp;
die "Normally terminated!\n";
