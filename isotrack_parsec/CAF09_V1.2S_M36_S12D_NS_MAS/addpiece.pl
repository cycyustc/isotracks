#!/usr/bin/perl
#apparently this script was done to add a missing HB track at the end of .HB files, creating the .HB.new 
@all1TPs = glob "*.INP";
@allINT = glob "*.INT";
foreach $INTfile (@allINT) {
    $HBfile = $INTfile; $HBfile =~ s/.INT/.HB/;
    foreach $file (@all1TPs) {
 	$ZY = $file; $ZY =~ s/S12_//; $ZY =~ s/_1TP.INP//; 
	$TPfile = $file;
	$_ = $INTfile; 
	last if (m/$ZY/);
    } 
    foreach $file ($INTfile,$HBfile) {
	$_ = $file; $thisisHB = 0; $thisisHB = 1 if (m/HB/); 
	$file2 = $file . ".new";
###############################################
	print "***** Doing $file $TPfile $ZY *****\n";
	open (GRID, "<$file");
	open (GRIDZ, ">$file2");
	#open (GRIDZ, "/dev/tty");
	undef $mass;
	print "$file is HB $thisisHB\n";
	if ($thisisHB==0) {
	    $line=<GRID>; ($n1,$n2) = split " ", $line; print GRIDZ "$line";
	}
	$line=<GRID>; ($ngroups) = split " ", $line; print GRIDZ "$line";
	if ($thisisHB==1) {
	    $line=<GRID>; print GRIDZ "$line";
	}
	foreach $n (1..$ngroups) {
	    $line=<GRID>; ($ng1,$ng2) = split " ", $line; print GRIDZ "$line";
	}
	$line=<GRID>; ($ntracks) = split " ", $line; print GRIDZ "$line";
	$count=0; 
	@npoints = undef; @masses = undef;
	print "number of tracks = $ntracks\n";
	while ($#npoints<$ntracks) {
	    $line=<GRID>; (@temp) = split " ", $line; print GRIDZ "$line";
	    push (@npoints, @temp); 
	}
	while ($#masses<$ntracks) {
	    $line=<GRID>; (@temp) = split " ", $line; print GRIDZ "$line";
	    push (@masses, @temp); 
	}
	print "masses @masses\n";	
	print "npoints @npoints\n";	
	foreach $m (1..$ntracks) {
	    $nstars = $npoints[$m]; 
	    print "number of stars in $m = $nstars\n";
	    foreach $i (0..$nstars-1) {
		#print " $i\n";
		$line=<GRID>; print GRIDZ "$line";
		($age,$logl,$logte,undef,$model,undef,$label,$mass) = split " ", $line;
		#print " $model $age $logl $logte \n";
################# start of track: look for its 1TP true value
		if ($age eq '0.000000000000E+00') {
		    print "mass in grid file = $mass\n";
		    open (TPFILE, "<$TPfile");
		    while ($_=<TPFILE>) {
			($TPmodel,undef,undef,undef,undef,$TPmass,
			 $TPlogl, $TPlogte) = split;
			last if ($TPmass==$mass);
		    }
		    close TPFILE;
		    if ($TPmodel>100) {
			print "1TP in $TPfile is at model $TPmodel, $TPlogl, $TPlogte \n";
		    } else {
			next;
		    }	
		}	
################ end of track: look for last tabulated point
		if ($i==$nstars-1) {
		    print "last tabulated point = $model $age $logl $logte $line \n";
		    #print "finds missing piece in TRACK file\n";
		    $dir = "S12D_NS_".$ZY;
		    if ($thisisHB==0) {
			@allpossiblefiles = glob "$dir/*.PMS";
		    } else {
			@allpossiblefiles = glob "$dir/*.PMS.HB";
		    }   
		    foreach $trackfile (@allpossiblefiles) {
			$masstrackfile = $trackfile;
			$masstrackfile =~ s/[^M]*M//; $masstrackfile =~ s/.PMS[\S]*//;
			if ($masstrackfile==$TPmass) {
			    print "trackfile might be $trackfile\n" ;
			    open (TRACK, "<$trackfile");
			    while ($linet = <TRACK>) {
				 ($temp) = split " ", $linet ;
				 if ($temp eq 'Starting') {
				     $deltal=0; $flag=0; $count=0;
				     while ($linet = <TRACK>) {
					 ($tmod,undef,$tage,undef,$tlogl,$tlogte) = 
					     split " ", $linet ;
					  if ($flag==0) {
					      $tlogl0 = $tlogl;
					      $flag=1;
					  }
					 $deltal = abs($tlogl-$tlogl0);
					 if ($tmod>$model && $tmod<$TPmodel && $deltal>0.02) {
					     $aget = 1.0*$tage;
					     print "$aget $tlogl $tlogte 0 $tmod\n"; 
					     print GRIDZ " $aget $tlogl $tlogte 0 $tmod\n";
					     $deltal=0; $flag=0; $count++;
					 }
				     }
				 }
			    }
			    close TRACK;
			    $npoints[$m] = $npoints[$m]+$count ;
			}
		    }
		}
##############	
	    }
	}
	while($line=<GRID>) {
	    print GRIDZ "$line";
	}
	close GRID;
	foreach $m (0..$ntracks) {
	    print GRIDZ " $npoints[$m]";
	} print GRIDZ "\n";
	close GRIDZ;
###############################################
    }
    die;
} 
