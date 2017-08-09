# HOW PROCESS AND REVISE PARSEC TRACKS 

## Introduction

PARSEC tracks are generally computed by Sandro and then processed by an IDL code that finds all equivalent evolutionary points (EEPs), and writes the entire grids into 3 files, hereafter dubbed `dbert` files:
* `*.LOW` - low mass stars up to the TRGB
* `*.HB` - low mass stars from HeB to 1st thermal pulse
* `*.INT` - intermediate- and high-mass stars until 1st TP or C Burning

The latest version of `dbert` files are those inside the `dbert_comp` subfolders, see for instance `CAF09_V1.2S_M36_S12D_NS_MAS3/dbert_comp/*.HB`. . Their structure is:

* *header* : number of tracks in grid, their masses, and number of points
* *track1* : track with L, Teff, mass, and chemical composition. the last column (if present) flags the EEPs
* *track2*, ..., *trackN* : the same
* *tail*: number of EEPs in every track, and their positions along the tracks

`dbert_comp` files are not ready to be used in TRILEGAL. They contain too many points, they do not have a suitable interface between LOW+HB and INT files, and their EAGB sections may have been recomputed by COLIBRI. Therefore, they have to be processed further:


## Using revisegrid_comp.pl

`revisegrid_comp.pl` will change the resolution of the tracks, revise/fix the interfaces between grids of tracks, and write warnings in case of some more evident problems.

First compile the C code:
```
cd revisegrid; make; cd ..
```

Then run the perl script. Usage: 
```
./revisegrid_comp.pl <revisedgrid_exe> <grids_dir> <tppattern> <newgridname.dat>
```
It will check files in `dbert_comp` format in folder `<grids_dir>/dbert_comp` and list correct files into `<newgridname.dat>`. 1TP files from Paola should have names `<grids_dir>/<tppattern>*_1TP.INP`: they will be used to check if the 1TP is the same between PARSEC and COLIBRI. <revisedgrid_exe> is inserted by Yang to speficy the location of the `revisedgrid/main`.

Example: 
```
./revisegrid_comp.pl ./revisegrid/main CAF09_V1.2S_M36_S12D_NS_MAS3 S12_ CAF09_V1.2S_M36_S12D_NS_MAS3_comp.dat
```

*WARNING*: `revisegrid_comp.pl` contains the resolution to be adopted in the new tracks. The one suitable for new PARSEC tracks seems to be:
```
$a=`revisegrid/main tmp $outfile 0.01 0.0025\n`; #ususally '0.04 0.01' is enough
```

Note: there is also a `revisegrid.pl` which was used in older version of dbert files.

