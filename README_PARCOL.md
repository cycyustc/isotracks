**What's this?**

*PARSEC tracks with EAGB replaced by that from COLIBRI*

`revisegrid/dbert_colibri_replace.py` is the essentiall cdoe for replacing PARSEC E-AGB with that from COLIBRI.


**To run**
* create the soft link 'parcol' (isotracks/isotrack/parcol) to the `isotrack_parcol` dir
* organize the COLIBRI tracks (and INP files) in the same way as those already included (S_001/, S_002/, ... and INP/)
* create the .list file, as those already included (S_001_colibri.list, ...): one line for one Z. Each line contains the following information: dbert_comp file from PARSEC (remove the extension .HB, .LOW or .INT), the COLIBRI track dir, COLIBRI INP file, and output dir. Lines start with '#' are ignored.
* Change the set name list in `parcol.py`, then run `python parcol.py`. This script integrate the following three steps:

> 
- ```python ./dbert_colibri_replace.py ./S_EAGB_colibri.list```.
- `revisegrid_comp.pl` (As described in README_PARSEC.md)
- replace `parsec` with `parcol` in the `CAF09_V1.2S_M36_S12D_NS_MAS3_parcol_comp.dat` file (`isotrack/parsec` is written by the revisedgrid C code).

**LFS support and uncompress tgz files**
The INP files, COLIBRI tracks and the resulting dbert files are now stored with git LFS scheme, and they are in tgz format to largely save space.
After cloning the repository, execute the following steps to uncompress the tgz files:

>
- ```cd isotrack_parcol```
- ```cat *.tgz | tar zxvf - -i```
- ```cd CAF09_V1.2S_M36_S12D_NS_MAS3```
- ```cat *.tgz | tar zxvf - -i```


**CAVEATS**
* If the logL is applied to replace EAGB in PAESEC (for preventing the jump between PARSEC and COLIBRI EAGB), there will be a mass jump. This is 'cured' by interpolating Teff between PARSEC and COLIBRI E-AGB within an adjustable logL interval (eg., 2.--2.5).
* Some COLIBRI tracks start before 'He0' (the label in PARSEC).
* In default COLIBRI tracks with less than two TPs are taken as no TP phase, to be consistent with the colibri2trilegal code. This can be switched off.