#PARSEC tracks with EAGB replaced by that from COLIBRI

##To run

`revisegrid/dbert_colibri_replace.py` replaces PARSEC E-AGB with that from COLIBRI, to run:
* create the `isotrack_parcol` dir, and soft link 'parcol' (isotracks/isotrack/parcol) to this `isotrack_parcol` dir.
* get COLIBRI tracks with early-AGB phase from Paola, place them into S_EAGB_colibri (or anywhere you like), extract tgz file into subfolders of different metallicities, e.g., S12_Z0.02_Y0.284_EAGB_1TP
* write an S_EAGB_colibri.list, one line for one Z. Each line contains the following information: dbert_comp file from PARSEC (remove the extension .HB, .LOW or .INT), the COLIBRI track dir, COLIBRI INP file, and output dir. Lines start with '#' are ignored.
* run ```python ./dbert_colibri_replace.py ./S_EAGB_colibri.list```.
* run `revisegrid_comp.pl` (As described in README_PARSEC.md)
* replace `isotrack/parsec` with `isotrack/parcol` in the `CAF09_V1.2S_M36_S12D_NS_MAS3_parcol_comp.dat` file (`isotrack/parsec` is written by the revisedgrid C code).
* the above three steps has been integrated into `parcol.py`. So, you simply have to create the list file and put the set name in `parcol.py`, then run with `python parcol.py`


##CAVEATS

* If the logL is applied to replace EAGB in PAESEC (for preventing the jump between PARSEC and COLIBRI EAGB), there will be a mass jump. This is 'cured' by interpolating Teff between PARSEC and COLIBRI E-AGB within an adjustable logL interval (eg., 2.--2.5).
* Some COLIBRI tracks start before 'He0' (the label in PARSEC).
* In default COLIBRI tracks with less than two TPs are taken as no TP phase, to be consistent with the colibri2trilegal code. This can be switched off.