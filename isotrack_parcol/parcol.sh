#python ../revisegrid/dbert_colibri_replace.py ./S_EAGB_colibri.list > log.txt
#../revisegrid/revisegrid_comp.pl    ../revisegrid/main    CAF09_V1.2S_M36_S12D_NS_MAS3    S12_    CAF09_V1.2S_M36_S12D_NS_MAS3_parcol_comp.dat > log1.txt
#sed -i 's/isotrack\/parsec/isotrack\/parcol/g' CAF09_V1.2S_M36_S12D_NS_MAS3_parcol_comp.dat > log0.txt

#../revisegrid/main tmp CAF09_V1.2S_M36_S12D_NS_MAS3_parcol_comp.dat 0.04 0.01 > log0.txt

#the isotrack_parsec/CAF09_V1.2S_M36_S12D_NS_MAS3/dbert_comp/ptcri_CAF09_V1.2S_M36_S12D_NS_MAS_Z0.017_Y0.279.dat.INT contains two exactly the same M=6.00 tracks.

#S_004
SET=004

python ../revisegrid/dbert_colibri_replace.py S_$SET_colibri.list > error_$SET.txt
../revisegrid/revisegrid_comp.pl    ../revisegrid/main    CAF09_V1.2S_M36_S12D_NS_MAS3  $SET  S12_    CAF09_V1.2S_M36_S12D_NS_MAS3_parcol_comp$SET.dat > log1.txt
sed -i 's/isotrack\/parsec/isotrack\/parcol/g' CAF09_V1.2S_M36_S12D_NS_MAS3_parcol_comp$SET.dat


