python ../revisegrid/dbert_colibri_replace.py ./S_EAGB_colibri.list > log.txt
../revisegrid/revisegrid_comp.pl    ../revisegrid/main    CAF09_V1.2S_M36_S12D_NS_MAS3    S12_    CAF09_V1.2S_M36_S12D_NS_MAS3_parcol_comp.dat
sed -i 's/isotrack\/parsec/isotrack\/parcol/g' CAF09_V1.2S_M36_S12D_NS_MAS3_parcol_comp.dat > log0.txt
