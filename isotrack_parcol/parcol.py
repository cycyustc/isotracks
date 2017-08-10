import os

for SET in ['001', '002', '003', '004', '005', '006', '007']:
    exe1='python ../revisegrid/dbert_colibri_replace.py S_'+SET+'_colibri.list > error_S_'+SET+'.txt'
    exe2='../revisegrid/revisegrid_comp.pl    ../revisegrid/main    CAF09_V1.2S_M36_S12D_NS_MAS3  '+SET+'  S12_    CAF09_V1.2S_M36_S12D_NS_MAS3_parcol_comp'+SET+'.dat > log_S_'+SET+'.txt'
    exe3="sed -i \'s/isotrack\/parsec/isotrack\/parcol/g\' CAF09_V1.2S_M36_S12D_NS_MAS3_parcol_comp"+SET+".dat"

    print exe1
    os.system(exe1)
    print exe2
    os.system(exe2)
    print exe3
    os.system(exe3)

exit(0)
