import csv
import os
import SYMMETRY_DETECTION_tolerances

#SYMMETRY_DETECTION.cif_dir = "OutputCIFs/Modified"

atom_tolerances = [x / 10.0 for x in range(1, 25, 1)] # [x / 10.0 for x in range(50, 51, 1)]
lat_tolerances  = [x / 100.0 for x in range(10, 100, 5)]#  [x / 10.0 for x in range(1, 12, 3)]
relative_tolerances = [0.999]



print(atom_tolerances)
print(relative_tolerances)
#print(relative_tolerances)
file = open('tolerances.csv', 'w')
file = csv.writer(file)
file.writerow(['atom', 'lattice', 'relative', 'found / 156'])
#os.mkdir('Tolerances')
for tol in atom_tolerances:
    print(tol)
    SYMMETRY_DETECTION_tolerances.atom_tol = tol
    for tol2 in lat_tolerances:
        print(tol2)
        SYMMETRY_DETECTION_tolerances.lat_tol = tol2
        for tol3 in relative_tolerances:
            print(tol3)
            SYMMETRY_DETECTION_tolerances.atom_rel_tol = tol3
            if os.path.exists("Tolerances/atom{}lat{}rel{}".format(tol,tol2,tol3)) is False:
                os.mkdir("Tolerances/atom{}lat{}rel{}".format(tol,tol2,tol3))
            hits = SYMMETRY_DETECTION_tolerances.main(tol, tol2, tol3)
            file.writerow([tol,tol2,tol3,hits])




