#fileroot="$1"
xdrlibdir="/home/dillion/pkg/surf-master/depends/lib"
xdrincldir="/home/dillion/pkg/surf-master/depends/xdrfile-1.1.4/include"
#/opt/intel/bin/icpc -DCPLUSPLUS -I $xdrincldir -L $xdrlibdir -lxdrfile -o $fileroot.out $fileroot.cpp -Wall
#/opt/intel/bin/icpc -DCPLUSPLUS -I $xdrincldir -L $xdrlibdir -lxdrfile -o calcVRprot.out calcVRprot.cpp -Wall
/opt/intel/bin/icpc -DCPLUSPLUS -I $xdrincldir -L $xdrlibdir -lxdrfile -o interface_potential.out interface_potential.cpp -Wall
#/opt/intel/bin/icpc -DCPLUSPLUS -I $xdrincldir -L $xdrlibdir -lxdrfile -o protein_proIIcal.out protein_proIIcal.cpp -Wall
