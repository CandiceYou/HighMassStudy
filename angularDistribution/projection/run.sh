source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.09/x86_64-slc6-gcc47-dbg/root/bin/thisroot.sh

#for i in {0..6}
for i in {0..1}
do
#root -q -b "angularDistributions_spin2_mh750_2e2mu.C(2,0,0)"
root -q -b "angularDistributions_spin2_mh750_2e2mu.C($i,0,0)"
done
