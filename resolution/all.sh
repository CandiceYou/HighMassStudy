for subdir in singleMassFit fixParamFit params 
do
mkdir $subdir
done

rm params/individual_fit_param.txt

#for ch in {2e2mu,4e,4mu,2l2q}
for ch in {2l2q,}
do
root -l -n -q -b readData.cc\(\"$ch\"\) simfit.cc\(\"$ch\"\)
bash readParam_single.sh $ch
root -l -n -b -q "plotVar_fitInd_$ch.C" 
bash readParam_ind.sh $ch
done
