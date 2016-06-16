 for j in 0 1 
do	 
 for k in 0 1
do
# for i in {0..6}
 for i in 5
do
 for decay in zz2l2nu #zz4q zz2l2nu zz2l2q  
do
root -q -b "angularDistributions_spin2_mh750_2e2mu.C($i,$j,$k,\"$decay\")"
#root -q -b "angularDistributions_spin2_mh750_4l.C($i,$j,$k)"
done
done
done
done
