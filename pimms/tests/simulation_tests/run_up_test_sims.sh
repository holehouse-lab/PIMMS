#!/sh
for i in $(seq 1 13)
do
	
	if [ -d test_${i} ]
	then
		cd test_${i}
		rm *pdb 2>/dev/null
		rm *xtc 2>/dev/null
		rm *dat 2>/dev/null
		rm restart.pimms 2>/dev/null
		rm parameters_used.prm 2>/dev/null
		rm absolute_energies_of_angles.txt 2>/dev/null
		rm log.txt 2>/dev/null
		cd ..
	fi

done