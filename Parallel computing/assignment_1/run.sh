make -f Makefile clean

make
rm -r test
mkdir test

for stencil in 5 9
do
	for N in 262144 4194304
       	do
		 	
        	for execution in 1 2 3
        	do	
			
                       	mpirun -np 12 ./prob1.x 4 $N 10 2 $stencil >> ./test/test_N${N}_stencil${stencil}
                       	
        	done
	done
done

