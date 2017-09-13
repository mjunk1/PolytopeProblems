#!/bin/bash

for l in {10..100}; do
	for cd in `seq 5 5 50`; do
		mkdir data/2Q/${l}_${cd}
		./polytope -q 2 -f polytope_twirled_constraints2Q.dat -o data/2Q/${l}_${cd}/np -N 0.3 -n 0.01 -i 1000 -l $l -p `echo "scale=2; ${cd}/100" | bc` &> ${l}_${cd}.log 
	done
done