#!/bin/bash --login

for i in `seq 4 4`
do
	N=$i && ./enum_edge.py $N > edge_$N.dat && ./plot_edge.py $N edge_$N.dat edge_$N.pdf
done
