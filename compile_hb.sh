#!/bin/bash
for i in {27..34}
do
	make clean
	./configure --with-hash-bits=$i
	make CFLAGS="-Ofast -mcmodel=large -g3"
	mv barriers barriers_hb$i
done
