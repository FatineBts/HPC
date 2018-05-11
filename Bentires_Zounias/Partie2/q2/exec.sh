#mkdir /tmp/3674377

make

export OMP_NUM_THREADS=16

mpirun -n 4 ./bin/shalw -x 8192 -y 8192 -t 20

#./visu.py /tmp/3674377/shalw_512x512_T40.sav