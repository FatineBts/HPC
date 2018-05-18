#mkdir /tmp/3302011

make

export OMP_NUM_THREADS=1

mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw -x 8192 -y 8192 -t 20

#./visu.py /tmp/3674377/shalw_256x256_T1000.sav