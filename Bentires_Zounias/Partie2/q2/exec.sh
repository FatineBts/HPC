mkdir /tmp/3302011

make

export OMP_NUM_THREADS=2

mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw -x 8192 -y 8192 -t 20 

./visu.py /tmp/3302011/shalw_8192x8192_T20.sav