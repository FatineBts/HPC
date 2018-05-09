mkdir /tmp/3302011

make

export OMP_NUM_THREADS=2

mpirun -n 2 ./bin/shalw --export --export-path /tmp/3302011

./visu.py /tmp/3302011/shalw_256x256_T1000.sav