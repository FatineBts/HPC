mkdir /tmp/3302011

make

./bin/shalw -x 512 -y 512 -t 1000 --export --export-path /tmp/3302011

./visu.py /tmp/3302011/shalw_512x512_T1000.sav