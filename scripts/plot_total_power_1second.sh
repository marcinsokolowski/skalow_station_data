#!/bin/bash

dada_file=channel_0_1_1681955312.352758.dada

# 28935*32*1.08 = 999993.60 usec ~= 1 second
echo "skalow_spectrometer $dada_file -f test -p 0 -C 1 -c 0 -s 4096 -Z  -m -1 -F 200 -N 32 -O dynspec -a 28935  -A total_power_1sec/ -v -v -v "
skalow_spectrometer $dada_file -f test -p 0 -C 1 -c 0 -s 4096 -Z  -m -1 -F 200 -N 32 -O dynspec -a 28935  -A total_power_1sec/ -v -v -v   

cd total_power_1sec/
awk '{print $1" "$2;}' total_power_out.txt > total_power_2col_X.txt
awk '{print $1" "$3;}' total_power_out.txt > total_power_2col_Y.txt
awk '{print $1" "$4;}' total_power_out.txt > total_power_2col_I.txt


echo "python ./plot_power_vs_time.py total_power_2col --comment="Total power on Sun" $plot_options"
python ./plot_power_vs_time.py total_power_2col --comment="Total power on Sun" $plot_options



