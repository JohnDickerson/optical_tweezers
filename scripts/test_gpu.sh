#!/bin/bash
##
## test_gpu.sh
## 
## Made by Robert Patro
## Login   <rob@gvilwks03>
## 
## Started on  Wed Dec  2 15:06:26 2009 Robert Patro
## Last update Wed Dec  2 15:06:26 2009 Robert Patro
##

start=$1
stop=$2
values=$3

i=$start
while [ $i -le $stop ]; do
    ./exp_tweezers $i >> $values 
    i=$(( $i*2 ));
done
