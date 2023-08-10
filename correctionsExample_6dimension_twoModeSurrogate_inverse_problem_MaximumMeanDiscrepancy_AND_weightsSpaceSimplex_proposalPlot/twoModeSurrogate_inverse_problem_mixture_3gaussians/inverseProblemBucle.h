#! /bin/bash

i=0
while [ 1 -le 3 ] 
do 
  mkdir $i
  inverse
  ((i++))
done

exit
