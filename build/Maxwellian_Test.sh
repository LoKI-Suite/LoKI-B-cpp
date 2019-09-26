#!/bin/bash
base="Model/Maxwellian"
for ((nCells=100;nCells<=10000;nCells+=100))
do
  path="$base/cellNumber_$nCells"
  for ((i=1;i<=51;i++))
  do
    setup="$path/setup_$i.in"
    echo $setup
    ./loki $setup
  done
done
