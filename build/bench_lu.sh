#!/bin/bash
base="Bench/lu"
for ((nCells=100;nCells<=8000;nCells+=100))
do
  setup="$base/cellNumber_$nCells.in"
  echo $setup
  ./loki $setup >> Bench/bench_lu.out
done
