#!/bin/bash
base="Bench/hessenberg"
for ((nCells=100;nCells<=8000;nCells+=100))
do
  setup="$base/cellNumber_$nCells.in"
  echo $setup
  ./loki $setup >> Bench/bench_hess.out
done
