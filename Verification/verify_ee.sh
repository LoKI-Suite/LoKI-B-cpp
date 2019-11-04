#!/bin/bash
base="../Verification/input/ee"
for ((nCells=1000;nCells<=10000;nCells+=1000))
do
  setup="$base/$nCells.in"
  echo $nCells
  ../build/loki $setup > /dev/null
done
