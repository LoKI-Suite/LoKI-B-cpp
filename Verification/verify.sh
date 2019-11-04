#!/bin/bash
base="../Verification/input/growth"
for ((nCells=6000;nCells<=10000;nCells+=100))
do
  setup="$base/$nCells.in"
  echo $nCells
  ../build/loki $setup > /dev/null
done
