#!/bin/bash

for i in {2..10}
do
  mpiexec -np "$i" ./main.exe
done
