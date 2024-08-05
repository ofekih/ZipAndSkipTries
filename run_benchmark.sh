#!/usr/bin/env bash

for i in {1..1000}
do
	echo "Iteration $i"
	./bin/benchmark
done
