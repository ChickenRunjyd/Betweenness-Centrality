#! /bin/bash

for index in `seq 1 10`;
do
    make test
    sleep 5
done
