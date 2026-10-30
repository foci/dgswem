#!/bin/bash
cd /home/wukenton/dgswem/tests/slopingbeach
for i in {1..10}; do
  ../../build-debug/dagswem > /dev/null 2>&1
  printf "post\nN\nquit\n" | ../../build-debug/adcpost > /dev/null 2>&1
  echo -n "Run $i: "
  md5sum fort.63
done
