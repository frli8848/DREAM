#!/bin/bash

rm test_dreamrect
rm test_dreamrect.o
make -f Makefile_p4 test_dreamrect
./test_dreamrect
