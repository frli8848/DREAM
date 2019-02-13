#!/bin/bash

# On supercalc
#export LD_LIBRARY_PATH=/disk2/users/staff/fl/usr/local/lib
#/disk2/users/staff/fl/usr/local/bin/runacovea -config gcc34_pentium4.acovea -bench ./test_dreamrect.c
#/disk2/users/staff/fl/usr/local/bin/runacovea -config gcc33_pentium4.acovea -bench ./test_dreamrect.c

#runacovea -config gcc34_pentiumM_dream.acovea4 -bench ./test_dreamrect.c
#runacovea -config gcc34_pentium_dream.acovea5 -input ./test_dreamrect.c
#runacovea -g 80 -config gcc34_sun32.acovea5 -input ./test_dreamrect.c
#runacovea -g 80 -config gcc34_fm.acovea5 -input ./test_dreamrect.c

# P4
#runacovea -g 80 -config gcc34_pentium4_dream.acovea5 -input ./test_dreamrect.c
#runacovea -g 80 -config gcc34_pentium4_dream_fm.acovea5 -input ./test_dreamrect.c

# Prescott
#runacovea -g 80 -config gcc34_prescott_dream.acovea5 -input ./test_dreamrect.c
#runacovea -g 80 -config gcc34_prescott_dream_fm.acovea5 -input ./test_dreamrect.c

# PM
#runacovea -g 80 -config gcc34_pentiumM_dream.acovea5 -input ./test_dreamrect.c
#runacovea -g 80 -config gcc34_pentiumM_dream_fm.acovea5 -input ./test_dreamrect.c

# P4-M
#runacovea -g 80 -config gcc34_pentium4m_dream.acovea5 -input ./test_dreamrect.c
#runacovea -g 80 -config gcc34_pentium4m_dream_fm.acovea5 -input ./test_dreamrect.c

# P-M
#runacovea -g 80 -config gcc34_pentium4m_dream.acovea5 -input ./test_dreamrect.c
runacovea -g 80 -config gcc411_pentiumM.acovea -input ./test_dreamrect.c




