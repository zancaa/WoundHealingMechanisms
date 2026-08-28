#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt ts=projects/WoundHealingMethods/test/TestWoundPurseCrawlingProlif.hpp
#

# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
# ">" directs std::cout to the file.
# "2>&1" directs std::cerr to the same place.
# "&" on the end lets the script carry on and not wait until this has finished.
nice -20 ../build/optimised/TestWoundPurseCrawlingProlifRunner -num_param_vals 1 -num_prolif_sims 3 > output/Wound_output.txt 2>&1 &

echo "Jobs submitted"