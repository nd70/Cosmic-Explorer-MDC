#!/bin/bash

HOME=/home/sgwynne.crowder
export HOME

#source /archive/home/shivaraj/.bashrc
#source /home/s.gwynne.crowder/.bashrc
#source /archive/home/shivaraj/sgwb/S5/matlab/matlab_script_64bit.sh
#source /home/sgwynne.crowder/checkOut/trunk/S5/matlab/matlab_script_64bit.sh
source /home/sgwynne.crowder/matlab_scripts/matlab_script_2013a.sh

hostname > /home/sgwynne.crowder/logs/condor_$3.info

#/archive/home/shivaraj/sgwb/S5/command/stochastic $1 $2 $3
/home/sgwynne.crowder/bin/matlab/stochastic $1 $2 $3