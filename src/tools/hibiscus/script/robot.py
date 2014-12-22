#!/usr/bin/env python

##
##
## Introduction
## ============
##
## It is a python script. The purpose of this script is to perform some
## automatic tests for the ctqmc quantum impurity solvers.
##
## This script should be used by the developer only.
##
## Usage
## =====
##
## ./robot.py robot_file
##
## Author
## ======
##
## This python script is designed, created, implemented, and maintained by
##
## Li Huang // email: huangli712@gmail.com
##
## History
## =======
##
## 12/20/2014 by li huang
##
##

import os
import sys
import json

def parse_robot_json(robot_json):
    json_data = {}
    cfg_robot = {}
    with open(robot_json) as json_file:
        json_data = json.load(json_file)
    cfg_robot = json_data.copy()
    return cfg_robot

def parse_jobs_json(jobs_json):
    json_data = {}
    cfg_jobs  = []
    with open(jobs_json) as json_file:
        json_data = json.load(json_file)
    cfg_robot = json_data["jobs"]
    return cfg_robot

if __name__ == '__main__':

# setup the robot_file we want to parse
    argu = sys.argv[1:]
    if ( len(argu) > 0 ):
        robot_file = argu[0]
    else:
        robot_file = './robot.json'

# parse the robot_file
    my_robot = parse_robot_json(robot_file)

# prepare the directory

# prepare the files

# main loop
    while True:

# step 1: parse the jobs_file
#        my_jobs = parse_jobs_json(my_robot["jobs"])

# step 2: generate the sub-jobs


