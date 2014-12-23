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
import time
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

def write_robot_info(cfg_robot):
    print "  Welcome to QUARK"
    print "  >>> Starting in " + time.asctime()
    print 
    print "  robot name        :", cfg_robot["robot"   ]
    print "  machine name      :", cfg_robot["machine" ]
    print "  queue system      :", cfg_robot["queue"   ]
    print "  scratch directory :", cfg_robot["scratch" ]
    print "  finish signal     :", cfg_robot["finish"  ]
    print "  log data          :", cfg_robot["logdata" ]
    print "  err data          :", cfg_robot["errdata" ]
    print "  time interval     :", cfg_robot["interval"]
    print "  number of slots   :", cfg_robot["slots"   ]
    print "  job lists         :", cfg_robot["jobs"    ]

if __name__ == '__main__':

# setup the robot_file we want to parse
    argu = sys.argv[1:]
    if ( len(argu) > 0 ):
        robot_file = argu[0]
    else:
        robot_file = './robot.json'

# parse the robot_file
    my_robot = parse_robot_json(robot_file)
    write_robot_info(my_robot)

# prepare the directory: scratch
    if not os.path.exists(my_robot["scratch"]):
        os.makedirs(my_robot["scratch"])

# prepare the directory: slots
    for i in range(my_robot["slots"]):
        if not os.path.exists(my_robot["scratch"] + "/slot" + str(i)):
            os.makedirs(my_robot["scratch"] + "/slot" + str(i))

# prepare the log files
    flog = open(my_robot["logdata"],"w")

# prepare the err files
    ferr = open(my_robot["errdata"],"w")

# main loop
#    while True:

# step 1: parse the jobs_file
#        my_jobs = parse_jobs_json(my_robot["jobs"])

# step 2: generate the sub-jobs
