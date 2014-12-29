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
import stat
import sys
import time
import json
import glob
import shutil
import subprocess

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

def get_empty_slot(slot_status):
    for i in range( len(slot_status) ):
        if slot_status[i]:
            return i
    return None

def submit_job(slot_id, my_robot, my_job):
    print 'submit job in', slot_id
    for file in glob.glob(str(my_job["job"])+"/*.in"):
        shutil.copy(file, my_robot["scratch"]+"/slot"+str(slot_id))
    submit_shell = open("submit.sh", "w")
    print >> submit_shell, "nohup " + my_job["exe"] + "> job.log &"
    submit_shell.close()
    st = os.stat('submit.sh')
    os.chmod('submit.sh', st.st_mode | stat.S_IEXEC)
    shutil.move('submit.sh', my_robot["scratch"]+"/slot"+str(slot_id))
    subprocess.Popen(["./submit.sh"], shell=True, cwd=my_robot["scratch"]+"/slot"+str(slot_id))
    my_job["status"] = 1

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

    rescan = True
    num_running_josb = 0
    slots_status = [True] * my_robot["slots"]

# main loop
    while True:

# step: parse the jobs_file
        if rescan:
            my_jobs = parse_jobs_json(my_robot["jobs"])
            for i in range( len(my_jobs) ):
                my_jobs[i]["status"] = 0 # 0: to do; 1: running; 99 : finished
            rescan = False

# step: try to submit jobs
        for i in range( len(my_jobs) ):
            if my_jobs[i]["status"] == 0:
                slot_id = get_empty_slot(slots_status)
                if slot_id is not None:
                    print slot_id
                    slots_status[slot_id] = False
                    submit_job(slot_id, my_robot, my_jobs[i])

# step: update the status of rescan
        rescan = True
        for i in range( len(my_jobs) ):
            if my_jobs[i]["status"] != 99:
                rescan = False

# step: determine whether we should exit
        if os.path.exists(my_robot["finish"]) and os.path.isfile(my_robot["finish"]):
            print "found", my_robot["finish"]
            break
