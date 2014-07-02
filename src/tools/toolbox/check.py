#!/usr/bin/env python

import sys

argu = sys.argv[1:]

i=0
f = file(argu[0])
while True:
    line = f.readline()
    i=i+1
    if len(line) == 0:
        break

#    if len(line) > 75:
#        print 'check line:'
#        print i, line

    if len(line) != len(line.rstrip()) + 1:
        print 'check space:'
        print i, line,
f.close()
