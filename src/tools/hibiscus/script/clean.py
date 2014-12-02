import os

for path, dirs, files in os.walk('.'):
  for f in files:
    if f.endswith('DS_Store'):
        print path
        print f
