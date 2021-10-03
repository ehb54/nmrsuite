#!/usr/bin/python

import os, sys

# Path to be created
myp = "/opt/genapp/rotdif/bin/test/hour"

if os.path.isdir(myp):
    print("Existed")
else:
    os.mkdir(path)
    print("Path is created")

print(os.getcwd())
def hi():
    os.chdir(myp)
    print(os.getcwd())
    with open('hello.txt','w') as in_f:
        in_f.write("yes")
hi()
print(os.getcwd())
sys.stderr.write("Hey here")
