import os
import sys
import subprocess

path = os.path.dirname(sys.argv[0])
mingw_make = "mingw32-make"

def get_each_file(files):
    for f in files:
        dirpath,_,names = f
        if len(names)==0: continue
        for name in names:
            yield (dirpath,name)

# make everything
for d,n in get_each_file(os.walk(path)):
    if n!="makefile": continue
    subprocess.run(args=[mingw_make],cwd=d)

# run it
for d,n in get_each_file(os.walk(path)):
    if not n.endswith(".exe"): continue
    file_path = os.path.join(d,n)
    args = [file_path]
    subprocess.run(args=args,cwd=d,creationflags=subprocess.ABOVE_NORMAL_PRIORITY_CLASS)