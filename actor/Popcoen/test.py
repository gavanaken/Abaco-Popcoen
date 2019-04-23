import subprocess, sys, time

stime = time.time()

pdb = open(sys.argv[1], 'r').read()
p = subprocess.Popen(['python2', 'popcoenMain.py'], stdin=subprocess.PIPE)
res = p.communicate(input=pdb)[0]
print(res)
print(time.time()-stime)
