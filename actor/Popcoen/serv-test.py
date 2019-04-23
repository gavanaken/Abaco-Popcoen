import subprocess, sys, time, os


dir = sys.argv[1]
total = 0
n = 0 
FNULL = open(os.devnull, 'w')
for filename in os.listdir(dir):
    if n > 100:
        break
    if filename.endswith('.pdb'):
        start_time = time.time()
        pdb = open(os.path.join(dir,filename), 'r').read()
        p = subprocess.Popen(['python2', 'client_for_popcoen.py'], stdin=subprocess.PIPE, stdout=FNULL)
        res = p.communicate(input=pdb)[0]
        n += 1
        total += (time.time() - start_time)
        print("--- %s seconds ---" % (time.time() - start_time))
        print("%s average" % (total/n))
