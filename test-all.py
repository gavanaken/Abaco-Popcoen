import os, sys, subprocess, time

execs = []

if len(sys.argv) > 1:
    l = sys.argv[1]
    json = '{ ' + '"debug":"true"' + ' }'
    os.system("abaco run -m \'{0}\' {1}".format(json, l))
else:
    rep=os.popen('abaco ls').read()
    for l in rep.split('\n'):
        if 'popcoen-stability2' in l.split(' '):
            json = '{ ' + '"debug":"true"' + ' }'
            a = l.split('  ')[-2]
            res = subprocess.check_output(["abaco", "run", "-m", json, a])
            execs.append((a, res.split('\n')[0]))

    time.sleep(20)
    for a, e in execs:
        print(a,e)
        os.system('abaco logs {0} {1}'.format(a, e))


