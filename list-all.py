import os

rep=os.popen('abaco ls').read()
res = []
for l in rep.split('\n'):
    if 'popcoen-stability2' in l.split(' '):
        res.append(l.split('  ')[-2])

print(' '.join(res))


