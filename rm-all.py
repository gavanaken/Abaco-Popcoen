import os

rep=os.popen('abaco ls').read()
for l in rep.split('\n'):
    if 'popcoen-stability2' in l.split(' '):
        os.system('abaco rm {0}'.format(l.split('  ')[-2]))


