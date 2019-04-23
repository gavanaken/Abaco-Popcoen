import os

execs = open('executions.tmp').read().split('\n')
for l in execs:
    if l != '':
        os.system('abaco workers {0}'.format(l))


