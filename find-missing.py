import sys
import os

datafile = sys.argv[1]
data = open(datafile, 'r').read()
tar = sys.argv[2]

output =[]

toc = os.popen('tar -tvf {0}'.format(tar)).read()
# Make sure tar exists and works
if 'tar: command not found' in  toc:
    print('Critical error. tar could not be located in current environment')
    sys.exit(0)
else:
    for l in toc.split('\n'):
        # perms owner  size date filename
        pdb = l.split(' ')[-1]
        # just in case there are weird files ...
        if pdb.endswith('.pdb'):
            if pdb.split('/')[-1].replace('.pdb','') not in data:
                output.append(pdb)
if os.path.exists('missing.tmp'):
    os.remove('missing.tmp')

with open('missing.tmp', 'w') as missing:
    missing.write('\n'.join(output))

