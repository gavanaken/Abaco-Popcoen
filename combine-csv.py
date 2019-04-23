import sys, os

parent = sys.argv[1]
child = sys.argv[2]
output = sys.argv[3]

data = open(parent, 'r').read()
add = '\n'.join(open(child, 'r').read().split('\n')[1:])

with open(output, 'w') as ofile:
    ofile.write(data + add)
