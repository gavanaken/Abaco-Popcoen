import sys
csv = open(sys.argv[1], 'r').read()
lines = csv.split('\n')[1:]
for i in range(len(lines)):
    for j in range(i+1,len(lines)):
        if lines[i] == lines[j]:
            print ("error! {0} in lines {1} and {2}".format(lines[i],str(i),str(j)))
