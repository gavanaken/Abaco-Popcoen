from agavepy.agave import Agave
import sys, csv, os

ag = Agave.restore()
file = sys.argv[1]
with open(file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            line_count += 1
        else:
            filename = row[0]
            lib = row[1]
            print(lib, filename)
            if lib == 'Eva1':
                cmd = 'files-upload -S data-tacc-work-gvanaken -F {0}.pdb pdb-data/eva1/'.format(filename)
                print(cmd)
                os.system('files-upload -S data-tacc-work-gvanaken -F {0}.pdb pdb-data/eva1/'.format(filename))
            else:
                break
            line_count += 1

