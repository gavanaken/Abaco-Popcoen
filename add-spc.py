import sys, os, time

def parse_args():
    actorIds = []
    filelist = []
    library = tarPath = tarSys = dsspPath = dsspSys = ''
    recurse = False
    # iterate over and look for flags
    for i in range (1,len(sys.argv)):
        if sys.argv[i] == "--tarPath":
            tarPath = sys.argv[i+1]
        if sys.argv[i] == "--tarSystemId":
            tarSys = sys.argv[i+1]
        if sys.argv[i] == "--dsspPath":
            dsspPath = sys.argv[i+1]
        if sys.argv[i] == "--dsspSystemId":
            dsspSys = sys.argv[i+1]
        if sys.argv[i] == "--library":
            library = sys.argv[i+1]
        if sys.argv[i] == "--actorIds":
            j = i + 1
            while j <= len(sys.argv)-1 and not sys.argv[j].startswith('-'):
                actorIds.append(sys.argv[j])
                j+= 1
        if sys.argv[i] == "-r":
            recurse = True
        if sys.argv[i] == "--filelist":
            filelist = open(sys.argv[i+1],'r').read().split('\n')
    return tarPath, tarSys, dsspPath, dsspSys, library, actorIds, recurse, filelist

def checkTar(tarPath):
    # Check local directory for the filename of the tar
    filename = os.path.split(tarPath)[-1]
    if os.path.exists(filename):
        print ('{0} found in local directory'.format(filename))
        return True
    else:
        print('{0} not found in local directory'.format(filename))
        return False

def getTar(tarPath, tarSys):
    # Use sd2e-cloud-cli to get the tar
    print ('Attempting to download {0} from {1}...'.format(tarPath, tarSys))
    os.system('files-get -S {0} {1}'.format(tarSys, tarPath))
    # Wait for it to download and ensure it did
    #os.wait()
    if not checkTar(tarPath):
        if raw_input('Failed to download. Try again? [y/n]').lower().startswith('y'):
            getTar(tarPath, tarSys)
        else:
            print('Critical error. Could not download pdbs.')
            sys.exit(0)
    else:
        return

def parseTar(tarPath):
    pdbFiles = []
    # Run table of contents on the tar
    filename = os.path.split(tarPath)[-1]
    toc = os.popen('tar -tvf {0}'.format(filename)).read()
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
                pdbFiles.append(pdb)
    return pdbFiles

def runActor(tarPath, tarSys, dsspPath, dsspSys, library, actor, pdbset):
    print("running actor {0} with {1} pdbs".format(actor, str(len(pdbset))))
    pdbSplits=list(pdbset[i:i+(100)] for i in xrange(0, len(pdbset), 100))
    resp = ''
    for subset in pdbSplits:
        fP = '"filePaths":"{0}"'.format(','.join(subset))
        l = '"library":"{0}"'.format(library)
        sP = '"tarPath":"{0}"'.format(tarPath)
        tS = '"tarSystemId":"{0}"'.format(tarSys)
        dP = '"dsspPath":"{0}"'.format(dsspPath)
        dS = '"dsspSystemId":"{0}"'.format(dsspSys)
        json ='{ ' + ','.join([fP,l,sP,tS,dP,dS]) + ' }'
        #print('abaco run -m \'{0}\' {1}'.format(json, actor))
        resp = resp + os.popen('abaco run -m \'{0}\' {1}'.format(json, actor)).read().split('\n')[0]
    return resp.split('\n')[0]

def collectOutput(executions, jobs):
    t = 0
    data = ''
    finished={}
    failed = []
    if os.path.exists('executions.tmp'):
        os.remove('executions.tmp')
    if os.path.exists('complete.tmp'):
        os.remove('complete.tmp')
    with open('executions.tmp', 'w') as log:
        for k in executions:
            for e in executions[k]:
                log.write('{0} {1}\n'.format(k, e))
                finished[e] = False
    outputfile = 'versioned_data-{0}.csv'.format(time.strftime("%Y%m%d-%H%M%S"))
    while False in finished.values():
        time.sleep(120)
        t += 120
        if t%12000==0:
            print("refreshing authentication tokens")
            os.system('auth-tokens-refresh')
        tick = 0
        for k in executions:
            for e in executions[k]:
                tick += 1
                log = os.popen('abaco logs {0} {1}'.format(k, e)).read()
                if 'Execution finished successfully' in log.split('\n')[1] and finished[e] == False:
                    if data == '':
                        data = log.split('\n')[2]
                        add = '\n'.join(log.split('\n')[2:-5])
                        with open(outputfile, 'a') as output:
                            output.write(add)
                    else:
                        add = '\n'.join(log.split('\n')[3:-5])
                        with open(outputfile, 'a') as output:
                            output.write(add)
                    finished[e] = True
                    print('{0} {1} complete\n'.format(k,e))
                    with open('complete.tmp', 'a') as complete:
                        complete.write('{0} {1}\n'.format(k, e))
                elif log.split('\n')[1] != "" and 'Execution finished successfully' not in log.split('\n')[1] and finished[e] == False:
                    print('execution failed: %s %s' % (k,e))
                    failed.append(','.join(jobs[e]))
                    finished[e] = True
                else:
                    sys.stdout.write('waiting for {0} executions to finish'.format(str(list(finished.values()).count(False))) + '..'*(tick//10) +'\r')

    return outputfile, failed


def runExecutions(pdbFiles, actorIds, tarPath, tarSys, dsspPath, dsspSys, library):
    k = len(actorIds)
    n = len(pdbFiles)
    bundle=list(pdbFiles[i:i+(n//k)] for i in xrange(0, n, n//k))
    executions={}
    jobs={}
    for a in actorIds:
        executions[a] = []
    for a,b in zip(actorIds, bundle):
        l = 100
        for j in list(b[i:i+l] for i in xrange(0,len(b),l)):
            exId = runActor(tarPath, tarSys, dsspPath, dsspSys, library, a, j)
            executions[a].append(exId)
            jobs[exId] = j
    outputfile, failed = collectOutput(executions, jobs)
    return outputfile, failed

def main():
    start = time.time()
    print("refreshing authentication tokens")
    os.system('auth-tokens-refresh')
    if os.path.exists("versioned_data.csv"):
        print("removing contents of 'versioned_data.csv'")
        os.remove("versioned_data.csv")
    tarPath, tarSys, dsspPath, dsspSys, library, actorIds, recurse, filelist = parse_args()
    localTar = checkTar(tarPath)
    if not localTar:
        getTar(tarPath, tarSys)
    if len(filelist) == 0:
        pdbFiles = parseTar(tarPath)
    else:
        pdbFiles = filelist
    output, failed = runExecutions(pdbFiles, actorIds, tarPath, tarSys, dsspPath, dsspSys, library)
    if recurse:
        while len(failed) != 0:
            print("re-running failed executions")
            print(failed)
            output, failed = runExecutions(failed, actorIds, tarPath, tarSys, dsspPath, dsspSys, library)
    print('output written to {0}'.format(output))
    stop = time.time()
    print("TIME ELAPSED:")
    print(stop - start)


if __name__ == '__main__':
    main()
