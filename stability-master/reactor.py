import sys, os, time, tarfile, json
from reactors.utils import Reactor
from agavepy.agave import Agave

class Args:

    def __init__(self):
        payload = r.context.message_dict
        # The assay library name
        self.library = payload["library"]
        # The compressed file of PDBs...
        self.tarPath = payload["tarPath"]
        # ...and the systemId on which it resides
        self.systemId = payload["tarSystemId"]
        # What workers do we have available to us?
        self.workers = payload["workers"]
        # The current uid is master
        self.master = r.uid
        # The current execid is masterExecId
        self.masterExecId = r.execid

def parse_opt_args():
    return {}

def pull_tar(tarPath, systemId):
    #   Download the tar and save it to a temporary file
    files_response = r.client.files.download(filePath=tarPath, systemId=systemId)
    filename = 'temp_pdbs.tar.gz'
    with open(filename, 'wb') as tempfile:
        for chunk in files_response.iter_content(chunk_size=1024):
            if chunk:
                tempfile.write(chunk)
    return filename

def parse_tar(temp_pdbs_tar_file):
    pdbFiles = []
    temp_pdbs_tar = tarfile.open(temp_pdbs_tar_file, 'r')
    for pdb in temp_pdbs_tar.getmembers():
        # In case there are weird file types
        if pdb.name.endswith('.pdb'):
            pdbFiles.append(pdb.name)
    return pdbFiles

def make_bundles(seq, num):
    avg = len(seq) / float(num)
    out = []
    last = 0.0
    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg
    return out

def run_actor(args, optArgs, w, j):
    message = {}
    message["type"] = 'request'
    message["filePaths"] = ','.join(j)
    message["library"] = args.library
    message["tarPath"] = args.tarPath
    message["tarSystemId"] = args.systemId
    message["masterId"] = args.master
    message["masterExecId"] = args.masterExecId
    if 'dsspPath' in optArgs:
        message["dsspPath"] = optArgs['dsspPath']
        message["dsspSystemId"] = optArgs['dsspSystemId']
    jsonData =json.dumps(message)
    resp = r.client.actors.sendMessage(actorId=w, body=jsonData)
    return resp["executionId"]

def write_manifest(executions):
    print("Writing manifest")
    man = ''
    for w, e in executions:
        man = man + '{0} {1}\n'.format(w,e)
    print(man)
    body = {'{0}-manifest'.format(r.execid): man}
    r.client.actors.updateState(actorId=r.uid, body=body)


def run_executions(pdbFiles, args, optArgs):
    print("Running executions")
    k = len(args.workers)
    n = len(pdbFiles)
    if "execSize" in optArgs:
        execSize = int(optArgs["execSize"])
    else:
        execSize = 100
    bundle=make_bundles(pdbFiles, k)
    execution_list = []
    for w,b in zip(args.workers, bundle):
        l = -(-len(bundle)//execSize)
        jobs = make_bundles(b, l)
        for j in jobs:
            exId = run_actor(args, optArgs, w, j)
            execution_list.append((w,exId))
    write_manifest(execution_list)

def request_main():
    args = Args()
    optArgs = parse_opt_args()
    temp_pdbs_tar = pull_tar(args.tarPath, args.systemId)
    # User-specified subset or entire library
    if "filePaths" not in optArgs:
        pdbFiles = parse_tar(temp_pdbs_tar)
    run_executions(pdbFiles, args, optArgs)
    os.remove(temp_pdbs_tar)

def success_main():
    payload = r.context.message_dict
    state1=r.client.actors.getState(actorId=r.uid)
    key = '{0}-output'.format(payload["masterExecId"])
    if key not in state1:
        body = {key: output}
    else:
        body = {key: state1[key] + '\n' + '\n'.join(output.split('\n')[1:])}
    r.client.actors.updateState(actorId=r.uid, body=body)
    print(r.client.actors.getState(actorId=r.uid))

def query_main():
    payload = r.context.message_dict
    execId = payload["execId"]
    state1=r.client.actors.getState(actorId=r.uid)
    print(state1['state']["{0}-output".format(execId)])

def main():
    global r
    r = Reactor()
    payload = r.context.message_dict
    if payload["type"] == "request":
        request_main()
    elif payload["type"] == "response":
        if payload["response"] == "success":
            success_main()
        else:
            failure_main()
    elif payload["type"] == "query":
        query_main()

if __name__ == '__main__':
    main()


