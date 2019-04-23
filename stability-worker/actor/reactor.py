from reactors.utils import Reactor
import subprocess, os, urllib2, time, re, tarfile, shutil, sys, json
import imp
import pandas as pd
import numpy as np
from StringIO import StringIO
from collections import OrderedDict
from agavepy.agave import Agave

def processResult(popcoenData):
    S_PC = 0
    reliability_n = 0
    Si_list = []
    for line in popcoenData.split('\n'):
        if line[:10]=='Si-predict':
            Si_value=line[22:29]
            new_Si=Si_value.strip()
            Si_list.append(new_Si)

        if line[:4] == 'S_PC':
            S_PC = line[15:-1].strip()

        if line[:25] == 'Reliability-number lambda':
            reliability_n=line[27:-1]
    return (S_PC, Si_list)

# def simpleOutput(entropyDict, library):
#     newOutput = "dataset,name,S_PC\n"
#     for e in entropyDict:
#         dir, name = e.split('/')
#         S_PC, SiVals = entropyDict[e]
#         newOutput = newOutput + ','.join([library,name,str(S_PC)]) + '\n'
#     return newOutput

def parseFiles(files_response):
    #this is a json object which houses the paths for all the files we want
    paths = []
    for d in files_response[1:]:
        paths.append(d['path'])
    return paths

def globals():
    AgaveAPI = Reactor().client

def add_zeroes(dssp_output, ss):
    dssp_output['Mean_{0}_entropy'.format(ss)].append(0)
    dssp_output['Sum{0}_entropies'.format(ss)].append(0)
    dssp_output['{0}_max_entropy'.format(ss)].append(0)
    dssp_output['{0}_min_entropy'.format(ss)].append(0)
    dssp_output['{0}_range_entropy'.format(ss)].append(0)
    return dssp_output

def add_values(dssp_output, ss, siVals):
    dssp_output['Mean_{0}_entropy'.format(ss)].append(np.mean(siVals))
    dssp_output['Sum{0}_entropies'.format(ss)].append(np.sum(siVals))
    dssp_output['{0}_max_entropy'.format(ss)].append(max(siVals))
    dssp_output['{0}_min_entropy'.format(ss)].append(min(siVals))
    dssp_output['{0}_range_entropy'.format(ss)].append(max(siVals) - min(siVals))
    return dssp_output

def addDssp(dssp_csv, entropyDict, library):
    dssp_output = OrderedDict()
    columns = ['dataset','name','S_PC','Mean_H_entropy','Mean_L_entropy',
                                'Mean_E_entropy','Mean_res_entropy','SumH_entropies','SumL_entropies',
                                'SumE_entropies','H_max_entropy','H_min_entropy','H_range_entropy',
                                'L_max_entropy','L_min_entropy','L_range_entropy','E_max_entropy',
                                'E_min_entropy','E_range_entropy']
    for c in columns:
        dssp_output[c] = [] # all init to empty string
#    i = 0
#    for l in dssp_csv.split('\n'):
#        if 'dataset' in l or 'name' in l or 'dssp' in l:
#            break
#        i+=1
    # Get rid of header
#    dssp_csv = '\n'.join(dssp_csv.split('\n')[i:])
#    dssp_csv = '\n'.join(dssp_csv.split('\n')[11:])
    dssp_df = pd.read_csv(StringIO(dssp_csv), sep=',', encoding="utf-8-sig")
    dssp_df.set_index('name',inplace=True)
    for k in entropyDict:
        try:
            dir, name = k.split('/')
            dssp_data = dssp_df.loc[name] # see if the name matches
            S_PC = entropyDict[k][0]
            si_values = [float(s) for s in entropyDict[k][1]]
            sequence = dssp_data['sequence']
            dssp = list(dssp_data['dssp']) # cast it to a list
            # ASSERT THAT len(si_values) == len(dssp)
        except: #if we couldn't find the key in dssp, just skip it
            continue
        L_s = []
        H_s = []
        E_s = []
        for ss, si in zip(dssp, si_values):
            #print(ss, si)
            if ss == 'L':
                L_s.append(float(si))
            elif ss == 'H':
                H_s.append(float(si))
            elif ss == 'E':
                E_s.append(float(si))
        dssp_output['dataset'].append(library)
        dssp_output['name'].append(name)
        dssp_output['S_PC'].append(S_PC)
        dssp_output['Mean_res_entropy'].append(np.mean(si_values))
        if len(L_s) == 0:
            dssp_output = add_zeroes(dssp_output, 'L')
        else:
            dssp_output = add_values(dssp_output, 'L', L_s)
        if len(H_s) == 0:
            dssp_output = add_zeroes(dssp_output, 'H')
        else:
            dssp_output = add_values(dssp_output, 'H', H_s)
        if len(E_s) == 0:
            dssp_output = add_zeroes(dssp_output, 'E')
        else:
            dssp_output = add_values(dssp_output, 'E', E_s)
    df = pd.DataFrame(data=dssp_output, columns = dssp_output.keys())
    return df.to_csv(index=False)

def send_response(message, actorId, r):
    data = {}
    data["type"] = "response"
    data["response"] = "success"
    data["masterExecId"] = r.context.message_dict["masterExecId"]
    data["sender_id"] = r.uid
    data["sender_execution"] = r.execid
    data["payload"] = message
    json_data = json.dumps(data)
    resp = r.client.actors.sendMessage(actorId=actorId, body=json_data)
    print('response sent to {0} {1}'.format(actorId, resp["executionId"]))


def parse_message(payload):
    # Who is our master?
    master = payload["masterId"]
    # What was the execId?
    masterExecId = payload["masterExecId"]
    # A comma separated string of filepaths: any/number/of/leading.../filename.pdb
    filePaths = payload["filePaths"]
    # The assay library name
    library = payload["library"]
    # The compressed file of PDBs...
    tarPath = payload["tarPath"]
    # ...and the systemId on which it resides
    systemId = payload["tarSystemId"]

    return master, tarPath, systemId, filePaths, library

def pull_dssp(payload, library, r):
    # User may have specified a path for DSSP...
    if 'dsspPath' in payload:
        dsspPath = payload["dsspPath"]
        dsspId = payload["dsspSystemId"]
    # Or they may not have, which means it is probably here...
    else:
        dsspPath = 'versioned-dataframes/protein-design/metadata/{0}.metadata.csv'.format(library)
        dsspId = 'data-sd2e-community'
    # Download the rosetta info
    dssp_csv = r.client.files.download(filePath=dsspPath, systemId=dsspId).text
    data_lines = dssp_csv.split('\n')
    i = 0
    while 'dataset' not in data_lines[i] and 'dssp' not in data_lines[i]:
        i+=1
    return '\n'.join(data_lines[i:])

def pull_tar(tarPath, systemId, r):
    #   Download the tar and save it to a temporary file
    files_response = r.client.files.download(filePath=tarPath, systemId=systemId)
    filename = 'temp_pdbs.tar.gz'
    with open(filename, 'wb') as tempfile:
        for chunk in files_response.iter_content(chunk_size=1024):
            if chunk:
                tempfile.write(chunk)
    return filename

def request_main(r):
    start = time.time()
    popPath=os.getcwd()
    popServer = os.path.join(popPath, 'popcoen_server.py')
    popClient = os.path.join(popPath, 'client_for_popcoen.py')

    devnull = open(os.devnull)
    # Start the server
    serv = subprocess.Popen(['python2', popServer], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=devnull, shell=False)

    # Parse the message elements
    master, tarPath, systemId, filePaths, library = parse_message(r.context.message_dict)

    # Get topology data from TACC
    dssp_csv = pull_dssp(r.context.message_dict, library, r)

    # Save the tar locally
    temp_pdbs_tar = pull_tar(tarPath, systemId, r)

    # Split out the files we actually want to extract, these come from the table of contents anyway so their paths are correct
    entropyDict = {}
    time.sleep(200)
    fileList = filePaths.split(',')
    temp_pdbs_dir = 'temp_pdbs'
    os.mkdir(temp_pdbs_dir)
    tmpTar = tarfile.open(temp_pdbs_tar, 'r')
    # Iterate over TOC and check for files specified for this execution
    for pdb in tmpTar.getmembers():
        if pdb.name in fileList:
            # Save to our temporary directory
            tmpTar.extract(pdb, temp_pdbs_dir)

    t = 0
    for path in fileList:
        pdb = open(os.path.join(temp_pdbs_dir, path), 'r').read()
        t += 1
        p = subprocess.Popen(['python2', popClient], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=devnull, shell=False)
        res = p.communicate(input=pdb)[0]
        key = '/'.join(path.split('/')[-2:]).replace('.pdb','') # dir/protein
        entropyDict[key] = processResult(res) # (S_PC, [si, si, si, si])
    # Remove the temporary dir of pdb's
    shutil.rmtree(temp_pdbs_dir)
    # Remove the tar so nothing is left on the reactor
    os.remove(temp_pdbs_tar)
    # Process Sconf values according to topology
    output = addDssp(dssp_csv, entropyDict, library)
    stop = time.time()
    send_response(output, master, r)

def main():
    r = Reactor()
    if r.context.message_dict["type"] == "request":
        request_main(r)

if __name__ == '__main__':
    main()
