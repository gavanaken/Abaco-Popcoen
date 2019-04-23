#from reactors.utils import Reactor
import subprocess, os, urllib2, time, re, tarfile, shutil, sys, json
import imp
import pandas as pd
import numpy as np
from StringIO import StringIO
from collections import OrderedDict
#from agavepy.agave import Agave

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
    data["sender_id"] = r.uid
    data["sender_execution"] = r.execid
    data["payload"] = message
    json_data = json.dumps(data)
    r.client.actors.sendMessage(actorId=actorId, body=json_data)


def parse_message(payload):
    # Who is our master?
    #master = payload["masterId"]
    master = "WIP"
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

def main():
    start = time.time()
    library = 'test'
    # Split out the files we actually want to extract, these come from the table of contents anyway so their paths are correct
    entropyDict = {}
    dir = sys.argv[1]
    dssp_csv = open(sys.argv[2], 'r').read()
    total = 0
    n = 0
    FNULL = open(os.devnull, 'w')
    for filename in os.listdir(dir):
        print(" %s: %s" % (n, filename))
        if n > 100:
            break
        if filename.endswith('.pdb'):
            path = os.path.join(dir,filename)
            pdb = open(os.path.join(dir,filename), 'r').read()
            p = subprocess.Popen(['python2', 'client_for_popcoen.py'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            res = p.communicate(input=pdb)[0]
#            print(res)
            n += 1
        # for testing
            key = '/'.join(path.split('/')[-2:]).replace('.pdb','') # dir/protein
            entropyDict[key] = processResult(res) # (S_PC, [si, si, si, si])
    # Remove the temporary dir of pdb's
    # Remove the tar so nothing is left on the reactor
    # Process Sconf values according to topology
    output = addDssp(dssp_csv, entropyDict, library)
    # FOR TESTING
    print("Execution finished successfully")
    print(n)
    print(output)
    stop = time.time()
    print(stop-start)


if __name__ == '__main__':
    main()
