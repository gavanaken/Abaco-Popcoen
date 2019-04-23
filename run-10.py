import os

run = 'python add-spc.py --tarPath protein-design/data_v1_April_2018/Eva1_Dec_2017/pdbs.tar.gz --tarSystemId data-sd2e-community --dsspPath protein-design/data_v1_April_2018/aggregated_v1_data_Hamed_May_23_2018/data_v1_aggregated.csv --dsspSystemId data-sd2e-community --library Eva1 --actorIds'
id = []
rep=os.popen('abaco ls').read()
for l in rep.split('\n'):
    if 'popcoen-stability2' in l.split(' '):
        id.append(l.split('  ')[-2])
print(id)
print(run + ' '.join(id))
os.popen(run + ' '.join(id))


