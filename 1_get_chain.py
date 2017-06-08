import sys
name = '1WZ2.pdb'

chain = 'C'
if len(sys.argv) > 1:
    name = sys.argv[1]
    chain = sys.argv[2]
f1 = open(name,'r')

f2=open('chain.pdb','w')
first_res=False
for line in f1:
 if len(line) >21:
    if line[21].upper() == chain and line[0:4] == 'ATOM':
        if first_res is False:
            first_res = int(line[22:26].strip())
        current_res = 1+int(line[22:26].strip()) - first_res
        line=line[0:22]+' '*(4-len(str(current_res)))+str(current_res)+line[26:]
        f2.write(line),

f2.write('TER\nEND')
f2.close()
