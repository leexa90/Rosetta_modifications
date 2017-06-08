import sys
sys.path.append ('/cluster/apps/x86_64/packages/lib/python2.4/site-packages')
import numpy as np
KK =  1000
sysargv =20
if len(sys.argv) == 3:
    KK =  int(sys.argv[2])
    sysargv =int(sys.argv[1])
#pi_backbone = ['O3',
dict_atom = {" O4'": 77, " O2'": 77,  " C2'": 77, 
 " C4'": 77, " O5'": 77, " O3'": 77, 
 " C1'": 77, " C3'": 77, ' P  ': 77, " C5'": 77}
atom_keys = [" C1'", " C2'", " C3'", " C4'", " C5'", 
             " O2'", " O3'", " O4'", " O5'", ' P  ']

atom_keys = [' P  ']

atom_keys = [" C1'", " C2'", " C3'", " C4'", " C5'", 
             " O2'", " O3'", " O4'", " O5'", ' P  ']

dict_atom = { 'CYT':  [" N1 ", " C2 ", " C5 "],
              'URA':  [" N1 ", " C2 ", " C5 "],
              'ADE':  [" C4 ", " N9 ", " N7 "],
              'GUA':  [" C4 ", " N9 ", " N7 "] }

# important residues to be used
dict_atom = { 'CYT':  ['P',"C4'","N1", "C2", "C5"],
              'URA':  ['P',"C4'","N1", "C2", "C5"],
              'ADE':  ['P',"C4'","C4", "N9", "N7"],
              'GUA':  ['P',"C4'","C4", "N9", "N7"],
              
              'C':  ['P',"C4'","C2", "N1", "C5"],
              'U':  ['P',"C4'","C2", "N1", "C5"],
              'A':  ['P',"C4'","C4", "N9", "N7"],
              'G':  ['P',"C4'","C4", "N9", "N7"]
              }

dict_atom = { 'CYT':  ['P',"C2", "N1", "C5"],
              'URA':  ['P',"C2", "N1", "C5"],
              'ADE':  ['P',"C4", "N9", "N7"],
              'GUA':  ['P',"C4", "N9", "N7"],
              
              'C':  ['P',"C2", "N1", "C5"],
              'U':  ['P',"C2", "N1", "C5"],
              'A':  ['P',"C4", "N9", "N7"],
              'G':  ['P',"C4", "N9", "N7"],
              
              'RC':  ['P',"C2", "N1", "C5"],
              'RU':  ['P',"C2", "N1", "C5"],
              'RA':  ['P',"C4", "N9", "N7"],
              'RG':  ['P',"C4", "N9", "N7"],

              'RC3':  ['P',"C2", "N1", "C5"],
              'RU3':  ['P',"C2", "N1", "C5"],
              'RA3':  ['P',"C4", "N9", "N7"],
              'RG3':  ['P',"C4", "N9", "N7"],

              'RC5':  ['P',"C2", "N1", "C5"],
              'RU5':  ['P',"C2", "N1", "C5"],
              'RA5':  ['P',"C4", "N9", "N7"],
              'RG5':  ['P',"C4", "N9", "N7"]
              }
num_atom_type=4
##dict_atom = { 'CYT':  ['   P',],
##              'URA':  ['   P',],
##              'ADE':  ['   P',],
##              'GUA':  ['   P',] }
#dict_atom = { 'CYT':  ["  C4", "  C2", "  C6"],
#              'URA':  ["  C4", "  C2", "  C6"],
#              'ADE':  ["  C4", "  C2", "  C6"],
#              'GUA':  ["  C4", "  C2", "  C6"] }



seq_template = {}
data = []
COUNTER =0
for i in open('chain.pdb', 'r'):
    if 'ATOM' in i and i[12:16].strip() in dict_atom[i[17:20].strip()] and int(i[22:26])%1==0:
        xyz = [i[30:38],i[38:46],i[46:54]]
        xyz = [float(j) for  j in xyz]
        temp = [i[6:11],i[12:16].strip(),np.array(xyz),int(i[22:26])]
        seq_template[int(i[22:26])]=i[17:20].strip()[0] # get sequence
        data += [temp,]
        COUNTER += 1

seq_target = {}
data_model= []
COUNTER =0
for i in open('chain.pdb', 'r'):
    if 'ATOM' in i and i[12:16].strip() in dict_atom[i[17:20].strip()] and int(i[22:26])%1==0:
        xyz = [i[30:38],i[38:46],i[46:54]]
        xyz = [float(j) for  j in xyz]
        temp = [i[6:11],i[12:16].strip(),np.array(xyz),int(i[22:26])]
        seq_target[int(i[22:26])]=i[17:20].strip()[0] # get sequence
        data_model += [temp,]
        COUNTER += 1
        

resi_num = [x[3] for x in data]

D_resnum = {}
for i in range(0,len(data)):
    resNum=data[i][3]
    if resNum not in D_resnum.keys():
        D_resnum[resNum] = [data[i],]
    else:
        D_resnum[resNum] += [data[i],]

target = ''
template = ''
for i in seq_target:
    target = target + seq_target[i]
for i in seq_template:
    template = template + seq_template[i]

result = []
f3=open('alignment.input','r')
for i in f3:
    if i[-1]=='\n':
        i=i[0:-1]
    if '>' not in i and i != '\n':
        result += [i,]
a=result[0]
b=result[1]

A,B='',''
counterA,counterB=0,0
Resi_map,template_resi ={},[] #target : template
Resi_map_inv={} # template: target
counter = 0
for i in open("PA.secstruct",'r'):
    if i[-1]=='\n':
        i=i[0:-1]
    if counter ==0:
        seq = i[0:len(i)]
    elif counter == 1:
        ss = i[0:len(i)]
    counter  += 1
# a= model b= #template , we need to add missing info the template
for i in range(len(b)):
    if a[i].lower() in ['a','u','g','c'] and b[i].lower() in ['a','u','g','c'] and counterA < len(ss) and ss[counterA].lower() !='m' :
        Resi_map[counterB+1]=counterA+1
        Resi_map_inv[counterA+1]=counterB+1
        
        template_resi += [counterB+1]
    if a[i].lower() in ['a','u','g','c','m'] :
        counterA+=1
        A = A + a[i].lower() 
    if b[i].lower() in ['a','u','g','c','m'] :
        counterB+=1
        B = B + b[i].lower()
        
##    if  ss[counterA].lower() == 'm':
##        B=B+seq[i]
##        counterB+=1
        
def simillar_atoms(i,j):
    if i[1] == j[1]:
        return True
    if ((j[1] == "N1" or j[1] == "N9") and
        (i[1] == "N1" or i[1] == "N9")) : #takes care of 
        return True
    if ((j[1] == "C2" or j[1] == "C4") and
        (i[1] == "C2" or i[1] == "C4")) :
        return True
    if ((j[1] == "C5" or j[1] == "N7") and
        (i[1] == "C5" or i[1] == "N7")) :
        return True
    

chk =[]
result = []
def got_close_by(data_i,data_j,sysargv=sysargv):
    ''' This finds out if the atoms in template are 3-20nm away. 
    '''
    for i in data_i:        
        for j in data_j :
            if simillar_atoms(i,j) is True :
                vec = i[2]-j[2]
                vec = round(sum(vec**2)**.5,2)
                if (vec > 3 and vec < sysargv) or vec > 190:
                        return True
           
list_of_intr_residues = []
list_of_intr_residues_dd = []
list_of_intr_residues_sssd = []
## Iterate all residue numbers of template, find residues close together ###
## Close is defined as having 1 atom < sysargv### 
for resnum_index in range(len(D_resnum.keys())-1):
    for resnum2_index in range(resnum_index+1, len(D_resnum.keys())):
        resnum=D_resnum.keys()[resnum_index]
        resnum2=D_resnum.keys()[resnum2_index]
        if resnum2 in D_resnum.keys():
            if got_close_by(D_resnum[resnum],D_resnum[resnum2],sysargv) is True:
                list_of_intr_residues += [[resnum,resnum2],]
                list_of_intr_residues_dd += [[resnum,resnum2],]


### Find equailvalent from model ###
pair_wise_constrains = []
pair_wise_constrains_rosetta = []
for i in list_of_intr_residues:
    resnum,resnum2 = i[0],i[1]
    for i in D_resnum[resnum]:        # iterate across template residue pairs
        for j in D_resnum[resnum2] :  # iterate across tempalte reidues pairs
            if simillar_atoms(i,j) is True :
                vec = i[2]-j[2]
                vec = round(sum(vec**2)**.5,2)
                if i[-1] in Resi_map and j[-1] in Resi_map: # resi-Map contains template_resi_num as index and target_resi_num as record
                    sort = sorted([int(i[0]),int(j[0])])
                    pair_wise_constrains += [[sort[0],sort[1],vec],]
                    atomTypei,atomTypej =B[i[-1]-1],B[j[-1]-1] #gets atom_type of i,j
                    #mapping atom type from template to target
                    # get atom index from dict_atom of atom_type
                    atomTypei=[w for w in range(num_atom_type) if dict_atom[atomTypei.upper()][w]==i[1].strip()][0]
                    atomTypej=[w for w in range(num_atom_type) if dict_atom[atomTypej.upper()][w]==j[1].strip()][0]
                    # get atom_number of template
                    atomNumi=Resi_map[i[-1]]
                    atomNumj=Resi_map[j[-1]]
                    atomTypei=dict_atom[A[atomNumi-1].upper()][atomTypei]
                    atomTypej=dict_atom[A[atomNumj-1].upper()][atomTypej]
                    
                    x = '%s %s %s %s HARMONIC %s 1' %(atomTypei,atomNumi,atomTypej,atomNumj,str(vec))
                    if vec < 180:
                        pair_wise_constrains_rosetta += [x]
                    else:
                        if i[1]=='P' and j[1]=='P':
                            pair_wise_constrains_rosetta += [x]
                            print x
                        

f1=open('Rosetta_spring_model_constraints','w')
f1.write('[ atompairs ]\n')
for i in pair_wise_constrains_rosetta[0::]:
    f1.write(i+'\n')
f1.close()
f1=open('constraints','w')
f1.write('[ atompairs ]\n')
for i in pair_wise_constrains_rosetta[0::]:
    f1.write(i+'\n')
f1.close()
###pair_wise_constrains = pair_wise_constrains[0:5] + pair_wise_constrains[0:5] 
##for i in range(0,len(pair_wise_constrains)):
##    if pair_wise_constrains[i][0:2] not in chk:
##        chk += [pair_wise_constrains[i][0:2],]
##        result += [pair_wise_constrains[i],]
##    elif pair_wise_constrains[i][0:2] in chk:
##        None
##        print 'double'
##    
###print len(result)
##for i in sorted(result):
##    None
##    print ('~constraint['+"\\"+'atom1{'+str(i[0])+'}'
##           +"\\"+'atom2{'+str(i[1])+'}'+"\\"+'fk{'+str(KK)
##           +'}'+"\\"+'r0{'+str(i[2])+'}]')



##### base pair constraints ###
##def mid(seq,ss,i): #returns area where there is ')' 
##    for i in range(i+1,len(ss)):
##        if ss[i] == '.':
##            None
##        elif ss[i] == '(':
##            return False,i
##        elif ss[i] == ')':
##            return True,i
###### MAIN(seq,ss) ####
#### CODE goes through all residues , find the inner hairpin loop "{}", exmaple  : ((((...({...})).... , 
#### and record base pair residues, and changes sec str to ((((...(.....)).... , 
##def main(seq,ss): 
##    if len(seq) != len(ss):
##        print 'diff length'
##    else:
##        base_pair = []
##        for i in range(0,len(seq)): 
##            if ss[i] == '(':
##                if mid(seq,ss,i)[0] is True: # found index of ')'
##                    i_1 = mid(seq,ss,i)[1]
##                    ss = ss[0:i]+'.'+ ss[i+1:i_1]+'.'+ss[i_1+1:len(ss)]
##                    base_pair += [[i,i_1,]]
##        return ss,base_pair                    
##                
##main(seq,ss)
###ans  =[[2, 22], [3, 21], [4, 7], [9, 20], [10, 19], [11, 18], [14, 17]]
##
##### recursive_solve(seq,ss,ans=[]) ####
##### Code recursively records and removes ss elements using main(seq,ss), and outputs connected helixes ####
##def recursive_solve(seq,ss,ans=[]):
##    if ('(' not in ss) and (')' not in ss):
##        ans = sorted(ans)
##        helix = []
##        i = 0
##        temp = []
##        for i in range(0,len(ans)-1):
##            ## THIS FINDS IF i and i +1 forms A CONTINOUS HELIX 
##            if (ans[i][0] - ans[i+1][0])**2 == 1 and (ans[i][1] - ans[i+1][1])**2 ==1 and i != len(ans) -2:
##                temp += [ans[i],]
##            ### FINAL END OF LIST "ans" 
##            elif i == len(ans) -2:
##                if (ans[i][0] - ans[i+1][0])**2 == 1 and (ans[i][1] - ans[i+1][1])**2 ==1 :
##                    temp += [ans[i],]
##                    temp += [ans[i+1],]
##                elif (ans[i][0] - ans[i-1][0])**2 == 1 and (ans[i][1] - ans[i-1][1])**2 ==1 :
##                    temp += [ans[i],]
##                    helix += [[ans[i+1]],]
##                else: ##if last one is a solo BP
##                    helix += [[ans[i]],]
##                    helix += [[ans[i+1]],]
##                    
##                if temp != []:
##                    helix += [temp,]
##            ### if not at final end of list, and i and i+1 does not form continous helix
##            else:
##                if (ans[i][0] - ans[i-1][0])**2 == 1 and (ans[i][1] - ans[i-1][1])**2 ==1 :
##                    temp += [ans[i],]
##                ### SOLO ONES ###    
##                elif  ((ans[i][0] - ans[i-1][0])**2 != 1 or (ans[i][1] - ans[i-1][1])**2 !=1 or
##                    (ans[i][0] - ans[i+1][0])**2 != 1 or (ans[i][1] - ans[i+1][1])**2 !=1 ):
##                            helix += [[ans[i]],]
##                if temp != []:
##                    helix += [temp,]
##                temp  = []
##        return helix
##
##
##    else:
##        y = main(seq,ss)
##        return recursive_solve(seq,y[0],ans+y[1])
##
##
##x=recursive_solve(seq,ss)
##
##dict_atom = { 'CYT':  ['N3',"O2","N4"],
##              'URA':  ['O4','N3'],
##              'ADE':  ['N6','N1'],
##              'GUA':  ['N2','O6','N1'] }
##
##bp=((('GUA', 'N2'), ('CYT', 'O2')), (('GUA', 'N1'), ('CYT', 'N3')), (('ADE', 'N6'), ('URA', 'O4')),
##    (('GUA', 'O6'), ('CYT', 'N4')), (('ADE', 'N1'), ('URA', 'N3')))
##
##ade= ((('ADE', 'N6'), ('URA', 'O4')),(('ADE', 'N1'), ('URA', 'N3')))
##cyt= ((('GUA', 'N2'), ('CYT', 'O2')), (('GUA', 'N1'), ('CYT', 'N3')),(('GUA', 'O6'), ('CYT', 'N4')))
###dict_atom = { 'CYT':  ["  C4", "  C2", "  C6"],
###              'URA':  ["  C4", "  C2", "  C6"],
###              'ADE':  ["  C4", "  C2", "  C6"],
###              'GUA':  ["  C4", "  C2", "  C6"] }
##
##data_model =[]
##for i in open('rna.pos_out', 'r'):
##    if 'ATOM' in i and i[12:16].strip() in dict_atom[i[17:20]] and int(i[22:26])%1==0:
##        xyz = [i[30:38],i[38:46],i[46:54]]
##        xyz = [float(j) for  j in xyz]
##        temp = [i[6:11],i[12:16].strip(),np.array(xyz),int(i[22:26])]
##        data_model += [temp,]
##        COUNTER += 1
##        
##ade= (( 'N6', 'O4'),( 'N1','N3'))
##cyt= (('N2', 'O2'), ( 'N1', 'N3'),( 'O6', 'N4'))
##
##for helix in x:
##    for basepair in helix:
##        atom1,atom2=basepair[0]+1,basepair[1]+1
##        atom1_list= [x for x in data_model if x[-1] == atom1]
##        atom2_list= [x for x in data_model if x[-1] == atom2]
##        if len(atom1_list) == 2:
##            for ii in ade:       
##                for atom1 in atom1_list:
##                    for atom2 in atom2_list:
##                        if (atom1[1],atom2[1]) == ii or (atom2[1],atom1[1]) == ii :
##                            i = sorted([int(atom1[0]),int(atom2[0])])
##                            print ('~constraint['+"\\"+'atom1{'+str(i[0])+'}'
##                           +"\\"+'atom2{'+str(i[1])+'}'+"\\"+'fk{'+str(KK*5)
##                           +'}'+"\\"+'r0{3.1}]')
##        else:
##            for ii in cyt:       
##                for atom1 in atom1_list:
##                    for atom2 in atom2_list:
##                        if (atom1[1],atom2[1]) == ii or (atom2[1],atom1[1]) == ii :
##                            i = sorted([int(atom1[0]),int(atom2[0])])
##                            print ('~constraint['+"\\"+'atom1{'+str(i[0])+'}'
##                           +"\\"+'atom2{'+str(i[1])+'}'+"\\"+'fk{'+str(KK*5)
##                           +'}'+"\\"+'r0{3.1}]')
                            
                        
                
        
        


