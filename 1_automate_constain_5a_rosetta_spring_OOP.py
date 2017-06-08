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


class PDB(object):
    def __init__(self):#,Resi_map=Resi_map,seqTemp=B,seqTar=A):
        None
##        self.Resi_map = Resi_map
##        self.seqTemp = seqTemp
##        self.seqTar = seqTar
    def pre_processing(self,pdb='chain.pdb'):
        seq_template = {}
        data = []
        COUNTER =0
        for i in open(pdb, 'r'):
            if 'ATOM' in i and i[12:16].strip() in dict_atom[i[17:20].strip()] and int(i[22:26])%1==0:
                xyz = [i[30:38],i[38:46],i[46:54]]
                xyz = [float(j) for  j in xyz]
                temp = [i[6:11],i[12:16].strip(),np.array(xyz),int(i[22:26])]
                seq_template[int(i[22:26])]=i[17:20].strip()[0] # get sequence
                data += [temp,]
                COUNTER += 1
        self.data=data
        resi_num = [x[3] for x in data]
        D_resnum = {}
        for i in range(0,len(data)):
            resNum=data[i][3]
            if resNum not in D_resnum.keys():
                D_resnum[resNum] = [data[i],]
            else:
                D_resnum[resNum] += [data[i],]

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
        self.Resi_map = Resi_map
        self.seqTemp = B
        self.seqTar = A
        self.D_resnum=D_resnum

    def got_close_by(self,data_i,data_j,sysargv=sysargv):
        ''' This finds out if the atoms in template are 3-20nm away. 
        '''
        for i in data_i:        
            for j in data_j :
                if self.simillar_atoms(i,j) is True :
                    vec = i[2]-j[2]
                    vec = round(sum(vec**2)**.5,2)
                    if vec > 3 and vec < sysargv:
                            return True
    def simillar_atoms(self,i,j):
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
    def get_interacting_residues(self):
        list_of_intr_residues = []
        list_of_intr_residues_dd = []
        list_of_intr_residues_sssd = []
        D_resnum=self.D_resnum
        ## Iterate all residue numbers of template, find residues close together ###
        ## Close is defined as having 1 atom < sysargv### 
        for resnum_index in range(len(D_resnum.keys())-1):
            for resnum2_index in range(resnum_index+1, len(D_resnum.keys())):
                resnum=D_resnum.keys()[resnum_index]
                resnum2=D_resnum.keys()[resnum2_index]
                if resnum2 in D_resnum.keys():
                    if self.got_close_by(D_resnum[resnum],D_resnum[resnum2],sysargv) is True:
                        list_of_intr_residues += [[resnum,resnum2],]
                        list_of_intr_residues_dd += [[resnum,resnum2],]
        self.list_of_intr_residues=list_of_intr_residues
    def get_pair_wise_constraints(self):
        D_resnum=self.D_resnum
        Resi_map = self.Resi_map
        ### Find equailvalent from model ###
        pair_wise_constrains = []
        pair_wise_constrains_rosetta = []
        for i in self.list_of_intr_residues:
            resnum,resnum2 = i[0],i[1]
            for i in D_resnum[resnum]:        # iterate across template residue pairs
                for j in D_resnum[resnum2] :  # iterate across tempalte reidues pairs
                    if self.simillar_atoms(i,j) is True :
                        vec = i[2]-j[2]
                        vec = round(sum(vec**2)**.5,2)
                        if i[-1] in Resi_map and j[-1] in Resi_map: # resi-Map contains template_resi_num as index and target_resi_num as record
                            sort = sorted([int(i[0]),int(j[0])])
                            pair_wise_constrains += [[sort[0],sort[1],vec],]
                            atomTypei,atomTypej = self.seqTemp[i[-1]-1], self.seqTemp[j[-1]-1] #gets atom_type of i,j
                            #mapping atom type from template to target
                            # get atom index from dict_atom of atom_type
                            atomTypei=[w for w in range(num_atom_type) if dict_atom[atomTypei.upper()][w]==i[1].strip()][0]
                            atomTypej=[w for w in range(num_atom_type) if dict_atom[atomTypej.upper()][w]==j[1].strip()][0]
                            # get atom_number of template
                            atomNumi=Resi_map[i[-1]]
                            atomNumj=Resi_map[j[-1]]
                            atomTypei=dict_atom[self.seqTar[atomNumi-1].upper()][atomTypei]
                            atomTypej=dict_atom[self.seqTar[atomNumj-1].upper()][atomTypej]
                            
                            x = '%s %s %s %s HARMONIC %s 1' %(atomTypei,atomNumi,atomTypej,atomNumj,str(vec))
                            pair_wise_constrains_rosetta += [x]
        self.pair_wise_constrains_rosetta = pair_wise_constrains_rosetta
    def write_constraint(self):
        pair_wise_constrains_rosetta =self.pair_wise_constrains_rosetta
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
a=PDB()
a.pre_processing()
a.get_interacting_residues()
a.get_pair_wise_constraints()
a.write_constraint()

