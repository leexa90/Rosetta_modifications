## rna_denovo_setup.py -fasta PA.fasta -secstruct_file PA.secstruct -working_res 6-125 -s *.pdb -fixed_stems -tag TPP_LONG_FL_sln 

'''
This code uses PA.secstruct and PA.fasta to write code to make RNA helices.
RNA helices include canonical and non-concical BP and NO pseudoknots
'''



pdb_string = ''
# list of differemnt base pairs 
pseudo_list=[ ['(',')'],['[',']'],['{','}'],['<','>'],['A', 'a'],
              ['B', 'b'], ['C', 'c'], ['D', 'd'], ['E', 'e'], ['F', 'f'],
              ['G', 'g'], ['H', 'h'], ['I', 'i'], ['J', 'j']][0:1]
def read_dbn_file(file_n ="PA.secstruct"): # read dotbracket notation file
    counter = 0
    for i in open(file_n,'r'):
        if i[-1]=='\n':
            i=i[0:-1]
        if counter ==0:
            seq = i.lower()
        elif counter == 1:
            ss = i
        counter  += 1
    return seq,ss
def mid(seq,ss,i): #returns areas where there is a hairpin 
    for i in range(i+1,len(ss)):
        if ss[i] == '.':
            None
        elif ss[i] == bp1:
            return False,i
        elif ss[i] == bp2:
            return True,i
        


def main(seq,ss): # removes a pair of base pairs from ss file, and updates base_pairs
    if len(seq) != len(ss):
        print 'diff length'
    else:
        base_pair = []
        for i in range(0,len(seq)):
            if ss[i] == bp1:
                if mid(seq,ss,i)[0] is True:
                    i_1 = mid(seq,ss,i)[1]
                    ss = ss[0:i]+'.'+ ss[i+1:i_1]+'.'+ss[i_1+1:len(ss)]
                    base_pair += [[i,i_1,]]
        #print ss,base_pair
        return ss,base_pair                    
                
def recursive_solve(seq,ss,ans=[]):
    ''' recursively finds base pairs, and outputs res numbers of RNA hairpins
    input ((..))..(.)
    output : [ [[1,6],[2,5]] , [[9,1]] ]
    '''
    if (bp1 not in ss) and (bp2 not in ss):
        #print ans
        ans = sorted(ans)
        helix = []
        j = 0
        while j in range(0,len(ans)-1):
            temp = [ans[j]]
            for i in range(j,len(ans)-1): 
                if (ans[i][0] - ans[i+1][0])**2 == 1 and (ans[i][1] - ans[i+1][1])**2 ==1 and i != len(ans) -2:
                    temp += [ans[i+1],]
                    
                elif i == len(ans) -2: # compare last pair of the list
                    #print ans[i],ans[i+1]
                    if (ans[i][0] - ans[i+1][0])**2 == 1 and (ans[i][1] - ans[i+1][1])**2 ==1 :
                        temp += [ans[i+1],]
                        helix += [temp,]
                    else: # if last element is a single basepair
                        helix += [temp,]
                        temp = [ans[i+1]]
                        helix += [temp,]
                    
                    j = i + 1
                    break                   
                else:
                    helix += [temp,]
                    j = i + 1
                    break
        if len(ans) ==1:
            helix=[ans]
        return helix
    else:
        y = main(seq,ss)
        return recursive_solve(seq,y[0],ans+y[1])
    
for k in range(len(pseudo_list)):
    seq,ss=read_dbn_file() #read dot bracket nottation file
    bp1,bp2=pseudo_list[k][0],pseudo_list[k][1] # get type of base pairs (eg.[],())
    ### print ###
    counter =0
    f1=open('commands','w')
    for i in recursive_solve(seq,ss): 
        x = [x[0] for x in i]
        y = [y[1] for y in i]
        x = sorted(x)
        y = sorted(y)
        seqx,seqy = '' ,''
        for j in range(0,len(x)):
            seqx += seq[x[j]]
            seqy += seq[y[j]]
        # printing output
        if len(x) >= 2:
            f1.write('rna_helix.py  -o '
                   + str(counter)+'_'+str(k)+'.pdb -seq '
                   +seqx+' '+seqy+' -resnum ' + str(x[0]+1)+'-'+str(x[len(x)-1]+1) + ' ' + str(y[0]+1)+'-'+str(y[len(y)-1]+1) +'\n' )
            pdb_string = pdb_string + ' '+ str(counter)+'_'+str(k)+'.pdb '
            counter += 1
        else:
            if x[0] -y[0] != -1:
                f1.write ('rna_helix.py  -o '+ str(counter)+'_'+str(k)+'.pdb -seq '
                       +seqx+' '+seqy+' -resnum ' + str(x[0]+1) + ' ' + str(y[0]+1) + '\n')
                pdb_string = pdb_string + ' '+ str(counter)+'_'+str(k)+'.pdb '
                counter += 1
    print '\n'

f1.write('rna_denovo_setup.py -fasta PA.fasta -secstruct_file PA.secstruct -working_res'+
       ' 1'+'-'+str(len(ss))+' -s '+pdb_string+' -fixed_stems -tag TPP_FL_sln  -nstruct 80 -j 16 -cycles 20000 \
-cst_file constraints -staged_constraints -minimize_rna\n')
f1.close()
#print ('#source README_FARFAR')
# rna_helix.py -o H2.pdb -seq cugug aggag -resnum 6-10 122-126 
##def recursive_find_helix(seq,ss,ans,helix=[],i=0):
##    if  len(ans) ==i-1:
##        return helix
##    elif (ans[i][0] - ans[1+i][0])**2 != 1 or (ans[i][1] - ans[1+i][1])**2 !=1:
##        print helix + ans[0]
##        return recursive_find_helix(seq,ss,ans,helix,i+1)
##    else:
##        print helix + ans[0]
##        return recursive_find_helix(seq,ss,ans,helix,i+1)
        

