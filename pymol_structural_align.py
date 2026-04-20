import Bio.SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO

sequence_list = []
for record in Bio.SeqIO.parse("muscle_seq.fna", "fasta"):
     #sequence order: ABLIM1, MYBPC, MYL2
    sequence_list.append(str(record.seq))

#saving muscle sequences as variables for later 
ablim1_seq = sequence_list[0]
mybpc_seq = sequence_list[1]
myl2_seq = sequence_list[2]


#calling the generated best motifs file: best_overall_matches 
matches_list = ""
with open ("best_overall_matches.tsv") as f: 
    var = f.read().split()
del var[0:8] #deleting titles 
print(var) 

count = len(var)
motif_dict = {}
list = []
protein_list = ["MYL2", "MYBPC", "ABLIM1"]

while count > 0:    
    for i in range(len(var)): 
        count -= 1 
        #checking if i is the name of a mucle protein 
        if var[i] in protein_list: 
            temp_list = var[i], var[i+6] #if yes then will shift 6 down to find motif sequence 
    
            list.append(temp_list) #the sequence is then added to a list to be further filtered 

            if var[i+7] not in protein_list: #if there is a second motif in chart 
                temp_list2 = var[i], var[i+7]
                list.append(temp_list2)

            if var[i+8] not in protein_list: # to account for any variability within chart 
                temp_list3 = var[i], var[i+8]
                list.append(temp_list3)

myl2_motifs = []
mybpc_motifs = []
ablim1_motifs = []

#identified terms (repeats/not motif) from matching motif to protein which need to be removed from list 
out_list = ["13.0", "12.0", "ID:50.0%", "ID:14.3%", "LSDIIQ", "LSSLIQ", "GLALSDLIQKYFF", "GTVLSDIIQKYFF","FKVGHGLAC"]

for i in list: 
    #extra filtering done this way to differenciate between muscle proteins + eliminate repeats 
    if i[0] == "MYL2" and i[1] not in out_list and i[1] not in myl2_motifs:  
        myl2_motifs.append(i[1])

    if i[0] == "MYBPC" and i[1] not in out_list and i[1] not in mybpc_motifs: 
        mybpc_motifs.append(i[1])

    if i[0] == "ABLIM1" and i[1] not in out_list and i[1] not in ablim1_motifs:
        ablim1_motifs.append(i[1])

print(myl2_motifs)
print(mybpc_motifs)
print(ablim1_motifs)
count = 0
count2 = 0 
count3 = 0 

#### Finding Motif Location in Muscle Sequences #### 
#ABLIM1 
for i in range (len(ablim1_seq)): 
        for x in ablim1_motifs: 
            temp_motif = str(x)

            if "".join(ablim1_seq[i:i+len(temp_motif)]) == temp_motif: #checking to see if motif is in sequence 
                start_site = i+1 # start location for motif  
                end_site = i+len(temp_motif) # end location for motif 
                count += 1
                
                #generating a Script for PyMOl ** the script made can be added to pyMOL 
                pymol_script_1 = ["load /Users/siennad0304/Documents/Senior Year/Comp Biology /ABLIM1.pdb", 
                                  #your own directory must be used due to visualize strucutres 
                                  # protein PDB files in GitHUb 
                                          "as cartoon", 
                                          (f"select i. {start_site}-{end_site}"), 
                                          "set_name sele, motif",
                                          "hide everything",
                                          "show sticks, motif",
                                          "fetch 6zu0", #this this the Nav1.5 sodium channel which we will align to  
                                          "align 6zu0, motif",
                                          "orient motif"]
                
                with open(f"ABLIM1_{count}.pml", "w") as l:
                    l.write("\n".join(pymol_script_1))
                
#MYBPC3 
for i in range (len(mybpc_seq)): 
        for x in mybpc_motifs: 
            temp_motif = str(x)

            if "".join(mybpc_seq[i:i+len(temp_motif)]) == temp_motif: 
                start_site = i+1
                end_site = i+len(temp_motif)
                count2 += 1

                #generating a Script for PyMOl ** the script made can be added to pyMOL 
                pymol_script_2 = ["load /Users/siennad0304/Documents/Senior Year/Comp Biology /MYBPC3.pdb",
                                  #your own directory must be used due to visualize strucutres 
                                  # protein PDB files in GitHUb 
                                          "as cartoon", 
                                          (f"select i. {start_site}-{end_site}"), 
                                          "set_name sele, motif",
                                          "hide everything",
                                          "show sticks, motif",
                                          "fetch 6zu0", #this this the Nav1.5 sodium channel which we will align to  
                                          "align 6zu0, motif",
                                          "orient motif"
                                            ]
                with open(f"MYBPC_{count2}.pml", "w") as x:
                    x.write("\n".join(pymol_script_2))

#MYL2 
for i in range (len(myl2_seq)): 
        for x in myl2_motifs: 
            temp_motif = str(x)


            if "".join(myl2_seq[i:i+len(temp_motif)]) == temp_motif: 
                start_site = i+1
                end_site = i+len(temp_motif)
                count3 += 1

                # you don't need to download anything PDB number can be used instead 
                pymol_script_3 = ["fetch 5TBY"
                                    "as cartoon", 
                                    (f"select i. {start_site}-{end_site}"), 
                                    "set_name sele, motif",
                                    "hide everything",
                                    "show sticks, motif",
                                    "fetch 6zu0", #this this the Nav1.5 sodium channel which we will align to  
                                    "align 6zu0, motif",
                                    "orient motif"
                                        ]
                        
                with open(f"MYL2_{count3}.pml", "w") as z:
                    z.write("\n".join(pymol_script_3))