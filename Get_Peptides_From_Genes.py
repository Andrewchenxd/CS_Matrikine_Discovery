#!/usr/bin/env python

#Peptide = "LLVLH";
#count_number_of_ocurances_in_MSP(Peptide)
#Main
#protein_gene_name=["EDN1","PLAT","SERPINF2","FBLN1","EFEMP2","DPT","FBN1","LRG1","FGG","BGN","PRELP","FMOD","DCN","LAMB3","FN1","EMILIN1","LTBP4","MFAP4","MFAP2","FBLN5","ELN","FBN1","COL7A1","COL6A3","COL4A3","COL4A1","COL3A1","COL1A1"];

#Writen by Matiss Ozols on Jauary 2018
#Code preforms peptide generation analysis based on Prosper outputs generated on stand by prosper run on Condor Manchester.
#the code also searches which other skin proteins according to MSP contain the particular peptide and in which domain.
#this is further used to derive 4mer and 5mer peptides used for testing in skin cosmetics by Sederma and WBA.


import pandas as pd
import re

#functions
def grep(x, y):
    # print 'grep'
    count = 0
    array = list()

    for line in x:

        if y == line:
            array.append(count)

        count = count + 1
    return array
def get_protein_protease_cleavage_sites(Gene):

    Position_of_Gene = grep(Protein_Protease_sites.iloc[:, 0], Gene)
    try:
        Position_of_Gene1 = Protein_Protease_sites.iloc[Position_of_Gene, :]
    except IndexError:
        Position_of_Gene1 =pd.DataFrame()
    return Position_of_Gene1

def process_peptides_matrix(input):
    
    protein_sequence = get_protein_sequence(The_protein_gene_name_to_search_for)
    protein_sequence = protein_sequence.replace(" ", "")

    peptide_info=pd.DataFrame() #a data frame to store the peptides.

    count=0
    for each in range(0,len(input)):

        try:

            each1 = input.iloc[each]
            pos = each1['Position'].split(':')[1]
            pos=int(pos.replace('"',''))
            cleavages_score=each1['Cleavage score'].split(':')[1]
            cleavages_score = float(cleavages_score.replace('"', ''))
            ID = each1['ID']
            name1 = each1['Enzyme Name']
            name1 = re.sub(r"\s+", " ", name1)
            for each_second_entry in range(each+1,len(input)):
                try:
                    each2 = input.iloc[each_second_entry]
                    pos2=each2['Position'].split(':')[1]
                    pos2 = int(pos2.replace('"', ''))
                    cleavages_score2=each2['Cleavage score'].split(':')[1]
                    cleavages_score2 = float(cleavages_score2.replace('"', ''))
                    ID2 = each2['ID']
                    name2 = each2['Enzyme Name']
                    name2 = re.sub(r"\s+", " ", name2)

                    average_cleavage_score=(cleavages_score2+cleavages_score)/2
                    #here we define which one is the C terminal cleavage position.
                    if pos>pos2:
                        peptide=protein_sequence[pos2:pos]
                        start_cleavage_site=pos2
                        end_cleavage_site=pos
                        Enzyme1=name2
                        Enzyme2=name1
                    else:
                        peptide = protein_sequence[pos:pos2]
                        start_cleavage_site=pos
                        end_cleavage_site = pos2
                        Enzyme1=name1
                        Enzyme2=name2
                    peptide_length = len(peptide)

                    #the cut of thereshold.
                    if peptide_length>2 and peptide_length<8:
                        peptide_info=peptide_info.append({'GN': The_protein_gene_name_to_search_for, 'Enzyme': str(Enzyme1)+' and '+ str(Enzyme2), 'start_cleavage_site': start_cleavage_site, 'end_cleavage_site': end_cleavage_site, 'average_cleavage_score': average_cleavage_score, 'length_of_peptide': peptide_length, 'peptide': peptide}, ignore_index=True)
                    else:
                        pass
                except:
                    pass

        except TypeError:
            print('No Peptides')
            peptide_info=pd.DataFrame()
    count=count+1
    return peptide_info

def get_protein_sequence(Gene):

    Position_of_Gene = grep(Protein_info.iloc[:, 0], Gene)

    # if the gene name is not found give error message
    if len(Position_of_Gene) == 0:
        # print 'no gene found'
        import Tkinter
        import tkMessageBox
        print("Gene Not Found in Your Dataset")
        tkMessageBox.showwarning('Error', 'Gene Name Not Found')
        return 'pass'


    else:
        Position_of_Gene1 = Protein_info.iloc[Position_of_Gene[0], 24]

        return Position_of_Gene1
def Get_protein_domains(The_protein_gene_name_to_search_for):
    Gene_Entries = grep(Protein_domains2.iloc[:, 0], The_protein_gene_name_to_search_for)
    Gene_Entry_lines = Protein_domains2.iloc[Gene_Entries,
                       :]  # the protein that we are willing to display domain information
    return Gene_Entry_lines
def find_str(s, char):
    #from http://stackoverflow.com/questions/21842885/python-find-a-substring-in-a-string-and-returning-the-index-of-the-substring
    index = 0

    if char in s:
        c = char[0]
        for ch in s:
            if ch == c:
                if s[index:index+len(char)] == char:
                    return index

            index += 1

    return -1
def count_number_of_ocurances_in_MSP(peptide):
    count = 0
    a_list_with_ids_of_peptide_cotaining_proteins = []
    count_occurances =0
    count_of_ecm =0;
    count_of_extracellular = 0;
    # loops through all the original sequences
    for line in Protein_info[' The_original_sequence']:

        try:
            if peptide in line:
                count_occurances=count_occurances+1
                Gene = str(Protein_info['Protein Gene Name'][count])
                Location = str(Protein_info['Subcellular Location'][count])
                if (Location=="ECM"):
                    count_of_ecm=count_of_ecm+1
                elif (Location=="extracellular"):
                    count_of_extracellular=count_of_extracellular+1;


                print (Gene + " " +Location)
                # print("peptide matched "+Gene)
                # check if the location is part of any domain
                Sequence = get_protein_sequence(Gene)
                Protein_domains = Get_protein_domains(
                    Gene)  # this finds the domains of the protein that the peptide has been matched with
                Start = find_str(Sequence, peptide)
                End = Start + len(peptide)

                for number2 in range(0, len(Protein_domains)):
                    # susceptibility_of_domain = 0
                    Domain_name = Protein_domains.iloc[number2, 3]

                    Domain_start = Protein_domains.iloc[number2, 4]
                    Domain_finish = Protein_domains.iloc[number2, 5]
                    z = Start
                    q = End
                    x = Domain_start
                    y = Domain_finish
                    gene_peptide_domain = []
                    if x <= z and z < y <= q or z <= x < q and y >= q or z <= x <= q and z <= y <= q or x <= z and y >= q and not x > q and y > q or not y < z and x < z:
                        # if this domain is fully within.
                        gene_peptide_domain = str(Gene + '(' + Domain_name + ')')

                if gene_peptide_domain == []:
                    a_list_with_ids_of_peptide_cotaining_proteins = a_list_with_ids_of_peptide_cotaining_proteins + [
                        str(Gene + '(' + ' ' + ')')]
                else:
                    a_list_with_ids_of_peptide_cotaining_proteins = a_list_with_ids_of_peptide_cotaining_proteins + [
                        str(gene_peptide_domain)]



        except:
            pass
        # print count
        count = count + 1
    peptide_proteins_domains = ';'.join(a_list_with_ids_of_peptide_cotaining_proteins)
    #print peptide_proteins_domains;
    print(count_occurances);

#files we are working with

if __name__ == '__main__':
    Protein_domains2 = pd.read_csv('Recourses/Domain_Info.csv.gz', index_col=False,compression='gzip')
    Protein_Protease_sites = pd.read_csv('Recourses/Prosper.csv.gz', index_col=False,compression='gzip')
    Protein_info = pd.read_csv('Recourses/Summary_MSP_E.csv', index_col=False)

    protein_gene_name=["HIST1H3A","HGFAC","HGF","HADHA","H3F3A","GSN","GRP","GPX3","GPI","GPC2","GPC1","GNLY","GLIPR2","GLG1","GLA","GHR","GGH","GDNF","GBP1","GBA","GAPDH","GANAB","FST","FRAS1","FOLR2","FN1","FMOD","FLT1","FLNB","FLNA","FKBP1A","FGG","FGFR4","FGFR3","FGFR2","FGFR1","FGF7","FGF3","FGF2","FGB","FGA","FBN2","FBN1","FBLN5","FBLN2","FBLN1","FAM3C","F2R","F2","ESR2","ERAP1","EPO","EPHB4","ENDOD1","EMILIN1","EMCN","ELN","EGFR","EFTUD2","EFNA3","EFNA1","EFEMP2","EFEMP1","EEF2","EDN1","ECM1","DYNC1H1","DSG1","DSC2","DPT","DPP7","DMKN","DMBT1","DEFB4A","DEFB1","DEFA3","DEFA1","DCN","DCD","DAG1","CXCL9","CXCL8","CXCL12","CXCL1","CXADR","CTSV","CTSK","CTSG","CTSD","CTHRC1","CTGF","CSTA","CST6","CST4","CST3","CSPG4","CRTAP","CRIP2","CRH","CREG1","CPQ","CPE","CPA3","CP","COPA","COMP","COLQ","COL7A1","COL6A6","COL6A5","COL6A3","COL6A2","COL6A1","COL5A2","COL5A1","COL4A6","COL4A5","COL4A4","COL4A3","COL4A2","COL4A1","COL3A1","COL2A1","COL1A2","COL1A1","COL18A1","COL17A1","COL16A1","COL15A1","COL14A1","COL12A1","CMA1","CLU","CLTC","CLCA4","CKLF","CKAP4","CHGB","CHGA","CFI","CFHR3","CFHR1","CFH","CFD","CFB","CELA1","CDSN","CDH1","CD8B","CD8A","CD5L","CD55","CD40","CD34","CD209","CD163","CD14","CD109","CCT6A","CCT2","CCL7","CCL5","CCL27","CCL26","CCL2","CCL11","CASP4","CASP1","CAPZA2","CAPZA1","CANX","CAMP","CALU","CALR","CALM1","C9","C8B","C7","C5","C4BPA","C4B","C4A","C3","C1S","C1R","C1QC","C1QBP","C1QB","BGN","BDNF","B2M","AZGP1","ATP5O","ATP5B","ATP5A1","ASPN","ARF4","APOH","APOE","APOD","APOB","APOA4","APOA2","APOA1","APCS","ANXA2P2","ANXA2","ANXA1","ANOS1","ANGPT2","ANGPT1","AMBP","ALDOA","ALCAM","ALB","AIMP1","AHSG","AGTR2","AGT","AFM","AEBP1","ADM2","ADM","ADIPOQ","ADAMTSL5","ADAMTS17","ACTN4","ACTN2","ACTN1","ACPP","ACE2","ACE","ABI3BP","A2ML1","A2M","A1BG"]
    all_the_protease_ids=['C01.036','M10.003','M10.004', 'M10.005','M10.008','S01.131','S01.133','S01.010'];

    for The_protein_gene_name_to_search_for in protein_gene_name :
        print(f"Analysing: {The_protein_gene_name_to_search_for}")

        cleavage_sites_for_gene = Protein_Protease_sites[Protein_Protease_sites.iloc[:, 0]==The_protein_gene_name_to_search_for]
        
        all_indexes=[];
        for g in all_the_protease_ids:
            indexes=grep(cleavage_sites_for_gene['ID'],g)
            all_indexes=all_indexes+indexes


        Protease_sites = cleavage_sites_for_gene.iloc[all_indexes, :]; #cleavage_sites_for_gene = all the cleavage sites for a specific gene name. all_indexes = all the columns that has the protease ID
        Peptides=process_peptides_matrix(Protease_sites)
        new_column_same_sequences_bio_f = pd.DataFrame()

        #this part finds all the other sequences that contain this pattern

        if (len(Peptides) < 1):
            #print "true, it is smaler then 1"
            pass
        else:


            for peptide in Peptides['peptide']:
                count = 0
                a_list_with_ids_of_peptide_cotaining_proteins = []

                #loops through all the original sequences
                for line in Protein_info[' The_original_sequence']:
                    try:
                        if peptide in line:

                            Gene = str(Protein_info['Protein Gene Name'][count])
                            #print("peptide matched "+Gene)
                            # check if the location is part of any domain
                            Sequence = get_protein_sequence(Gene)
                            Protein_domains = Get_protein_domains(Gene) #this finds the domains of the protein that the peptide has been matched with
                            Start = find_str(Sequence, peptide)
                            End = Start + len(peptide)

                            for number2 in range(0, len(Protein_domains)):
                                # susceptibility_of_domain = 0
                                Domain_name = Protein_domains.iloc[number2, 3]

                                Domain_start = Protein_domains.iloc[number2, 4]
                                Domain_finish = Protein_domains.iloc[number2, 5]
                                z = Start
                                q = End
                                x = Domain_start
                                y = Domain_finish
                                gene_peptide_domain = []
                                if x <= z and z < y <= q or z <= x < q and y >= q or z <= x <= q and z <= y <= q or x <= z and y >= q and not x > q and y > q or not y < z and x < z:
                                    # if this domain is fully within.
                                    gene_peptide_domain = str(Gene + '(' + Domain_name + ')')

                            if gene_peptide_domain == []:
                                a_list_with_ids_of_peptide_cotaining_proteins = a_list_with_ids_of_peptide_cotaining_proteins + [
                                    str(Gene + '(' + ' ' + ')')]
                            else:
                                a_list_with_ids_of_peptide_cotaining_proteins = a_list_with_ids_of_peptide_cotaining_proteins + [
                                    str(gene_peptide_domain)]

                    except:
                        pass
                    # print count
                    count = count + 1


                #here we prepear for the export
                peptide_proteins_domains = ';'.join(a_list_with_ids_of_peptide_cotaining_proteins)
                new_column_same_sequences_bio_f = new_column_same_sequences_bio_f.append(
                    {'peptide_proteins_domains': peptide_proteins_domains}, ignore_index=True)
            Peptides = pd.concat([Peptides.reset_index(), new_column_same_sequences_bio_f], axis=1)
            Peptides.to_csv("Results/"+The_protein_gene_name_to_search_for+".csv", index=False)
            #Export_function(Peptides)


