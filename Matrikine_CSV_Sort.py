from os import listdir
from os.path import isfile, join
import pandas as pd
import re

F = open('Candidate_Matrikines_extra_cleavage_score.csv','w');


def record_file_Excel_Format():
    F.write(a + "\n")

    F.write("4mer\n")
    F.write("Tetramer,Enzyme,Average Cleavage Score,Nr.O.,Proteins Domains\n")

    for index, row in df_4mer.iterrows():
        F.write(str(row['peptide']) + "," + str(row['Enzyme']) + "," + str(row['average_cleavage_score']) + "," + str(
            row['How many other protein has this peptide:']) + "," + str(row['peptide_proteins_domains:']) + "\n")
    F.write("\n")
    F.write("5mer\n")
    F.write("Pentamer,Enzyme,Average Cleavage Score,Nr.O.,Peptide,Proteins Domains\n")

    for index, row in df_5mer.iterrows():
        F.write(str(row['peptide']) + "," + str(row['Enzyme']) + "," + str(row['average_cleavage_score']) + "," + str(
            row['How many other protein has this peptide:']) + "," + str(row['peptide_proteins_domains:']) + "\n")
    F.write("\n")
def record_4mers_5mers():
    df_4mer.to_csv("Processed/4mer/4mer"+a+".csv", columns=None)
    df_5mer.to_csv("Processed/5mer/5mer"+a+".csv", columns=None)
def determine_if_peptide_is_cleaved(entry, peptide,a):
    list = re.split(r'(?<! );', entry);
    list_of_elements=[]
    if (peptide == 'FLID'):
        print (peptide)

    for element in list:

        if __name__ == '__main__':
            if(element.split('(')[0]==a):
                #here the peptide represents itself, so we ignore this.
                record = '';
                #pass
            else:

                filename3 = element.split('(')[0]
                try:
                    ending = element.split('(')[1]
                except:
                    print (filename3)
                try:
                    df_filename = pd.read_csv('All_The_ECM_Proteins/' + filename3+'.csv')

                    frame = df_filename['peptide']
                    s =frame[frame == peptide]
                    ##Add a debug if here is no peptide

                    if len(s)>0:
                        combined_cleavage_score = []
                        index4 = s.index.tolist()
                        ##here loop through all the indexes.

                        for ind in index4:
                            cleavage_score = df_filename['average_cleavage_score'][ind]
                            combined_cleavage_score.append(str(cleavage_score))
                        combine_cleavage_score_string =  ';'.join(combined_cleavage_score)
                        #record the cleavage score and assemble the string back together
                        record = filename3 + ' ([' + combine_cleavage_score_string + '] ' + ending;
                    else:
                        #here we record the output when there is no matched cleaved peptides.
                        record = element
                    list_of_elements.append(record)
                except IOError:
                    record = element
    total_string = ';'.join(list_of_elements)
   #print total_string
    return total_string

onlyfiles = [f for f in listdir("Exported_Peptides_grB/") if isfile(join("Exported_Peptides_grB/", f))]
for filename in onlyfiles:
    a,b= (filename.split("."))
    df = pd.read_csv('Exported_Peptides_grB/'+filename)

    df2 = df[df['average_cleavage_score'] >= 0.9] #only values with cleavage score higher then 0.9
    #we reset index to be able to iterate through each of the peptides
    df2 = df2.reset_index(drop=True)
    #count the number of unique entries
    df_number_occurances = pd.DataFrame(columns=['How many other protein has this peptide:'])
    df_cleavage_score_in_other_proteins = pd.DataFrame(columns=['peptide_proteins_domains:'])
    count =0;
    for entry in df2['peptide_proteins_domains']:

        #get the peptide sequence
        peptide = df2['peptide'][count]
        count=count+1


        peptide_cleavage_score_string = determine_if_peptide_is_cleaved(entry, peptide, a)
        number_of_occurances = peptide_cleavage_score_string.count('(')
        #here we add a function that searches for the peptide in each of the protein peptide files and outputs the confidence score if peptide is found.

        df_number_occurances=df_number_occurances.append({'How many other protein has this peptide:':number_of_occurances}, ignore_index=True)
        df_cleavage_score_in_other_proteins=df_cleavage_score_in_other_proteins.append({'peptide_proteins_domains:':peptide_cleavage_score_string}, ignore_index=True)
    df2 = df2.reset_index(drop=True)

    df_number_occurances=df_number_occurances.reset_index(drop=True)
    #print df_cleavage_score_in_other_proteins
    df2 = pd.concat([df2, df_number_occurances, df_cleavage_score_in_other_proteins], axis=1)
    del df2['peptide_proteins_domains']

    df_4mer = df2[df2['length_of_peptide'] == 4]
    df_5mer = df2[df2['length_of_peptide'] == 5]


    df_4mer=df_4mer.sort_values('peptide')
    df_5mer=df_5mer.sort_values('peptide')

    #activite foloving to record excel file
    record_file_Excel_Format()

    #activate folowing to record 4mers and 5mers in gene individual files
    record_4mers_5mers()
