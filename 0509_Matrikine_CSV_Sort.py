import os
import pandas as pd
import re

OUTPUT_FILE = 'Candidate_Matrikines_extra_cleavage_score.csv'
PROCESSED_4MER_DIR = "Processed/4mer"
PROCESSED_5MER_DIR = "Processed/5mer"
EXPORTED_PEPTIDES_DIR = "Results"
ALL_ECM_PROTEINS_DIR = "All_The_ECM_Proteins"
os.makedirs(PROCESSED_4MER_DIR, exist_ok=True)
os.makedirs(PROCESSED_5MER_DIR, exist_ok=True)
os.makedirs(EXPORTED_PEPTIDES_DIR, exist_ok=True)
os.makedirs(ALL_ECM_PROTEINS_DIR, exist_ok=True)

def record_file_excel_format(file, a, df_4mer, df_5mer):
    """Records 4mer and 5mer data in a specified Excel format."""
    file.write(f"{a}\n")
    file.write("4mer\n")
    file.write("Tetramer,Enzyme,Average Cleavage Score,Nr.O.,Proteins Domains\n")
    for _, row in df_4mer.iterrows():
        file.write(
            f"{row['peptide']},{row['Enzyme']},{row['average_cleavage_score']},{row['How many other protein has this peptide:']},{row['peptide_proteins_domains:']}\n")
    file.write("\n")
    file.write("5mer\n")
    file.write("Pentamer,Enzyme,Average Cleavage Score,Nr.O.,Peptide,Proteins Domains\n")
    for _, row in df_5mer.iterrows():
        file.write(
            f"{row['peptide']},{row['Enzyme']},{row['average_cleavage_score']},{row['How many other protein has this peptide:']},{row['peptide_proteins_domains:']}\n")
    file.write("\n")


def record_4mers_5mers(a, df_4mer, df_5mer):
    """Records 4mer and 5mer data to separate CSV files."""
    df_4mer.to_csv(os.path.join(PROCESSED_4MER_DIR, f"4mer{a}.csv"), index=False)
    df_5mer.to_csv(os.path.join(PROCESSED_5MER_DIR, f"5mer{a}.csv"), index=False)


def determine_if_peptide_is_cleaved(entry, peptide, a):
    """Determines if a peptide is cleaved and records cleavage scores."""
    elements = re.split(r'(?<! );', entry)
    list_of_elements = []

    if peptide == 'FLID':
        print(peptide)

    for element in elements:
        filename3 = element.split('(')[0]
        if filename3 == a:
            # Here the peptide represents itself, so we ignore this.
            continue
        try:
            ending = element.split('(')[1]
        except IndexError:
            print(filename3)
            ending = ''  # Handle cases where ending is missing

        try:
            df_filename = pd.read_csv(os.path.join(ALL_ECM_PROTEINS_DIR, f"{filename3}.csv"))
            peptide_series = df_filename['peptide']
            matching_peptides = peptide_series[peptide_series == peptide]

            if not matching_peptides.empty:
                combined_cleavage_scores = []
                for ind in matching_peptides.index:
                    cleavage_score = df_filename['average_cleavage_score'][ind]
                    combined_cleavage_scores.append(str(cleavage_score))
                combined_cleavage_score_string = ';'.join(combined_cleavage_scores)
                record = f"{filename3} ([{combined_cleavage_score_string}] {ending}"
            else:
                record = element
            list_of_elements.append(record)
        except FileNotFoundError:
            record = element
            list_of_elements.append(record)

    total_string = ';'.join(list_of_elements)
    return total_string


def main():
    """Main function to process peptide data."""
    onlyfiles = [f for f in os.listdir(EXPORTED_PEPTIDES_DIR) if
                 os.path.isfile(os.path.join(EXPORTED_PEPTIDES_DIR, f))]

    with open(OUTPUT_FILE, 'w') as outfile:
        for filename in onlyfiles:
            a, _ = os.path.splitext(filename)
            df = pd.read_csv(os.path.join(EXPORTED_PEPTIDES_DIR, filename))

            # Filter values with cleavage score higher than 0.9
            df2 = df[df['average_cleavage_score'] >= 0.9].copy()

            # Reset index to be able to iterate through each of the peptides
            df2 = df2.reset_index(drop=True)

            # DataFrames for storing additional information
            df_number_occurrences = pd.DataFrame(columns=['How many other protein has this peptide:'])
            df_cleavage_score_in_other_proteins = pd.DataFrame(columns=['peptide_proteins_domains:'])

            # Iterate through entries and populate the dataframes
            for peptide, entry in zip(df2['peptide'], df2['peptide_proteins_domains']):
                peptide_cleavage_score_string = determine_if_peptide_is_cleaved(entry, peptide, a)
                number_of_occurrences = peptide_cleavage_score_string.count('(')

                new_row = pd.DataFrame(
                    {'How many other protein has this peptide:': [number_of_occurrences]})  # Create DataFrame
                df_number_occurrences = pd.concat([df_number_occurrences, new_row], ignore_index=True)

                new_row = pd.DataFrame({'peptide_proteins_domains:': [peptide_cleavage_score_string]})
                df_cleavage_score_in_other_proteins = pd.concat([df_cleavage_score_in_other_proteins, new_row],
                                                                 ignore_index=True)

            # Concatenate the new columns to df2
            df2 = pd.concat([df2, df_number_occurrences, df_cleavage_score_in_other_proteins], axis=1)

            # Delete the original 'peptide_proteins_domains' column
            df2 = df2.drop('peptide_proteins_domains', axis=1)

            # Separate 4mers and 5mers
            df_4mer = df2[df2['length_of_peptide'] == 4].copy()
            df_5mer = df2[df2['length_of_peptide'] == 5].copy()

            # Sort values by peptide
            df_4mer = df_4mer.sort_values('peptide')
            df_5mer = df_5mer.sort_values('peptide')

            # Record data to files
            record_file_excel_format(outfile, a, df_4mer, df_5mer)
            record_4mers_5mers(a, df_4mer, df_5mer)


if __name__ == "__main__":
    main()
