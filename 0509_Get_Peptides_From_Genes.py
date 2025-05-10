#!/usr/bin/env python3

import pandas as pd
import re
import concurrent.futures
import os
import logging
from tqdm import tqdm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Constants
PEPTIDE_LENGTH_MIN = 3
PEPTIDE_LENGTH_MAX = 7
MAX_WORKERS = 25  # Number of threads for parallel processing
OUTPUT_DIR = "Results"

# Data loading function
def load_data():
    """Loads data from CSV files."""
    try:
        protein_domains = pd.read_csv('Recourses/Domain_Info.csv.gz', index_col=False, compression='gzip')
        protein_protease_sites = pd.read_csv('Recourses/Prosper.csv.gz', index_col=False, compression='gzip')
        protein_info = pd.read_csv('Recourses/Summary_MSP_E.csv', index_col=False)
        return protein_domains, protein_protease_sites, protein_info
    except FileNotFoundError as e:
        logging.error(f"Error loading data: {e}")
        raise

def find_indices(data_list, target_value):
    """Finds indices where elements in data_list are equal to target_value."""
    return [i for i, value in enumerate(data_list) if value == target_value]

def get_protein_protease_cleavage_sites(protein_protease_sites, gene_name):
    """Gets cleavage sites for a given gene."""
    positions = find_indices(protein_protease_sites.iloc[:, 0].tolist(), gene_name)
    if not positions:
        return pd.DataFrame()
    return protein_protease_sites.iloc[positions, :]

def get_protein_sequence(protein_info, gene_name):
    """Gets protein sequence from Summary_MSP_E.csv."""
    positions = find_indices(protein_info.iloc[:, 0].tolist(), gene_name)
    if not positions:
        logging.warning(f"⚠️ Gene Not Found: {gene_name}")
        return ''
    return protein_info.iloc[positions[0], 24]

def get_protein_domains(protein_domains, gene_name):
    """Gets domain information for a given gene."""
    return protein_domains[protein_domains.iloc[:, 0] == gene_name]

def find_substring_index(sequence, substring):
    """Finds substring index in string."""
    try:
        return sequence.index(substring)
    except ValueError:
        return -1

def process_peptides_matrix(input_df, gene_name, protein_seq):
    """Generates peptides from cleavage sites and filter by length (3-7)."""
    peptide_list = []
    seq_len = len(protein_seq)

    for i in range(len(input_df)):
        try:
            row_i = input_df.iloc[i]
            pos_str = row_i['Position'].split(':')[1].strip()
            pos = int(round(float(pos_str.replace('"', ''))))
            score_i = float(row_i['Cleavage score'].split(':')[1].replace('"', ''))
            enzyme_i = re.sub(r"\s+", " ", row_i['Enzyme Name'])

            for j in range(i + 1, len(input_df)):
                row_j = input_df.iloc[j]
                pos_str_j = row_j['Position'].split(':')[1].strip()
                pos_j = int(round(float(pos_str_j.replace('"', ''))))
                score_j = float(row_j['Cleavage score'].split(':')[1].replace('"', ''))
                enzyme_j = re.sub(r"\s+", " ", row_j['Enzyme Name'])

                if pos > pos_j:
                    start, end = pos_j, pos
                    enzyme1, enzyme2 = enzyme_j, enzyme_i
                else:
                    start, end = pos, pos_j
                    enzyme1, enzyme2 = enzyme_i, enzyme_j

                if start < 0 or end > seq_len:
                    continue

                peptide = protein_seq[start:end]
                length = len(peptide)

                if PEPTIDE_LENGTH_MIN <= length <= PEPTIDE_LENGTH_MAX:
                    peptide_list.append({
                        'GN': gene_name,
                        'Enzyme': f"{enzyme1} and {enzyme2}",
                        'start_cleavage_site': start,
                        'end_cleavage_site': end,
                        'average_cleavage_score': (score_i + score_j) / 2,
                        'length_of_peptide': length,
                        'peptide': peptide
                    })

        except Exception as e:
            logging.error(f"Error processing row {i}: {e}")

    return pd.DataFrame(peptide_list)

def find_matching_proteins(peptide, protein_info, protein_domains):
    """Finds other proteins containing this peptide and check domains."""
    matches = []

    for idx, line in enumerate(protein_info[' The_original_sequence']):
        try:
            if isinstance(line, str) and peptide in line:
                gene = protein_info['Protein Gene Name'][idx]
                seq = get_protein_sequence(protein_info, gene)
                if not seq:
                    continue

                start = find_substring_index(seq, peptide)
                end = start + len(peptide)

                domains = get_protein_domains(protein_domains, gene)
                matched_domain = ""

                for _, domain_row in domains.iterrows():
                    d_start, d_end = domain_row.iloc[4], domain_row.iloc[5]
                    if d_start <= start < d_end <= end or start <= d_start < end:
                        matched_domain = f"{gene}({domain_row.iloc[3]})"
                        break

                if not matched_domain:
                    matches.append(f"{gene}( )")
                else:
                    matches.append(matched_domain)

        except Exception as e:
            logging.error(f"Error checking match at index {idx}: {e}")

    return ';'.join(matches)

if __name__ == '__main__':
    # Load data
    try:
        protein_domains, protein_protease_sites, protein_info = load_data()
    except Exception as e:
        logging.error(f"Failed to load data: {e}")
        exit(1)

    protease_ids = ['C01.036', 'M10.003', 'M10.004', 'M10.005', 'M10.008',
                    'S01.131', 'S01.133', 'S01.010']

    # List of genes to analyze
    genes_to_analyze = [
        "HIST1H3A", "HGFAC", "HGF", "HADHA", "H3F3A", "GSN", "GRP", "GPX3", "GPI", "GPC2", "GPC1",
        "GNLY", "GLIPR2", "GLG1", "GLA", "GHR", "GGH", "GDNF", "GBP1", "GBA", "GAPDH", "GANAB",
        "FST", "FRAS1", "FOLR2", "FN1", "FMOD", "FLT1", "FLNB", "FLNA", "FKBP1A", "FGG", "FGFR4",
        "FGFR3", "FGFR2", "FGFR1", "FGF7", "FGF3", "FGF2", "FGB", "FGA", "FBN2", "FBN1", "FBLN5",
        "FBLN2", "FBLN1", "FAM3C", "F2R", "F2", "ESR2", "ERAP1", "EPO", "EPHB4", "ENDOD1", "EMILIN1",
        "EMCN", "ELN", "EGFR", "EFTUD2", "EFNA3", "EFNA1", "EFEMP2", "EFEMP1", "EEF2", "EDN1",
        "ECM1", "DYNC1H1", "DSG1", "DSC2", "DPT", "DPP7", "DMKN", "DMBT1", "DEFB4A", "DEFB1",
        "DEFA3", "DEFA1", "DCN", "DCD", "DAG1", "CXCL9", "CXCL8", "CXCL12", "CXCL1", "CXADR",
        "CTSV", "CTSK", "CTSG", "CTSD", "CTHRC1", "CTGF", "CSTA", "CST6", "CST4", "CST3", "CSPG4",
        "CRTAP", "CRIP2", "CRH", "CREG1", "CPQ", "CPE", "CPA3", "CP", "COPA", "COMP", "COLQ",
        "COL7A1", "COL6A6", "COL6A5", "COL6A3", "COL6A2", "COL6A1", "COL5A2", "COL5A1", "COL4A6",
        "COL4A5", "COL4A4", "COL4A3", "COL4A2", "COL4A1", "COL3A1", "COL2A1", "COL1A2", "COL1A1",
        "COL18A1", "COL17A1", "COL16A1", "COL15A1", "COL14A1", "COL12A1", "CMA1", "CLU", "CLTC",
        "CLCA4", "CKLF", "CKAP4", "CHGB", "CHGA", "CFI", "CFHR3", "CFHR1", "CFH", "CFD", "CFB",
        "CELA1", "CDSN", "CDH1", "CD8B", "CD8A", "CD5L", "CD55", "CD40", "CD34", "CD209", "CD163",
        "CD14", "CD109", "CCT6A", "CCT2", "CCL7", "CCL5", "CCL27", "CCL26", "CCL2", "CCL11",
        "CASP4", "CASP1", "CAPZA2", "CAPZA1", "CANX", "CAMP", "CALU", "CALR", "CALM1", "C9",
        "C8B", "C7", "C5", "C4BPA", "C4B", "C4A", "C3", "C1S", "C1R", "C1QC", "C1QBP", "C1QB",
        "BGN", "BDNF", "B2M", "AZGP1", "ATP5O", "ATP5B", "ATP5A1", "ASPN", "ARF4", "APOH",
        "APOE", "APOD", "APOB", "APOA4", "APOA2", "APOA1", "APCS", "ANXA2P2", "ANXA2", "ANXA1",
        "ANOS1", "ANGPT2", "ANGPT1", "AMBP", "ALDOA", "ALCAM", "ALB", "AIMP1", "AHSG", "AGTR2",
        "AGT", "AFM", "AEBP1", "ADM2", "ADM", "ADIPOQ", "ADAMTSL5", "ADAMTS17", "ACTN4", "ACTN2",
        "ACTN1", "ACPP", "ACE2", "ACE", "ABI3BP", "A2ML1", "A2M", "A1BG"
    ]

    def process_gene(gene):
        output_path = os.path.join(OUTPUT_DIR, f"{gene}.csv")
        if os.path.exists(output_path):
            logging.info(f"Skipping {gene}: Result file already exists.")
            return

        logging.info(f"Analysing: {gene}")

        # Get cleavage sites
        cleavage_sites = get_protein_protease_cleavage_sites(protein_protease_sites, gene)
        # Debug: print column names
        logging.info(f"Columns in cleavage_sites: {cleavage_sites.columns.tolist()}")
        # Find the correct column name for protease IDs
        id_column = [col for col in cleavage_sites.columns if 'ID' in col or 'id' in col or 'Id' in col]
        if not id_column:
            logging.error(f"No ID column found in cleavage_sites for {gene}")
            return
        filtered_sites = cleavage_sites[cleavage_sites[id_column[0]].isin(protease_ids)]

        if filtered_sites.empty:
            logging.info(f"No matching protease sites for {gene}.")
            return

        # Get protein sequence
        protein_seq = get_protein_sequence(protein_info, gene)
        if not protein_seq:
            logging.info(f"No sequence found for {gene}. Skipping.")
            return

        # Generate peptides
        peptides = process_peptides_matrix(filtered_sites, gene, protein_seq)

        if peptides.empty:
            logging.info(f"No peptides generated for {gene}.")
            return

        # Match peptides in other proteins with progress bar
        domain_annotations = []
        for peptide in tqdm(peptides['peptide'],
                          total=len(peptides),
                          desc=f"Matching peptides for {gene}",
                          leave=False):
            result = find_matching_proteins(peptide, protein_info, protein_domains)
            domain_annotations.append({'peptide_proteins_domains': result})

        peptides = pd.concat([peptides.reset_index(drop=True),
                              pd.DataFrame(domain_annotations)], axis=1)

        # Save to CSV
        try:
            peptides.to_csv(output_path, index=False)
            logging.info(f"Saved results to: {output_path}")
        except Exception as e:
            logging.error(f"Error saving results for {gene}: {e}")

    # Use a thread pool to process genes in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # Wrap executor.map with tqdm for progress bar
        list(tqdm(executor.map(process_gene, genes_to_analyze),
             total=len(genes_to_analyze),
             desc="Processing genes"))