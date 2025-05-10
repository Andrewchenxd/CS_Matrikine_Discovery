# Matrikine Discovery Pipeline

<p align="center">
  <img src="./images/Ozols_FINAL.jpg" width="60%"/>
</p>


### Repository Overview

This repository hosts a collection of code designed for the comprehensive analysis of matrikines—specifically, peptides with potential applications in skin health and cosmetics. Utilizing a rich dataset that includes protein domain information, cleavage site data, and UniProt entries, our scripts perform detailed peptide generation and analysis. This process leverages output from the Prosper algorithm. The analysis not only identifies peptides derived from skin proteins but also determines their presence within specific protein domains as cataloged in the MSP database. The ultimate goal is to isolate 4mer and 5mer peptides for empirical testing as discussed in [manuscript](https://academic.oup.com/bjd/advance-article/doi/10.1093/bjd/ljae061/7610994?login=true).

#### Key Components of the Analysis

- **Peptide Generation**: Automated analysis based on Prosper algorithm outputs, focusing on peptides that may play a role in skin health.
- **Protein Domain Search**: Identification of skin proteins and the domains that contain specific peptides of interest.
- **Peptide Application**: Isolation of short peptides (4mer and 5mer) for subsequent testing in skincare products.

### Manuscript Reference

The methodologies and predictions utilized in this project are detailed in a manuscript that leverages Prosper prediction data. This foundational work is accessible [here](https://academic.oup.com/bjd/advance-article/doi/10.1093/bjd/ljae061/7610994?login=true).

### Data Preparation

Prepared files and datasets are integral to our analysis:

- **Skin Proteome Analysis**: The skin proteome, encompassing 2,859 proteins listed in the [Manchester Proteome Database](https://www.manchesterproteome.manchester.ac.uk/#/Proteome), was analyzed using the Prosper algorithm. The resultant data are committed to this repository (`./Resources/Prosper.csv.gz`).
- **Protein Domain Information**: Data on protein domains, including their positions, were meticulously extracted to assess the occurrence of specific peptides within active domains (`./Resources/Domains_Info.csv.gz`).

This repository serves as a comprehensive resource for researchers and industry professionals interested in the intersection of bioinformatics and cosmetic science, providing tools and data for the advanced analysis of peptides in skin health applications.
