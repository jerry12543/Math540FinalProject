import pandas as pd
import config


def get_metadata():
    sample_metadata = {}
    for file_path in config.SERIESMATRIXPATHS:
        gsm_ids = []
        conditions = []
        celltype = []

        with open(file_path, "r") as file:
            for line in file:
                if line.startswith("!Sample_geo_accession"):
                    gsm_ids = line.strip().split("\t")[1:]  # extract all GSM IDs
                    gsm_ids = [entry.replace('"', "") for entry in gsm_ids]
                if line.startswith("!Sample_characteristics_ch1") and "cell type" in line:
                    celltype = [entry.split(": ")[1] for entry in line.strip().split("\t")[1:]] # extract celltype
                    celltype = [entry.replace('"', "") for entry in celltype]
                if line.startswith("!Sample_characteristics_ch1") and "diagnosis" in line:
                    conditions = [entry.split(": ")[1] for entry in line.strip().split("\t")[1:]]  # extract diagnosis
                    conditions = [entry.replace('"', "") for entry in conditions]
            if gsm_ids and conditions and celltype:
                for gsm, condition, celltype in zip(gsm_ids, conditions, celltype):
                    sample_metadata[gsm] = {"diagnosis": condition, "cell_type": celltype}
    return sample_metadata


# Load the raw counts matrix
raw_counts = pd.read_csv(config.RAWCOUNTSPATH, sep='\t', index_col=0)

# Load the normalized TPM matrix
tpm_counts = pd.read_csv(config.TPMCOUNTSPATH, sep='\t', index_col=0)

# Load the gene annotation table
annotations = pd.read_csv(config.ANNOTATIONSPATH, sep='\t', index_col=0)

metadata = get_metadata()

# Display data summaries
# print(raw_counts.head())
# print(tpm_counts.head())
# print(annotations.head())
