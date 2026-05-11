import pandas as pd

RODENT_TAXA = [
    '10088|Mus',
    '10090|Mus musculus',
    '10105|Mus minutoides',
    '10114|Rattus',
    '10116|Rattus norvegicus'
]


def generate_file_counts():
    df = pd.read_csv('fig/large_heatmap/normalized/data/all_sampleinformation.tsv', sep='\t', low_memory=False)

    # Exclude rows with missing body site
    df = df[df['UBERONBodyPartName'].notna() & (df['UBERONBodyPartName'] != 'missing value')]

    # Human: count unique LCMS files per body site
    human_df = df[df['NCBITaxonomy'].str.contains('Homo sapiens', na=False)]
    human_counts = human_df.groupby('UBERONBodyPartName')['filename'].nunique().reset_index()
    human_counts.columns = ['UBERONBodyPartName', 'file_count']
    human_counts = human_counts.sort_values('file_count', ascending=False).reset_index(drop=True)
    human_counts.to_csv('fig/large_heatmap/normalized/data/human_bodysite_file_count.tsv', sep='\t', index=False)
    print(f"Human body site file counts ({len(human_counts)} sites):")
    print(human_counts.to_string(index=False))

    # Rodent: count unique LCMS files per body site
    rodent_pattern = '|'.join([taxa.split('|')[1].strip() for taxa in RODENT_TAXA])
    rodent_df = df[df['NCBITaxonomy'].str.contains(rodent_pattern, na=False, regex=True)]
    rodent_counts = rodent_df.groupby('UBERONBodyPartName')['filename'].nunique().reset_index()
    rodent_counts.columns = ['UBERONBodyPartName', 'file_count']
    rodent_counts = rodent_counts.sort_values('file_count', ascending=False).reset_index(drop=True)
    rodent_counts.to_csv('fig/large_heatmap/normalized/data/rodent_bodysite_file_count.tsv', sep='\t', index=False)
    print(f"\nRodent body site file counts ({len(rodent_counts)} sites):")
    print(rodent_counts.to_string(index=False))

    return human_counts, rodent_counts


if __name__ == "__main__":
    generate_file_counts()
    print('\nFile count generation complete!')
