# alkamine_analysis

Code and data for the analyses and figure panels of the alkamine paper.

## Contents

| Folder | Analysis |
|---|---|
| `fig/large_heatmap/` | Alkamines detected in public LC-MS/MS data (MASST): amine × ketone/aldehyde heatmap, human and rodent body-site plots, steroid alkamines |
| `fig/large_heatmap/normalized/` | Body-site plots normalized by the number of public files per body site |
| `fig/overall_dist/` | UpSet plot of alkamines across sample categories (food, humans, microbes, plants, rodents) |
| `fig/structure_upset/` | UpSet plot of alkamine structures found in other databases |
| `culture_data/` | Phenylalanine–3-oxo-LCA in bacterial cultures |
| `dose_response/` | Dose–response curves of N-acetylcadaverine–LCA and its isomers |
| `immune_data/` | Immune-cell assay of alkamines |
| `T_cell/` | IL-10 and IL-17A in CD4⁺ T cells treated with 5-ASA, sulfapyridine and their conjugates |

Each script saves its figures as SVG files in its own folder.

## Data

The MASST tables list matches between the alkamine reference MS/MS spectra and public LC-MS/MS data in GNPS/MassIVE. Each match is a public spectrum identified by its [USI](https://www.psidev.info/usi), with sample metadata (organism, body site, disease) from [ReDU](https://redu.gnps2.org). Alkamine names have the form `<amine>_<ketone or aldehyde>`.

| File | Contents |
|---|---|
| `fig/large_heatmap/data/all_alkamines_new.csv.gz` | MASST matches of amine–ketone/aldehyde alkamines |
| `fig/large_heatmap/data/steroid_alkamines.csv` | MASST matches of steroid alkamines |
| `fig/large_heatmap/normalized/data/*_bodysite_file_count.tsv` | Number of public LC-MS files per body site (human, rodent) |
| `fig/overall_dist/ALL_MASST_category.csv.gz` | MASST matches with the sample category of each file |
| `fig/structure_upset/structure_overlap.tsv` | Alkamine structures and the databases that contain them |
| `culture_data/data/` | LC-MS peak areas from bacterial cultures |
| `dose_response/data/` | Fluorescence at each concentration, three replicates |
| `immune_data/data/` | Fluorescence for each metabolite, three replicates; DMSO is the control |
| `T_cell/data/` | % IL-10⁺ or IL-17A⁺ among CD4⁺ cells; one column per treatment, one row per replicate |

The files-per-body-site tables were made with `fig/large_heatmap/normalized/generate_organ_file_count.py` from the full ReDU metadata table, which is too large to include (~790 MB) and can be downloaded from ReDU.

## License

[Apache License 2.0](LICENSE)
