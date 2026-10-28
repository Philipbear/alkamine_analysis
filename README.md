# alkamine_analysis

This repository contains the code and input data for the data analysis and figure panels of the alkamine paper.

Each folder holds one analysis: a Python script, the input table(s) it reads, and the figure panels it writes as SVG files.

| Folder | Analysis | Input data |
|---|---|---|
| [`fig/large_heatmap/`](fig/large_heatmap) | Alkamine detections in public LC-MS/MS data (MASST): amine × ketone/aldehyde heatmap, human and rodent body-site distributions, steroid alkamines | MASST matches with ReDU sample metadata |
| [`fig/large_heatmap/normalized/`](fig/large_heatmap/normalized) | Body-site distributions normalized by the number of public LC-MS files per body site | Outputs of `fig/large_heatmap/` + file counts per body site |
| [`fig/overall_dist/`](fig/overall_dist) | UpSet plot of alkamines detected across sample categories (food, humans, microbes, plants, rodents) | MASST matches with sample category |
| [`fig/structure_upset/`](fig/structure_upset) | UpSet plot of alkamine structures found in other structure databases | Structure–database overlap table |
| [`culture_data/`](culture_data) | Phenylalanine–3-oxo-LCA in bacterial cultures | LC-MS peak areas |
| [`dose_response/`](dose_response) | Dose–response curves of N-acetylcadaverine–LCA and its isomers | Fluorescence readout |
| [`immune_data/`](immune_data) | Immune-cell assay of alkamines and N-acetylcadaverine–LCA isomers | Fluorescence readout |
| [`T_cell/`](T_cell) | IL-10 and IL-17A in CD4⁺ T cells treated with 5-ASA, sulfapyridine and their conjugates | Percentage of cytokine-positive CD4⁺ cells |

## Requirements

- Python 3.10 (tested on macOS). Each script takes under 15 seconds on a laptop.
- The packages in [`requirements.txt`](requirements.txt): numpy, pandas, matplotlib, scipy, seaborn, upsetplot, openpyxl.
- The Arial font, which the figures use. Without it matplotlib falls back to its default font and prints a warning; the plots are still produced.

## Quick start

```bash
git clone https://github.com/Philipbear/alkamine_analysis.git
cd alkamine_analysis

python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# The two largest input tables are stored gzipped; decompress them once
gunzip -k fig/large_heatmap/data/all_alkamines_new.csv.gz fig/overall_dist/ALL_MASST_category.csv.gz
```

All file paths in the scripts are relative to the repository root, so **run every script from the repository root** (`python <folder>/<script>.py`, not `cd <folder>` first).

To regenerate every figure panel in order:

```bash
# MASST analyses
python fig/large_heatmap/overall_MASST_heatmap.py
python fig/large_heatmap/human_bodysite_scatterplot_amine.py
python fig/large_heatmap/human_bodysite_scatterplot_ketone_aldehyde.py
python fig/large_heatmap/rodent_bodysite_scatterplot_amine.py
python fig/large_heatmap/rodent_bodysite_scatterplot_ketone_aldehyde.py
python fig/large_heatmap/normalized/human_bodysite_scatterplot_amine.py
python fig/large_heatmap/normalized/human_bodysite_scatterplot_ketone_aldehyde.py
python fig/large_heatmap/normalized/rodent_bodysite_scatterplot_amine.py
python fig/large_heatmap/normalized/rodent_bodysite_scatterplot_ketone_aldehyde.py
python fig/large_heatmap/steroid_MASST_heatmap.py
python fig/large_heatmap/bodysite_scatterplot_steroid.py
python fig/overall_dist/overall_dist_upset.py
python fig/structure_upset/structure_overlap_upset.py

# Experimental data
python culture_data/heatmap_culture.py
python dose_response/dose_curve_plot.py
python immune_data/barplot.py
python T_cell/T_cell_barplots.py
```

Order matters within `fig/large_heatmap/`: each `*_amine.py` script writes the per-species table that the matching `*_ketone_aldehyde.py` script reads, and the `normalized/` scripts read the grouped tables written by the scripts one level up. The other analyses are independent of each other.

## Data

### MASST data (`fig/`)

Alkamine detections in public data come from searching the reference MS/MS spectra of the alkamine library against public LC-MS/MS data in GNPS/MassIVE with MASST. Each match is one public MS/MS spectrum, identified by its [Universal Spectrum Identifier](https://www.psidev.info/usi) (USI), joined with the sample metadata from [ReDU](https://redu.gnps2.org). Any USI can be viewed in the [GNPS2 Metabolomics USI resolver](https://metabolomics-usi.gnps2.org).

Alkamine names have the form `<amine>_<ketone or aldehyde>`. In the tables, `Part1` is the amine building block and `Part2` is the ketone or aldehyde building block.

| File | Contents |
|---|---|
| `fig/large_heatmap/data/all_alkamines_new.csv.gz` | MASST matches of the amine–ketone/aldehyde alkamines, one row per (USI, alkamine). Columns: `USI`; `name` (alkamine); `Part1`, `Part2` (building blocks); `Part1class`, `Part2class` (building-block classes, e.g. *Amino acids (Proteinogenic)*, *Polyamines*, *ketone*, *aldehyde*); `NCBITaxonomy`, `UBERONBodyPartName`, `DOIDCommonName` (ReDU metadata: organism, body site, disease); `Taxaname_file`. |
| `fig/large_heatmap/data/steroid_alkamines.csv` | MASST matches of the steroid alkamines (amines conjugated to 3-oxo bile acids and cholestanones). Columns: `usi`; `NAME` (`<amine>_<steroid>`); `cosine`, `matching_peaks`, `delta_mass` (spectral match scores); `MassIVE` (dataset); `Taxaname_file`, `UBERONBodyPartName`, `DOIDCommonName` (ReDU metadata). |
| `fig/large_heatmap/normalized/data/human_bodysite_file_count.tsv`, `rodent_bodysite_file_count.tsv` | Number of unique public LC-MS files per body site (`UBERONBodyPartName`, `file_count`), for human and rodent (*Mus*, *Rattus*) samples. Used to normalize the body-site plots. |
| `fig/overall_dist/ALL_MASST_category.csv.gz` | MASST matches with the sample category of each matched file. Columns: `USI`, `NAME` (alkamine), `Category_group` (`food`, `human`, `microbe`, `plant`, `rodent`, `personalCareProduct`, `other`). |
| `fig/structure_upset/structure_overlap.tsv` | One row per alkamine structure (`inchi`). For each database, `dbid_<database>` lists the matching database IDs and `dbid_<database>_bool` says whether the structure is present. The plot uses the alkamine library, GNPS, PubChem, PubChemLite and HMDB. |

The body-site file counts were generated from the full ReDU metadata table (`all_sampleinformation.tsv`, about 790 MB), which is not stored here. To regenerate them, download the full ReDU metadata table from [ReDU](https://redu.gnps2.org), save it as `fig/large_heatmap/normalized/data/all_sampleinformation.tsv`, and run:

```bash
python fig/large_heatmap/normalized/generate_organ_file_count.py
```

### Experimental data

| File | Contents |
|---|---|
| `culture_data/data/Peak Area extraction_bacteria culture.csv` | LC-MS peak areas from bacterial cultures. Columns: `Molecule`, `Replicate Name`, `Product Mz`, `Product Adduct`, `Retention Time`, `Area`. `Replicate Name` encodes the strain or medium, whether the bile acid was added (`_drug`), the time point (`T0` or `T72`) and the replicate; pooled QC (`poolqc`), `ID` and `Blank` injections are excluded by the script. |
| `dose_response/data/N-acetylcadaverine_LCA_dose response curve.xlsx` | Relative fluorescence units at increasing concentrations of the N-acetylcadaverine–LCA mixture and its two isomers. Each compound takes five columns: concentration (µM), log10 concentration (M), and three replicates. |
| `immune_data/data/Test_immune_cells_nonspecific.tsv` | Relative fluorescence units for each tested metabolite (`Metabolite`, one per row; `DMSO ctrl` is the vehicle control) in three replicates (`Rep1`–`Rep3`). |
| `immune_data/data/Test_immune_cells_nonspecific_Isomer_LCA_Nacetylcadaverine_DMSO.csv` | Same assay for (S)- and (R)-N-acetylcadaverine–LCA against the DMSO control. |
| `T_cell/data/*.csv` (4 files) | Percentage of IL-10⁺ or IL-17A⁺ cells among CD4⁺ T cells, for Th17 and Treg cells (one file per cell type × cytokine). One column per treatment (`Contr`, `5-ASA`, `5-ASA_pyruvate`, `5-ASA_ribose`, `Sulfapyridine`, `Sulfapyridine_glyceraldehyde`, `Sulfapyridine_AKG`), one row per replicate. |

## Scripts and outputs

Each script prints a short summary as it runs. Figures are written inside the analysis folder and intermediate tables go to its `data/` folder; paths in the tables below are relative to the analysis folder unless they start at the repository root.

### `fig/large_heatmap/`

| Script | What it does | Writes |
|---|---|---|
| `overall_MASST_heatmap.py` | Cleans building-block names, removes metformin matches, corrects three building-block class labels, keeps one row per (alkamine, USI), and plots the number of matched spectra for every amine × ketone/aldehyde pair, with amines grouped by class. | `overall_MASST_heatmap.svg`, `.png`; `data/all_refined.tsv`, `data/heatmap_data.tsv`, `data/pivot_data.tsv`, `data/pivot_data_sorted.tsv` |
| `human_bodysite_scatterplot_amine.py` | Keeps matches from *Homo sapiens* samples and counts unique USIs per amine × body site. | `human_bodysite_scatterplot_amine.svg`; `data/human_data.tsv`, `data/human_bodysite_amine_grouped.tsv` |
| `human_bodysite_scatterplot_ketone_aldehyde.py` | Same for ketones/aldehydes. Reads `data/human_data.tsv`. | `human_bodysite_scatterplot_ketone_aldehyde.svg`; `data/human_bodysite_part2_grouped.tsv` |
| `rodent_bodysite_scatterplot_amine.py` | Keeps matches from *Mus* and *Rattus* samples and counts unique USIs per amine × body site. | `rodent_bodysite_scatterplot_amine.svg`; `data/rodent_data.tsv`, `data/rodent_bodysite_amine_grouped.tsv` |
| `rodent_bodysite_scatterplot_ketone_aldehyde.py` | Same for ketones/aldehydes. Reads `data/rodent_data.tsv`. | `rodent_bodysite_scatterplot_ketone_aldehyde.svg`; `data/rodent_bodysite_part2_grouped.tsv` |
| `steroid_MASST_heatmap.py` | Splits steroid alkamine names into amine and steroid, keeps one row per (alkamine, USI), and plots matched spectra per steroid × amine. | `steroid_MASST_heatmap.svg`; `data/steroid_refined.tsv`, `data/steroid_heatmap_data.tsv`, `data/steroid_pivot_data.tsv` |
| `bodysite_scatterplot_steroid.py` | Counts unique USIs per body site for each amine and each steroid of the steroid alkamines. | `bodysite_scatterplot_steroid.svg`; `data/steroid_data.tsv`, `data/bodysite_steroid_grouped.tsv` |

In the body-site scatter plots, circle size is the number of unique matched USIs; the side bar charts show totals per building block and per body site.

### `fig/large_heatmap/normalized/`

The four `*_bodysite_scatterplot_*.py` scripts redraw the human and rodent body-site plots with circle size = matched USIs / number of public LC-MS files from that body site. They read `fig/large_heatmap/data/*_bodysite_*_grouped.tsv` (written by the scripts above) and `data/*_bodysite_file_count.tsv`, and write `*_normalized.svg` in this folder. `generate_organ_file_count.py` regenerates the file counts (see [Data](#masst-data-fig)).

### `fig/overall_dist/` and `fig/structure_upset/`

| Script | What it does | Writes |
|---|---|---|
| `fig/overall_dist/overall_dist_upset.py` | Marks each alkamine as detected or not in food, human, microbe, plant and rodent samples and draws an UpSet plot of the overlaps; prints summary counts. | `fig/overall_dist/overall_dist_upset.svg` |
| `fig/structure_upset/structure_overlap_upset.py` | Draws an UpSet plot of alkamine structures present in the alkamine library, GNPS, PubChem, PubChemLite and HMDB; prints summary counts. | `fig/structure_upset/structure_overlap_upset.svg` |

### Experimental data

| Script | What it does | Writes |
|---|---|---|
| `culture_data/heatmap_culture.py` | Takes the Phenylalanine–3-oxo-LCA (`Phenylalanine_Dehydrolithocholicacid`) peak areas, removes QC and blank injections, and plots the mean peak area per strain or medium at 0 h and 72 h, with and without added bile acid. | `culture_data/heatmap_culture.svg`; `culture_data/data/refined_data.csv`, `culture_data/data/heatmap_matrix_avg_intensity.csv` |
| `dose_response/dose_curve_plot.py` | Fits a four-parameter logistic (Hill) curve on log10 concentration for each compound and plots replicates, mean ± SD and the fit within the measured range. | `dose_response/plots/dose_response_mix.svg`, `_isomer1.svg`, `_isomer2.svg`, `_combined.svg` |
| `immune_data/barplot.py` | Plots replicates with mean ± SEM; p-values from Welch's t-test against the DMSO control (N-acetylcadaverine–LCA in the first plot, each isomer in the second). | `immune_data/plots/immune_cells_nonspecific.svg`, `immune_data/plots/nacetylcadaverine_isomers.svg` |
| `T_cell/T_cell_barplots.py` | For each cell type × cytokine and each conjugate, plots control, parent drug (5-ASA or sulfapyridine) and conjugate as replicates with mean ± SEM, with pairwise Welch's t-test p-values. | 16 SVGs in `T_cell/plots/`, named `<th17 or treg>_<IL10pos or IL17Apos>_<conjugate>.svg` |

The p-values are shown unadjusted, as computed by `scipy.stats.ttest_ind(..., equal_var=False)`.

## License

[Apache License 2.0](LICENSE)
