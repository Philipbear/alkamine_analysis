import pandas as pd
import numpy as np


def create_scatter_plot():
    """
    Normalized version: circle size = USI matches / number of LCMS files for that body site.
    """
    grouped = pd.read_csv('fig/large_heatmap/data/human_bodysite_amine_grouped.tsv', sep='\t')
    file_counts = pd.read_csv('fig/large_heatmap/normalized/data/human_bodysite_file_count.tsv', sep='\t')

    # Normalize counts by file count per body site
    grouped = grouped.merge(file_counts, on='UBERONBodyPartName', how='left')
    missing = grouped['file_count'].isna().sum()
    if missing > 0:
        print(f"Warning: {missing} rows have no file count — dropping them")
        grouped = grouped.dropna(subset=['file_count'])
    grouped['count_norm'] = grouped['count'] / grouped['file_count']

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 6

    part1_class_order = [
        'Amino acids (Proteinogenic)',
        'Amino acids (Proteinogenic; Modified)',
        'Amino acids (Non-proteinogenic)',
        'Polyamines',
        'Peptides and peptide derivatives',
        'Neurotransmitters',
        'Biogenic amines',
        'Drugs and drug metabolites',
        'Nucleic',
        'Other'
    ]

    class_display_names = {
        'Amino acids (Proteinogenic; Modified)': 'Amino acids (Proteinogenic; modified)',
        'Peptides and peptide derivatives': 'Peptides',
        'Drugs and drug metabolites': 'Drugs',
        'Nucleic': 'Nucleic acid derivatives',
        'Other': 'Others'
    }

    unique_body_sites = sorted(grouped['UBERONBodyPartName'].unique(), reverse=True)

    class_order_map = {cls: i for i, cls in enumerate(part1_class_order)}
    grouped['class_order'] = grouped['Part1class'].map(lambda x: class_order_map.get(x, 999))

    amine_info = grouped[['Part1', 'Part1class', 'class_order']].drop_duplicates()
    amine_info = amine_info.sort_values(['class_order', 'Part1'])
    unique_amines = amine_info['Part1'].tolist()

    amine_to_pos = {amine: i for i, amine in enumerate(unique_amines)}
    bodysite_to_pos = {site: i for i, site in enumerate(unique_body_sites)}

    x = grouped['Part1'].map(amine_to_pos)
    y = grouped['UBERONBodyPartName'].map(bodysite_to_pos)
    sizes = grouped['count_norm']

    classes_in_data = [cls for cls in part1_class_order if cls in grouped['Part1class'].values]

    palette1 = plt.cm.Set2.colors
    palette2 = plt.cm.Paired.colors
    combined_palette = list(palette1) + list(palette2)
    class_to_color = {cls: combined_palette[i] for i, cls in enumerate(classes_in_data)}
    colors = grouped['Part1class'].map(class_to_color)

    # Normalized totals for bar plots
    amine_totals = grouped.groupby('Part1')['count_norm'].sum()
    amine_totals = amine_totals.reindex(unique_amines, fill_value=0)

    bodysite_totals = grouped.groupby('UBERONBodyPartName')['count_norm'].sum()
    bodysite_totals = bodysite_totals.reindex(unique_body_sites, fill_value=0)

    fig = plt.figure(figsize=(8, 3.))
    gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[12, 1], height_ratios=[1, 5],
                           hspace=0.025, wspace=0.025)

    ax_main = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)

    size_scale = 0.0025
    scatter = ax_main.scatter(x, y, s=np.sqrt(sizes) / size_scale, alpha=0.6, c=colors, linewidth=0)

    amine_colors = [class_to_color[amine_info[amine_info['Part1'] == amine]['Part1class'].values[0]]
                    for amine in unique_amines]
    ax_top.bar(range(len(unique_amines)), amine_totals, color=amine_colors, alpha=0.6, width=0.8)
    ax_top.set_ylabel('Normalized\ncount', fontsize=6)
    ax_top.tick_params(axis='y', labelsize=5, length=1.5, width=0.5, color='0.2', pad=1.5)
    ax_top.tick_params(axis='x', labelbottom=False, length=0)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['bottom'].set_visible(False)
    for spine in ['left']:
        ax_top.spines[spine].set_linewidth(0.5)
        ax_top.spines[spine].set_color('0.4')
    ax_top.set_xlim(-0.5, len(unique_amines) - 0.5)

    ax_right.barh(range(len(unique_body_sites)), bodysite_totals, color='0.6', alpha=0.6, height=0.8)
    ax_right.set_xlabel('Normalized\ncount', fontsize=6)
    ax_right.tick_params(axis='x', labelsize=5, length=1.5, width=0.5, color='0.2', pad=1.5)
    ax_right.tick_params(axis='y', labelleft=False, length=0)
    ax_right.spines['top'].set_visible(False)
    ax_right.spines['right'].set_visible(False)
    ax_right.spines['left'].set_visible(False)
    for spine in ['bottom']:
        ax_right.spines[spine].set_linewidth(0.5)
        ax_right.spines[spine].set_color('0.4')

    ax_main.set_xticks(range(len(unique_amines)))
    ax_main.set_xticklabels(unique_amines, rotation=90, ha='center', fontsize=6)

    ax_main.set_yticks(range(len(unique_body_sites)))
    ax_main.set_yticklabels(unique_body_sites, fontsize=6)

    ax_main.tick_params(axis='both', length=1.5, width=0.5, color='0.2', pad=1.5)
    ax_main.set_xlim(-0.5, len(unique_amines) - 0.5)
    ax_main.set_ylim(-0.5, len(unique_body_sites) - 0.5)

    ax_main.set_xlabel('Amines', fontsize=7, weight='semibold', labelpad=4.5)
    ax_main.set_ylabel('', fontsize=7, weight='semibold')

    ax_main.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
    ax_main.set_axisbelow(True)

    for spine in ax_main.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color('0.4')
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    class_legend_elements = [Patch(facecolor=class_to_color[cls], alpha=0.6,
                                   label=class_display_names.get(cls, cls))
                             for cls in classes_in_data]

    class_legend = ax_main.legend(handles=class_legend_elements,
                                  loc='upper left',
                                  bbox_to_anchor=(-0.14, -0.6),
                                  ncol=5,
                                  frameon=False,
                                  fontsize=5.5,
                                  title=None,
                                  title_fontsize=6.5,
                                  columnspacing=.8,
                                  handletextpad=0.25,
                                  handlelength=1.5)

    size_values = [0.05, 0.2, 1.0]
    size_legend_elements = [Line2D([0], [0], marker='o', color='w',
                                   markerfacecolor='0.5', markersize=np.sqrt(np.sqrt(s) / size_scale),
                                   alpha=0.6, label=f'{s}',
                                   markeredgewidth=0,
                                   linestyle='None')
                            for s in size_values]

    size_legend = ax_main.legend(handles=size_legend_elements,
                                 loc='upper left',
                                 bbox_to_anchor=(0.85, -0.555),
                                 ncol=3,
                                 frameon=False,
                                 fontsize=5.5,
                                 title='Matches per LCMS file',
                                 title_fontsize=6.5,
                                 columnspacing=1.5,
                                 handletextpad=0.8)
    size_legend._legend_box.sep = 6

    ax_main.add_artist(class_legend)

    plt.savefig('fig/large_heatmap/normalized/human_bodysite_scatterplot_amine_normalized.svg',
                format='svg', bbox_inches='tight', transparent=True)

    print('\nScatter plot saved to fig/large_heatmap/normalized/')
    print('  - human_bodysite_scatterplot_amine_normalized.svg')

    return fig, class_to_color, classes_in_data


if __name__ == "__main__":
    create_scatter_plot()
    print('\nAnalysis complete!')
