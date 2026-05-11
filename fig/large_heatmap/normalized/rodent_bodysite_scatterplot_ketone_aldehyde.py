import pandas as pd
import numpy as np


def create_scatter_plot():
    """
    Normalized version: circle size = USI matches / number of LCMS files for that body site.
    """
    grouped = pd.read_csv('fig/large_heatmap/data/rodent_bodysite_part2_grouped.tsv', sep='\t')
    file_counts = pd.read_csv('fig/large_heatmap/normalized/data/rodent_bodysite_file_count.tsv', sep='\t')

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

    part2_class_order = ['ketone', 'aldehyde']

    class_display_names = {
        'ketone': 'Ketone',
        'aldehyde': 'Aldehyde'
    }

    unique_body_sites = sorted(grouped['UBERONBodyPartName'].unique(), reverse=True)

    class_order_map = {cls: i for i, cls in enumerate(part2_class_order)}
    grouped['class_order'] = grouped['Part2class'].map(lambda x: class_order_map.get(x, 999))

    part2_info = grouped[['Part2', 'Part2class', 'class_order']].drop_duplicates()
    part2_info = part2_info.sort_values(['class_order', 'Part2'])
    unique_part2 = part2_info['Part2'].tolist()

    part2_to_pos = {part2: i for i, part2 in enumerate(unique_part2)}
    bodysite_to_pos = {site: i for i, site in enumerate(unique_body_sites)}

    x = grouped['Part2'].map(part2_to_pos)
    y = grouped['UBERONBodyPartName'].map(bodysite_to_pos)
    sizes = grouped['count_norm']

    classes_in_data = [cls for cls in part2_class_order if cls in grouped['Part2class'].values]

    colors_map = {
        'ketone': '#8dd3c7',
        'aldehyde': '#fb8072'
    }
    class_to_color = {cls: colors_map[cls] for cls in classes_in_data}
    colors = grouped['Part2class'].map(class_to_color)

    # Normalized totals for bar plots
    part2_totals = grouped.groupby('Part2')['count_norm'].sum()
    part2_totals = part2_totals.reindex(unique_part2, fill_value=0)

    bodysite_totals = grouped.groupby('UBERONBodyPartName')['count_norm'].sum()
    bodysite_totals = bodysite_totals.reindex(unique_body_sites, fill_value=0)

    fig = plt.figure(figsize=(8, 3.4))
    gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[12, 1], height_ratios=[1, 5],
                           hspace=0.025, wspace=0.025)

    ax_main = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)

    size_scale = 0.004
    scatter = ax_main.scatter(x, y, s=sizes / size_scale, alpha=0.6, c=colors, linewidth=0)

    part2_colors = [class_to_color[part2_info[part2_info['Part2'] == part2]['Part2class'].values[0]]
                    for part2 in unique_part2]
    ax_top.bar(range(len(unique_part2)), part2_totals, color=part2_colors, alpha=0.6, width=0.8)
    ax_top.set_ylabel('Normalized\ncount', fontsize=6)
    ax_top.tick_params(axis='y', labelsize=5, length=1.5, width=0.5, color='0.2', pad=1.5)
    ax_top.tick_params(axis='x', labelbottom=False, length=0)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['bottom'].set_visible(False)
    for spine in ['left']:
        ax_top.spines[spine].set_linewidth(0.5)
        ax_top.spines[spine].set_color('0.4')
    ax_top.set_xlim(-0.5, len(unique_part2) - 0.5)

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

    ax_main.set_xticks(range(len(unique_part2)))
    ax_main.set_xticklabels(unique_part2, rotation=90, ha='center', fontsize=6)

    ax_main.set_yticks(range(len(unique_body_sites)))
    ax_main.set_yticklabels(unique_body_sites, fontsize=6)

    ax_main.tick_params(axis='both', length=1.5, width=0.5, color='0.2', pad=1.5)
    ax_main.set_xlim(-0.5, len(unique_part2) - 0.1)
    ax_main.set_ylim(-0.5, len(unique_body_sites) - 0.2)

    ax_main.set_xlabel('Ketones and Aldehydes', fontsize=7, weight='semibold')
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
                                  bbox_to_anchor=(-0.12, -0.41),
                                  ncol=2,
                                  frameon=False,
                                  fontsize=5.5,
                                  title=None,
                                  title_fontsize=6.5,
                                  columnspacing=1,
                                  handletextpad=0.25,
                                  handlelength=1.5)

    size_values = [0.05, 0.2, 1.0]
    size_legend_elements = [Line2D([0], [0], marker='o', color='w',
                                   markerfacecolor='0.5', markersize=np.sqrt(s / size_scale),
                                   alpha=0.6, label=f'{s}',
                                   markeredgewidth=0,
                                   linestyle='None')
                            for s in size_values]

    size_legend = ax_main.legend(handles=size_legend_elements,
                                 loc='upper left',
                                 bbox_to_anchor=(-0.13, -0.5),
                                 ncol=3,
                                 frameon=False,
                                 fontsize=5.5,
                                 title='Matches per LCMS file',
                                 title_fontsize=6.5,
                                 columnspacing=1.5,
                                 handletextpad=0.8)
    size_legend._legend_box.sep = 6

    ax_main.add_artist(class_legend)

    plt.savefig('fig/large_heatmap/normalized/rodent_bodysite_scatterplot_ketone_aldehyde_normalized.svg',
                format='svg', bbox_inches='tight', transparent=True)

    print('\nScatter plot saved to fig/large_heatmap/normalized/')
    print('  - rodent_bodysite_scatterplot_ketone_aldehyde_normalized.svg')

    return fig, class_to_color, classes_in_data


if __name__ == "__main__":
    create_scatter_plot()
    print('\nAnalysis complete!')
