import pandas as pd
import numpy as np


def group_data_for_plot(df):
    """
    Group by Part2 (ketone/aldehyde) and UBERONBodyPartName (body site)
    Count unique USIs for each combination
    """
    # Remove rows with missing body part names
    df_filtered = df[df['UBERONBodyPartName'].notna()].copy()
    
    print(f"Filtered to {len(df_filtered)} matches with valid body part names")
    
    # Group by Part2, Part2class and UBERONBodyPartName, count unique USIs
    grouped = df_filtered.groupby(['Part2', 'Part2class', 'UBERONBodyPartName'])['USI'].nunique().reset_index()
    grouped.columns = ['Part2', 'Part2class', 'UBERONBodyPartName', 'count']
    
    # Remove Part2 and body sites with zero total counts
    # Calculate sum of counts for each Part2
    part2_totals = grouped.groupby('Part2')['count'].sum()
    print(f"Total unique USI counts per Part2:\n{part2_totals}")
    valid_part2 = part2_totals[part2_totals > 0].index
    
    # Calculate sum of counts for each body site
    bodysite_totals = grouped.groupby('UBERONBodyPartName')['count'].sum()
    print(f"Total unique USI counts per body site:\n{bodysite_totals}")
    valid_bodysites = bodysite_totals[bodysite_totals > 0].index
    
    # Filter to keep only valid Part2 and body sites
    grouped = grouped[grouped['Part2'].isin(valid_part2) & 
                      grouped['UBERONBodyPartName'].isin(valid_bodysites)]
    
    # Sort by Part2 and count for better visualization
    grouped = grouped.sort_values(['Part2', 'count'], ascending=[True, False])
    
    print(f"Grouped data: {len(grouped)} unique Part2-BodySite combinations")
    print(f"Number of unique Part2: {grouped['Part2'].nunique()}")
    print(f"Number of unique body sites: {grouped['UBERONBodyPartName'].nunique()}")
    
    # Save grouped data
    grouped.to_csv('fig/large_heatmap/data/human_bodysite_part2_grouped.tsv', sep='\t', index=False)
    
    return grouped


def create_scatter_plot():
    """
    Create scatter plot with ketones/aldehydes on x-axis, body sites on y-axis,
    and circle size representing count
    """
    
    grouped = pd.read_csv('fig/large_heatmap/data/human_bodysite_part2_grouped.tsv', sep='\t')
    
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 6
    
    # Define Part2 class order (ketone first, then aldehyde)
    part2_class_order = ['ketone', 'aldehyde']
    
    # Define display name mapping
    class_display_names = {
        'ketone': 'Ketone',
        'aldehyde': 'Aldehyde'
    }
    
    # Get unique body sites (reverse order so A starts at top)
    unique_body_sites = sorted(grouped['UBERONBodyPartName'].unique(), reverse=True)
    
    # Sort Part2 by class (ketone first, then aldehyde), then by name
    # Create a mapping for class order
    class_order_map = {cls: i for i, cls in enumerate(part2_class_order)}
    
    # Add sort keys to grouped data
    grouped['class_order'] = grouped['Part2class'].map(lambda x: class_order_map.get(x, 999))
    
    # Get unique Part2 and sort by class order, then by name
    part2_info = grouped[['Part2', 'Part2class', 'class_order']].drop_duplicates()
    part2_info = part2_info.sort_values(['class_order', 'Part2'])
    unique_part2 = part2_info['Part2'].tolist()
    
    # Create mapping for positions
    part2_to_pos = {part2: i for i, part2 in enumerate(unique_part2)}
    bodysite_to_pos = {site: i for i, site in enumerate(unique_body_sites)}
    
    # Convert to positions
    x = grouped['Part2'].map(part2_to_pos)
    y = grouped['UBERONBodyPartName'].map(bodysite_to_pos)
    sizes = grouped['count']
    
    # Create color mapping for Part2class
    classes_in_data = [cls for cls in part2_class_order if cls in grouped['Part2class'].values]
    
    # Use distinct colors for ketone and aldehyde
    colors_map = {
        'ketone': '#8dd3c7',  # Teal/cyan
        'aldehyde': '#fb8072'  # Salmon/red
    }
    
    class_to_color = {cls: colors_map[cls] for cls in classes_in_data}
    colors = grouped['Part2class'].map(class_to_color)
    
    # Calculate totals for bar plots
    part2_totals = grouped.groupby('Part2')['count'].sum()
    part2_totals = part2_totals.reindex(unique_part2, fill_value=0)
    
    bodysite_totals = grouped.groupby('UBERONBodyPartName')['count'].sum()
    bodysite_totals = bodysite_totals.reindex(unique_body_sites, fill_value=0)
    
    # Create figure with gridspec for marginal bar plots
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(8, 3.))
    gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[12, 1], height_ratios=[1, 5],
                          hspace=0.025, wspace=0.025)
    
    # Main scatter plot (bottom left)
    ax_main = fig.add_subplot(gs[1, 0])
    
    # Top bar plot (top left)
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    
    # Right bar plot (bottom right)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)
    
    # Create scatter plot without edges, colored by Part2class
    # Note: scatter `s` is marker area in pt^2
    size_scale = 2
    scatter = ax_main.scatter(x, y, s=sizes / size_scale, alpha=0.6, c=colors, linewidth=0)
    
    # Top bar plot: total count per Part2
    part2_colors = [class_to_color[part2_info[part2_info['Part2'] == part2]['Part2class'].values[0]] 
                    for part2 in unique_part2]
    ax_top.bar(range(len(unique_part2)), part2_totals, color=part2_colors, alpha=0.6, width=0.8)
    ax_top.set_ylabel('Total count', fontsize=6)
    ax_top.tick_params(axis='y', labelsize=5, length=1.5, width=0.5, color='0.2', pad=1.5)
    ax_top.tick_params(axis='x', labelbottom=False, length=0)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['bottom'].set_visible(False)
    for spine in ['left']:
        ax_top.spines[spine].set_linewidth(0.5)
        ax_top.spines[spine].set_color('0.4')
    ax_top.set_xlim(-0.5, len(unique_part2) - 0.5)
    
    # Right bar plot: total count per body site
    ax_right.barh(range(len(unique_body_sites)), bodysite_totals, color='0.6', alpha=0.6, height=0.8)
    ax_right.set_xlabel('Total count', fontsize=6)
    ax_right.tick_params(axis='x', labelsize=5, length=1.5, width=0.5, color='0.2', pad=1.5)
    ax_right.tick_params(axis='y', labelleft=False, length=0)
    ax_right.spines['top'].set_visible(False)
    ax_right.spines['right'].set_visible(False)
    ax_right.spines['left'].set_visible(False)
    for spine in ['bottom']:
        ax_right.spines[spine].set_linewidth(0.5)
        ax_right.spines[spine].set_color('0.4')
    
    # Set ticks and labels for main scatter plot
    ax_main.set_xticks(range(len(unique_part2)))
    ax_main.set_xticklabels(unique_part2, rotation=90, ha='center', fontsize=6)
    
    ax_main.set_yticks(range(len(unique_body_sites)))
    ax_main.set_yticklabels(unique_body_sites, fontsize=6)
    
    # Adjust tick parameters: shorter length and lighter color
    ax_main.tick_params(axis='both', length=1.5, width=0.5, color='0.2', pad=1.5)
    
    # Set x-axis limits to reduce space
    ax_main.set_xlim(-0.5, len(unique_part2) - 0.1)
    
    # Set y-axis limits to prevent cutting scatter points at top (add extra margin)
    ax_main.set_ylim(-0.5, len(unique_body_sites) - 0.2)
    
    # Set labels
    ax_main.set_xlabel('Ketones and Aldehydes', fontsize=7, weight='semibold')
    ax_main.set_ylabel('', fontsize=7, weight='semibold')
    
    # Add grid
    ax_main.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
    ax_main.set_axisbelow(True)
    
    # Style spines, remove top and right spines
    for spine in ax_main.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color('0.4')
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)
    
    # Add legends below the main plot
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    
    # Part2 class legend (2 items - ketone and aldehyde)
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
    
    # Scatter size legend (horizontal, to the right of class legend)
    size_values = [50, 100, 500]
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
                                title='MS/MS matches', 
                                title_fontsize=6.5,
                                columnspacing=1.5,
                                handletextpad=0.8)
    size_legend._legend_box.sep = 6

    # Add back the class legend (since creating size_legend removes it)
    ax_main.add_artist(class_legend)
    
    # Save figure
    plt.savefig('fig/large_heatmap/human_bodysite_scatterplot_ketone_aldehyde.svg', 
                format='svg', bbox_inches='tight', transparent=True)
    
    print('\nScatter plot saved to fig/large_heatmap/')
    print(f'  - human_bodysite_scatterplot_ketone_aldehyde.svg')
    
    return fig, class_to_color, classes_in_data


if __name__ == "__main__":
    df = pd.read_csv('fig/large_heatmap/data/human_data.tsv', sep='\t')
    
    # Group data for plotting
    grouped = group_data_for_plot(df)
    
    # Create scatter plot
    create_scatter_plot()
    
    print('\nAnalysis complete!')
