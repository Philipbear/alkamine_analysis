import pandas as pd
import numpy as np


# Rodent taxa to filter for
RODENT_TAXA = [
    '10088|Mus',
    '10090|Mus musculus',
    '10105|Mus minutoides',
    '10114|Rattus',
    '10116|Rattus norvegicus'
]


def preprocess_data():
    df = pd.read_csv('fig/large_heatmap/data/all_alkamines_new.csv', low_memory=False)

    # Filter for rodent species
    rodent_pattern = '|'.join([taxa.split('|')[1] for taxa in RODENT_TAXA])  # Extract species names
    df = df[(df['NCBITaxonomy'].str.contains(rodent_pattern, na=False, regex=True))].reset_index(drop=True)
    
    # Capitalize the first letter found, keeping the rest as is
    def capitalize_first_letter(s):
        if s.startswith('p-') or s.startswith('cis-'):
            return s
        for i, char in enumerate(s):
            if char.isalpha():
                return s[:i] + char.upper() + s[i+1:]
        return s
    
    df['Part1'] = df['Part1'].apply(capitalize_first_letter)
    df['Part2'] = df['Part2'].apply(capitalize_first_letter)
    df['name'] = df['Part1'] + '_' + df['Part2']
    
    df = df[['USI', 'name', 'Part1', 'Part2', 'Part1class', 'Part2class', 'UBERONBodyPartName']]
    
    print(f"Loaded {len(df)} matches, before deduplication")
    
    # dereplicate by name and USI, keeping the first occurrence
    df = df.drop_duplicates(subset=['name', 'USI'], keep='first')
        
    print(f"Loaded {len(df)} unique matches")
    
    # save preprocessed data
    df.to_csv('fig/large_heatmap/data/rodent_data.tsv', sep='\t', index=False)
    
    return df


def group_data_for_plot(df):
    """
    Group by Part1 (amine) and UBERONBodyPartName (body site)
    Count unique USIs for each combination
    """
    # Remove rows with missing body part names
    df_filtered = df[df['UBERONBodyPartName'].notna()].copy()
    
    print(f"Filtered to {len(df_filtered)} matches with valid body part names")
    
    # Group by Part1, Part1class and UBERONBodyPartName, count unique USIs
    grouped = df_filtered.groupby(['Part1', 'Part1class', 'UBERONBodyPartName'])['USI'].nunique().reset_index()
    grouped.columns = ['Part1', 'Part1class', 'UBERONBodyPartName', 'count']
    
    # Remove amines and body sites with zero total counts
    # Calculate sum of counts for each amine
    amine_totals = grouped.groupby('Part1')['count'].sum()
    print(f"Total unique USI counts per amine:\n{amine_totals}")
    valid_amines = amine_totals[amine_totals > 0].index
    
    # Calculate sum of counts for each body site
    bodysite_totals = grouped.groupby('UBERONBodyPartName')['count'].sum()
    print(f"Total unique USI counts per body site:\n{bodysite_totals}")
    valid_bodysites = bodysite_totals[bodysite_totals > 0].index
    
    # Filter to keep only valid amines and body sites
    grouped = grouped[grouped['Part1'].isin(valid_amines) & 
                      grouped['UBERONBodyPartName'].isin(valid_bodysites)]
    
    # Sort by Part1 and count for better visualization
    grouped = grouped.sort_values(['Part1', 'count'], ascending=[True, False])
    
    print(f"Grouped data: {len(grouped)} unique Part1-BodySite combinations")
    print(f"Number of unique amines: {grouped['Part1'].nunique()}")
    print(f"Number of unique body sites: {grouped['UBERONBodyPartName'].nunique()}")
    
    # Save grouped data
    grouped.to_csv('fig/large_heatmap/data/rodent_bodysite_amine_grouped.tsv', sep='\t', index=False)
    
    return grouped


def create_scatter_plot():
    """
    Create scatter plot with amines on x-axis, body sites on y-axis,
    and circle size representing count
    """
    
    grouped = pd.read_csv('fig/large_heatmap/data/rodent_bodysite_amine_grouped.tsv', sep='\t')
    
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 6
    
    # Define Part1 class order
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
    
    # Define display name mapping
    class_display_names = {
        'Amino acids (Proteinogenic; Modified)': 'Amino acids (Proteinogenic; modified)',
        'Peptides and peptide derivatives': 'Peptides',
        'Drugs and drug metabolites': 'Drugs',
        'Nucleic': 'Nucleic acid derivatives',
        'Other': 'Others'
    }
    
    # Get unique body sites (reverse order so A starts at top)
    unique_body_sites = sorted(grouped['UBERONBodyPartName'].unique(), reverse=True)
    
    # Sort amines by class (using class order), then by name
    # Create a mapping for class order
    class_order_map = {cls: i for i, cls in enumerate(part1_class_order)}
    
    # Add sort keys to grouped data
    grouped['class_order'] = grouped['Part1class'].map(lambda x: class_order_map.get(x, 999))
    
    # Get unique amines and sort by class order, then by name
    amine_info = grouped[['Part1', 'Part1class', 'class_order']].drop_duplicates()
    amine_info = amine_info.sort_values(['class_order', 'Part1'])
    unique_amines = amine_info['Part1'].tolist()
    
    # Create mapping for positions
    amine_to_pos = {amine: i for i, amine in enumerate(unique_amines)}
    bodysite_to_pos = {site: i for i, site in enumerate(unique_body_sites)}
    
    # Convert to positions
    x = grouped['Part1'].map(amine_to_pos)
    y = grouped['UBERONBodyPartName'].map(bodysite_to_pos)
    sizes = grouped['count']
    
    # Create color mapping for Part1class using class order
    # Combine two color palettes to have enough colors
    classes_in_data = [cls for cls in part1_class_order if cls in grouped['Part1class'].values]
    
    # Get colors from Set2 and Paired palettes
    palette1 = plt.cm.Set2.colors  # 8 colors
    palette2 = plt.cm.Paired.colors  # 12 colors
    combined_palette = list(palette1) + list(palette2)
    
    class_to_color = {cls: combined_palette[i] for i, cls in enumerate(classes_in_data)}
    colors = grouped['Part1class'].map(class_to_color)
    
    # Calculate totals for bar plots
    amine_totals = grouped.groupby('Part1')['count'].sum()
    amine_totals = amine_totals.reindex(unique_amines, fill_value=0)
    
    bodysite_totals = grouped.groupby('UBERONBodyPartName')['count'].sum()
    bodysite_totals = bodysite_totals.reindex(unique_body_sites, fill_value=0)
    
    # Create figure with gridspec for marginal bar plots
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(8, 3.4))
    gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[12, 1], height_ratios=[1, 5],
                          hspace=0.025, wspace=0.025)
    
    # Main scatter plot (bottom left)
    ax_main = fig.add_subplot(gs[1, 0])
    
    # Top bar plot (top left)
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    
    # Right bar plot (bottom right)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)
    
    # Create scatter plot without edges, colored by Part1class
    scatter = ax_main.scatter(x, y, s=sizes/2, alpha=0.6, c=colors, linewidth=0)
    
    # Top bar plot: total count per amine
    amine_colors = [class_to_color[amine_info[amine_info['Part1'] == amine]['Part1class'].values[0]] 
                    for amine in unique_amines]
    ax_top.bar(range(len(unique_amines)), amine_totals, color=amine_colors, alpha=0.6, width=0.8)
    ax_top.set_ylabel('Total count', fontsize=6)
    ax_top.tick_params(axis='y', labelsize=5, length=1.5, width=0.5, color='0.2', pad=1.5)
    ax_top.tick_params(axis='x', labelbottom=False, length=0)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['bottom'].set_visible(False)
    for spine in ['left']:
        ax_top.spines[spine].set_linewidth(0.5)
        ax_top.spines[spine].set_color('0.4')
    ax_top.set_xlim(-0.5, len(unique_amines) - 0.5)
    
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
    ax_main.set_xticks(range(len(unique_amines)))
    ax_main.set_xticklabels(unique_amines, rotation=90, ha='center', fontsize=6)
    
    ax_main.set_yticks(range(len(unique_body_sites)))
    ax_main.set_yticklabels(unique_body_sites, fontsize=6)
    
    # Adjust tick parameters: shorter length and lighter color
    ax_main.tick_params(axis='both', length=1.5, width=0.5, color='0.2', pad=1.5)
    
    # Set x-axis limits to reduce space
    ax_main.set_xlim(-0.5, len(unique_amines) - 0.5)
    
    # Set y-axis limits to prevent cutting scatter points at top
    ax_main.set_ylim(-0.5, len(unique_body_sites) - 0.5)
    
    # Set labels
    ax_main.set_xlabel('Amines', fontsize=7, weight='semibold', labelpad=4.5)
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
    
    # Amine class legend (5 columns)
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
    
    # Scatter size legend (horizontal, to the right of amine class legend)
    size_values = [50, 100, 500]
    size_legend_elements = [Line2D([0], [0], marker='o', color='w', 
                                   markerfacecolor='0.5', markersize=np.sqrt(s / 2), 
                                   alpha=0.6, label=f'{s}',
                                   markeredgewidth=0) 
                           for s in size_values]
    
    size_legend = ax_main.legend(handles=size_legend_elements, 
                                loc='upper left', 
                                bbox_to_anchor=(0.85, -0.555), 
                                ncol=3, 
                                frameon=False, 
                                fontsize=5.5, 
                                title='MS/MS matches', 
                                title_fontsize=6.5,
                                columnspacing=1.5,
                                handletextpad=0.8)
    
    # Add back the class legend (since creating size_legend removes it)
    ax_main.add_artist(class_legend)
    
    # Save figure
    plt.savefig('fig/large_heatmap/rodent_bodysite_scatterplot_amine.svg', 
                format='svg', bbox_inches='tight', transparent=True)
    
    print('\nScatter plot saved to fig/large_heatmap/')
    print(f'  - rodent_bodysite_scatterplot_amine.svg')
    
    return fig, class_to_color, classes_in_data


if __name__ == "__main__":
    # Preprocess data
    df = preprocess_data()
    
    # Group data for plotting
    group_data_for_plot(df)
    
    # Create scatter plot
    create_scatter_plot()
    
    print('\nAnalysis complete!')
