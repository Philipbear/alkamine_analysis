import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def preprocess_data():
    """
    Preprocess steroid_alkamines.csv for body site scatter plot
    """
    df = pd.read_csv('fig/large_heatmap/data/steroid_alkamines.csv', low_memory=False)
    
    # # Filter for human data
    # df = df[df['Taxaname_file'].str.contains('Homo sapiens', na=False)].reset_index(drop=True)
    
    # Parse NAME into Part1 (amine) and Part2 (steroid)
    df['Part1'] = df['NAME'].apply(lambda x: x.split('_')[0] if '_' in str(x) else x)
    df['Part2'] = df['NAME'].apply(lambda x: '_'.join(x.split('_')[1:]) if '_' in str(x) else '')
    
    # Rename columns for clarity
    df = df.rename(columns={'usi': 'USI'})
    df = df[['USI', 'NAME', 'Part1', 'Part2', 'UBERONBodyPartName']]
    
    # Dereplicate by NAME and USI
    df = df.drop_duplicates(subset=['NAME', 'USI'], keep='first')
    
    print(f"After deduplication: {len(df)} unique matches")
    
    # Capitalize first letter
    def capitalize_first_letter(s):
        if pd.isna(s):
            return s
        s = str(s)
        if s.startswith('p-') or s.startswith('cis-'):
            return s
        for i, char in enumerate(s):
            if char.isalpha():
                return s[:i] + char.upper() + s[i+1:]
        return s
    
    df['Part1'] = df['Part1'].apply(capitalize_first_letter)
    # df['Part2'] = df['Part2'].apply(capitalize_first_letter)
    
    # Apply name mapping for steroids (from steroid_MASST_heatmap.py)
    name_mapping = {
        '3-keto-cholic acid': '3-oxo-CA',
        '3-keto-deoxycholic acid': '3-oxo-DCA',
        '3-oxo-glyco_DCA': '3-oxo-glycoDCA',
        '3-oxo-glycolithocholate': '3-oxo-glycoLCA',
        '5-Cholestan-3-one': 'Cholestan-3-one',
        '5-Cholesten-3-one': '5-cholesten-3-one',
        'Dehydrolithocholicacid': '3-oxo-LCA'
    }
    df['Part2'] = df['Part2'].replace(name_mapping)
    
    # Label missing body sites as 'Missing value' instead of filtering them out
    df['UBERONBodyPartName'] = df['UBERONBodyPartName'].fillna('missing value')
    
    print(f"Total matches: {len(df)}")
    
    # Save preprocessed data
    df.to_csv('fig/large_heatmap/data/steroid_data.tsv', sep='\t', index=False)
    
    return df


def group_data_for_plot(df):
    """
    Group by Part1 (amine) and Part2 (steroid) separately with UBERONBodyPartName (body site)
    Count unique USIs for each combination
    Returns two dataframes: one for amines, one for steroids
    """
    # Group by Part1 (amine) and UBERONBodyPartName
    grouped_amine = df.groupby(['Part1', 'UBERONBodyPartName'])['USI'].nunique().reset_index()
    grouped_amine.columns = ['molecule', 'UBERONBodyPartName', 'count']
    grouped_amine['type'] = 'amine'
    
    # Group by Part2 (steroid) and UBERONBodyPartName
    grouped_steroid = df.groupby(['Part2', 'UBERONBodyPartName'])['USI'].nunique().reset_index()
    grouped_steroid.columns = ['molecule', 'UBERONBodyPartName', 'count']
    grouped_steroid['type'] = 'steroid'
    
    # Remove molecules with zero total counts
    for grouped in [grouped_amine, grouped_steroid]:
        mol_totals = grouped.groupby('molecule')['count'].sum()
        valid_mols = mol_totals[mol_totals > 0].index
        grouped = grouped[grouped['molecule'].isin(valid_mols)]
    
    # Combine for body site filtering
    combined = pd.concat([grouped_amine, grouped_steroid], ignore_index=True)
    
    # Sort each group alphabetically
    grouped_amine = grouped_amine.sort_values(['molecule', 'count'], ascending=[True, False])
    grouped_steroid = grouped_steroid.sort_values(['molecule', 'count'], ascending=[True, False])
    
    print(f"Amine grouped data: {len(grouped_amine)} unique Amine-BodySite combinations")
    print(f"Number of unique amines: {grouped_amine['molecule'].nunique()}")
    print(f"Steroid grouped data: {len(grouped_steroid)} unique Steroid-BodySite combinations")
    print(f"Number of unique steroids: {grouped_steroid['molecule'].nunique()}")
    print(f"Number of unique body sites: {combined['UBERONBodyPartName'].nunique()}")
    
    # Save grouped data
    combined.to_csv('fig/large_heatmap/data/bodysite_steroid_grouped.tsv', sep='\t', index=False)
    
    return grouped_amine, grouped_steroid


def create_scatter_plot(grouped_amine, grouped_steroid):
    """
    Create scatter plot with amines and steroids on x-axis (concatenated), 
    body sites on y-axis, and circle size representing count
    """
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 5
    
    # Sort amines and steroids alphabetically
    unique_amines = sorted(grouped_amine['molecule'].unique())
    unique_steroids = sorted(grouped_steroid['molecule'].unique())
    
    # Combine all molecules: steroids first, then amines
    all_molecules = unique_steroids + unique_amines
    n_steroids = len(unique_steroids)
    n_amines = len(unique_amines)
    
    # Get unique body sites from both datasets (reverse order so A starts at top)
    all_body_sites = set(grouped_amine['UBERONBodyPartName'].unique()) | set(grouped_steroid['UBERONBodyPartName'].unique())
    # Sort with 'Missing value' at the bottom (position 0)
    unique_body_sites = sorted([s for s in all_body_sites if s != 'missing value'], reverse=True)
    if 'missing value' in all_body_sites:
        unique_body_sites.insert(0, 'missing value')
    
    # Create mapping for positions
    mol_to_pos = {mol: i for i, mol in enumerate(all_molecules)}
    bodysite_to_pos = {site: i for i, site in enumerate(unique_body_sites)}
    
    # Combine grouped data
    combined = pd.concat([grouped_amine, grouped_steroid], ignore_index=True)
    
    # Convert to positions
    x = combined['molecule'].map(mol_to_pos)
    y = combined['UBERONBodyPartName'].map(bodysite_to_pos)
    sizes = combined['count']
    
    # Assign colors by type
    amine_color = '#5A9BD5'  # Blue for amines
    steroid_color = '#ED7D31'  # Orange for steroids
    colors = combined['type'].map({'amine': amine_color, 'steroid': steroid_color})
    
    # Calculate totals for bar plots
    mol_totals = combined.groupby('molecule')['count'].sum()
    mol_totals = mol_totals.reindex(all_molecules, fill_value=0)
    
    bodysite_totals = combined.groupby('UBERONBodyPartName')['count'].sum()
    bodysite_totals = bodysite_totals.reindex(unique_body_sites, fill_value=0)
    
    # Create figure with gridspec for marginal bar plots
    fig = plt.figure(figsize=(4.9, 1.8))
    gs = gridspec.GridSpec(2, 2, figure=fig, width_ratios=[16, 1], height_ratios=[2, 9],
                          hspace=0.025, wspace=0.025)
    
    # Main scatter plot (bottom left)
    ax_main = fig.add_subplot(gs[1, 0])
    
    # Top bar plot (top left)
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    
    # Right bar plot (bottom right)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)
    
    # Create scatter plot colored by type (sqrt scale for sizes)
    size_scale = 5  # scaling factor for sqrt sizes
    scatter = ax_main.scatter(x, y, s=np.sqrt(sizes) * size_scale, alpha=0.6, c=colors, linewidth=0)
    
    # Add vertical line separator between steroids and amines
    separator_x = n_steroids - 0.5
    ax_main.axvline(x=separator_x, color='0.5', linestyle='--', linewidth=0.8, alpha=0.6)
    ax_top.axvline(x=separator_x, color='0.5', linestyle='--', linewidth=0.8, alpha=0.6)
    
    # Top bar plot: total count per molecule, colored by type
    bar_colors = [steroid_color] * n_steroids + [amine_color] * n_amines
    ax_top.bar(range(len(all_molecules)), mol_totals, color=bar_colors, alpha=0.6, width=0.8)
    ax_top.set_ylabel('Total\ncount', fontsize=6)
    ax_top.tick_params(axis='y', labelsize=5, length=1.5, width=0.5, color='0.2', pad=1.5)
    ax_top.tick_params(axis='x', labelbottom=False, length=0)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['bottom'].set_visible(False)
    for spine in ['left']:
        ax_top.spines[spine].set_linewidth(0.5)
        ax_top.spines[spine].set_color('0.4')
    ax_top.set_xlim(-0.5, len(all_molecules) - 0.5)
    
    # Right bar plot: total count per body site
    ax_right.barh(range(len(unique_body_sites)), bodysite_totals, color='0.6', alpha=0.6, height=0.8)
    ax_right.set_xlabel('Total\ncount', fontsize=6)
    ax_right.tick_params(axis='x', labelsize=5, length=1.5, width=0.5, color='0.2', pad=1.5)
    ax_right.tick_params(axis='y', labelleft=False, length=0)
    ax_right.spines['top'].set_visible(False)
    ax_right.spines['right'].set_visible(False)
    ax_right.spines['left'].set_visible(False)
    for spine in ['bottom']:
        ax_right.spines[spine].set_linewidth(0.5)
        ax_right.spines[spine].set_color('0.4')
    
    # Set ticks and labels for main scatter plot
    ax_main.set_xticks(range(len(all_molecules)))
    ax_main.set_xticklabels(all_molecules, rotation=90, ha='center', fontsize=6)
    
    ax_main.set_yticks(range(len(unique_body_sites)))
    ax_main.set_yticklabels(unique_body_sites, fontsize=6)
    
    # Adjust tick parameters
    ax_main.tick_params(axis='both', length=1.5, width=0.5, color='0.2', pad=1.5)
    
    # Set axis limits
    ax_main.set_xlim(-0.5, len(all_molecules) - 0.5)
    ax_main.set_ylim(-0.5, len(unique_body_sites) - 0.5)
    
    # Set labels
    ax_main.set_xlabel('', fontsize=7, weight='semibold', labelpad=4.5)
    ax_main.set_ylabel('', fontsize=7, weight='semibold')
    
    # Add grid
    ax_main.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
    ax_main.set_axisbelow(True)
    
    # Style spines
    for spine in ax_main.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color('0.4')
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)
    
    # Add legends
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    
    # Type legend (steroid vs amine)
    type_legend_elements = [
        Patch(facecolor=steroid_color, alpha=0.6, label='Steroids'),
        Patch(facecolor=amine_color, alpha=0.6, label='Amines')
    ]
    
    type_legend = ax_main.legend(handles=type_legend_elements, 
                                loc='upper left', 
                                bbox_to_anchor=(0.25, -0.72), 
                                ncol=2, 
                                frameon=False, 
                                fontsize=6,
                                columnspacing=0.8,
                                handletextpad=0.25,
                                handlelength=1.5)
    
    # Size legend (sqrt scale)
    size_values = [10, 50, 200]
    size_legend_elements = [Line2D([0], [0], marker='o', color='w', 
                                   markerfacecolor='0.5', markersize=np.sqrt(np.sqrt(s) * size_scale), 
                                   alpha=0.6, label=f'{s}',
                                   markeredgewidth=0) 
                           for s in size_values]
    
    size_legend = ax_main.legend(handles=size_legend_elements, 
                                loc='upper left', 
                                bbox_to_anchor=(0.75, -0.6), 
                                ncol=3, 
                                frameon=False, 
                                fontsize=5.5, 
                                title='MS/MS matches', 
                                title_fontsize=6.,
                                columnspacing=0.9,
                                handletextpad=0.5)
    
    # Add back the type legend
    ax_main.add_artist(type_legend)
    
    # Save figure
    plt.savefig('fig/large_heatmap/bodysite_scatterplot_steroid.svg', 
                format='svg', bbox_inches='tight', transparent=True)
    
    return fig


if __name__ == "__main__":
    # Preprocess data
    df = preprocess_data()
    
    # Group data for plotting
    grouped_amine, grouped_steroid = group_data_for_plot(df)
    
    # Create scatter plot
    fig = create_scatter_plot(grouped_amine, grouped_steroid)
    
    print('\nAnalysis complete!')
