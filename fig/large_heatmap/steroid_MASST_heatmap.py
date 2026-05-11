import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


def preprocess_data():
    df = pd.read_csv('fig/large_heatmap/data/steroid_alkamines.csv', low_memory=False)
    df = df[['usi', 'NAME']]
    df = df.rename(columns={'usi': 'USI', 'NAME': 'name'})
    
    # Split name by '_' to get Part1 and Part2
    df['Part1'] = df['name'].apply(lambda x: x.split('_')[0] if '_' in str(x) else x)
    df['Part2'] = df['name'].apply(lambda x: '_'.join(x.split('_')[1:]) if '_' in str(x) else '')
    
    df = df[['USI', 'name', 'Part1', 'Part2']]
    
    print(f"Loaded {len(df)} matches, before deduplication")
    
    # dereplicate by name and USI, keeping the first occurrence
    df = df.drop_duplicates(subset=['name', 'USI'], keep='first')
    
    # Capitalize the first letter found, keeping the rest as is
    def capitalize_first_letter(s):
        for i, char in enumerate(s):
            if char.isalpha():
                return s[:i] + char.upper() + s[i+1:]
        return s
    
    df['Part1'] = df['Part1'].apply(capitalize_first_letter)
    # df['Part2'] = df['Part2'].apply(capitalize_first_letter)
    
    # name mapping
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
    
    print(f"Loaded {len(df)} unique matches")
    
    # save
    df.to_csv('fig/large_heatmap/data/steroid_refined.tsv', sep='\t', index=False)
    
    return df


def prepare_heatmap_data():
    """
    Prepare data for heatmap by creating count matrix
    """
    
    df = pd.read_csv('fig/large_heatmap/data/steroid_refined.tsv', sep='\t')
    
    # Count occurrences of each Part1-Part2 combination
    heatmap_df = df.groupby(['Part1', 'Part2']).size().reset_index(name='count')
    
    # Create pivot table
    pivot_df = heatmap_df.pivot_table(
        index='Part2',  # steroids as rows (y-axis)
        columns='Part1',  # amines as columns (x-axis)
        values='count', 
        fill_value=0
    )
    
    print(f"Heatmap dimensions: {pivot_df.shape}")
    
    # save
    heatmap_df.to_csv('fig/large_heatmap/data/steroid_heatmap_data.tsv', sep='\t', index=False)
    pivot_df.to_csv('fig/large_heatmap/data/steroid_pivot_data.tsv', sep='\t', index=False)
    
    return pivot_df, heatmap_df


def sort_molecules(pivot_df):
    """
    Sort molecules alphabetically (no class ordering)
    """
    # Add three additional steroids with all zeros
    additional_steroids = ['Estrone', 'Tetrahydrocortisol', 'Progesterone']
    for steroid in additional_steroids:
        if steroid not in pivot_df.index:
            pivot_df.loc[steroid] = 0
    
    # Sort rows (steroids) alphabetically
    pivot_df = pivot_df.reindex(sorted(pivot_df.index))
    
    # Sort columns (amines) alphabetically
    pivot_df = pivot_df.reindex(columns=sorted(pivot_df.columns), fill_value=0)
    
    print(f"Sorted heatmap dimensions: {pivot_df.shape}")
    print(f"Steroids (rows): {len(pivot_df.index)}")
    print(f"Amines (columns): {len(pivot_df.columns)}")
    
    return pivot_df


def create_heatmap(pivot_df, figsize=(4.45, 1.9)):
    """
    Create and save the heatmap
    """
    print('Creating heatmap...')
    
    # Apply log transformation (log10(x + 1) to handle zeros)
    log_data = np.log10(pivot_df.values + 1)
    
    # Set up the figure
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 5
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create the colormap
    colors = [(1, 1, 1), (0.78, 0.84, 0.94), (0.94, 0.69, 0.65)]
    positions = [0, 0.25, 1.0]
    cmap_wbr = LinearSegmentedColormap.from_list('white_blue_red', list(zip(positions, colors)), N=30)
    
    # Create heatmap
    im = ax.imshow(log_data, cmap=cmap_wbr, aspect='auto', interpolation='nearest')
    
    # Define compounds to mark with stars (lowercase for matching)
    # User will provide the list later
    starred_compounds = [
            'glycine_3-oxo-lca',
            'alanine_3-oxo-lca',
            'valine_3-oxo-lca',
            'valine_3-oxo-dca',
            'leucine_3-oxo-dca',
            'isoleucine_3-oxo-dca',
        
            'phenylalanine_3-oxo-lca',
            'isoleucine_3-oxo-lca',
            'leucine_3-oxo-lca',
            'n-acetylcadaverine_3-oxo-lca',
            
            'n-acetylcadaverine_3-oxo-dca',            
            'phenylalanine_3-oxo-dca',
            'phenethylamine_3-oxo-dca',
            'tyramine_3-oxo-dca',
            'tryptophan_3-oxo-dca',
            'tyrosine_3-oxo-dca',
            'lysine_3-oxo-dca',
            'putrescine_3-oxo-dca',
            
            'valine_cholestan-3-one',
            'leucine_cholestan-3-one',
            # 'phenylalanine_cholestan-3-one',
            'tyramine_cholestan-3-one',
            'n-acetylcadaverine_cholestan-3-one',
            'n-acetylputrescine_cholestan-3-one',
            'alanine_cholestan-3-one',
            'lysine_cholestan-3-one',
            'proline_cholestan-3-one',
            'glycine_cholestan-3-one',
            'tyrosine_cholestan-3-one'
        ]
    
    # Add stars to specified compounds
    found_compounds = set()
    for row_idx, steroid in enumerate(pivot_df.index):  # rows = steroids
        for col_idx, amine in enumerate(pivot_df.columns):  # columns = amines
            compound_key = f"{amine.lower()}_{steroid.lower()}"
            if compound_key in starred_compounds:
                # Add star marker
                ax.text(col_idx, row_idx + 0.5, '*', fontsize=9, ha='center', va='center', color='black')
                found_compounds.add(compound_key)
    
    # Print compounds that were not found
    not_found = set(starred_compounds) - found_compounds
    if not_found:
        print(f"\nWarning: {len(not_found)} starred compounds not found in heatmap:")
        for compound in sorted(not_found):
            print(f"  - {compound}")
    elif starred_compounds:
        print(f"\nAll {len(starred_compounds)} starred compounds found in heatmap!")
    
    # Define compounds to highlight with frames
    highlighted_compounds = [
        'phenylalanine_3-oxo-lca',
        # 'isoleucine_3-oxo-lca',
        # 'leucine_3-oxo-lca',
        'n-acetylcadaverine_3-oxo-lca',
        
        'n-acetylcadaverine_3-oxo-dca',            
        'phenylalanine_3-oxo-dca',
        'phenethylamine_3-oxo-dca',
        'tyramine_3-oxo-dca',
        'tryptophan_3-oxo-dca',
        # 'tyrosine_3-oxo-dca',
        'lysine_3-oxo-dca',
        'putrescine_3-oxo-dca',
        
        'valine_cholestan-3-one',
        'leucine_cholestan-3-one',
        # 'phenylalanine_cholestan-3-one',
        'tyramine_cholestan-3-one',
        # 'n-acetylcadaverine_cholestan-3-one',
        # 'n-acetylputrescine_cholestan-3-one',
        'alanine_cholestan-3-one',
        # 'lysine_cholestan-3-one',
        # 'proline_cholestan-3-one',
        'glycine_cholestan-3-one',
        # 'tyrosine_cholestan-3-one'
    ]
    
    # Add frames around highlighted compounds
    from matplotlib.patches import Rectangle
    highlighted_found = set()
    for row_idx, steroid in enumerate(pivot_df.index):  # rows = steroids
        for col_idx, amine in enumerate(pivot_df.columns):  # columns = amines
            compound_key = f"{amine.lower()}_{steroid.lower()}"
            if compound_key in highlighted_compounds:
                # Add rectangle frame around the cell
                rect = Rectangle((col_idx - 0.5, row_idx - 0.5), 1, 1, 
                               linewidth=.8, edgecolor='#FF6B00', facecolor='none', zorder=10)
                ax.add_patch(rect)
                highlighted_found.add(compound_key)
    
    # Print highlighted compounds status
    not_found_highlighted = set(highlighted_compounds) - highlighted_found
    if not_found_highlighted:
        print(f"\nWarning: {len(not_found_highlighted)} highlighted compounds not found:")
        for compound in sorted(not_found_highlighted):
            print(f"  - {compound}")
    elif highlighted_compounds:
        print(f"\nAll {len(highlighted_compounds)} highlighted compounds found and framed!")
    
    # Set ticks and labels
    ax.set_xticks(range(len(pivot_df.columns)))
    ax.set_xticklabels(pivot_df.columns, rotation=90, ha='center', fontsize=6, color='0.2')
    
    # Move x-axis labels to top
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    
    ax.set_yticks(range(len(pivot_df.index)))
    ax.set_yticklabels(pivot_df.index, fontsize=6, color='0.2')
    
    # Remove bottom ticks, keep top and left ticks
    ax.tick_params(left=True, bottom=False, top=True, pad=1.5, width=0.5, color='0.4', length=1)
    
    # Style the heatmap frame (spines)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_edgecolor('0.2')
    
    # Add small horizontal colorbar at bottom left
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(ax, width="18%", height="8%", loc='lower left', 
                      bbox_to_anchor=(-0.16, 1.6, 0.6, 0.6), bbox_transform=ax.transAxes)
    cbar = plt.colorbar(im, cax=axins, orientation='horizontal')
    cbar.set_label('Log(Count + 1)', fontsize=5., color='0.2', labelpad=1)
    
    # Set colorbar ticks to show 0, middle, and max values
    vmin, vmax = log_data.min(), log_data.max()
    cbar.set_ticks([vmin, vmax/2, vmax])
    cbar.set_ticklabels([f'{vmin:.1f}', f'{vmax/2:.1f}', f'{vmax:.1f}'])
    
    cbar.ax.tick_params(labelsize=4.75, width=0.25, color='0.2', length=1.2, pad=1)
    cbar.outline.set_linewidth(0.25)
    cbar.outline.set_edgecolor('0.2')
    
    # Set labels
    ax.set_xlabel('', fontsize=7, labelpad=0)
    ax.set_ylabel('', fontsize=7, labelpad=0)
    
    # Add legend for star markers and highlighted compounds (only if lists are not empty)
    if starred_compounds or highlighted_compounds:
        from matplotlib.patches import Patch
        from matplotlib.lines import Line2D
        
        legend_elements = []
        if starred_compounds:
            legend_elements.append(
                Line2D([0], [0], marker='*', color='w', markerfacecolor='black', 
                       markersize=9, label='Validated by RT')
            )
        if highlighted_compounds:
            legend_elements.append(
                Patch(facecolor='none', edgecolor='#FF6B00', linewidth=0.8, 
                      label='Validated by CCS')
            )
        
        if legend_elements:
            legend = ax.legend(handles=legend_elements, loc='lower left', 
                              bbox_to_anchor=(-0.21, 1.02), frameon=False, 
                              fontsize=5., handlelength=1., handleheight=1., 
                              handletextpad=0.5, borderpad=0.2)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    plt.savefig('fig/large_heatmap/steroid_MASST_heatmap.svg', format='svg', bbox_inches='tight', transparent=True)
    
    print('Heatmap saved to fig/large_heatmap/')
    
    return fig


if __name__ == "__main__":
    # Preprocess data
    preprocess_data()
    
    # Prepare heatmap data
    pivot_df, heatmap_df = prepare_heatmap_data()
    
    # Sort molecules alphabetically
    pivot_df = sort_molecules(pivot_df)
    
    # Create heatmap
    fig = create_heatmap(pivot_df)
    
    print('\nAnalysis complete!')
