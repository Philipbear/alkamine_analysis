import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import LinearSegmentedColormap


def preprocess_data():
    df = pd.read_csv('fig/large_heatmap/data/all_alkamines_new.csv', low_memory=False)
    df = df[['USI', 'Part1', 'Part2', 'Part1class', 'Part2class', 'NCBITaxonomy', 'Taxaname_file', 
             'UBERONBodyPartName', 'name', 'DOIDCommonName']]
    
    # Part 1 and Part, strip
    df['Part1'] = df['Part1'].str.strip()
    df['Part2'] = df['Part2'].str.strip()
    
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
    df = df[['USI', 'name', 'Part1', 'Part2', 'Part1class', 'Part2class']]
    
    print(f"Loaded {len(df)} matches, before deduplication")
    
    # dereplicate by name and USI, keeping the first occurrence
    df = df.drop_duplicates(subset=['name', 'USI'], keep='first').reset_index(drop=True)
    
    # remove metformin matches
    df = df[~df['Part1'].str.contains('Metformin', case=False)].reset_index(drop=True)
    
    print(f"Loaded {len(df)} unique matches")
    
    # print nunique names
    print(f"Number of unique Part1 molecules: {df['Part1'].nunique()}")
    print(f"Number of unique Part2 molecules: {df['Part2'].nunique()}")
    print(f"Number of unique names: {df['name'].nunique()}")
    
    # # remove all bile acids
    # bile_acids = [
    #     '12-Oxo-glyco-cholic acid', '3-Keto-deoxycholic acid', '3-Keto-cholic acid', 
    #     '3-Oxo-glycolithocholate', '3-Oxo-glyco-deoxycholic acid', 'Dehydrolithocholicacid'
    # ]
    # df = df[~df['Part2'].isin(bile_acids)].reset_index(drop=True)
    
    # save
    df.to_csv('fig/large_heatmap/data/all_refined.tsv', sep='\t', index=False)
    
    return df


def prepare_heatmap_data():
    """
    Prepare data for heatmap by creating count matrix
    """
    
    df = pd.read_csv('fig/large_heatmap/data/all_refined.tsv', sep='\t')
    
    # Count occurrences of each Part1-Part2 combination
    heatmap_df = df.groupby(['Part1', 'Part2', 'Part1class', 'Part2class']).size().reset_index(name='count')
    
    # Create pivot table
    pivot_df = heatmap_df.pivot_table(
        index=['Part1', 'Part1class'], 
        columns='Part2', 
        values='count', 
        fill_value=0
    )
    
    print(f"Heatmap dimensions: {pivot_df.shape}")
    
    # save
    heatmap_df.to_csv('fig/large_heatmap/data/heatmap_data.tsv', sep='\t', index=False)
    pivot_df.to_csv('fig/large_heatmap/data/pivot_data.tsv', sep='\t', index=False)
    
    return pivot_df, heatmap_df


def sort_molecules(pivot_df, heatmap_df):
    """
    Sort molecules according to specified order
    """
    # Define Part1 class order (using actual dataframe values)
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
    
    # Sort Part1 (rows) by class order
    row_order = []
    for class_name in part1_class_order:
        # Get all Part1 molecules in this class
        molecules_in_class = pivot_df[pivot_df.index.get_level_values('Part1class') == class_name].index.get_level_values('Part1').unique()
        row_order.extend([(mol, class_name) for mol in sorted(molecules_in_class)])
    
    # Reindex rows according to the order
    pivot_df = pivot_df.reindex(row_order)
    
    # Sort Part2 (columns) by class: ketones first, then aldehydes
    part2_info = heatmap_df[['Part2', 'Part2class']].drop_duplicates()
    part2_ketones = sorted(part2_info[part2_info['Part2class'] == 'ketone']['Part2'].tolist())
    part2_aldehydes = sorted(part2_info[part2_info['Part2class'] == 'aldehyde']['Part2'].tolist())
    col_order = part2_ketones + part2_aldehydes
    
    # Reindex columns
    pivot_df = pivot_df.reindex(columns=col_order, fill_value=0)
    
    print(f"Sorted heatmap dimensions: {pivot_df.shape}")
    print(f"Part1 classes order: {part1_class_order}")
    print(f"Ketones: {len(part2_ketones)}, Aldehydes: {len(part2_aldehydes)}")
    
    # save sorted pivot
    pivot_df.to_csv('fig/large_heatmap/data/pivot_data_sorted.tsv', sep='\t', index=False)
    
    return pivot_df, part1_class_order, class_display_names


def create_heatmap(pivot_df, part1_class_order, class_display_names, figsize=(8, 9.5)):
    """
    Create and save the heatmap
    """
    print('Creating heatmap...')
    
    # for plotting only
    heatmap_values = np.where(pivot_df.values > 0, np.maximum(pivot_df.values, 1), 0)

    # Apply log transformation
    log_data = np.log10(heatmap_values + 1)
    
    # Set up the figure
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 5
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Create the colormap
    colors = [(1, 1, 1), (0.78, 0.84, 0.94), (0.94, 0.69, 0.65)]
    positions = [0, 0.20, 1.0]
    cmap_wbr = LinearSegmentedColormap.from_list('white_blue_red', list(zip(positions, colors)), N=40)
    
    # Create heatmap
    im = ax.imshow(log_data, cmap=cmap_wbr, aspect='auto', interpolation='nearest')
    
    # Define compounds to mark with stars (lowercase for matching)
    starred_compounds = [
        'glycine_3-oxo-lithocholic acid',
        'alanine_3-oxo-lithocholic acid',
        'valine_3-oxo-lithocholic acid',
        'valine_3-oxo-deoxycholic acid',
        'leucine_3-oxo-deoxycholic acid',
        'isoleucine_3-oxo-deoxycholic acid',
        
        'phenylalanine_3-oxo-lithocholic acid',
        'isoleucine_3-oxo-lithocholic acid',
        'leucine_3-oxo-lithocholic acid',
        'phenylalanine_3-oxo-deoxycholic acid',
        'phenethylamine_3-oxo-deoxycholic acid',
        'tyramine_3-oxo-deoxycholic acid',
        'tryptophan_3-oxo-deoxycholic acid',
        'tyrosine_3-oxo-deoxycholic acid',
        'lysine_3-oxo-deoxycholic acid',
        'putrescine_3-oxo-deoxycholic acid',
        
        'n-acetylcadaverine_3-oxo-lithocholic acid',
        'n-acetylcadaverine_3-oxo-deoxycholic acid',
        '5-asa_ribose',
        '5-asa_glyceraldehyde',
        '5-asa_propanal',
        '5-asa_decanal',
        '5-asa_pyruvate',
        '5-asa_2-hydroxypropanal',
        '5-asa_deoxyribose',
        '5-asa_glycoaldehyde',
        'tryptophan_indole-3-carbaldehyde',
        'tryptophan_salicylaldehyde',
        '2-phenylglycine_octadecanal',
        '6-aminocaproic acid_hexanal',
        'ac-lys-oh_vanillin',
        'aspartic acid_cis-9-octadecenal',
        'aspartic acid_hexadecanal',
        'taurine_alpha-ketoisocaproic acid',
        'aspartic acid_decanal',
        'aspartic acid_tetradecanal',
        'adenosine_furfural',
        '5-hydroxytryptophan_dihydroxyacetone',
        'phosphoserine_galacturonic acid',
        'proline_citral',
        'glutamine_vanillin',
        'asparagine_salicylaldehyde',
        'gmp_salicylaldehyde',
        'tyrosine_vanillin',
        'sulfapyridine_glyceraldehyde',
        'sulfapyridine_2-oxoglutaric acid',
        'valine_cholestan-3-one',
        'leucine_cholestan-3-one',
        # 'phenylalanine_cholestan-3-one',
        'tyramine_cholestan-3-one',
        'n-acetylcadaverine_cholestan-3-one',
        'n-acetylputrescine_cholestan-3-one',
        'alanine_cholestan-3-one',
        'lysine_cholestan-3-one',
        'proline_cholestan-3-one',
        '5-asa_nonanal',
        '5-asa_octadecanal',
        'glycine_cholestan-3-one',
        'tyrosine_cholestan-3-one'
    ]
    
    # Add stars to specified compounds
    found_compounds = set()
    for row_idx, (part1, _) in enumerate(pivot_df.index):
        for col_idx, part2 in enumerate(pivot_df.columns):
            compound_key = f"{part1.lower()}_{part2.lower()}"
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
    else:
        print(f"\nAll {len(starred_compounds)} starred compounds found in heatmap!")
    
    highlighted_compounds = [
        '5-asa_ribose',
        '5-asa_glyceraldehyde',
        '5-asa_propanal',
        '5-asa_pyruvate',
        '5-asa_2-hydroxypropanal',
        '5-asa_deoxyribose',
        
        'n-acetylcadaverine_3-oxo-lithocholic acid',
        'phenylalanine_3-oxo-lithocholic acid',
        
        'n-acetylcadaverine_3-oxo-deoxycholic acid',
        'phenylalanine_3-oxo-deoxycholic acid',
        'phenethylamine_3-oxo-deoxycholic acid',
        'tyramine_3-oxo-deoxycholic acid',
        'tryptophan_3-oxo-deoxycholic acid',
        # 'tyrosine_3-oxo-deoxycholic acid',
        'lysine_3-oxo-deoxycholic acid',
        'putrescine_3-oxo-deoxycholic acid',
        
        'valine_cholestan-3-one',
        'leucine_cholestan-3-one',
        'tyramine_cholestan-3-one',
        'alanine_cholestan-3-one',
        'glycine_cholestan-3-one'
    ]
    
    # Add frames around highlighted compounds
    from matplotlib.patches import Rectangle
    highlighted_found = set()
    for row_idx, (part1, _) in enumerate(pivot_df.index):
        for col_idx, part2 in enumerate(pivot_df.columns):
            compound_key = f"{part1.lower()}_{part2.lower()}"
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
    else:
        print(f"\nAll {len(highlighted_compounds)} highlighted compounds found and framed!")
    
    # Set ticks and labels
    ax.set_xticks(range(len(pivot_df.columns)))
    ax.set_xticklabels(pivot_df.columns, rotation=90, ha='center', fontsize=6, color='0.2')
    
    # Move x-axis labels to top
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    
    ax.set_yticks(range(len(pivot_df.index)))
    ax.set_yticklabels([idx[0] for idx in pivot_df.index], fontsize=6, color='0.2')
    
    # Remove bottom ticks, keep top and left ticks
    ax.tick_params(left=True, bottom=False, top=True, pad=1.5, width=0.5, color='0.4', length=1)
    
    # Style the heatmap frame (spines)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_edgecolor('0.2')
    
    # Add small horizontal colorbar at bottom left
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    axins = inset_axes(ax, width="15%", height="2%", loc='lower left', 
                      bbox_to_anchor=(-0.13, -0.04, 0.6, 0.6), bbox_transform=ax.transAxes)
    cbar = plt.colorbar(im, cax=axins, orientation='horizontal')
    cbar.set_label('Log(Count + 1)', fontsize=5.5, color='0.2', labelpad=1)
    
    # Set colorbar ticks to show 0, middle, and max values
    vmin, vmax = log_data.min(), log_data.max()
    cbar.set_ticks([vmin, vmax/2, vmax])
    cbar.set_ticklabels([f'{vmin:.1f}', f'{vmax/2:.1f}', f'{vmax:.1f}'])
    
    cbar.ax.tick_params(labelsize=5, width=0.25, color='0.2', length=1.25, pad=1)
    cbar.outline.set_linewidth(0.25)
    cbar.outline.set_edgecolor('0.2')
    
    # Add class separators for Part1
    current_pos = 0
    for class_name in part1_class_order:
        n_molecules = sum(1 for idx in pivot_df.index if idx[1] == class_name)
        if n_molecules > 0:
            # Add horizontal line
            if current_pos > 0:
                ax.axhline(y=current_pos - 0.5, color='0.4', linewidth=.4, linestyle='-')
            # Add class label on the right
            # Use display name if available, otherwise use original name
            display_name = class_display_names.get(class_name, class_name)
            # Add line breaks for better readability
            label_text = display_name.replace(' (', '\n(')
            if display_name == 'Nucleic acid derivatives':
                label_text = 'Nucleic acid\nderivatives'
            elif display_name == 'Biogenic amines':
                label_text = 'Biogenic\namines'
            elif display_name == 'Neurotransmitters':
                label_text = 'Neuro-\ntransmitters'
            elif display_name == 'Polyamines':
                label_text = 'Poly-\namines'
            ax.text(len(pivot_df.columns) + 0.72, current_pos + n_molecules/2 - 0.5, 
                   label_text, fontsize=6.2, va='center', rotation=270, ha='center', weight='semibold', color='0.2')
            current_pos += n_molecules
    
    # Add separator between ketones and aldehydes
    n_ketones = sum(1 for col in pivot_df.columns if any(heatmap_df[heatmap_df['Part2'] == col]['Part2class'] == 'ketone'))
    if n_ketones > 0 and n_ketones < len(pivot_df.columns):
        ax.axvline(x=n_ketones - 0.5, color='0.4', linewidth=.4, linestyle='-')
        # Add labels at bottom
        ax.text(n_ketones/2 - 0.5, len(pivot_df.index) + 0.9, 'Ketones', fontsize=6.5, ha='center', weight='semibold', color='0.2')
        ax.text(n_ketones + (len(pivot_df.columns) - n_ketones)/2 - 0.5, len(pivot_df.index) + 0.9, 
               'Aldehydes', fontsize=6.5, ha='center', weight='semibold', color='0.2')
    
    # Set labels
    ax.set_xlabel('', fontsize=7, labelpad=0)
    ax.set_ylabel('', fontsize=7, labelpad=0)
    
    # Add legend for star markers and highlighted compounds
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    
    legend_elements = [
        Line2D([0], [0], marker='*', color='w', markerfacecolor='black', 
               markersize=9, label='Further validated by retention time'),
        Patch(facecolor='none', edgecolor='#FF6B00', linewidth=0.8, 
              label='Further validated by ion mobility')
    ]
    
    legend = ax.legend(handles=legend_elements, loc='lower left', 
                      bbox_to_anchor=(-0., -0.065), frameon=False, 
                      fontsize=5.5, handlelength=1.2, handleheight=1.2, 
                      handletextpad=0.5, borderpad=0.2)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save figure
    plt.savefig('fig/large_heatmap/overall_MASST_heatmap.svg', format='svg', bbox_inches='tight', transparent=True)
    plt.savefig('fig/large_heatmap/overall_MASST_heatmap.png', format='png', dpi=600, bbox_inches='tight', transparent=True)
    
    print('Heatmap saved to fig/large_heatmap/')
    
    return fig


if __name__ == "__main__":
    preprocess_data()
    
    pivot_df, heatmap_df = prepare_heatmap_data()
    
    # Sort molecules
    pivot_df, part1_class_order, class_display_names = sort_molecules(pivot_df, heatmap_df)
    
    # Create heatmap
    fig = create_heatmap(pivot_df, part1_class_order, class_display_names)
    
    print('\nAnalysis complete!')
    