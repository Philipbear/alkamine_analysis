import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators
import numpy as np
import os


def load_data():
    """
    Load overall distribution data
    """
    print('Loading overall distribution data...')
    # df = pd.read_csv('fig/overall_dist/overall_dist.tsv', sep='\t')
    
    df = pd.read_csv('fig/overall_dist/ALL_MASST_category.csv')

    df['Category_group'] = df['Category_group'].str.strip()

    # Deduplicate so each (USI, NAME, category) combination is unique
    df = df.drop_duplicates(subset=['USI', 'NAME', 'Category_group']).reset_index(drop=True)

    print(f"Loaded {len(df)} rows (one per unique USI/NAME/category)")
    print(f"Columns: {list(df.columns)}")
    
    return df


def prepare_upset_data(df):
    """
    Convert category data to format suitable for UpSet plot
    
    Args:
        df: DataFrame with Compound and Category_group columns
    
    Returns:
        pandas Series suitable for UpSet plotting
    """
    print('Preparing data for UpSet plot...')
    
    # Categories to include (excluding 'other' and 'personalCareProduct')
    categories_to_include = ['food', 'human', 'microbe', 'plant', 'rodent']

    # Filter to relevant categories (each row is already a single category)
    df_filtered = df[df['Category_group'].isin(categories_to_include)].copy()

    # Pivot: one row per unique compound name, one column per category
    binary_df = df_filtered.assign(_val=True).pivot_table(
        index='NAME',
        columns='Category_group',
        values='_val',
        aggfunc='any',
        fill_value=False
    )

    # Ensure all expected categories are present as columns
    for cat in categories_to_include:
        if cat not in binary_df.columns:
            binary_df[cat] = False
    binary_df = binary_df[categories_to_include]

    # Rename columns to display names
    binary_df = binary_df.rename(columns={
        'food': 'Food',
        'human': 'Humans',
        'microbe': 'Microbes',
        'plant': 'Plants',
        'rodent': 'Rodents'
    })

    # Remove rows where all values are False (no category match)
    binary_df = binary_df[binary_df.any(axis=1)]
    
    # Convert to format expected by upsetplot
    upset_data = from_indicators(binary_df)
    
    print(f"Found {len(binary_df)} compounds with at least one category match")
    print(f"Categories: {list(binary_df.columns)}")
    
    return upset_data, binary_df


def create_upset_plot(upset_data, binary_df, figsize=(4.8, 1.8), min_subset_size=1):
    """
    Create and save UpSet plot
    """
    print('Creating UpSet plot...')
    
    # Set Arial font and figure size
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['figure.figsize'] = figsize
    plt.rcParams['font.size'] = 7
    
    # Create UpSet plot with visible elements and custom styling
    upset = UpSet(
        upset_data,
        subset_size='count',
        intersection_plot_elements=5,  # Show intersection bars
        totals_plot_elements=3,        # Show totals bars
        min_subset_size=min_subset_size,
        show_counts=True,
        sort_by='cardinality',         # Sort by number of sets in intersection
        sort_categories_by='-cardinality',  # Sort categories by total count
        element_size=10,               # Size of dots in matrix
        facecolor='0.4',               # Dot color
        connecting_line_width=1,     # Width of lines connecting dots
    )
    
    # Generate the plot - this returns a dict of axes, not a figure
    axes_dict = upset.plot()
    
    # Get the figure from one of the axes
    fig = list(axes_dict.values())[0].figure
    
    # Now set the figure size
    fig.set_size_inches(figsize)
    
    # Adjust tick parameters and spine line width for bar plots
    for ax_name, ax in axes_dict.items():
        if ax_name == 'intersections':            
            # Intersection bar plot
            ax.tick_params(
                axis='both',
                which='major',
                length=1.5,      # Tick length
                width=0.8,     # Tick width
                labelsize=5.5,   # Label size
                color='0.3',   # Tick color
                pad=1.5          # Distance between tick and label
            )
            # Remove minor ticks if present
            ax.tick_params(axis='both', which='minor', length=0, width=0)
            
            # Add grid
            ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, axis='y')
            
            # Adjust spine line width for intersection bars
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)
                spine.set_color('0.4')
            
            # Color bars
            for patch in ax.patches:
                patch.set_facecolor('0.4')      # Bar fill color
                patch.set_edgecolor('0.4')      # Bar edge color
            
            # Shift x-axis left to avoid overlap with first bar count label
            xlim = ax.get_xlim()
            ax.set_xlim(xlim[0] - 0.3, xlim[1])
            
            ax.set_ylabel('Intersection\nsize', fontsize=7, labelpad=3.5)
            
        elif ax_name == 'totals':
            # Totals bar plot
            ax.tick_params(
                axis='both',
                which='major',
                length=1.5,
                width=0.8,
                labelsize=5.5,
                color='0.3',
                pad=1.5
            )
            ax.tick_params(axis='both', which='minor', length=0, width=0)
            
            # add x labels
            ax.set_xlabel('Set size', fontsize=6, labelpad=2)
            
            # Add grid
            ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5, axis='x')
            
            # Adjust spine line width for totals bars
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)
                spine.set_color('0.4')
            
            # Color bars with different colors for each category
            # Define color mapping for categories
            category_colors = {
                'Humans': '#E64B35',     # Red/coral (blood/life)
                'Rodents': '#8B7355',    # Brown (rodent color)
                'Plants': '#4DAF4A',     # Green (plant color)
                'Food': '#FF9933',       # Orange (food color)
                'Microbes': '#9467BD'    # Purple (microbes)
            }
            
            # Get category labels from y-axis
            category_labels = [label.get_text() for label in ax.get_yticklabels()]
            
            # Color each bar according to its category
            for i, patch in enumerate(ax.patches):
                if i < len(category_labels):
                    category = category_labels[i]
                    color = category_colors.get(category, '#476f95')  # Default color if not found
                    patch.set_facecolor(color)      # Bar fill color
                    patch.set_edgecolor(color)      # Bar edge color
                    patch.set_alpha(0.75)           # Add transparency
    
    # Format numbers with commas for intersection bars
    if 'intersections' in axes_dict:
        ax = axes_dict['intersections']
        for text in ax.texts:
            current_text = text.get_text()
            if current_text.isdigit():
                formatted_text = f"{int(current_text):,}"
                text.set_text(formatted_text)
                text.set_fontsize(5.5)
    
    # Format numbers with commas for totals bars
    if 'totals' in axes_dict:
        ax = axes_dict['totals']
        for text in ax.texts:
            current_text = text.get_text()
            if current_text.isdigit():
                formatted_text = f"{int(current_text):,}"
                text.set_text(formatted_text)
                text.set_fontsize(5.5)
                # Shift text position to the left
                x, y = text.get_position()
                text.set_position((x + 20, y))
    
    # Save the plot
    plt.tight_layout()
    plt.savefig('fig/overall_dist/overall_dist_upset.svg', 
                format='svg', bbox_inches='tight', transparent=True)
    
    print('UpSet plot saved to fig/overall_dist/')
    
    return fig


def create_summary_statistics(binary_df):
    """
    Create summary statistics for the dataset
    """
    print('\n=== Summary Statistics ===')
    
    total_compounds = len(binary_df)
    total_categories = len(binary_df.columns)
    
    print(f'Total compounds: {total_compounds:,}')
    print(f'Total categories: {total_categories}')
    print(f'Categories: {", ".join(binary_df.columns)}')
    
    # Compounds found in multiple categories
    compounds_per_category = binary_df.sum(axis=1)
    multi_category_compounds = (compounds_per_category > 1).sum()
    single_category_compounds = (compounds_per_category == 1).sum()
    
    print(f'\nCompounds found in multiple categories: {multi_category_compounds:,} ({multi_category_compounds/total_compounds*100:.1f}%)')
    print(f'Compounds found in single category: {single_category_compounds:,} ({single_category_compounds/total_compounds*100:.1f}%)')
    
    # Most ubiquitous compounds
    max_categories = compounds_per_category.max()
    most_ubiquitous = compounds_per_category[compounds_per_category == max_categories].index
    print(f'\nMost ubiquitous compound(s) found in {max_categories} categories: {len(most_ubiquitous)} compound(s)')
    
    # Category counts
    print('\nCompounds per category:')
    for category in binary_df.columns:
        count = binary_df[category].sum()
        print(f'  {category}: {count:,}')


def main():
    """
    Main function to run the complete analysis
    """
    # Load data
    df = load_data()
    
    # Prepare data for UpSet plot
    upset_data, binary_df = prepare_upset_data(df)
    
    # Create summary statistics
    create_summary_statistics(binary_df)
    
    # Create UpSet plot
    fig = create_upset_plot(upset_data, binary_df)
    
    print('\nAnalysis complete! Plot saved to fig/overall_dist/')
    
    return df, binary_df, upset_data


if __name__ == '__main__':
    # Run the complete analysis
    df, binary_df, upset_data = main()
