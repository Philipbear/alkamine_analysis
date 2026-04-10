import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def load_data():
    df = pd.read_csv("culture_data/data/Peak Area extraction_bacteria culture.csv")
    
    # for this molecule
    df = df[df['Molecule'] == 'Phenylalanine_Dehydrolithocholicacid'].reset_index(drop=True)
    
    # remove qc
    df = df[~df['Replicate Name'].str.contains('poolqc')].reset_index(drop=True)
    df = df[~df['Replicate Name'].str.contains('ID')].reset_index(drop=True)
    df = df[~df['Replicate Name'].str.contains('Blank')].reset_index(drop=True)
    
    df = df[['Replicate Name', 'Area']]
    df.columns = ['file', 'peak_area']
    
    df['BA_added'] = df['file'].apply(lambda x: True if '_drug' in x else False)
    df['time'] = df['file'].apply(lambda x: 'T0' if 'T0' in x else ('T72' if 'T72' in x else 'Unknown'))
    
    def get_microbe(file_name):
        if '_drug' in file_name:
            a = file_name.split('_drug')[0]
        elif '_T0' in file_name:
            a = file_name.split('_T0')[0]
        else:
            a = file_name.split('_T72')[0]
        
        a = a.replace('_', ' ').replace('110723 9', '110723-9').replace('T0', '').replace('T72', '').strip()
        return a
    def get_short_name(a):
        # count how many times ' ' appears in the name
        count = a.count(' ')
        if count > 1:
            # Basaras XXXX YYY  => B. XXXX YYY
            return a.split(' ')[0][0] + '. ' + ' '.join(a.split(' ')[1:])
        else:
            return a
    
    df['microbe'] = df['file'].apply(get_microbe)
    df['short_name'] = df['microbe'].apply(get_short_name)

    # save
    df.to_csv("culture_data/data/refined_data.csv", index=False)
    return df


def build_heatmap_matrix_from_refined(refined_path="culture_data/data/refined_data.csv"):
    """
    Build a condition x microbe matrix from refined_data.csv.
    Each cell is the average intensity across replicates.
    """
    df = pd.read_csv(refined_path)

    condition_order = [
        '0 h (w/o BA)',
        '72 h (w/o BA)',
        '0 h (w/ BA)',
        '72 h (w/ BA)'
    ]

    df['condition'] = np.where(
        (df['time'] == 'T0') & (~df['BA_added']),
        '0 h (w/o BA)',
        np.where(
            (df['time'] == 'T72') & (~df['BA_added']),
            '72 h (w/o BA)',
            np.where(
                (df['time'] == 'T0') & (df['BA_added']),
                '0 h (w/ BA)',
                np.where((df['time'] == 'T72') & (df['BA_added']), '72 h (w/ BA)', 'Unknown')
            )
        )
    )

    df = df[df['condition'] != 'Unknown'].copy()
    df['condition'] = pd.Categorical(df['condition'], categories=condition_order, ordered=True)

    # Average intensity over replicates for each condition and microbe.
    mean_df = (
        df.groupby(['condition', 'short_name'], observed=True)['peak_area']
        .mean()
        .reset_index(name='avg_peak_area')
    )

    matrix = (
        mean_df.pivot(index='condition', columns='short_name', values='avg_peak_area')
        .reindex(condition_order)
        .fillna(0)
    )

    matrix = matrix.reindex(sorted(matrix.columns), axis=1)

    # Keep selected media controls at the end of x-axis order.
    trailing_cols = ['cCOM', 'BHIS', 'BruPT-FE', 'RCM']
    trailing_present = [col for col in trailing_cols if col in matrix.columns]
    leading_cols = [col for col in matrix.columns if col not in trailing_present]
    matrix = matrix[leading_cols + trailing_present]

    matrix.to_csv("culture_data/data/heatmap_matrix_avg_intensity.csv")
    return matrix


def plot_heatmap(matrix):
    """
    Plot heatmap from the average-intensity matrix.
    """
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 5.5

    display_values = matrix.values

    fig, ax = plt.subplots(figsize=(2.4, 0.63))

    # Same palette as steroid heatmap.
    colors = [(1, 1, 1), (0.78, 0.84, 0.94), (0.94, 0.69, 0.65)]
    positions = [0, 0.75, 1.0]
    cmap_wbr = LinearSegmentedColormap.from_list('white_blue_red', list(zip(positions, colors)), N=30)

    im = ax.imshow(display_values, cmap=cmap_wbr, aspect='auto', interpolation='nearest')

    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation=50, ha='left')
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(matrix.index)

    # Match frame and tick style used in steroid heatmap.
    ax.tick_params(left=True, bottom=False, top=True, pad=1., width=0.5, color='0.4', length=1)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_edgecolor('0.2')

    # Vertical colorbar at the left-top corner.
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cax = inset_axes(
        ax,
        width="4%",
        height="55%",
        loc='upper left',
        bbox_to_anchor=(-0.17, 0.93, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=0
    )
    cbar = plt.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label('Avg. peak area', fontsize=5.5, labelpad=1.5)
    cbar.ax.yaxis.set_label_position('left')

    vmin, vmax = display_values.min(), display_values.max()
    vmid = (vmin + vmax) / 2
    cbar.set_ticks([vmin, vmid, vmax])
    cbar.set_ticklabels([f'{vmin:.1e}', f'{vmid:.1e}', f'{vmax:.1e}'])
    cbar.ax.tick_params(labelsize=5., width=0.25, color='0.2', length=1.2, pad=1)
    cbar.outline.set_linewidth(0.25)
    cbar.outline.set_edgecolor('0.2')

    ax.set_xlabel('')
    ax.set_ylabel('')

    fig.subplots_adjust(left=0.16, right=0.98, top=0.90, bottom=0.16)
    plt.savefig("culture_data/heatmap_culture.svg", format='svg', bbox_inches='tight', transparent=True)
    return fig


if __name__ == "__main__":
    # Keep your original preprocessing step.
    load_data()

    # Build and plot heatmap from refined_data.csv.
    heatmap_matrix = build_heatmap_matrix_from_refined("culture_data/data/refined_data.csv")
    plot_heatmap(heatmap_matrix)

    print("Heatmap generated: culture_data/heatmap_culture.svg")
