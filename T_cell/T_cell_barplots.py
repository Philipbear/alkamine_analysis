import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
from itertools import combinations
import os

matplotlib.rcParams['font.family'] = 'Arial'

# Font sizes
FS_TICK  = 4.5
FS_LABEL = 5.5
FS_STAT  = 4.25

# Load data
dataset_configs = [
    {
        'cell_type': 'Th17',
        'marker_name': 'IL-10+',
        'path': 'T_cell/data/CD4+IL10+_Th17_Name_Plot for all Th17_IL10A_second from right.csv',
    },
    {
        'cell_type': 'Treg',
        'marker_name': 'IL-10+',
        'path': 'T_cell/data/CD4+IL10+_Treg_Name_Plot for all Treg_IL10A_fouth from right.csv',
    },
    {
        'cell_type': 'Th17',
        'marker_name': 'IL-17A+',
        'path': 'T_cell/data/CD4+IL17A+_Th17_Name_Plot for all Th17_IL17A_first from right.csv',
    },
    {
        'cell_type': 'Treg',
        'marker_name': 'IL-17A+',
        'path': 'T_cell/data/CD4+IL17A+_Treg_Name_Plot for all Treg_IL17A_third from right.csv',
    },
]

# 4 comparison pairs × 2 markers = 8 plots
comparison_pairs = [
    ('5-ASA_pyruvate',              '5-ASA'),
    ('5-ASA_ribose',                '5-ASA'),
    ('Sulfapyridine_glyceraldehyde', 'Sulfapyridine'),
    ('Sulfapyridine_AKG',           'Sulfapyridine'),
]

label_map = {
    'Contr':                        'Ctrl',
    '5-ASA_pyruvate':               '5-ASA_\npyruvate',
    '5-ASA_ribose':                 '5-ASA_\nribose',
    '5-ASA':                        '5-ASA',
    'Sulfapyridine_glyceraldehyde': 'Sulfapyridine_\nglyceraldehyde',
    'Sulfapyridine_AKG':            'Sulfapyridine_\nAKG',
    'Sulfapyridine':                'Sulfapyridine',
}

# Colors: [control, parent, conjugate]
colors_5asa  = ['#000000', '#C0392B', '#E07B39']
colors_sulfa = ['#000000', '#1A4F82', '#3A7EBF']
summary_color = '0.45'

# Output directory
out_dir = 'T_cell/plots'
os.makedirs(out_dir, exist_ok=True)


def fmt_pval(p):
    return f'{p:.2e}'


def add_stat_brackets(ax, x_positions, group_vals, y_start_frac=0.2):
    """Draw all pairwise Welch's t-test p-value brackets above the data."""
    all_vals = np.concatenate(group_vals)
    y_top = float(np.max(all_vals))
    y_bot = float(np.min(all_vals))
    span = max(y_top - y_bot, 1e-6)

    # Keep a minimum spacing so p-value labels do not overlap in low-range panels.
    step = max(span * y_start_frac, max(abs(y_top), 1.0) * 0.18)
    bh = step * 0.12 # bracket height
    text_pad = step * 0.1 # space between bracket and p-value text

    pairs = list(combinations(range(len(x_positions)), 2))
    for k, (i, j) in enumerate(pairs):
        _, p = stats.ttest_ind(group_vals[i], group_vals[j],
                               equal_var=False, nan_policy='omit')
        y_br = y_top + step * (k + 1)
        xi, xj = x_positions[i], x_positions[j]

        ax.plot([xi, xi, xj, xj],
                [y_br - bh, y_br, y_br, y_br - bh],
                color='black', linewidth=0.5, clip_on=False)
        ax.text((xi + xj) / 2, y_br + text_pad,
                fmt_pval(p),
                ha='center', va='bottom',
                fontsize=FS_STAT, color='black')

    return y_top + step * (len(pairs) + 1.4)


def make_plot(df, groups, colors, marker_name, cell_type, out_path):
    fig, ax = plt.subplots(figsize=(0.95, 1.1))

    x_pos = np.arange(len(groups))
    rng   = np.random.default_rng(42)

    group_vals = []
    for i, (col, color) in enumerate(zip(groups, colors)):
        vals = df[col].values.astype(float)
        group_vals.append(vals)
        mean = vals.mean()
        sem  = stats.sem(vals)

        jitter = rng.uniform(-0.12, 0.12, len(vals))
        ax.plot([x_pos[i] - 0.17, x_pos[i] + 0.17], [mean, mean],
            color=summary_color, linewidth=1., zorder=2, solid_capstyle='round')

        ax.errorbar(x_pos[i], mean, yerr=sem,
                fmt='none', color=summary_color,
                elinewidth=0.8, capsize=0, zorder=2)

        ax.scatter(x_pos[i] + jitter, vals,
               color=color, s=8, zorder=4,
             linewidths=0.2, edgecolors='white')

    y_ceil = add_stat_brackets(ax, x_pos, group_vals)

    ax.set_xlim(-0.5, len(groups) - 0.5)
    ax.set_ylim(bottom=0, top=y_ceil)

    ax.set_xticks(x_pos)
    ax.set_xticklabels([label_map.get(g, g) for g in groups],
                       fontsize=FS_TICK, rotation=30, ha='right')

    ax.set_ylabel(f'{marker_name}\namong CD4+ (%)', fontsize=FS_LABEL)
    ax.tick_params(axis='y', labelsize=FS_TICK, length=1., width=0.5, pad=1)
    ax.tick_params(axis='x', length=0, pad=1.)

    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(0.5)
        ax.spines[sp].set_color('0.3')

    plt.tight_layout(pad=0.3)
    plt.savefig(out_path, format='svg', bbox_inches='tight', transparent=True)
    plt.close(fig)
    print(f'Saved: {out_path}')


# ── Generate individual SVGs for each cell type/marker and comparison ───────
for cfg in dataset_configs:
    df = pd.read_csv(cfg['path'])
    marker_name = cfg['marker_name']
    cell_type = cfg['cell_type']

    for pair_idx, (conjugate, parent) in enumerate(comparison_pairs):
        colors = colors_5asa if pair_idx < 2 else colors_sulfa
        groups = ['Contr', parent, conjugate]

        safe_conj = conjugate.replace('-', '').replace('_', '_')
        safe_marker = marker_name.replace('+', 'pos').replace('-', '')
        safe_cell = cell_type.lower()
        fname = f'{out_dir}/{safe_cell}_{safe_marker}_{safe_conj}.svg'

        make_plot(df, groups, colors, marker_name, cell_type, fname)
