import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
import os

matplotlib.rcParams['font.family'] = 'Arial'

FS_TICK     = 4.5
FS_LABEL    = 5.5
FS_STAT     = 4.25
SUMMARY_CLR = '0.45'
REP_COLS    = ['Rep1', 'Rep2', 'Rep3']
OUT_DIR     = 'immune_data/plots'

os.makedirs(OUT_DIR, exist_ok=True)


def _draw_dots_and_mean(ax, x_pos, vals, color, rng):
    mean = vals.mean()
    sem  = stats.sem(vals)
    jitter = rng.uniform(-0.12, 0.12, len(vals))
    ax.plot([x_pos - 0.17, x_pos + 0.17], [mean, mean],
            color=SUMMARY_CLR, linewidth=1.0, zorder=2, solid_capstyle='round')
    ax.errorbar(x_pos, mean, yerr=sem,
                fmt='none', color=SUMMARY_CLR, elinewidth=0.8, capsize=0, zorder=2)
    ax.scatter(x_pos + jitter, vals,
               color=color, s=8, zorder=4, linewidths=0.2, edgecolors='white')


def _pval_bracket(ax, xi, xj, vals_i, vals_j, y_top, step):
    _, p = stats.ttest_ind(vals_i, vals_j, equal_var=False, nan_policy='omit')
    bh   = step * 0.12
    tpad = step * 0.10
    ax.plot([xi, xi, xj, xj], [y_top - bh, y_top, y_top, y_top - bh],
            color='black', linewidth=0.5, clip_on=False)
    ax.text((xi + xj) / 2, y_top + tpad, f'{p:.2e}',
            ha='center', va='bottom', fontsize=FS_STAT, color='black')


def _style_ax(ax, x_pos, tick_labels):
    ax.set_xticks(x_pos)
    ax.set_xticklabels(tick_labels, fontsize=FS_TICK, rotation=45, ha='right')
    ax.set_ylabel('Relative fluorescence units', fontsize=FS_LABEL)
    ax.tick_params(axis='y', labelsize=FS_TICK, length=1.0, width=0.5, pad=1)
    ax.tick_params(axis='x', length=0, pad=1.0)
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(0.5)
        ax.spines[sp].set_color('0.3')


def plot_nonspecific():
    df = pd.read_csv('immune_data/data/Test_immune_cells_nonspecific.tsv', sep='\t')
    df = df.set_index('Metabolite')
    metabolites = df.index.tolist()

    palette = [
        '#C0392B', '#1A4F82', '#1A7A6E', '#8E44AD', '#E07B39',
        '#2ECC71', '#E74C3C', '#3498DB', '#F39C12', '#16A085',
    ]
    colors = ['#000000'] + [palette[i % len(palette)] for i in range(len(metabolites) - 1)]

    fig, ax = plt.subplots(figsize=(2.7, 1.3))
    x_pos = np.arange(len(metabolites))
    rng = np.random.default_rng(42)

    for i, (met, color) in enumerate(zip(metabolites, colors)):
        vals = df.loc[met, REP_COLS].values.astype(float)
        _draw_dots_and_mean(ax, x_pos[i], vals, color, rng)

    # P-value bracket: DMSO ctrl vs N-acetylcadaverine_LCA
    i_ctrl = metabolites.index('DMSO ctrl')
    i_nac  = metabolites.index('N-acetylcadaverine_LCA')
    vals_ctrl = df.loc['DMSO ctrl', REP_COLS].values.astype(float)
    vals_nac  = df.loc['N-acetylcadaverine_LCA', REP_COLS].values.astype(float)

    y_top_data = df[REP_COLS].values.max()
    span  = max(y_top_data - df[REP_COLS].values.min(), 1e-6)
    step  = max(span * 0.15, abs(y_top_data) * 0.18)
    y_br  = max(vals_ctrl.max(), vals_nac.max()) + step
    _pval_bracket(ax, x_pos[i_ctrl], x_pos[i_nac], vals_ctrl, vals_nac, y_br, step)

    ax.set_xlim(-0.5, len(metabolites) - 0.5)
    ax.set_ylim(bottom=0, top=y_br + step * 0.4)
    _style_ax(ax, x_pos, [m.replace('_', '_\n') for m in metabolites])

    plt.tight_layout(pad=0.3)
    out_path = f'{OUT_DIR}/immune_cells_nonspecific.svg'
    plt.savefig(out_path, format='svg', bbox_inches='tight', transparent=True)
    plt.close(fig)
    print(f'Saved: {out_path}')


def plot_nacetylcadaverine_isomers():
    df = pd.read_csv(
        'immune_data/data/Test_immune_cells_nonspecific_Isomer_LCA_Nacetylcadaverine_DMSO.csv'
    )
    df = df.set_index('Metabolite')
    metabolites = df.index.tolist()
    colors = ['#000000', '#8E44AD', '#8E44AD']

    fig, ax = plt.subplots(figsize=(1, 0.65))
    x_pos = np.arange(len(metabolites))
    rng = np.random.default_rng(42)

    group_vals = []
    for i, (met, color) in enumerate(zip(metabolites, colors)):
        vals = df.loc[met, REP_COLS].values.astype(float)
        group_vals.append(vals)
        _draw_dots_and_mean(ax, x_pos[i], vals, color, rng)

    y_top_data = df[REP_COLS].values.max()
    span = max(y_top_data - df[REP_COLS].values.min(), 1e-6)
    step = max(span * 0.35, abs(y_top_data) * 0.25)

    y_base = max(v.max() for v in group_vals)
    for k, j in enumerate([1, 2]):
        y_br = y_base + step * (k + 1)
        _pval_bracket(ax, x_pos[0], x_pos[j], group_vals[0], group_vals[j], y_br, step)

    y_ceil = y_base + step * 3
    label_map = {
        'DMSO contr':              'DMSO',
        '(S) N-acetylcadaverine_LCA': '(S)-isomer',
        '(R) N-acetylcadaverine_LCA': '(R)-isomer',
    }
    tick_labels = [label_map.get(m, m) for m in metabolites]

    ax.set_xlim(-0.5, len(metabolites) - 0.5)
    ax.set_ylim(bottom=0, top=y_ceil)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(tick_labels, fontsize=FS_TICK, rotation=15, ha='center')
    ax.tick_params(axis='y', labelsize=FS_TICK, length=1.0, width=0.5, pad=1)
    ax.tick_params(axis='x', length=1.0, width=0.5, pad=1.0)
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_linewidth(0.5)
        ax.spines[sp].set_color('0.3')

    plt.tight_layout(pad=0.3)
    out_path = f'{OUT_DIR}/nacetylcadaverine_isomers.svg'
    plt.savefig(out_path, format='svg', bbox_inches='tight', transparent=True)
    plt.close(fig)
    print(f'Saved: {out_path}')


plot_nonspecific()
plot_nacetylcadaverine_isomers()
