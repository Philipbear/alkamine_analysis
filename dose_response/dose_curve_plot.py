import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


DATA_PATH = 'dose_response/data/N-acetylcadaverine_LCA_dose response curve.xlsx'
OUT_DIR = 'dose_response/plots'

# Sheet column layout: each compound spans 5 columns
#   [Concentration (uM), log10(M), rep1, rep2, rep3]
COMPOUNDS = [
    {
        'name': 'N-acetylcadaverine_LCA (mix)',
        'short': 'mix',
        'conc_col': 'Concentration (uM)',
        'rep_cols': ['Unnamed: 2', 'Unnamed: 3', 'Unnamed: 4'],
        'color': '#e54616',
    },
    {
        'name': 'N-acetylcadaverine_LCA isomer 1',
        'short': 'isomer1',
        'conc_col': 'Concentration (uM).1',
        'rep_cols': ['Unnamed: 7', 'Unnamed: 8', 'Unnamed: 9'],
        'color': '#4e639e',
    },
    {
        'name': 'N-acetylcadaverine_LCA isomer 2',
        'short': 'isomer2',
        'conc_col': 'Concentration (uM).2',
        'rep_cols': ['Unnamed: 12', 'Unnamed: 13', 'Unnamed: 14'],
        'color': '#dba053',
    },
]


def four_pl(x, bottom, top, log_ec50, hill):
    """4-parameter logistic (Hill) model on log10(concentration in uM)."""
    return bottom + (top - bottom) / (1.0 + 10.0 ** ((log_ec50 - x) * hill))


def fit_curve(log_x, y):
    bottom_0 = float(np.min(y))
    top_0 = float(np.max(y))
    log_ec50_0 = float(np.median(log_x))
    p0 = [bottom_0, top_0, log_ec50_0, 1.0]
    bounds = (
        [-np.inf, -np.inf, float(np.min(log_x)) - 2, 0.1],
        [np.inf, np.inf, float(np.max(log_x)) + 2, 10.0],
    )
    try:
        popt, _ = curve_fit(four_pl, log_x, y, p0=p0, bounds=bounds, maxfev=20000)
        return popt
    except Exception as exc:
        print(f'  curve_fit failed: {exc}')
        return None


def load_compound(df, comp):
    conc = df[comp['conc_col']].to_numpy(dtype=float)
    reps = df[comp['rep_cols']].to_numpy(dtype=float)
    long_x, long_y = [], []
    for i in range(conc.shape[0]):
        for j in range(reps.shape[1]):
            if np.isfinite(conc[i]) and np.isfinite(reps[i, j]):
                long_x.append(conc[i])
                long_y.append(reps[i, j])
    long_x = np.array(long_x)
    long_y = np.array(long_y)
    log_x = np.log10(long_x)
    return long_x, long_y, log_x


def style_axes(ax):
    ax.tick_params(axis='both', pad=1.2, width=0.5, color='0.4', length=1.5,
                   labelsize=5.5)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_edgecolor('0.2')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def plot_single(comp, output_path):
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 6
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['mathtext.bf'] = 'Arial:bold'
    plt.rcParams['mathtext.default'] = 'regular'

    df = pd.read_excel(DATA_PATH)
    x_uM, y, log_x = load_compound(df, comp)

    # Per-concentration mean / std for error bars.
    unique_conc = np.unique(x_uM)
    means = np.array([y[x_uM == c].mean() for c in unique_conc])
    stds = np.array([y[x_uM == c].std(ddof=1) for c in unique_conc])
    log_unique = np.log10(unique_conc)

    popt = fit_curve(log_x, y)

    fig, ax = plt.subplots(figsize=(2.2, 1.7))

    color = comp['color']

    # Replicate scatter.
    ax.scatter(log_x, y, s=10, color=color, alpha=0.55,
               edgecolors='white', linewidths=0.3, zorder=3)

    # Mean + std error bars (mean shown as a short horizontal line to avoid
    # being confused with the replicate scatter dots).
    ax.errorbar(log_unique, means, yerr=stds, fmt='_',
                color=color, ecolor=color, elinewidth=0.5, capsize=1.5,
                capthick=0.5, markersize=6, markeredgewidth=0.9,
                zorder=4, label='Mean ± SD')

    # Fitted curve — clipped to the measured concentration range (no extrapolation).
    if popt is not None:
        xx = np.linspace(log_x.min(), log_x.max(), 400)
        yy = four_pl(xx, *popt)
        ax.plot(xx, yy, color=color, linewidth=1.0, zorder=2)
        ec50_uM = 10.0 ** popt[2]
        print(f"  {comp['short']}: EC50 = {ec50_uM:.3f} uM, "
              f"Hill = {popt[3]:.3f}, bottom = {popt[0]:.3f}, top = {popt[1]:.3f}")

    # Log-scale ticks on x-axis (label as concentration in uM).
    tick_logs = np.arange(int(np.floor(log_x.min())), int(np.ceil(log_x.max())) + 1)
    ax.set_xticks(tick_logs)
    ax.set_xticklabels([f'$10^{{{int(t)}}}$' for t in tick_logs])

    ax.set_xlabel('Concentration (µM)', labelpad=1.5)
    ax.set_ylabel('Relative fluorescence units', labelpad=1.5)
    ax.set_title(comp['name'], fontsize=6, pad=2)

    style_axes(ax)
    plt.tight_layout()
    plt.savefig(output_path, format='svg', bbox_inches='tight', transparent=True)
    plt.close(fig)
    print(f'  saved: {output_path}')


def plot_combined(output_path):
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 6
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['mathtext.bf'] = 'Arial:bold'
    plt.rcParams['mathtext.default'] = 'regular'

    df = pd.read_excel(DATA_PATH)
    fig, ax = plt.subplots(figsize=(3.2, 2.0))

    all_log_x = []
    for comp in COMPOUNDS:
        x_uM, y, log_x = load_compound(df, comp)
        all_log_x.append(log_x)
        unique_conc = np.unique(x_uM)
        means = np.array([y[x_uM == c].mean() for c in unique_conc])
        stds = np.array([y[x_uM == c].std(ddof=1) for c in unique_conc])
        log_unique = np.log10(unique_conc)

        color = comp['color']
        popt = fit_curve(log_x, y)

        ax.errorbar(log_unique, means, yerr=stds, fmt='o',
                    color=color, ecolor=color, elinewidth=0.5, capsize=1.5,
                    capthick=0.5, markersize=3.2, markeredgecolor='white',
                    markeredgewidth=0.3, zorder=4, label=comp['name'])

        if popt is not None:
            xx = np.linspace(log_x.min(), log_x.max(), 400)
            yy = four_pl(xx, *popt)
            ax.plot(xx, yy, color=color, linewidth=1.0, zorder=2)

    all_log_x = np.concatenate(all_log_x)
    tick_logs = np.arange(int(np.floor(all_log_x.min())),
                          int(np.ceil(all_log_x.max())) + 1)
    ax.set_xticks(tick_logs)
    ax.set_xticklabels([f'$10^{{{int(t)}}}$' for t in tick_logs])

    ax.set_xlabel('Concentration (µM)', labelpad=1.5)
    ax.set_ylabel('Relative fluorescence units', labelpad=1.5)

    leg = ax.legend(frameon=False, fontsize=5.5, loc='upper left',
                    handlelength=1.2, handletextpad=0.5, borderaxespad=0.3)
    for txt in leg.get_texts():
        txt.set_color('0.2')

    style_axes(ax)
    plt.tight_layout()
    plt.savefig(output_path, format='svg', bbox_inches='tight', transparent=True)
    plt.close(fig)
    print(f'  saved: {output_path}')


if __name__ == '__main__':
    os.makedirs(OUT_DIR, exist_ok=True)

    print('Individual dose-response curves:')
    for comp in COMPOUNDS:
        out = os.path.join(OUT_DIR, f"dose_response_{comp['short']}.svg")
        plot_single(comp, out)

    print('Combined dose-response curve:')
    plot_combined(os.path.join(OUT_DIR, 'dose_response_combined.svg'))
