import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import vonmises, pearsonr

# --- CONFIG ---
file_path = "phases.xlsx"
increase_sheet = "Increase"
decrease_sheet = "Decrease"
cols = ["Post_optimal"]            # column(s) to use
bins = 24                          # 15° bins
kappa = 4                          # von Mises concentration
cosine_ref_amplitude = 1.0         # scale of red reference cosine
figsize = (8, 8)
bar_alpha = 0.6
inc_bar_color = "#FFD580"          # light orange
dec_bar_color = "#ADD8E6"          # light blue
inc_line_color = "#CC7700"         # dark orange
dec_line_color = "#1E3F66"         # dark blue
ref_line_color = "red"             # red dotted reference
ring_step = 2                      # radial rings every 2 sessions
theta_zero = "E"
theta_direction_clockwise = True

def load_angles(xlsx_path, sheet_name, columns):
    df = pd.read_excel(xlsx_path, sheet_name=sheet_name)
    cols_map = {c.lower(): c for c in df.columns}
    arr_list = []
    for c in columns:
        key = c.lower()
        if key not in cols_map:
            raise ValueError(f"Column '{c}' not found in sheet '{sheet_name}'.")
        arr_list.append(df[cols_map[key]].to_numpy())
    if not arr_list:
        return np.array([])
    return np.concatenate(arr_list, axis=0)

def circ_hist(theta, nbins=24):
    edges = np.linspace(0, 2*np.pi, nbins+1)
    counts, _ = np.histogram(np.mod(theta, 2*np.pi), bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2.0
    width = edges[1] - edges[0]
    return centers, counts, width

def vonmises_kde(theta, grid, kappa=4):
    n = theta.size
    if n == 0:
        return np.zeros_like(grid)
    c = 1.0 / (2*np.pi*np.i0(kappa))
    dens = np.zeros_like(grid, dtype=float)
    for t in theta:
        dens += c * np.exp(kappa * np.cos(grid - t))
    dens /= n
    return dens

def scale_kde_to_counts(density, n_samples, bin_width_rad):
    # Convert density (area=1) to approx counts/bin, so it overlays with bars
    return density * n_samples * bin_width_rad

def main():
    # 1) load raw codes from both sheets/columns
    inc_df = pd.read_excel(file_path, sheet_name=increase_sheet)
    dec_df = pd.read_excel(file_path, sheet_name=decrease_sheet)

    inc_codes = inc_df[cols].to_numpy().flatten()
    dec_codes = dec_df[cols].to_numpy().flatten()

    # mapping (1..8) in radians -> degrees, wrapped to [0, 360)
    phase_mapping_degrees = {
        1: np.degrees(-2.7489) % 360,
        2: np.degrees(-1.9635) % 360,
        3: np.degrees(-1.1781) % 360,
        4: np.degrees(-0.3927) % 360,
        5: np.degrees(0.3927)  % 360,
        6: np.degrees(1.1781)  % 360,
        7: np.degrees(1.9635)  % 360,
        8: np.degrees(2.7489)  % 360,
    }

    # map codes -> degrees
    inc_deg = [phase_mapping_degrees[int(c)] for c in inc_codes if pd.notna(c) and int(c) in phase_mapping_degrees]
    dec_deg = [phase_mapping_degrees[int(c)] for c in dec_codes if pd.notna(c) and int(c) in phase_mapping_degrees]

    # shift all phases by +30° and wrap to [0, 360)
    inc_deg = [(d + 30) % 360 for d in inc_deg]
    dec_deg = [(d + 30) % 360 for d in dec_deg]

    # convert to radians (unchanged)
    inc_rad = np.radians(inc_deg)
    dec_rad = np.radians(dec_deg)

    # 2) histogram + KDE for each condition
    inc_centers, inc_counts, inc_width = circ_hist(inc_rad, nbins=bins)
    dec_centers, dec_counts, dec_width = circ_hist(dec_rad, nbins=bins)

    grid = np.linspace(0, 2*np.pi, 720, endpoint=False)
    inc_dens = vonmises_kde(inc_rad, grid, kappa=kappa)
    dec_dens = vonmises_kde(dec_rad, grid, kappa=kappa)

    inc_kde_counts = scale_kde_to_counts(inc_dens, len(inc_rad), inc_width)
    dec_kde_counts = scale_kde_to_counts(dec_dens, len(dec_rad), dec_width)

    # Shared radial calibration (sessions)
    max_count = max(int(inc_counts.max()) if inc_counts.size else 1,
                    int(dec_counts.max()) if dec_counts.size else 1)
    if ring_step > 0 and max_count % ring_step:
        max_count += (ring_step - (max_count % ring_step))

    # 3) single polar overlay
    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(1, 1, 1, projection="polar")
    ax.set_theta_zero_location(theta_zero)
    ax.set_theta_direction(-1 if theta_direction_clockwise else 1)
    ax.set_thetagrids(range(0, 360, 45))

    # Bars
    if inc_counts.size:
        ax.bar(inc_centers, inc_counts, width=inc_width, bottom=0,
               color=inc_bar_color, edgecolor="none", alpha=bar_alpha, align='center', label="Increase")
    if dec_counts.size:
        ax.bar(dec_centers, dec_counts, width=dec_width, bottom=0,
               color=dec_bar_color, edgecolor="none", alpha=bar_alpha, align='center', label="Decrease")

    # KDE lines
    if inc_kde_counts.size:
        ax.plot(grid, inc_kde_counts, color=inc_line_color, lw=2)
    if dec_kde_counts.size:
        ax.plot(grid, dec_kde_counts, color=dec_line_color, lw=2)

    # μ reference (cosine), dotted red
    mu = (np.cos(grid) + 1.0)/2.0 * max_count * cosine_ref_amplitude
    ax.plot(grid, mu, ls=":", color=ref_line_color, lw=2, label="μ reference (cosine)")

    # Radial ticks/range (sessions)
    rticks = np.arange(ring_step, max_count+1, ring_step)
    ax.set_rticks(rticks)
    ax.set_yticklabels([str(int(v)) for v in rticks])
    ax.set_rlim(0, max_count)

    ax.legend(loc="upper right", frameon=True)
    fig.tight_layout()
    plt.savefig("circular_kde_post_optimal_SHIFT30.png", dpi=300, bbox_inches="tight")
    print("Saved: circular_kde_post_optimal_SHIFT30.png")

# (optional) circular correlation helpers
def js_circ_corr(alpha, beta):
    alpha = np.asarray(alpha); beta = np.asarray(beta)
    alpha = np.mod(alpha, 2*np.pi); beta = np.mod(beta, 2*np.pi)
    alpha = alpha - np.mean(alpha)
    beta = beta - np.mean(beta)
    sin_a = np.sin(alpha); sin_b = np.sin(beta)
    num = np.sum(sin_a * sin_b)
    den = np.sqrt(np.sum(sin_a**2) * np.sum(sin_b**2))
    return float(num/den) if den != 0 else np.nan

def js_circ_corr_p(alpha, beta, n_perm=2000, seed=42):
    rng = np.random.default_rng(seed)
    r_obs = js_circ_corr(alpha, beta)
    if not np.isfinite(r_obs):
        return np.nan
    count = 0
    for _ in range(n_perm):
        r_perm = js_circ_corr(alpha, rng.permutation(beta))
        if np.abs(r_perm) >= np.abs(r_obs):
            count += 1
    return (count + 1) / (n_perm + 1)

if __name__ == "__main__":
    main()
