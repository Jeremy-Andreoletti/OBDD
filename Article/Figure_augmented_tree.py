"""
Generate Figure: The OBDD model with data augmentation.

Panels:
  A) Empirical (observed) OBD tree with coded fossils and uncoded occurrences
  B) One complete (data-augmented) OBDD tree
  C) Latent speciation rate λ(t) on the augmented tree
  D) Latent extinction rate μ(t) on the augmented tree

Usage:
  python Figure_augmented_tree.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
import matplotlib.gridspec as gridspec

np.random.seed(42)

# ── Colour palette (matching paper conventions) ─────────────────────
COL_OBS       = "black"          # observed lineages
COL_AUG       = "#AAAAAA"        # augmented (gray) lineages
COL_FOSSIL    = "#8B5E3C"        # coded fossils (brown, paper convention)
COL_OCC       = "#E8833A"        # uncoded occurrences (orange, paper convention)

# Muted colormap for rates: steel blue (low) → warm tan → muted brick (high)
RATE_CMAP = LinearSegmentedColormap.from_list(
    'muted_rates',
    ['#2B6A99', '#7AABC5', '#D4CDB0', '#C8916E', '#A0452A'],
    N=256)

# ── Time axis ───────────────────────────────────────────────────────
T_ROOT = 0.0
T_PRESENT = 10.0

# ── Define the empirical (observed) tree topology ───────────────────

class Branch:
    def __init__(self, y, t0, t1, extinct=False, parent_y=None,
                 fossil_times=None):
        self.y = y
        self.t0 = t0
        self.t1 = t1
        self.extinct = extinct
        self.parent_y = parent_y
        self.fossil_times = fossil_times or []


def build_observed_tree():
    """Build an empirical FBD-like tree with coded fossils."""
    branches = []
    # Root stem
    branches.append(Branch(y=5.0, t0=0.0, t1=1.5, parent_y=None))

    # First split at t=1.5: y=5 → y=7.5 (upper) and y=2.5 (lower)
    # Upper clade
    branches.append(Branch(y=7.5, t0=1.5, t1=3.5, parent_y=5.0))
    # Upper-upper split at t=3.5
    branches.append(Branch(y=8.5, t0=3.5, t1=10.0, parent_y=7.5,
                           fossil_times=[5.0]))  # sampled ancestor
    branches.append(Branch(y=6.5, t0=3.5, t1=7.0, parent_y=7.5,
                           extinct=True))  # tip fossil (extinct)

    # Lower clade: split at t=1.5
    branches.append(Branch(y=2.5, t0=1.5, t1=3.0, parent_y=5.0))
    # Lower split at t=3.0
    branches.append(Branch(y=4.0, t0=3.0, t1=5.5, parent_y=2.5))
    # Sub-split at t=5.5
    branches.append(Branch(y=4.5, t0=5.5, t1=10.0, parent_y=4.0,
                           fossil_times=[7.5]))  # sampled ancestor
    branches.append(Branch(y=3.5, t0=5.5, t1=10.0, parent_y=4.0))

    branches.append(Branch(y=1.5, t0=3.0, t1=6.0, parent_y=2.5,
                           extinct=True,
                           fossil_times=[4.5]))  # tip fossil

    return branches


def build_augmented_tree(obs_branches):
    """Add augmented (unobserved) lineages to the empirical tree."""
    aug_branches = []

    # Augmented lineage 1: early divergence from stem, goes extinct
    aug_branches.append(Branch(y=9.2, t0=0.8, t1=2.8, extinct=True, parent_y=5.0))

    # Augmented lineage 2: diverges from upper clade branch (y=7.5, t=1.5–3.5)
    aug_branches.append(Branch(y=5.7, t0=2.2, t1=4.5, extinct=True, parent_y=7.5))

    # Augmented lineage 3: diverges from extinct lower branch (y=1.5, t=3.0–6.0)
    aug_branches.append(Branch(y=0.8, t0=4.0, t1=7.5, extinct=True, parent_y=1.5))

    # Augmented lineage 4: late divergence from branch y=4.5 (t=5.5–10.0)
    aug_branches.append(Branch(y=5.2, t0=6.5, t1=9.0, extinct=True, parent_y=4.5))

    # Augmented continuation after fossil tip at (7.0, 6.5) → hidden extinction
    aug_branches.append(Branch(y=6.5, t0=7.0, t1=8.5, extinct=True, parent_y=6.5))

    # Augmented continuation after fossil tip at (6.0, 1.5) → unsampled extant tip
    aug_branches.append(Branch(y=1.5, t0=6.0, t1=10.0, extinct=False, parent_y=1.5))

    return aug_branches


def generate_occurrences():
    """Generate uncoded occurrence times (not placed on tree)."""
    return [1.2, 2.5, 3.8, 4.2, 5.0, 5.8, 6.3, 7.0, 7.8, 8.5, 9.2]


def draw_tree(ax, branches, color=COL_OBS, lw=1.8, draw_fossils=True,
              show_tips=True, alpha=1.0):
    """Draw a phylogenetic tree on the given axes."""
    for br in branches:
        # Horizontal branch
        ax.plot([br.t0, br.t1], [br.y, br.y], color=color, lw=lw,
                alpha=alpha, solid_capstyle='round', zorder=2)

        # Vertical connector to parent
        if br.parent_y is not None:
            ax.plot([br.t0, br.t0], [br.parent_y, br.y], color=color,
                    lw=lw, alpha=alpha, solid_capstyle='round', zorder=2)

        if draw_fossils:
            # Coded fossils (brown diamonds)
            for ft in br.fossil_times:
                ax.plot(ft, br.y, 'D', color=COL_FOSSIL, markersize=7,
                        markeredgecolor='white', markeredgewidth=0.5, zorder=4)

            # Extinct tip markers (brown diamond)
            if br.extinct:
                ax.plot(br.t1, br.y, 'D', color=COL_FOSSIL, markersize=7,
                        markeredgecolor='white', markeredgewidth=0.5, zorder=4)

        if show_tips and not br.extinct:
            ax.plot(br.t1, br.y, 'o', color=color, markersize=3,
                    zorder=3, alpha=alpha)


def draw_occurrences(ax, occ_times, y_base):
    """Draw uncoded occurrences as orange circles on a timeline below the tree."""
    ax.plot([T_ROOT, T_PRESENT], [y_base, y_base], color='gray',
            lw=0.5, zorder=1)
    for t in occ_times:
        ax.plot(t, y_base, 'o', color=COL_OCC, markersize=7,
                markeredgecolor='white', markeredgewidth=0.5, zorder=4)


def generate_rates_along_branch(t0, t1, rate0, drift, sigma, dt=0.02):
    """Simulate GBM rates along a branch."""
    times = np.arange(t0, t1, dt)
    if len(times) == 0:
        return np.array([t0, t1]), np.array([rate0, rate0])
    log_rate = np.log(rate0)
    log_rates = [log_rate]
    for i in range(1, len(times)):
        log_rate += drift * dt + sigma * np.sqrt(dt) * np.random.randn()
        log_rates.append(log_rate)
    times = np.append(times, t1)
    log_rates.append(log_rates[-1] + drift * dt + sigma * np.sqrt(dt) * np.random.randn())
    return times, np.exp(np.array(log_rates))


def draw_rate_tree(ax, obs_branches, aug_branches, drift, sigma,
                   rate0=0.3):
    """Draw the augmented tree colored by rates using the muted colormap.

    Rates are inherited at speciation: daughter branches start with the
    parent's rate at the branching time, then evolve independently via GBM.
    """
    all_rates = []
    branch_data = []

    all_branches = obs_branches + aug_branches
    is_aug = [False]*len(obs_branches) + [True]*len(aug_branches)

    # Process branches in order of start time so parents are computed first
    order = sorted(range(len(all_branches)), key=lambda i: all_branches[i].t0)

    # Store (times, rates) indexed by branch position in all_branches
    branch_rate_lookup = {}

    for idx in order:
        br = all_branches[idx]

        # Determine starting rate: inherit from parent branch at branching time
        start_rate = rate0  # default for root
        if br.parent_y is not None:
            # Find parent branch: same y as parent_y, time range covers br.t0
            for pidx, pbr in enumerate(all_branches):
                if pbr.y == br.parent_y and pbr.t0 <= br.t0 <= pbr.t1:
                    if pidx in branch_rate_lookup:
                        ptimes, prates = branch_rate_lookup[pidx]
                        # Interpolate parent rate at branching time
                        start_rate = np.interp(br.t0, ptimes, prates)
                    break

        times, rates = generate_rates_along_branch(
            br.t0, br.t1, start_rate, drift, sigma)
        branch_rate_lookup[idx] = (times, rates)
        all_rates.extend(rates)
        branch_data.append((br, times, rates, is_aug[idx]))

    vmin = np.percentile(all_rates, 5)
    vmax = np.percentile(all_rates, 95)
    norm = Normalize(vmin=vmin, vmax=vmax)

    for br, times, rates, aug in branch_data:
        alpha_val = 0.5 if aug else 1.0
        lw = 2.5 if aug else 3.5

        # Colored horizontal segments
        for i in range(len(times)-1):
            col = RATE_CMAP(norm((rates[i] + rates[i+1])/2))
            col = (*col[:3], alpha_val)
            ax.plot([times[i], times[i+1]], [br.y, br.y],
                    color=col, lw=lw, solid_capstyle='butt', zorder=2)

        # Vertical connector
        if br.parent_y is not None:
            col = RATE_CMAP(norm(rates[0]))
            col_a = (*col[:3], alpha_val)
            ax.plot([br.t0, br.t0], [br.parent_y, br.y],
                    color=col_a, lw=max(lw*0.6, 1.0),
                    solid_capstyle='round', zorder=1)

        # Fossils (brown diamonds, smaller on rate panels)
        # Only draw fossils on observed branches, not augmented ones
        if not aug:
            for ft in br.fossil_times:
                ax.plot(ft, br.y, 'D', color=COL_FOSSIL, markersize=5,
                        markeredgecolor='white', markeredgewidth=0.3, zorder=4)
            if br.extinct:
                ax.plot(br.t1, br.y, 'D', color=COL_FOSSIL, markersize=5,
                        markeredgecolor='white', markeredgewidth=0.3, zorder=4)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=RATE_CMAP, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.02, aspect=15)
    cbar.set_label('Rate value', fontsize=13)
    cbar.ax.tick_params(labelsize=11)

    return branch_data


def setup_ax(ax, title, y_lim=(-2.2, 11.2), show_xlabel=True):
    """Common axis setup."""
    ax.set_xlim(-0.5, T_PRESENT + 0.5)
    ax.set_ylim(y_lim)
    ax.set_title(title, fontsize=16, fontweight='bold', loc='left', pad=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([])
    if show_xlabel:
        ax.set_xlabel('Time', fontsize=14)
        ax.tick_params(axis='x', labelsize=12)
    else:
        ax.tick_params(axis='x', labelsize=12, labelbottom=False)


# ═══════════════════════════════════════════════════════════════════
# Build the figure
# ═══════════════════════════════════════════════════════════════════

fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 2, hspace=0.30, wspace=0.30,
                       left=0.04, right=0.96, top=0.96, bottom=0.10)

obs_branches = build_observed_tree()
aug_branches = build_augmented_tree(obs_branches)
occ_times = generate_occurrences()

# ── Panel A: Empirical tree ─────────────────────────────────────────
ax_a = fig.add_subplot(gs[0, 0])
setup_ax(ax_a, 'A   Empirical (observed) OBD tree', show_xlabel=False)

draw_tree(ax_a, obs_branches, color=COL_OBS, lw=2.2)
draw_occurrences(ax_a, occ_times, y_base=-0.7)

# Annotations — positioned to avoid title and other elements
ax_a.annotate('coded fossil', xy=(4.5, 1.5), xytext=(2.0, 0.3),
              fontsize=12, ha='center', color=COL_FOSSIL, fontweight='bold',
              arrowprops=dict(arrowstyle='->', color=COL_FOSSIL, lw=1.0))

ax_a.annotate('uncoded\noccurrences', xy=(9.2, -0.7), xytext=(8.0, 0.7),
              fontsize=12, ha='center', color=COL_OCC, fontweight='bold',
              arrowprops=dict(arrowstyle='->', color=COL_OCC, lw=1.0))

ax_a.plot([T_PRESENT, T_PRESENT], [-2.2, 10.0], color='gray', ls=':', lw=0.5, zorder=0)
ax_a.text(T_PRESENT, 10.8, 'present', ha='center', fontsize=12, color='gray')

# ── Panel B: Augmented tree ─────────────────────────────────────────
ax_b = fig.add_subplot(gs[1, 0])
setup_ax(ax_b, 'B   One complete (data-augmented) OBDD tree')

# Augmented lineages first (gray), then observed tree on top (black)
draw_tree(ax_b, aug_branches, color=COL_AUG, lw=1.5, draw_fossils=False,
          show_tips=False, alpha=0.6)
draw_tree(ax_b, obs_branches, color=COL_OBS, lw=2.2)
draw_occurrences(ax_b, occ_times, y_base=-0.7)

# Extinction crosses and extant tips for augmented branches
for br in aug_branches:
    if br.extinct:
        ax_b.plot(br.t1, br.y, 'x', color=COL_AUG, markersize=6,
                  markeredgewidth=1.8, zorder=4)
    else:
        # Unsampled extant tip (open circle)
        ax_b.plot(br.t1, br.y, 'o', color=COL_AUG, markersize=5,
                  markeredgecolor=COL_AUG, markerfacecolor='white',
                  markeredgewidth=1.5, zorder=4)

ax_b.annotate('data-augmented\nlineage', xy=(1.8, 9.2),
              xytext=(3.2, 9.8),
              fontsize=12, ha='center', color='#777777', fontweight='bold',
              arrowprops=dict(arrowstyle='->', color='#777777', lw=1.0))

ax_b.text(10.2, 1.5, 'unsampled\nextant lineage', fontsize=12, color='#777777',
          fontweight='bold', va='center')

ax_b.annotate('extinction', xy=(9.0, 5.2), xytext=(9.0, 7.2),
              fontsize=12, ha='center', color='#777777', fontweight='bold',
              arrowprops=dict(arrowstyle='->', color='#777777', lw=1.0))

ax_b.annotate('hidden\nspeciation', xy=(0.8, 5.0), xytext=(0.3, 3.2),
              fontsize=12, ha='center', color='#777777', fontweight='bold',
              arrowprops=dict(arrowstyle='->', color='#777777', lw=1.0))

ax_b.plot([T_PRESENT, T_PRESENT], [-2.2, 10.0], color='gray', ls=':', lw=0.5, zorder=0)
ax_b.text(T_PRESENT, 10.8, 'present', ha='center', fontsize=12, color='gray')

# ── Panel C: Speciation rate λ(t) ──────────────────────────────────
ax_c = fig.add_subplot(gs[0, 1])
setup_ax(ax_c, r'C   Latent speciation rate $\lambda(t)$', show_xlabel=False)

draw_rate_tree(ax_c, obs_branches, aug_branches,
               drift=-0.02, sigma=0.35, rate0=0.35)
draw_occurrences(ax_c, occ_times, y_base=-0.7)
ax_c.axvline(x=T_PRESENT, color='gray', ls=':', lw=0.5, zorder=0)

# ── Panel D: Extinction rate μ(t) ──────────────────────────────────
ax_d = fig.add_subplot(gs[1, 1])
setup_ax(ax_d, r'D   Latent extinction rate $\mu(t)$')

draw_rate_tree(ax_d, obs_branches, aug_branches,
               drift=0.02, sigma=0.4, rate0=0.15)
draw_occurrences(ax_d, occ_times, y_base=-0.7)
ax_d.axvline(x=T_PRESENT, color='gray', ls=':', lw=0.5, zorder=0)

# ── Save ────────────────────────────────────────────────────────────
fig.savefig('Images/Figure_augmented_tree_OBDD.pdf',
            bbox_inches='tight', dpi=300)
fig.savefig('Images/Figure_augmented_tree_OBDD.png',
            bbox_inches='tight', dpi=300)
print("Figure saved to Images/Figure_augmented_tree_OBDD.pdf")
plt.close()
