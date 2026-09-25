"""Shared figure style. Categorical palette validated for CVD (Okabe-Ito derived)."""
import matplotlib as mpl, matplotlib.pyplot as plt
mpl.use("Agg")
CAT = ["#0072B2", "#D55E00", "#009E73", "#7A5195"]      # validated: all six checks pass
UP, DOWN, NEUTRAL = "#B2182B", "#2166AC", "#BFBFBF"      # diverging + neutral midpoint
LOW, HIGH = "#0072B2", "#D55E00"
SEQ = "Blues"
INK, INK2, MUTED = "#1a1a1a", "#444444", "#7a7a7a"

mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7.5,
    "axes.titleweight": "bold", "axes.labelcolor": INK, "text.color": INK,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.edgecolor": "#999999", "axes.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "grid.color": "#E4E4E4", "grid.linewidth": 0.5, "legend.frameon": False,
    "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.facecolor": "white",
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.4,
})

def panel(ax, letter, dx=-0.16, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=10, fontweight="bold", va="top", ha="left")

def save(fig, name, outdir="/home/claude/an/figs"):
    import os; os.makedirs(outdir, exist_ok=True)
    fig.savefig(f"{outdir}/{name}.png"); fig.savefig(f"{outdir}/{name}.pdf")
    plt.close(fig); print("wrote", name)
