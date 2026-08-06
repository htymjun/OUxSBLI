"""Shared matplotlib styling for the user-facing plotting scripts."""


def apply(fontsize=20):
    """Serif/STIX styling with inward ticks, matching the figure conventions
    used across the case plotting scripts."""
    import matplotlib.pyplot as plt

    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "stix"
    plt.rcParams["xtick.direction"] = "in"
    plt.rcParams["ytick.direction"] = "in"
    plt.rcParams["font.size"] = fontsize
