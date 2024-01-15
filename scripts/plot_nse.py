import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def load_species(species_file):
    species, a, z = np.genfromtxt(
        species_file,
        dtype=str,
        unpack=True,
        skip_header=1,
    )
    a = a.astype(int)
    z = z.astype(int)

    arrays = [a, z]
    tuples = list(zip(*arrays))
    index = pd.MultiIndex.from_tuples(tuples, names=["A", "Z"])
    s = pd.DataFrame({"Name": species, "A": a, "Z": z}, index=index)
    return s


def load_abundances(abundances_file):
    name, a, z, abundance = np.genfromtxt(
        abundances_file,
        dtype=str,
        unpack=True,
        skip_header=1,
        delimiter=",",
    )
    a = a.astype(int)
    z = z.astype(int)
    abundance = abundance.astype(float)

    abundances = pd.DataFrame(
        {"Name": name, "A": a, "Z": z, "Abundance": abundance},
    )
    return abundances


def _annotate(ax, x, y):
    # this all gets repeated below:
    X, Y = np.meshgrid(x, y)
    ax.plot(X.flat, Y.flat, "o", color="m")


def plot_magic_numbers(ax):
    ax.axhline(7.5, color="k", linestyle=":", alpha=0.5)
    ax.axhline(8.5, color="k", linestyle=":", alpha=0.5)
    ax.axhline(19.5, color="k", linestyle=":", alpha=0.5)
    ax.axhline(20.5, color="k", linestyle=":", alpha=0.5)
    ax.axhline(27.5, color="k", linestyle=":", alpha=0.5)
    ax.axhline(28.5, color="k", linestyle=":", alpha=0.5)
    ax.axhline(49.5, color="k", linestyle=":", alpha=0.5)
    ax.axhline(50.5, color="k", linestyle=":", alpha=0.5)

    ax.axvline(7.5, color="k", linestyle=":", alpha=0.5)
    ax.axvline(8.5, color="k", linestyle=":", alpha=0.5)
    ax.axvline(19.5, color="k", linestyle=":", alpha=0.5)
    ax.axvline(20.5, color="k", linestyle=":", alpha=0.5)
    ax.axvline(27.5, color="k", linestyle=":", alpha=0.5)
    ax.axvline(28.5, color="k", linestyle=":", alpha=0.5)
    ax.axvline(49.5, color="k", linestyle=":", alpha=0.5)
    ax.axvline(50.5, color="k", linestyle=":", alpha=0.5)


def main():
    species = load_species("data/species.txt")
    species["N"] = species["A"] - species["Z"]

    abundances = load_abundances("output.txt")
    abundances["N"] = abundances["A"] - abundances["Z"]
    abundances = abundances[abundances["Abundance"] > 1e-15]
    abundances["Abundance"] = np.log10(abundances["Abundance"])

    stable = pd.read_csv("data/stable_isotopes.csv")
    stable["A"] = stable["N"] + stable["Z"]
    stable = stable[stable["A"].isin(species["A"])]

    fig, ax = plt.subplots(1, 1, figsize=(6.4 * 2, 4.8 * 2))

    ax.scatter(
        species["N"],
        species["Z"],
        s=100,
        alpha=0.5,
        marker="s",
        facecolors="none",
        edgecolors="k",
        label="Species included in network",
    )
    ax.scatter(
        stable["N"],
        stable["Z"],
        s=100,
        marker="s",
        facecolors="none",
        edgecolors="k",
        label="Stable isotopes",
        zorder=10,
    )

    ax.scatter(
        abundances["N"],
        abundances["Z"],
        s=100,
        c=abundances["Abundance"],
        marker="s",
        cmap="viridis",
        label="NSE abundances",
    )

    plot_magic_numbers(ax)

    cbar = fig.colorbar(ax.collections[2], ax=ax)
    cbar.set_label("log(X)")

    ax.set_xlabel("Number of Neutrons")
    ax.set_ylabel("Number of Protons")
    plt.legend()

    plt.tight_layout()

    fig.savefig("test.png", bbox_inches="tight")

    return


if __name__ == "__main__":
    main()
