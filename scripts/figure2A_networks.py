import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

from utils_python.connect_db import read_translation_data


def get_giant_component(g):
    gcc = sorted(nx.connected_components(g), key=len, reverse=True)[0]
    gc = g.subgraph(gcc)

    return gc


coexp = read_translation_data()
seed = 10001
canonical_orfs = coexp.loc[coexp.is_canonical == "canonical"].transcript.values

canonical_only_file = (
    "/home/aar75/coexpression/20221110_rho/networkChanges/networks/canonical_999.csv"
)
canonical_noncanonical_file = (
    "/home/aar75/coexpression/20221110_rho/networkChanges/networks/whole_999.csv"
)

canonical_only_edgelist = pd.read_csv(canonical_only_file)
canonical_noncanonical_edgelist = pd.read_csv(canonical_noncanonical_file)

g_canonical_only = nx.from_pandas_edgelist(
    canonical_only_edgelist, source="row", target="column"
)
g_canonical_noncanonical = nx.from_pandas_edgelist(
    canonical_noncanonical_edgelist, source="row", target="column"
)

random_cn_nn_edges_file = (
    "/home/oma21/coexpression/data/interim/randomized_networks/0.csv"
)
random_cn_nn_edges_edgelist = pd.read_csv(random_cn_nn_edges_file, header=None)
g_random_cn_nn_edges_edgelist = nx.from_pandas_edgelist(
    random_cn_nn_edges_edgelist, source=1, target=2
)

gc_canonical_only = get_giant_component(g_canonical_only)
gc_canonical_noncanonical = get_giant_component(g_canonical_noncanonical)
gc_random_cn_nn_edges = get_giant_component(g_random_cn_nn_edges_edgelist)

pos_canonical_only = nx.spring_layout(gc_canonical_only, k=0.06, seed=seed)

canonical_color = "#7570b3"
noncanonical_color = "#1b9e77"

canonical_color_talk = "#CC6677"
noncanonical_color_talk = "#89CCED"

# Plot canonical only
fig, ax = plt.subplots(figsize=(6, 6), facecolor="white")

nx.draw_networkx_nodes(
    gc_canonical_only,
    pos=pos_canonical_only,
    node_size=5,
    ax=ax,
    node_color=[
        canonical_color if i in canonical_orfs else noncanonical_color
        for i in gc_canonical_only.nodes
    ],
    alpha=0.5,
)
nx.draw_networkx_edges(gc_canonical_only, pos=pos_canonical_only, width=0.1, alpha=0.5)
ax.axis("off")
plt.savefig(
    "../reports/figures/paper_figures_01062023/networkx_plot_canonical_only.png",
    dpi=150,
    facecolor="white",
)

# plot canonical only with talk colors
fig, ax = plt.subplots(figsize=(10, 10), facecolor="white")

nx.draw_networkx_nodes(
    gc_canonical_only,
    pos=pos_canonical_only,
    node_size=20,
    ax=ax,
    node_color=[
        canonical_color_talk if i in canonical_orfs else noncanonical_color_talk
        for i in gc_canonical_only.nodes
    ],
    alpha=0.5,
)
nx.draw_networkx_edges(gc_canonical_only, pos=pos_canonical_only, width=0.2, alpha=0.5)
ax.axis("off")
plt.savefig(
    "../reports/figures/paper_figures_01062023/networkx_plot_canonical_only_talk.png",
    dpi=150,
    facecolor="white",
)


pos_canonical_noncanonical = nx.spring_layout(
    g_canonical_noncanonical, k=0.06, seed=seed
)

fig, ax = plt.subplots(figsize=(6, 6), facecolor="white")

nx.draw_networkx_nodes(
    g_canonical_noncanonical,
    pos=pos_canonical_noncanonical,
    node_size=5,
    ax=ax,
    node_color=[
        canonical_color if i in canonical_orfs else noncanonical_color
        for i in g_canonical_noncanonical.nodes
    ],
    alpha=0.5,
)
nx.draw_networkx_edges(
    g_canonical_noncanonical, pos=pos_canonical_noncanonical, width=0.1, alpha=0.5
)
ax.axis("off")
plt.savefig(
    "../reports/figures/paper_figures_01062023/networkx_plot_expanded.png",
    dpi=150,
    facecolor="white",
)

fig, ax = plt.subplots(figsize=(10, 10), facecolor="white")

nx.draw_networkx_nodes(
    g_canonical_noncanonical,
    pos=pos_canonical_noncanonical,
    node_size=20,
    ax=ax,
    node_color=[
        canonical_color_talk if i in canonical_orfs else noncanonical_color_talk
        for i in g_canonical_noncanonical.nodes
    ],
    alpha=0.8,
)
nx.draw_networkx_edges(
    g_canonical_noncanonical, pos=pos_canonical_noncanonical, width=0.2, alpha=0.5
)
ax.axis("off")
plt.savefig(
    "../reports/figures/paper_figures_01062023/networkx_plot_expanded_talk.png",
    dpi=150,
    facecolor="white",
)

pos_random_cn_nn_06 = nx.spring_layout(gc_random_cn_nn_edges, k=0.06, seed=seed)
fig, ax = plt.subplots(figsize=(6, 6), facecolor="white")

nx.draw_networkx_nodes(
    gc_random_cn_nn_edges,
    pos=pos_random_cn_nn_06,
    node_size=5,
    ax=ax,
    node_color=[
        canonical_color if i in canonical_orfs else noncanonical_color
        for i in gc_random_cn_nn_edges.nodes
    ],
    alpha=0.5,
)
nx.draw_networkx_edges(
    gc_random_cn_nn_edges, pos=pos_random_cn_nn_06, width=0.1, alpha=0.5
)
ax.axis("off")
plt.savefig(
    "../reports/figures/paper_figures_01062023/networkx_plot_random.png",
    dpi=150,
    facecolor="white",
)

fig, ax = plt.subplots(figsize=(10, 10), facecolor="white")

nx.draw_networkx_nodes(
    gc_random_cn_nn_edges,
    pos=pos_random_cn_nn_06,
    node_size=20,
    ax=ax,
    node_color=[
        canonical_color_talk if i in canonical_orfs else noncanonical_color_talk
        for i in gc_random_cn_nn_edges.nodes
    ],
    alpha=0.5,
)
nx.draw_networkx_edges(
    gc_random_cn_nn_edges, pos=pos_random_cn_nn_06, width=0.2, alpha=0.5
)
ax.axis("off")
plt.savefig(
    "../reports/figures/paper_figures_01062023/networkx_plot_random_talk.png",
    dpi=150,
    facecolor="white",
)
