import pathlib
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product

_SIG_TYPE = "signature_type"
_META_SIG = "meta_signature"
_GROUND_TRUTH = "ground_truth_signature"


def get_contribution(cancer_type, scenario, signature):
    results = []
    n_sigs = int(signature.split("_")[1])
    for i in range(n_sigs, n_sigs+3):
        scores = None
        for method_path in pathlib.Path(f"../results/{cancer_type}/corrs/scrna").iterdir():
            method = method_path.stem
            corr_path = method_path.joinpath(f"{scenario}_all_0/{signature}/{i}/corrs.csv")
            if not corr_path.is_file():
                continue
            corr = pd.read_csv(corr_path, index_col=0)

            sig_type = corr.pop(_SIG_TYPE)
            sig_gts = corr.index[sig_type == _GROUND_TRUTH].tolist()
            if scores is None:
                scores = pd.DataFrame(columns=sig_gts)
            sig_names = corr.index[sig_type == _META_SIG].tolist()

            corr[corr < 0.] = 0.
            contribution = pd.DataFrame(index=sig_names, columns=sig_gts)
            for sig_name in sig_names:
                for sig_gt in sig_gts:
                    other_gts = sig_gts.copy()
                    other_gts.remove(sig_gt)
                    postive_score = corr.loc[sig_name, sig_gt]
                    negative_score = np.max(np.maximum(0, corr.loc[[sig_name], other_gts].values -
                                                       corr.loc[[sig_gt], other_gts].values))
                    score = postive_score - negative_score
                    score = np.maximum(0, score)
                    contribution.loc[sig_name, sig_gt] = score
            score = {}
            for _ in sig_gts:
                if contribution.shape[0] == 0:
                    break
                max_idx = np.unravel_index(np.argmax(contribution.values, axis=None), contribution.values.shape)
                score[contribution.columns[max_idx[1]]] = contribution.iloc[max_idx].item()
                contribution = contribution.drop(index=contribution.index[max_idx[0]],
                                                 columns=contribution.columns[max_idx[1]])
            scores.loc[method] = score
        scores["#Cluster"] = i
        results.append(scores)
    results = pd.concat(results, axis=0)

    return results


def plot_heatmap_with_bar_plot(df, method_colormap,figsize=(8, 12), xlabel="Dataset (#Signatures)", x_tick_top=True):
    _offset = -0.3
    n_ticks = df.shape[1]

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 2, width_ratios=[4, 1], hspace=0.0, wspace=0.0)
    main_ax = fig.add_subplot(gs[0, 0])
    bar_ax = fig.add_subplot(gs[0, 1], sharey=main_ax)

    bar_color = []
    for method, _ in product(df.index.get_level_values(level=0).unique(), range(3)):
            bar_color.append(method_colormap[method])


    hmap = sns.heatmap(df, vmax=1., vmin=0., cmap="OrRd", ax=main_ax, cbar=False, annot=True)
    means = df.T.mean(axis=1)
    bar_ax.barh(range(len(means)), means, height=1, align="edge", color=bar_color)
    #bar_ax.set_ylim(main_ax.get_ylim())
    bar_ax.set_xlim((means.min()-0.05, means.max()+0.05))
    bar_ax.set_title('Average\nScore')
    bar_ax.grid(True, axis='x')
    bar_ax.set_yticks([])

    for i, (method, cluster) in enumerate(df.index):
        if not i % 3:
            rectangle = Rectangle((_offset, 1 - (i+3)/n_ticks), -_offset, 3/n_ticks, color=method_colormap[method], transform=main_ax.transAxes, clip_on=False)
            main_ax.add_patch(rectangle)
            main_ax.text(_offset+0.015, 1 - (i+3)/n_ticks + 3/(2*n_ticks), method, ha="left", va="center", transform=main_ax.transAxes, rotation=0)
            current_method = method
        else:
            if method != current_method:
                raise ValueError("Dataframe was not sorted!")

        main_ax.text(-0.05, 1 - (i + 1)/(n_ticks) + 1/(2*(n_ticks)), i % 3, ha="center", va="center", transform=main_ax.transAxes)

    main_ax.hlines(y=[1 - i/(n_ticks) for i in range(1, (n_ticks))], xmin=[-0.1 - 0.2 * ((i % 3)==0) for i in range(1, (n_ticks))], xmax=0.0, colors="black", transform=main_ax.transAxes, clip_on=False)

    main_ax.yaxis.set_tick_params(labelleft=True)
    main_ax.set_yticks([])
    main_ax.set_yticklabels([])
    main_ax.set_ylabel("")
    main_ax.set_xlabel(xlabel=xlabel)
    if x_tick_top:
        main_ax.xaxis.tick_top()
        main_ax.set_xticks(main_ax.get_xticks(), main_ax.get_xticklabels(), rotation=45, ha='left')

    cbar = fig.colorbar(hmap.get_children()[0], ax=[main_ax, bar_ax], orientation='horizontal', label='Score', 
                    fraction=0.025, pad=0.05)
    return fig