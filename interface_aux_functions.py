
# Needed Imports
import itertools
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hier
from scipy.stats import norm, chi2
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import panel as pn
import venn
import upsetplot

import holoviews as hv

# metanalysis_standard.py file
import metanalysis_standard as metsta

# The initial pages, especially the read file one does not have the nomenclature that I started using later on
# for the different widgets as well as organization
pn.extension('plotly', 'floatpanel', notifications=True)
hv.extension('plotly')
pn.config.sizing_mode="stretch_width"


### SVG Icons for buttons in the interface

home_icon = '''<svg xmlns="http://www.w3.org/2000/svg" class="icon icon-tabler icon-tabler-home" width="24" height="24" viewBox="0 0 24 24" stroke-width="2" stroke="currentColor" fill="none" stroke-linecap="round" stroke-linejoin="round">
   <path stroke="none" d="M0 0h24v24H0z" fill="none"></path>
   <path d="M5 12l-2 0l9 -9l9 9l-2 0"></path>
   <path d="M5 12v7a2 2 0 0 0 2 2h10a2 2 0 0 0 2 -2v-7"></path>
   <path d="M9 21v-6a2 2 0 0 1 2 -2h2a2 2 0 0 1 2 2v6"></path>
</svg>'''

img_confirm_button = '''<svg xmlns="http://www.w3.org/2000/svg"
    class="icon icon-tabler icon-tabler-check" width="24" height="24" viewBox="0 0 24 24" stroke-width="2"
    stroke="currentColor" fill="none" stroke-linecap="round" stroke-linejoin="round">
    <path stroke="none" d="M0 0h24v24H0z" fill="none"></path>
    <path d="M5 12l5 5l10 -10"></path>
    </svg>'''

download_icon = '''<svg xmlns="http://www.w3.org/2000/svg" class="icon icon-tabler icon-tabler-download" width="24" height="24" viewBox="0 0 24 24" stroke-width="2" stroke="currentColor" fill="none" stroke-linecap="round" stroke-linejoin="round">
   <path stroke="none" d="M0 0h24v24H0z" fill="none"></path>
   <path d="M4 17v2a2 2 0 0 0 2 2h12a2 2 0 0 0 2 -2v-2"></path>
   <path d="M7 11l5 5l5 -5"></path>
   <path d="M12 4l0 12"></path>
</svg>'''


### Functions related to the common and exclusive compound page of the graphical interface

def _group_compounds_per_class(com_exc_compounds, target_list, DataFrame_Store):
    "Creates a DataFrame per biological class with the features that appear in samples of said class (and provides a description of each DataFrame)."

    # Creating the dictionary with information about which samples belong to which classes
    for cl in target_list.color_classes.keys(): # Setting up the keys (classes)
        com_exc_compounds.groups[cl] = []

    for c, t in zip(target_list.sample_cols, target_list.target): # Setting up the values
        for g in com_exc_compounds.groups:
            if g == t:
                com_exc_compounds.groups[g].append(c)

    # Creating the DataFrame for each class by dropping features of the processed complete DataFrame if they do not appear in any sample of the class
    # One DataFrame is made considering every feature in the dataset, and another considering only annotated features
    for g in com_exc_compounds.groups:
        for c, t in zip(target_list.sample_cols, target_list.target):
            if g == t:
                com_exc_compounds.group_dfs[g] = DataFrame_Store.processed_df.dropna(
                    subset=com_exc_compounds.groups[g], thresh=1)
                com_exc_compounds.group_dfs_ids[g] = com_exc_compounds.group_dfs[g].iloc[[
                    i for i in range(len(com_exc_compounds.group_dfs[
                    g]['Has Match?'])) if com_exc_compounds.group_dfs[g]['Has Match?'][i]]]

    # Description of the number of metabolites (and annotated metabolites) that appear in at least one sample of each class
    desc_string = ['**Nº of peaks per class:**', '', ]
    for g in com_exc_compounds.groups:
        group_string = f'**{g}**: **{len(com_exc_compounds.group_dfs[g])}** metabolites, from which **{len(com_exc_compounds.group_dfs_ids[g])}** have matches.'
        desc_string.append(group_string)

    com_exc_compounds.groups_description = '<br />'.join(desc_string) # <br /> leads to line breaks in Markdown


def _compute_com_exc_compounds(com_exc_compounds):
    "Compute the number of common compounds to all classes and exclusive of each class (providing a small description)"
    # Common to all classes
    com_exc_compounds.common_all = metsta.common(com_exc_compounds.group_dfs.values())
    com_exc_compounds.common_all_id = metsta.common(com_exc_compounds.group_dfs_ids.values())

    # Exclusive of each class
    exc = metsta.exclusive(com_exc_compounds.group_dfs.values())
    exc_id = metsta.exclusive(com_exc_compounds.group_dfs_ids.values())
    com_exc_compounds.exclusives = dict(zip(com_exc_compounds.group_dfs, exc))
    com_exc_compounds.exclusives_id = dict(zip(com_exc_compounds.group_dfs, exc_id))

    # Description of the number of metabolites (and annotated metabolites) that appear in all classes and are exclusive to each class
    desc_string = [f'**Nº of Compounds Common to All Classes**: **{len(com_exc_compounds.common_all.index)}** metabolites, **{len(com_exc_compounds.common_all_id.index)}** of which annotated.',
                   '', '**Nº of Compounds Exclusive to each Class:**', '',]

    for g in com_exc_compounds.exclusives:
        group_string = f'**{g}**: **{len(com_exc_compounds.exclusives[g].index)}** exclusive metabolites (**{len(com_exc_compounds.exclusives_id[g].index)}** with matches).'
        desc_string.append(group_string)

    com_exc_compounds.com_exc_desc = '<br />'.join(desc_string)


# From Stack Overflow (https://stackoverflow.com/questions/29643352/converting-hex-to-rgb-value-in-python)
def RGB(col): return '#%02x%02x%02x' % (int(col[0]), int(col[1]), int(col[2]))
def hex_to_rgb(value):
    value = value.lstrip('#')
    size = len(value)
    return tuple(int(value[i:i + size // 3], 16) for i in range(0, size, size // 3))


def _plot_Venn_diagram(com_exc_compounds, target_list):
    "Plot Venn diagrams using the venn package."

    # Calculate different intersections on all features and considering only annotated features
    labels = venn.get_labels([com_exc_compounds.group_dfs[i].index for i in com_exc_compounds.venn_class_subset],
                             fill=['number'])
    labels_ids = venn.get_labels([com_exc_compounds.group_dfs_ids[i].index for i in com_exc_compounds.venn_class_subset],
                                 fill=['number'])

    if com_exc_compounds.type_of_venn == 'All Metabolites (Annotated)': # Show the metabolites (and annotated) numbers
        # Join the labels of all features and considering only annotated features
        labels_all = {}
        for i, j in labels.items():
            labels_all[i] = j + f' ({labels_ids[i]})'
        x_label = 'Nº of peaks (Nº of matched compounds)'
    elif com_exc_compounds.type_of_venn == 'All Metabolites': # Only show the metabolite numbers
        labels_all = labels
        x_label = 'Nº of peaks'
    else: # Only show the annotated numbers
        labels_all = labels_ids
        x_label = 'Nº of matched compounds'

    # Get the colours in RGB format to draw the Venn
    colours = []
    for c in com_exc_compounds.venn_class_subset:
        rgb_c = hex_to_rgb(target_list.color_classes[c])
        colours.append((rgb_c[0]/255, rgb_c[1]/255, rgb_c[2]/255, com_exc_compounds.venn_alpha))

    n_class = len(com_exc_compounds.venn_class_subset)

    # Draw the Venn based on the number of classes you have
    if n_class == 2:
        fig, ax = venn.venn2(labels_all, names=com_exc_compounds.venn_class_subset, figsize=(8,8), fontsize=11,
                             colors=colours, dpi=com_exc_compounds.dpi_venn, constrained_layout=True) # 2 Classes
        plt.text(0.5,0, x_label, fontsize=12, horizontalalignment='center')
        leg = ax.legend(com_exc_compounds.venn_class_subset, loc='upper center', bbox_to_anchor=(0,1.1),
                        fancybox=True, framealpha=0.5)
        return fig
    elif n_class == 3:
        fig, ax = venn.venn3(labels_all, names=com_exc_compounds.venn_class_subset, figsize=(8,8), fontsize=11,
                             colors=colours, dpi=com_exc_compounds.dpi_venn, constrained_layout=True) # 3 Classes
        leg = ax.legend(com_exc_compounds.venn_class_subset, loc='upper center', bbox_to_anchor=(0,1.1),
                        fancybox=True, framealpha=0.5)
        plt.text(0.5,-0.05, x_label, fontsize=12, horizontalalignment='center')
        return fig
    elif n_class == 4:
        fig, ax = venn.venn4(labels_all, names=com_exc_compounds.venn_class_subset, figsize=(8,8), fontsize=11,
                             colors=colours, dpi=com_exc_compounds.dpi_venn, constrained_layout=True) # 4 Classes
        leg = ax.legend(com_exc_compounds.venn_class_subset, loc='upper center', bbox_to_anchor=(0,1.1),
                        fancybox=True, framealpha=0.5)
        plt.text(0.5,0.05, x_label, fontsize=12, horizontalalignment='center')
        return fig
    elif n_class == 5:
        fig, ax = venn.venn5(labels_all, names=com_exc_compounds.venn_class_subset, figsize=(8,8), fontsize=11,
                             colors=colours, dpi=com_exc_compounds.dpi_venn, constrained_layout=True) # 5 Classes
        leg = ax.legend(com_exc_compounds.venn_class_subset, loc='upper center', bbox_to_anchor=(0,1.1),
                        fancybox=True, framealpha=0.5)
        plt.text(0.5,0, x_label, fontsize=12, horizontalalignment='center')
        return fig
    elif n_class == 6:
        fig, ax = venn.venn6(labels_all, names=com_exc_compounds.venn_class_subset, figsize=(8,8), fontsize=11,
                             colors=colours, dpi=com_exc_compounds.dpi_venn, constrained_layout=True) # 6 Classes
        plt.text(0.5,0.2, x_label, fontsize=12, horizontalalignment='center')
        leg = ax.legend(com_exc_compounds.venn_class_subset, loc='upper center', bbox_to_anchor=(0,1.1),
                        fancybox=True, framealpha=0.5)
        return fig
    else:
        pn.state.notifications.info(f'Venn Diagram can only be made with 2 to 6 different classes. You currently have {n_class} classes.',
                                    duration=2000)


def _plot_upsetplots(com_exc_compounds, groups_dict, ups):
    "Plot an Upsetplot"
    # Plotting UpSetPlot
    f,ax = plt.subplots(1,1, constrained_layout=True, dpi=400)
    if com_exc_compounds.upset_include_counts_percentages == 'Show Nº and % of metabolites':
        include_counts, include_percentages = True, True
    elif com_exc_compounds.upset_include_counts_percentages == 'Show Nº of metabolites':
        include_counts, include_percentages = True, False
    else:
        include_counts, include_percentages = False, False
    ax.axis('Off')
    upsetplot.plot(ups, f, subset_size='count',
                   show_counts=include_counts,
                   show_percentages=include_percentages,
                   sort_categories_by='input', include_empty_subsets=False)
    return f


### Functions related to the PCA section of the unsupervised analysis page of the graphical interface
# Functions slightly changed from the elips.py file to be ran with plotly

def plot_confidence_ellipse(points, q=None, nstd=2):
    """
    Plots an `nstd` sigma ellipse based on the mean and covariance of a point
    "cloud" (points, an Nx2 array).

    Parameters
    ----------
        points : An Nx2 array of the data points.
        q : float, optional
            Confidence level, should be in (0, 1)
        nstd : int, optional
            Confidence level in unit of standard deviations.
        E.g. 1 stands for 68.3% and 2 stands for 95.4%.
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    if q != 0:
        q = np.asarray(q)
    elif nstd is not None:
        q = 2 * norm.cdf(nstd) - 1
    else:
        raise ValueError('One of `q` and `nsig` should be specified.')

    pos = points.mean(axis=0)
    cov = np.cov(points, rowvar=False)
    return plot_cov_ellipse(cov, pos, q)

def plot_cov_ellipse(cov, pos, q):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        q : float, optional
            Confidence level, should be in (0, 1)
        ax : The axis that the ellipse will be plotted on. Defaults to the
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    r2 = chi2.ppf(q, 2)

    vals, vecs = eigsorted(cov)
    # Width and height are "full" widths, not radius
    # width, height = 2 * nstd * np.sqrt(vals)
    width, height = 2 * np.sqrt(vals * r2)
    theta = np.arctan2(*vecs[:,0][::-1])

    ellip = hv.Ellipse(pos[0], pos[1], (width,height), orientation=theta)

    return ellip



### Functions related to the HCA section of the unsupervised analysis page of the graphical interface

def color_list_to_matrix_and_cmap(colors, ind, axis=0):
        if any(issubclass(type(x), list) for x in colors):
            all_colors = set(itertools.chain(*colors))
            n = len(colors)
            m = len(colors[0])
        else:
            all_colors = set(colors)
            n = 1
            m = len(colors)
            colors = [colors]
        color_to_value = dict((col, i) for i, col in enumerate(all_colors))

        matrix = np.array([color_to_value[c]
                           for color in colors for c in color])

        matrix = matrix.reshape((n, m))
        matrix = matrix[:, ind]
        if axis == 0:
            # row-side:
            matrix = matrix.T

        cmap = mpl.colors.ListedColormap(all_colors)
        return matrix, cmap

def plot_dendogram(Z, leaf_names, label_colors, title='', ax=None, x_axis_len=4, no_labels=False, labelsize=12, **kwargs):
    if ax is None:
        ax = plt.gca()
    hier.dendrogram(Z, labels=leaf_names, leaf_font_size=10, above_threshold_color='0.2', orientation='left',
                    ax=ax, **kwargs)
    #Coloring labels
    #ax.set_ylabel('Distance (AU)')
    ax.set_xlabel('Distance (AU)')
    ax.set_title(title, fontsize = 15)

    #ax.tick_params(axis='x', which='major', pad=12)
    ax.tick_params(axis='y', which='major', labelsize=labelsize, pad=x_axis_len*3)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    #xlbls = ax.get_xmajorticklabels()
    xlbls = ax.get_ymajorticklabels()
    rectimage = []
    for lbl in xlbls:
        col = label_colors[lbl.get_text()]
        lbl.set_color(col)
        #lbl.set_fontweight('bold')
        if no_labels:
            lbl.set_color('w')
        rectimage.append(col)

    cols, cmap = color_list_to_matrix_and_cmap(rectimage, range(len(rectimage)), axis=0)

    axins = inset_axes(ax, width="5%", height="100%",
                   bbox_to_anchor=(1, 0, 1, 1),
                   bbox_transform=ax.transAxes, loc=3, borderpad=0)

    axins.pcolor(cols, cmap=cmap, edgecolors='w', linewidths=1)
    axins.axis('off')