
# Needed Imports
import itertools
import pandas as pd
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hier
import random
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


# Functions related to Data pre-processing and pre-treatment that needed to be adapted for the graphical interface

def initial_filtering(df, sample_cols, target=None, filt_method='total_samples', filt_kw=2):
    "Performs feature filtering"
    meta_cols = [i for i in df.columns if i not in sample_cols]
    temp = metsta.basic_feat_filtering(df[sample_cols].T, target=target,
                                filt_method=filt_method, # Method
                                filt_kw=filt_kw) # Make a sample have to appear
    df = pd.concat((df[meta_cols].reindex(temp.columns), temp.T), axis=1)
    data_characteristics = metsta.characterize_data(df[sample_cols].T, target=target)

    return df, pd.DataFrame(pd.Series(data_characteristics)).iloc[1:]


# Read the database function
def read_database(filename, abv, ID_col, name_col, formula_col):
    # If the file has the right terminations
    if filename.endswith('.csv'):
        try:
            db = pd.read_csv(filename)
        except:
            pn.state.notifications.error(f'Database {abv} could not be read. File could not be found. Check if all parameters given are correct.',
                                     duration=5000)
            raise ValueError(f'Database {abv} could not be read. File could not be found. Check if all parameters given are correct.')
    elif filename.endswith('.xlsx'):
        try:
            db = pd.read_excel(filename)
        except:
            pn.state.notifications.error(f'Database {abv} could not be read. File could not be found. Check if all parameters given are correct.',
                                     duration=5000)
            raise ValueError(f'Database {abv} could not be read. File could not be found. Check if all parameters given are correct.')
    else:
        pn.state.notifications.error(f'Database {abv} could not be read. File Format not accepted. Only csv and xlsx files are accepted.',
                                    duration=5000)
        raise ValueError(f'Database {abv} could not be read. File Format not accepted. Only csv and xlsx files are accepted.')

    if abv == '': # If abbreviation was set
        pn.state.notifications.error(f'Database could not be read. Abbreviation was not set.',
                                    duration=5000)
        raise ValueError(f'Database could not be read. Abbreviation was not set.')

    if ID_col in db.columns: # If ID col is found
        db = db.set_index(ID_col)
        if abv == 'HMDB':
            db[name_col] = db[name_col].str.replace("b'", "")
            db[name_col] = db[name_col].str.replace("'", "")
    else:
        pn.state.notifications.error(f'Database {abv} could not be read. DB ID - Index Column given was not found.',
                                    duration=5000)
        f'Database {abv} could not be read. DB ID - Index Column given was not found.'

    if name_col not in db.columns: # If DB name col is not found
        pn.state.notifications.error(f'Database {abv} could not be read. DB Annotation Column given was not found.',
                                    duration=5000)
        raise ValueError(f'Database {abv} could not be read. DB Annotation Column given was not found.')
    elif formula_col not in db.columns: # If formula col is not found
        pn.state.notifications.error(f'Database {abv} could not be read. DB Formula Column given was not found.',
                                    duration=5000)
        raise ValueError(f'Database {abv} could not be read. DB Formula Column given was not found.')
    ##
    # Calculate Mass based on Formula column
    db['Mass'] = db[formula_col].dropna().apply(metsta.calculate_monoisotopic_mass)
    return db


def creating_has_match_column(DataFrame_Store, n_databases, checkbox_annotation):
    "Creating the has match column to be used later (common/exclusive compounds)."

    # If annotation was not performed in this software
    if n_databases.value == 0:
        for i in DataFrame_Store.original_df.index:
            df = DataFrame_Store.original_df.loc[[i]]
            hasmatch = df[checkbox_annotation.value].notnull().values.any()
            DataFrame_Store.original_df.at[i, 'Has Match?'] = hasmatch

    # If annotation was performed in this software
    else:
        # Join columns with annotations not made by this software with the ones made by this software
        cols_to_see = checkbox_annotation.value + list(DataFrame_Store.original_df.columns[-n_databases.value*4:])

        DataFrame_Store.original_df['Has Match?'] = np.nan # Creating the column
        for i in DataFrame_Store.original_df.index:
            df = DataFrame_Store.original_df.loc[[i]]
            hasmatch = df[cols_to_see].notnull().values.any()
            DataFrame_Store.original_df.at[i, 'Has Match?'] = hasmatch


def performing_pretreatment(PreTreatment_Method, DataFrame_Store, target, sample_cols):
    "Putting keywords to pass to filtering_pretreatment function and performing pre-treatment."

    # Missing Value Imputation
    mvi_translation = {'Minimum of Sample':'min_sample', 'Minimum of Feature':'min_feat', 'Minimum of Data':'min_data',
                      'Zero':'zero'}
    mvi = mvi_translation[PreTreatment_Method.mvi_method]
    mvi_kw = PreTreatment_Method.mvi_kw

    # Normalization
    norm_translation = {'Reference Feature':'ref_feat', 'Total Intensity Sum':'total_sum', 'PQN':'PQN',
                      'Quantile':'Quantile', 'None': None}
    norm = norm_translation[PreTreatment_Method.norm_method]

    norm_kw = PreTreatment_Method.norm_kw
    if norm in ['PQN', 'Quantile']:
        norm_kw = norm_kw[0]
    if norm_kw == 'Mean':
        norm_kw = 'mean'
    elif norm_kw == 'Median':
        norm_kw = 'median'

    # Transformation
    tf_translation = {'Generalized Logarithmic Transformation (glog)':'glog', None: None}
    tf = tf_translation[PreTreatment_Method.tf_method]

    tf_kw = PreTreatment_Method.tf_kw

    # Scaling
    scaling_translation = {'Pareto Scaling':'pareto', 'Auto Scaling':'auto', 'Mean Centering':'mean_center',
                      'Range Scaling':'range', 'Vast Scaling': 'vast', 'Level Scaling': 'level', None: None}
    scaling = scaling_translation[PreTreatment_Method.scaling_method]

    if PreTreatment_Method.scaling_kw == 'Average':
        scaling_kw = True
    else:
        scaling_kw = False

    # Data to pass
    data = DataFrame_Store.original_df

    a,b,c,d,e = metsta.filtering_pretreatment(data, target, sample_cols,
                      filt_method=None, filt_kw=2, # Filtering based on number of times features appear
                      extra_filt=None, # Filtering based on annotation of features ('Formula' or 'Name')
                      mvi=mvi, mvi_kw=mvi_kw, # Missing value imputation
                      norm=norm, norm_kw=norm_kw, # Normalization
                      tf=tf, tf_kw=tf_kw, # Transformation
                      scaling=scaling, scaling_kw=scaling_kw) # Scaling

    return a, b, c, d, e


# Functions to convert RGB colours to hex code and vice-versa
# From Stack Overflow (https://stackoverflow.com/questions/29643352/converting-hex-to-rgb-value-in-python)
def RGB(col): return '#%02x%02x%02x' % (int(col[0]), int(col[1]), int(col[2]))

def hex_to_rgb(value):
    value = value.lstrip('#')
    size = len(value)
    return tuple(int(value[i:i + size // 3], 16) for i in range(0, size, size // 3))


# Function to fill the target and colours of the above class
def TargetStorage_filling(target_list, colours):
    "Create target and assign default colours to the classes."

    temp_dict = {}
    target = target_list
    classes = pd.unique(target)
    for cl in range(len(classes)):

        # The first 10 different targets will follow the tab10 default colours
        if cl < len(colours):
            temp_dict[classes[cl]] = RGB(np.array(colours[cl])*255)

        # Generate random colours after the first 10
        else:
            # From https://www.geeksforgeeks.org/create-random-hex-color-code-using-python/
            # Generating a random number in between 0 and 2^24
            color = random.randrange(0, 2**24)
            # Converting that number from base-10 (decimal) to base-16 (hexadecimal)
            hex_color = hex(color)
            color_in_hex = "#" + hex_color[2:]
            temp_dict[classes[cl]] = color_in_hex
    return temp_dict



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