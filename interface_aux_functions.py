
# Needed Imports
import itertools
import pandas as pd
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hier
import seaborn as sns
import random
from scipy.stats import norm, chi2
import sklearn.ensemble as skensemble
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import panel as pn
import venn
import upsetplot

import holoviews as hv
import plotly.express as px
import plotly.graph_objects as go

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

hourglass_icon = '''<svg xmlns="http://www.w3.org/2000/svg" class="icon icon-tabler icon-tabler-hourglass-high" width="24" height="24" viewBox="0 0 24 24" stroke-width="2" stroke="currentColor" fill="none" stroke-linecap="round" stroke-linejoin="round">
    <path stroke="none" d="M0 0h24v24H0z" fill="none"/><path d="M6.5 7h11" />
    <path d="M6 20v-2a6 6 0 1 1 12 0v2a1 1 0 0 1 -1 1h-10a1 1 0 0 1 -1 -1z" />
    <path d="M6 4v2a6 6 0 1 0 12 0v-2a1 1 0 0 0 -1 -1h-10a1 1 0 0 0 -1 1z" />
</svg>'''


# Function related to File data reading
def read_file(filename, target_in_file):
    "Function to read the file given."

    # Samples names frequently have 00000.
    def renamer(colname):
        # Util to optionally remove all those 00000 from sample names
        return ''.join(colname.split('00000'))
    target_file = {}

    # No File Inputted
    if filename == '':
        file = pd.DataFrame()

    # If it is a .csv file
    elif filename.endswith('.csv'):

        if target_in_file: # If you have the target in the file
            file = pd.read_csv(filename, header=[0,1])
            colnames = [renamer(i) for i in file.columns.get_level_values(0)]
            target_file = dict(zip(colnames, file.columns.get_level_values(1)))
            file.columns = colnames

        else: # If you do not have the target in the file
            file = pd.read_csv(filename)
            file.columns = [renamer(i) for i in file.columns]
            target_file = {}

    # If it is a .xlsx file
    elif filename.endswith('.xlsx') or filename.endswith('.xls'):

        if target_in_file: # If you have the target in the file
            file = pd.read_excel(filename, header=[0,1])
            colnames = [renamer(i) for i in file.columns.get_level_values(0)]
            target_file = dict(zip(colnames, file.columns.get_level_values(1)))
            file.columns = colnames

        else: # If you do not have the target in the file
            file = pd.read_excel(filename)
            file.columns = [renamer(i) for i in file.columns]
            target_file = {}

    else:
        pn.state.notifications.error('Provided file is not an Excel or a csv file.')

    # Treated the read file to put them as we want it - # Important for database match
    try:
        file.insert(1, 'Neutral Mass', file['Bucket label'].str.replace('Da', '').astype('float'))
    except:
        pn.state.notifications.warning('Neutral Mass could not be inferred from Bucket Label. No annotation can be performed.')

    file = file.set_index('Bucket label')
    # Replaces zeros with numpy nans. Essential for data processing
    file = file.replace({0:np.nan})

    return file, target_file


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


def performing_pretreatment(PreTreatment_Method, original_df, target, sample_cols):
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
    data = original_df

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


def build_common_exclusive_dfs_to_save(com_exc_compounds, target_list, checkbox_annotation, checkbox_formula,
                                      radiobox_neutral_mass, checkbox_others):
    "Builds the common to all and exclusive to each dataframe annotated compounds of all classes to save as excel."

    # Setting the excel files for common compounds
    # Start creating the DataFrame
    common_df = pd.DataFrame(index=com_exc_compounds.common_all_id.index)
    # In how many samples and what percentage of samples of the dataset each annotated compound appears
    common_df['Appear in Samples'] = com_exc_compounds.common_all_id.loc[:, target_list.sample_cols].notnull().sum(axis=1)
    common_df['% of Samples'] = (com_exc_compounds.common_all_id.loc[
        :, target_list.sample_cols].notnull().sum(axis=1)) / len(target_list.target) * 100

    # Annotation information
    for col in com_exc_compounds.common_all_id.columns:
        if col not in target_list.sample_cols:
            # Annotations made in the dataset previous to using this software
            if col in checkbox_annotation.value:
                common_df[f'Prev. Annotation Match - {col}'] = com_exc_compounds.common_all_id[col]
            elif col in checkbox_formula.value:
                common_df[f'Prev. Annotation Formula - {col}'] = com_exc_compounds.common_all_id[col]

            # Annotations made with this software
            elif col not in ['Has Match?', radiobox_neutral_mass.value,] + checkbox_others.value:
                common_df[col] = com_exc_compounds.common_all_id[col]

    #common_df.index = com_exc_compounds.common_all_id['Neutral Mass']

    # Setting the excel files for exclusive compounds
    exclusive_dfs = {}
    for i in com_exc_compounds.exclusives_id.keys(): # For each class
        # Start creating the DataFrame
        current_exc_df = com_exc_compounds.exclusives_id[i]
        df_temp = pd.DataFrame(index=current_exc_df.index)
        # In how many samples and what percentage of samples of the exclusive class each annotated compound appears
        df_temp['Appear in Class Samples'] = current_exc_df.loc[:, target_list.sample_cols].notnull().sum(axis=1)
        df_temp['% of Class Samples'] = (
            current_exc_df.loc[:, target_list.sample_cols].notnull().sum(axis=1)) / target_list.target.count(i) * 100

        # Annotation information
        for col in current_exc_df.columns:
            if col not in target_list.sample_cols:
                # Annotations made in the dataset previous to using this software
                if col in checkbox_annotation.value:
                    df_temp[f'Prev. Annotation Match - {col}'] = current_exc_df[col]
                elif col in checkbox_formula.value:
                    df_temp[f'Prev. Annotation Formula - {col}'] = current_exc_df[col]

                # Annotations made with this software
                elif col not in ['Has Match?', radiobox_neutral_mass.value,] + checkbox_others.value:
                    df_temp[col] = current_exc_df[col]

        #df_temp.index = current_exc_df[radiobox_neutral_mass.value]
        exclusive_dfs[i] = df_temp # Storing the data

    return common_df, exclusive_dfs


def common_exclusive_compound_excel_writer(common_df, exclusive_dfs):
    "Writes to an excel the common and exclusive dataframes made with `build_common_exclusive_dfs_to_save`."

    # Initiate Excel File
    writer = pd.ExcelWriter('Common_Exclusive_Compounds.xlsx', engine='xlsxwriter')

    # For the common compounds DataFrame
    common_df.to_excel(writer, sheet_name='Common')

    # Format the columns based on the type of columns
    text_format = writer.book.add_format({'text_wrap' : True, 'valign': 'top'})
    for i in range(1, len(common_df.columns)+1):
        width=18
        if i in [1,2]:
            width=8
        elif common_df.columns[i-1].endswith('IDs'):
            width=15
        elif common_df.columns[i-1].endswith('count'):
            width=8
        elif common_df.columns[i-1].endswith('names') or common_df.columns[i-1].startswith('Prev. Annotation Match'):
            width=40
        writer.sheets['Common'].set_column(i,i,width,text_format)

    # Header formatting
    header_format = writer.book.add_format({'bold': True, 'text_wrap': True, 'valign': 'top'})
    # Overwrite both the value and the format of each header cell
    for col_num, value in enumerate(common_df.columns.values):
        writer.sheets['Common'].write(0, col_num + 1, value, header_format)

    # Same process repeated for each class for exclusive compounds DataFrame
    for a in exclusive_dfs.keys():
        exclusive_dfs[a].to_excel(writer, sheet_name=a+' Exclusive')

        # Format the columns based on the type of columns
        text_format = writer.book.add_format({'text_wrap' : True, 'valign': 'top'})
        for i in range(1, len(exclusive_dfs[a].columns)+1):
            width=18
            if i in [1,2]:
                width=8
            elif exclusive_dfs[a].columns[i-1].endswith('IDs'):
                width=15
            elif exclusive_dfs[a].columns[i-1].endswith('count'):
                width=8
            elif exclusive_dfs[a].columns[i-1].endswith('names') or exclusive_dfs[a].columns[i-1].startswith('Prev. Annotation Match'):
                width=40
            writer.sheets[a+' Exclusive'].set_column(i,i,width,text_format)

        # Header formatting
        header_format = writer.book.add_format({'bold': True, 'text_wrap': True, 'valign': 'top'})
        # Overwrite both the value and the format of each header cell
        for col_num, value in enumerate(exclusive_dfs[a].columns.values):
            writer.sheets[a+' Exclusive'].write(0, col_num + 1, value, header_format)

    writer.close()


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


def _plot_PCA(PCA_params, target_list):
    "Function to plot 2- or 3-D PCA projections"

    # 2 dimensions
    if PCA_params.n_dimensions == '2 Components':
        PCx_var_explained = PCA_params.explained_variance[int(PCA_params.PCx.split(" ")[-1])-1]
        PCy_var_explained = PCA_params.explained_variance[int(PCA_params.PCy.split(" ")[-1])-1]

        # Plot PCA
        PCA_plot = px.scatter(
            PCA_params.pca_scores, x=PCA_params.PCx, y=PCA_params.PCy, color=PCA_params.pca_scores['Label'],
            color_discrete_map=target_list.color_classes,
            title=f'''Total Explained Variance: {(PCx_var_explained + PCy_var_explained)*100:.2f}%''',
            labels={PCA_params.PCx: PCA_params.PCx + f' ({PCx_var_explained*100:.2f}%)',
                   PCA_params.PCy: PCA_params.PCy + f' ({PCy_var_explained*100:.2f}%)'})
        PCA_plot.update_traces(marker={'size': PCA_params.dot_size})

        # Draw ellipses if ellipses wanted
        if PCA_params.ellipse_draw:
            ellipses_df = pd.DataFrame()
            for lbl in target_list.color_classes.keys():
                subset_points = PCA_params.pca_scores[PCA_params.pca_scores['Label']==lbl]
                subset_points = subset_points[[PCA_params.PCx, PCA_params.PCy]]
                temp = plot_confidence_ellipse(
                    subset_points, q=PCA_params.confidence, nstd=PCA_params.confidence_std).array()
                temp = pd.concat((pd.DataFrame(temp), pd.DataFrame([lbl,]*len(temp), columns=['label'])), axis=1)
                ellipses_df = pd.concat((ellipses_df, temp))
            ellipses = px.line(ellipses_df, x=0, y=1, color='label', color_discrete_map=target_list.color_classes)

            # Final Plot joining the 2
            final_PCA_plot = go.Figure(data=PCA_plot.data + ellipses.data)
            final_PCA_plot.update_layout(
                title=f'''Total Explained Variance: {(PCx_var_explained + PCy_var_explained)*100:.2f}%''')
            final_PCA_plot.update_xaxes(title=PCA_params.PCx + f' ({PCx_var_explained*100:.2f}%)')
            final_PCA_plot.update_yaxes(title=PCA_params.PCy + f' ({PCy_var_explained*100:.2f}%)')

        else:
            final_PCA_plot = PCA_plot

    # 3 dimensions
    else:
        PCx_var_explained = PCA_params.explained_variance[int(PCA_params.PCx.split(" ")[-1])-1]
        PCy_var_explained = PCA_params.explained_variance[int(PCA_params.PCy.split(" ")[-1])-1]
        PCz_var_explained = PCA_params.explained_variance[int(PCA_params.PCz.split(" ")[-1])-1]

        final_PCA_plot = px.scatter_3d(
            PCA_params.pca_scores, x=PCA_params.PCx, y=PCA_params.PCy, z=PCA_params.PCz,
            color=PCA_params.pca_scores['Label'],
            color_discrete_map=target_list.color_classes,
            title=f'''Total Explained Variance: {(PCx_var_explained + PCy_var_explained + PCz_var_explained)*100:.2f}%''',
            labels={PCA_params.PCx: PCA_params.PCx + f' ({PCx_var_explained*100:.2f}%)',
                   PCA_params.PCy: PCA_params.PCy + f' ({PCy_var_explained*100:.2f}%)',
                   PCA_params.PCz: PCA_params.PCz + f' ({PCz_var_explained*100:.2f}%)'})
        final_PCA_plot.update_traces(marker={'size': PCA_params.dot_size})

    return final_PCA_plot

def _plot_PCA_explained_variance(PCA_params):
    "Function to plot cumulative explained variance"

    exp_var_cumul = np.cumsum(PCA_params.explained_variance)*100
    exp_var_df = pd.DataFrame((range(1, exp_var_cumul.shape[0] + 1), exp_var_cumul, PCA_params.explained_variance*100),
                             index=["Number of Components", "Cumulative Explained Variance (%)",
                                   "Explained Variance of This Component (%)"]).T
    exp_var_fig = px.area(exp_var_df,
        x="Number of Components",
        y="Cumulative Explained Variance (%)",
                          custom_data=["Explained Variance of This Component (%)"]
    )
    # Make the hover list contained the variance explained by the corresponding component
    exp_var_fig.update_traces(
        hovertemplate="<br>".join([
            "Number of Components: %{x}",
            "Cumulative Explained Variance (%): %{y}",
            "Explained Variance of This Component (%): %{customdata[0]:.5f}"]))
    exp_var_fig.update_yaxes(range=(0,100))
    return exp_var_fig

def _scatter_PCA_plot(PCA_params, target_list):
    "Function plotting a full matrix of PCA projections of the 6 (maximum) first Principal Components."

    columns = [PCA_params.pca_scores.columns[i] for i in range(6) if i < PCA_params.n_components]

    scatter_PCA_plot = px.scatter_matrix(
        PCA_params.pca_scores,
        dimensions=columns,
        color="Label",
        color_discrete_map=target_list.color_classes,
        labels={PCA_params.pca_scores.columns[i]: PCA_params.pca_scores.columns[i] +f" ({var:.2f}%)"
            for i, var in enumerate(PCA_params.explained_variance * 100)},
        height=800
    )

    scatter_PCA_plot.update_traces(diagonal_visible=False)
    return scatter_PCA_plot



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



### Functions related to the supervised analysis page of the graphical interface

## Functions for PLS-DA and RF section

def creating_importance_feat_table(imp_feat_metric, DataFrame_Store, model_feat_results):
    "Creates and organizes the DataFrame with the feature importance for supervised models."

    imp_feat_colname = imp_feat_metric + ' Score' # Name of Importance Column
    imp_feats = DataFrame_Store.metadata_df.copy() # Starting df
    imp_feats.insert(0, DataFrame_Store.metadata_df.index.name, imp_feats.index) # Include Column
    imp_feats.insert(1, imp_feat_colname, '') # Include Column for importance
    # Fill importance column
    for n in range(len(model_feat_results)):
        imp_feats[imp_feat_colname][model_feat_results[n][0]] = model_feat_results[n][1]
    # Sort from highest to lowest importance
    imp_feats = imp_feats.sort_values(by=imp_feat_colname, ascending=False)
    # Make Index be the place/position of each feature in the importance list
    imp_feats.index = range(1, len(imp_feats)+1)
    imp_feats.index.name = 'Place'

    return imp_feats


def _plot_permutation_test(perm_results, DataFrame_Store, n_fold, metric, title='Permutation Test'):
    "Plots the permutation test results with matplotlib."

    with plt.style.context('seaborn-whitegrid'):
        fig, ax = plt.subplots(1,1, figsize=(6,6))

        n_labels = len(DataFrame_Store.treated_df.index)
        tab20bcols = sns.color_palette('tab20b', 20)

        # Histogram with performance of permutated values
        hist_res = ax.hist(np.array(perm_results[1]), n_labels, range=(0, 1.00001),
                     edgecolor='black', color=tab20bcols[1], alpha = 1)

        # Plot the non-permutated model performance
        ylim = [0, hist_res[0].max()*1.2]
        ax.plot(2 * [perm_results[0]], ylim, '-', linewidth=3, color='darkred', #alpha = 0.5,
                     label='_p_-value %.5f)' % perm_results[2], solid_capstyle='round')
        ax.tick_params(labelsize=13)
        ax.set_xlabel(f'{n_fold}-fold Stratified CV Model {metric}', fontsize=14)
        ax.set_ylabel('Nº of occurrences', fontsize=14)
        if perm_results[0] >= 0.5:
            ax.text(perm_results[0]-0.45, hist_res[0].max()*1.1, 'p-value = %.3f' % perm_results[2], fontsize = 15)
        else:
            ax.text(perm_results[0]+0.05, hist_res[0].max()*1.1, 'p-value = %.3f' % perm_results[2], fontsize = 15)
        ax.set_title(title, size = 15)
        ax.set_axisbelow(True)

    return fig


## Functions for PLS-DA section

def _plot_PLS_optimization_components_fig(pls_scores):
    "Plot the figure with the optimization of the number of components."

    # Set up the initial line figure
    fig = px.line(pls_scores, x="Nº of Components", y='Scores', color='Class', title='PLS Optimization Plot',
             range_y=(0, 1.05), range_x=(0, max(list(pls_scores.index)) + 1))

    # Single out the q2 values and add annotation and marker to the maximum Q2 value found
    q2_val = pls_scores[pls_scores['Class'] == 'Q2']['Scores']
    rec_comp = q2_val.idxmax()
    fig.add_annotation(x=rec_comp,
                       y=q2_val.max(), text='Max.')
    fig.add_trace(go.Scatter(x=[rec_comp], y=[q2_val.max()],
                             marker=dict(color='Black', size=8), hoverinfo='skip', showlegend=False))

    return rec_comp, fig


def _plot_PLS(PLSDA_store, target_list):
    "Function to plot 2- or 3-D PLS projections. Variation of the `_plot_PCA` function."

    # 2 dimensions
    if PLSDA_store.n_dimensions == '2 Components':
        # Plot PLS
        PLS_plot = px.scatter(
            PLSDA_store.x_scores, x=PLSDA_store.LVx, y=PLSDA_store.LVy, color=PLSDA_store.x_scores['Label'],
            color_discrete_map=target_list.color_classes, title=f'''PLS Projection''')
        PLS_plot.update_traces(marker={'size': PLSDA_store.dot_size})

        # Draw ellipses if ellipses wanted
        if PLSDA_store.ellipse_draw:
            ellipses_df = pd.DataFrame()
            for lbl in target_list.color_classes.keys():
                subset_points = PLSDA_store.x_scores[PLSDA_store.x_scores['Label']==lbl]
                subset_points = subset_points[[PLSDA_store.LVx, PLSDA_store.LVy]]
                temp = plot_confidence_ellipse(
                    subset_points, q=PLSDA_store.confidence, nstd=PLSDA_store.confidence_std).array()
                temp = pd.concat((pd.DataFrame(temp), pd.DataFrame([lbl,]*len(temp), columns=['label'])), axis=1)
                ellipses_df = pd.concat((ellipses_df, temp))
            ellipses = px.line(ellipses_df, x=0, y=1, color='label', color_discrete_map=target_list.color_classes)

            # Final Plot joining the 2
            final_PLS_plot = go.Figure(data=PLS_plot.data + ellipses.data)
            final_PLS_plot.update_layout(
                title=f'''PLS Projection''')
            final_PLS_plot.update_xaxes(title=PLSDA_store.LVx)
            final_PLS_plot.update_yaxes(title=PLSDA_store.LVy)

        else:
            final_PLS_plot = PLS_plot

    # 3 dimensions
    else:
        final_PLS_plot = px.scatter_3d(
            PLSDA_store.x_scores, x=PLSDA_store.LVx, y=PLSDA_store.LVy, z=PLSDA_store.LVz,
            color=PLSDA_store.x_scores['Label'], color_discrete_map=target_list.color_classes,
            title=f'''PLS Projection''')
        final_PLS_plot.update_traces(marker={'size': PLSDA_store.dot_size})

    return final_PLS_plot


## Functions for RF section

def _optimization_n_trees_rf(RF_store, data, target):
    "Performs optimization of number of trees for Random Forest."

    # Vector with values for the parameter n_estimators
    values = {'n_estimators': range(RF_store.n_min_max_trees[0], RF_store.n_min_max_trees[1]+1, RF_store.n_interval)}

    # Setting up the RF model and Grid Search
    rf = skensemble.RandomForestClassifier(n_estimators=200)
    clf = GridSearchCV(rf, values, cv=RF_store.n_fold) # Change cv to change cross-validation

    # Fitting RF models
    clf.fit(data, target)

    return clf.cv_results_


### Functions related to the Univariate analysis page of the graphical interface

def _perform_univariate_analysis(UnivarA_Store, DataFrame_Store, target_list, filt_method, filt_kw):
    "Performs Univariate Analysis."

    # See if a T-Test or Mann-Whitney test will be made
    if UnivarA_Store.univariate_test == 'T-Test (Non-Parametric)':
        use_MW = False
    else:
        use_MW = True

    # Transform fold change (FC) threshold into logarithm
    abs_log2FC_threshold = np.log2(UnivarA_Store.fold_change_threshold)

    # If you have 2 classes
    if len(target_list.classes) == 2:
        univariate_df = DataFrame_Store.univariate_df

        # Perform Univariate Analysis
        univariate_results = metsta.compute_FC_pvalues_2groups(normalized=univariate_df, # Used for Fold-Change Computation
                                      processed=DataFrame_Store.treated_df, # Used for p-value computation
                                      labels=target_list.target, # Labels of the samples
                                      control_class=UnivarA_Store.control_class, # Control class
                                      test_class=UnivarA_Store.test_class, # Non-control class
                                      equal_var=UnivarA_Store.expected_equal_variance, # Consider variance between groups as equal
                                      useMW=use_MW) # Use Mann-Whitney Test or standard T-test

        # Select only Features considered significant (below a certain p-value threshold)
        filt_uni_results = univariate_results[univariate_results['FDR adjusted p-value'] < UnivarA_Store.p_value_threshold].copy()

        # Select features that have an absolute fold change (in log2) greater than abs_log2FC_threshold
        # Calculate absolute Log2 Fold-Change
        filt_uni_results['abs_log2FC'] = abs(filt_uni_results['log2FC'])
        # Select
        filt_uni_results = filt_uni_results[filt_uni_results['abs_log2FC'] > abs_log2FC_threshold]
        filt_uni_results = filt_uni_results.drop(columns='abs_log2FC')

        # Returns results (for 2 classes, there is no comprehensive dictionary of results)
        return univariate_df, {}, filt_uni_results, univariate_results, {}


    # More than 2 classes
    else:
        # To store results
        test_classes = [cl for cl in target_list.classes if cl != UnivarA_Store.control_class]
        univariate_results = {}
        filt_uni_results = {}
        univariate_df = {}

        for test_class in test_classes: # For each non-control class
            # Select only the samples of the control and current test class
            selection = [i in [UnivarA_Store.control_class, test_class] for i in target_list.target]
            target_temp = list(np.array(target_list.target)[selection])
            # Select only samples in the control and test classes
            file_temp = DataFrame_Store.original_df[target_list.sample_cols].copy()
            file_temp = file_temp.loc[:, selection]

            # Perform the same filtering and pre-treatments steps but using only the control and test class samples
            if filt_method.value == 'Total Samples':
                f_meth = 'total_samples'
            elif filt_method.value == 'Class Samples':
                f_meth = 'class_samples'
            filt_df, _ = initial_filtering(file_temp, file_temp.columns, target=target_list.target,
                                               filt_method=f_meth, filt_kw=filt_kw.value)
            t_data,_,filt_data,_,_  = performing_pretreatment(UnivarA_Store, filt_df,
                                                                 target_temp, file_temp.columns)

            univariate_df[test_class] = t_data

            # Perform Univariate Analysis on this newly acquired data
            univariate_results[test_class] = metsta.compute_FC_pvalues_2groups(
                                      normalized=filt_data, # Used for Fold-Change Computation
                                      processed=t_data, # Used for p-value computation
                                      labels=target_temp, # Labels of the samples
                                      control_class=UnivarA_Store.control_class, # Control class
                                      test_class=test_class, # Non-control class
                                      equal_var=UnivarA_Store.expected_equal_variance, # Consider variance between groups as equal
                                      useMW=use_MW) # Use Mann-Whitney Test if True or standard T-test if False

            # Select only Features considered significant
            filt_uni_results[test_class] = univariate_results[test_class][univariate_results[test_class][
                    'FDR adjusted p-value'] < UnivarA_Store.p_value_threshold].copy()

            # Select features that have an absolute fold change (in log2) greater than abs_log2FC_threshold
            # Calculate absolute Log2 Fold-Change
            filt_uni_results[test_class]['abs_log2FC'] = abs(filt_uni_results[test_class]['log2FC'])
            # Select
            filt_uni_results[test_class] = filt_uni_results[test_class][
                filt_uni_results[test_class]['abs_log2FC'] > abs_log2FC_threshold]
            filt_uni_results[test_class] = filt_uni_results[test_class].drop(columns='abs_log2FC')

        # Returns results
        return (univariate_df[UnivarA_Store.test_class], univariate_df,
                filt_uni_results[UnivarA_Store.test_class], univariate_results[UnivarA_Store.test_class], filt_uni_results)


def _plot_Volcano_plot(results_df, UnivarA_Store):
    "Plots Volcano Plots with plotly express."

    # Main Volcano Plot
    fig = px.scatter(results_df, x=results_df.log2FC, y=results_df['-log10(Adj. p-value)'],
          color=results_df.Expression, opacity=0.6,
          color_discrete_map={'Non-Significant': UnivarA_Store.color_non_sig,
                              'Downregulated': UnivarA_Store.color_down_sig,
                              'Upregulated': UnivarA_Store.color_up_sig},
          width=600, height=600,
          labels={'log2FC': f'log2 (Fold Change))',
               '-log10(Adj. p-value)': '- log10 (Adjusted (Benjamini-Hochberg) p-value)'})

    # Set Threshold lines
    fig.add_vline(x=np.log2(UnivarA_Store.fold_change_threshold), line_width=2, line_dash="dash", line_color="black")
    fig.add_vline(x=-np.log2(UnivarA_Store.fold_change_threshold), line_width=2, line_dash="dash", line_color="black")
    fig.add_hline(y=-np.log10(UnivarA_Store.p_value_threshold), line_width=2, line_dash="dash", line_color="black")

    fig.update_layout(title_text=f'Volcano Plot - {UnivarA_Store.test_class}/{UnivarA_Store.control_class}', title_x=0.45)

    return fig


def _univariate_intersections(UnivarA_Store, DataFrame_Store):
    "Filter metadata DF to include only features of the control class significant to all test classes chosen."

    idxs = DataFrame_Store.metadata_df.index # Start with all metabolites
    for cl in UnivarA_Store.test_class_subset: # For each test class chosen
        idxs = np.intersect1d(idxs, UnivarA_Store.univariate_results_set[cl].index) # Keep only the common features

    UnivarA_Store.specific_cl_df = DataFrame_Store.metadata_df.loc[idxs]
    n_inter_sig_feats = len(idxs) # Nº of features that are significant versus all selected test classes
    sig_annots_df = UnivarA_Store.specific_cl_df[UnivarA_Store.specific_cl_df['Has Match?']]
    n_inter_sig_annots = len(sig_annots_df) # Nº of annotated features that are significant versus all selected test classes
    # Select DataFrame with only annotated metabolites if this was chosen
    if UnivarA_Store.show_annots_only:
        UnivarA_Store.specific_cl_df = sig_annots_df

    # Building the string description to put in the page
    string = f'**{n_inter_sig_feats}** metabolites are significant against all the chosen test classes ('
    for cl in UnivarA_Store.test_class_subset:
        string = string+cl+', '
    string = string[:-2] + f'), **{n_inter_sig_annots}** of which are annotated.'
    UnivarA_Store.inter_description = string