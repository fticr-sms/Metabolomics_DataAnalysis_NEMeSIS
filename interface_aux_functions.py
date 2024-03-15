
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
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import panel as pn
import math
import upsetplot
import venn

import holoviews as hv
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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

instruction_icon = '''<svg xmlns="http://www.w3.org/2000/svg" class="icon icon-tabler icon-tabler-info-octagon" width="24" height="24" viewBox="0 0 24 24" stroke-width="2" stroke="currentColor" fill="none" stroke-linecap="round" stroke-linejoin="round">
    <path stroke="none" d="M0 0h24v24H0z" fill="none"/>
    <path d="M12.802 2.165l5.575 2.389c.48 .206 .863 .589 1.07 1.07l2.388 5.574c.22 .512 .22 1.092 0 1.604l-2.389 5.575c-.206 .48 -.589 .863 -1.07 1.07l-5.574 2.388c-.512 .22 -1.092 .22 -1.604 0l-5.575 -2.389a2.036 2.036 0 0 1 -1.07 -1.07l-2.388 -5.574a2.036 2.036 0 0 1 0 -1.604l2.389 -5.575c.206 -.48 .589 -.863 1.07 -1.07l5.574 -2.388a2.036 2.036 0 0 1 1.604 0z" />
    <path d="M12 9h.01" />
    <path d="M11 12h1v4h1" />
</svg>'''

img_confirm_button = '''<svg xmlns="http://www.w3.org/2000/svg"
    class="icon icon-tabler icon-tabler-check" width="24" height="24" viewBox="0 0 24 24" stroke-width="2"
    stroke="currentColor" fill="none" stroke-linecap="round" stroke-linejoin="round">
    <path stroke="none" d="M0 0h24v24H0z" fill="none"></path>
    <path d="M5 12l5 5l10 -10"></path>
</svg>'''

upload_icon = '''<svg xmlns="http://www.w3.org/2000/svg" class="icon icon-tabler icon-tabler-upload" width="24" height="24" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" fill="none" stroke-linecap="round" stroke-linejoin="round">
   <path stroke="none" d="M0 0h24v24H0z" fill="none"/>
   <path d="M4 17v2a2 2 0 0 0 2 2h12a2 2 0 0 0 2 -2v-2" />
   <path d="M7 9l5 -5l5 5" />
   <path d="M12 4l0 12" />
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
def read_file(filename, target_in_file, type_of_mass_values):
    "Function to read the file given."

    # Samples names frequently have 00000.
    def renamer(colname):
        # Util to optionally remove all those 00000 from sample names
        return ''.join(colname.split('00000'))
    target_file = {}
    nm_column = False

    # No File Inputted
    if filename == '':
        file = pd.DataFrame()

    # If it is a .csv file
    elif filename.endswith('.csv'):

        if target_in_file: # If you have the target in the file
            file = pd.read_csv(filename, header=[0,1])
            colnames = [renamer(i) for i in file.columns.get_level_values(0)]
            target_file = dict(zip(colnames, file.columns.get_level_values(1).astype(str)))
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
            target_file = dict(zip(colnames, file.columns.get_level_values(1).astype(str)))
            file.columns = colnames

        else: # If you do not have the target in the file
            file = pd.read_excel(filename)
            file.columns = [renamer(i) for i in file.columns]
            target_file = {}

    else:
        pn.state.notifications.error('Provided file is not an Excel or a csv file.')

    # Treated the read file to put them as we want it - # Important for database match
    try:
        # If the masses in the index are Neutral
        file.insert(1, 'Neutral Mass', file[file.columns[0]].astype('str').str.replace('Da', '').astype('float'))
        # If the masses are m/z values obtained in Positive Ionization Mode
        if type_of_mass_values == 'm/z (Positive)':
            electron_mass = 0.0005485799090654
            H_mass = 1.007825031898
            file['Neutral Mass'] = file['Neutral Mass'] - H_mass + electron_mass
        # If the masses are m/z values obtained in Negative Ionization Mode
        elif type_of_mass_values == 'm/z (Negative)':
            electron_mass = 0.0005485799090654
            H_mass = 1.007825031898
            file['Neutral Mass'] = file['Neutral Mass'] + H_mass + electron_mass
        nm_column = True
    except:
        pn.state.notifications.warning('Neutral Mass could not be inferred from the 1st column of your data. No annotation can be performed.')

    file = file.set_index(file.columns[0])
    file.index.name = 'Bucket label'
    # Replaces zeros with numpy nans. Essential for data processing
    file = file.replace({0:np.nan})

    return file, target_file, nm_column



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


# TODO: Observe this function better and update and improve it for the graphical interface
def duplicate_disambiguator(data_ann_deduplicator, annotated_data, sample_cols, neutral_mass_col, mz_col=True):
    "Adapted function from metanalysis_standard.py"
    """Attempts to remove duplicate (or more) annotations of peaks by merging them when possible.

       1) See peaks that have the same metabolite annotation by other databases.
       2) See if the other compound annotations do not have different annotations for those peaks.
       3) If not, save the meta data of the compound and formula annotations by the different databases.
       4) Situation Trouble - If yes, then we may have a problem. If for example, HMDB puts two different compounds for the
        2 m/z peaks and LOTUS puts the same compound, it is fair to treat them as different peaks. HOWEVER, if there are
        more than two peaks assigned with the same formula, the following can happen. Let's imagine a scenario where HMDB
        puts the same compound for 4 m/z peaks and LOTUS assigns to one of them one compound, to a second one a different
        compound and the last two ones does not assign a compound. What is the correct course of action? Right now, it just
        does not merge any of these peaks, but we could merge the two peaks that do not have an annotation by LOTUS. Would
        that be correct? Or should we merge with one of the two other peaks which have annotations by LOTUS. After all,
        they would normally be merged if not for the existence of two different LOTUS annotations. Hence, the problem.
       5) Then create the new peak, by keeping the highest intensity value in each sample from the different peaks (our
        intensity values come from the maximum value in the peak and not peak area.
       6) Situation 1: If all the highest intensity values come from one m/z peak, then that peak becomes the 'de facto'
        peak and all others are erased.
       7) Situation 2: The highest intensity comes from at least two different m/z peaks and ALL peaks come from the same
        adduct (including ones that are not used for the merge). Then, the peak 'Bucket Label' and 'Neutral Mass' columns
        become the weighted average (based on the average intensity of the peaks) of all the peaks with the same annotation.
        If there is no m/z column, this is the situation used. If there is no m/z column in the data, this situation is
        used instead of Situation 3.
       8) Situation 3: Identical to Situation 2 but there is at least one peak that comes from a different adduct based on
        m/z column. Then, the peak 'Bucket Label' and 'Neutral Mass' columns become identical to the peak which has the
        highest average intensity of all the peaks with the same annotation.
       9) This process is repeated for MetaboScape annotations and all databases first. And is then used for Smart Formula.
        Usually the number of de-duplications made by each database should decrease since when you de-duplciate duplicate
        assignments by one database, you are usually de-duplicating in others.

       The problem mentioned should be rare. However, check the merge_problems variable. If it is NOT empty, then those
        merge problems issues might exist. We do not currently have an automatic answer for them. They should be seen on a
        case by case basis.

       returns: annotated_data (pandas DataFrame with data after merging),
                mergings_performed (dictionary with summary of mergings made),
                merging_situations (dictionary with counts of times each situation of merging was used),
                merge_description (dictionary with descriptions of mergings made),
                merge_problems (dictionary with descriptions of possible problems to merge)
    """

    #
    mcid = data_ann_deduplicator.mcid

    # To store results
    merge_problems = {}
    merging_situations = {'Overwrite': 0, 'Merge same adducts (no m/z col)':0,
                          'Merge different adducts (Imp.)':0, 'Possible Problem Cases': 0}
    mergings_performed = dict(zip(data_ann_deduplicator.mcid_alt.values(), [0,]*len(data_ann_deduplicator.mcid_alt.values())))
    merge_description = {}

    # For the different databases used
    for col in mcid:
        new_idxs = {} # To save the idx where the merged peaks will be
        lost_idxs = [] # To save the idxs which will have to be removed

        # See repeating annotations
        repeating_names = annotated_data[col].value_counts()[annotated_data[col].value_counts()>1].index

        for name in repeating_names: # For each repeating annotation
            # Grab the part of the dataframe with the repeats
            subset_df = annotated_data.loc[[i for i in annotated_data[col].index if annotated_data.loc[i, col] == name]]

            # See if merging is possible and storing metadata
            merge=True
            saving_annotations = {} # To store the annotations to keep in the merged line
            for col_alt in mcid: # All other databases
                if col_alt != col:
                    # See if there are annotations in the other databases
                    subset_notnull = subset_df[col_alt][subset_df[col_alt].notnull()]
                    if len(subset_notnull) < 2:
                        if len(subset_notnull) == 1: # If there is only one, save it
                            saving_annotations[col_alt] = subset_notnull.value_counts().index[0]
                        continue
                        # If there is not, nothing to save
                    else:
                        n_annotations = len(subset_notnull.value_counts())
                        # If there is more than one but they are all the same, save it
                        if n_annotations == 1:
                            saving_annotations[col_alt] = subset_notnull.value_counts().index[0]
                            #print(subset_notnull.value_counts().index)
                            continue
                        # If there are multiple different annotations, PROBLEM
                        # Merging will not happen
                        else:
                            #if len(subset_df[col_alt]) > n_annotations:
                            merge_problems[
                                merging_situations['Possible Problem Cases']] = {'DB': data_ann_deduplicator.mcid_alt[col],
                                                                              'Annotation': name,
                                                                              'Indexes': list(subset_df.index),
                                                                              'Nº of peaks':len(subset_df),
                                                                              'Reason': data_ann_deduplicator.mcid_alt[col_alt]}
                            merging_situations['Possible Problem Cases'] = merging_situations['Possible Problem Cases'] + 1
                            merge=False

            # Starting the merging process if possible
            if merge:
                n_masses = subset_df[neutral_mass_col] # Get the neutral mass if it exists

                # Get the idxs positions of the repeated annotations in the dataframe
                idxs = [annotated_data.index.get_loc(subset_df.index[i]) for i in range(len(subset_df))]

                min_idx = min(idxs) # Grab the minimum that will become the merged peak

                idxs.remove(min_idx)
                lost_idxs.extend(idxs) # Add the other to lost_idxs that will be removed when all this is over

                overwrite = False #

                new_line = subset_df[sample_cols].max() # Get the intensity values for the new merged line
                # For each annotation, see if the highest intensity values all come from one line
                # Probably a better way to do this
                for i in range(len(subset_df)):
                    if new_line.notnull().sum() - (new_line == subset_df[sample_cols].iloc[i]).sum() == 0:
                        # If yes, then situation 1: Overwrite is the correct option here
                        overwrite = True

                        # The new id (bucket label) will be the same as in the line with all the higher intensities
                        # Store it
                        keep_id = subset_df.index[i]
                        new_idxs[annotated_data.index[min_idx]] = keep_id

                        # The new line will be identical to that line in terms of intensity
                        temp_full_new_line = subset_df.iloc[i].copy()

                        # Filling the meta data of the new line with the data stored in saving annotations
                        for key in saving_annotations:
                            if key.startswith('Matched') and key.endswith(' IDs'):
                                # Grab the rest of the columns with meta_data: IDs, names, formulas and match count
                                temp = annotated_data.loc[[
                                    i for i in annotated_data.index if annotated_data.loc[i,
                                                 key] == saving_annotations[key]]]
                                rel_cols = [key, 'Matched '+ key[8:-4] +' names', 'Matched '+ key[8:-4] +' formulas',
                                                    key[8:-4] +' match count']
                                temp_full_new_line[rel_cols] = temp.iloc[0][rel_cols]

                            else:
                                temp_full_new_line.loc[key] = saving_annotations[key]

                        merging_situations['Overwrite'] = merging_situations['Overwrite'] + 1
                        mergings_performed[data_ann_deduplicator.mcid_alt[col]] = mergings_performed[
                            data_ann_deduplicator.mcid_alt[col]] + 1
                        merge_description[keep_id] = {'DB': data_ann_deduplicator.mcid_alt[col],
                                                      'Repeating annotation': name,
                                                      'Nº merged peaks': len(subset_df),
                                                      'Merged peaks': list(subset_df.index),
                                                      'Situation': 'Overwrite'}

                        # Putting the merged line in the DataFrame
                        annotated_data.loc[annotated_data.index[min_idx]] = temp_full_new_line.copy()

                        continue

                # If Situation 1, we end here
                if overwrite:
                    continue

                # If not Situation 1
                temp_full_new_line = annotated_data.iloc[min_idx].copy()
                # The new line will have the higher intensities for each sample
                temp_full_new_line.loc[sample_cols] = new_line

                # Filling the meta data of the new line with the data stored in saving annotations
                for key in saving_annotations:
                    if key.startswith('Matched') and key.endswith(' IDs'):
                        # Grab the rest of the columns with meta_data: IDs, names, formulas and match count
                        temp = annotated_data.loc[[
                            i for i in annotated_data.index if annotated_data.loc[i,
                                         key] == saving_annotations[key]]]
                        rel_cols = [key, 'Matched '+ key[8:-4] +' names', 'Matched '+ key[8:-4] +' formulas',
                                            key[8:-4] +' match count']
                        temp_full_new_line[rel_cols] = temp.iloc[0][rel_cols]

                    else:
                        temp_full_new_line.loc[key] = saving_annotations[key]

                # All that's left is the Bucket Label and Neutral Mass
                if mz_col:
                    mz_col
                    # If an mz_col exists see if the distance between the maximum and minimum mass is low
                    # If yes, Situation 2: Bucket Label and Neutral Mass will be the weighted averages of all the
                    # possible peaks

                    # Removed since this does not happen in the interface

                else:
                    # Get the new bucket label
                    new_bucket = str(np.average(n_masses, weights=subset_df[sample_cols].mean(axis=1))) + ' Da'
                    new_idxs[annotated_data.index[min_idx]] = new_bucket
                    # Get the Neutral Mass
                    temp_full_new_line[neutral_mass_col] = np.average(
                        n_masses, weights=subset_df.loc[:,sample_cols].mean(axis=1))

                    # Storing info
                    #situation = 'merging'
                    merging_situations['Merge same adducts (no m/z col)'] = merging_situations['Merge same adducts (no m/z col)'] + 1
                    mergings_performed[data_ann_deduplicator.mcid_alt[col]] = mergings_performed[
                        data_ann_deduplicator.mcid_alt[col]] + 1
                    merge_description[new_bucket] = {'DB': data_ann_deduplicator.mcid_alt[col],
                                                     'Repeating annotation': name,
                                                     'Nº merged peaks': len(subset_df),
                                                     'Merged peaks': list(subset_df.index),
                                                     'Situation': 'Merging same adducts (no m/z col)'}

                # Putting the merged line in the DataFrame
                annotated_data.loc[annotated_data.index[min_idx]] = temp_full_new_line.copy()

        # Removing lost idxs peaks
        named_idxs = []
        for idx in lost_idxs:
            named_idxs.append(annotated_data.index[idx])
        for idx in named_idxs:
            annotated_data = annotated_data.drop(index=idx)

        # Assigning the new bucket labels
        for old_idx, new_idx in new_idxs.items():
            annotated_data = annotated_data.rename(index= {old_idx : new_idx})

    return annotated_data, mergings_performed, merging_situations, merge_description, merge_problems


def _merge_problems_creation(mp):
    "Creates a DataFrame with the problems in merging from the dictionary passed."

    mp = pd.DataFrame(mp).T
    new_merge_problems = pd.DataFrame(columns=['DBs with same annotation', 'Annotation',
                                            'Nº of Peaks', 'Indexes', 'Contradiction With'])
    idx_counts = mp['Indexes'].value_counts(sort=False)
    a = 0
    for case in idx_counts.index:
        subset = mp.loc[[i for i in mp.index if mp.loc[i, 'Indexes'] == case]]
        ann_val = {subset['DB'].iloc[i]: subset['Annotation'].iloc[i] for i in range(len(subset))}
        new_merge_problems.loc[a] =[' | '.join(pd.unique(subset['DB'].values)), ann_val,
                                    subset['Nº of peaks'].iloc[0], subset['Indexes'].iloc[0],
                                    ', '.join(pd.unique(subset['Reason'].values))]
        a+=1

    new_merge_problems = new_merge_problems.set_index('Indexes')

    return new_merge_problems


def individually_merging(data_ann_deduplicator, given_idxs, sample_cols, annotated_cols, neutral_mass_col, mz_col=False):
    """Attempts to merge a set of peaks given in the DataFrame.

       returns: annotated_data (pandas DataFrame with data after merging),
                merge_description (dictionary with description of the merging)
    """

    # To store results
    merge_description = {}

    # Grab the part of the dataframe you want to merge
    annotated_data = data_ann_deduplicator.annotated_df.copy()
    subset_df = annotated_data.loc[given_idxs]

    # Split it into the different parts
    intensities = subset_df[sample_cols]
    metadata_ann = subset_df[annotated_cols]
    other_metadata = subset_df[neutral_mass_col]

    # See if merging is possible and get the annotation metadata that will be in the merged line
    meta_vals = []
    for col in metadata_ann.columns:
        n_entries = metadata_ann[col].value_counts().index
        if len(n_entries) > 1:
            pn.state.notifications.error(f'Merging of {given_idxs} was not possible: multiple incompatible values in {col}.')
            raise ValueError(f'Merging of {given_idxs} was not possible: multiple incompatible values in {col}.')
        # No annotation found
        if len(n_entries) == 0:
            meta_vals.append(np.nan)
        # One annotation found
        else:
            meta_vals.append(n_entries[0])

    # For the different databases used
    new_idxs = {} # To save the idx where the merged peaks will be
    lost_idxs = [] # To save the idxs which will have to be removed

    # Get the idxs positions of the repeated annotations in the dataframe
    idxs = [annotated_data.index.get_loc(subset_df.index[i]) for i in range(len(subset_df))]

    min_idx = min(idxs) # Grab the minimum that will become the merged peak

    idxs.remove(min_idx)
    lost_idxs.extend(idxs) # Add the other to lost_idxs that will be removed when all this is over

    overwrite = False #

    # Set up new line
    temp_new_line = annotated_data.loc[annotated_data.index[min_idx]].copy()
    # Metadata
    temp_new_line[annotated_cols] = meta_vals
    # Intensities
    new_line = intensities.max() # Get the intensity values for the new merged line
    temp_new_line[sample_cols] = new_line
    # Neutral Mass
    # For each annotation, see if the highest intensity values all come from one line
    # Probably a better way to do this
    for i in range(len(intensities)):
        if new_line.notnull().sum() - (new_line == intensities.iloc[i]).sum() == 0:
            # If yes, then situation 1: Overwrite is the correct option here
            overwrite = True

            # The new id (bucket label) will be the same as in the line with all the higher intensities
            # Store it
            keep_id = intensities.index[i]
            new_idxs[annotated_data.index[min_idx]] = keep_id

            # Neutral Mass Column
            temp_new_line[neutral_mass_col] = other_metadata.iloc[i]

            merge_description[keep_id] = {'Nº merged peaks': len(subset_df),
                                          'Merged peaks': list(subset_df.index),
                                          'Situation': 'Individual - Overwrite'}

            # Putting the merged line in the DataFrame
            annotated_data.loc[annotated_data.index[min_idx]] = temp_new_line.copy()

            continue

    # If not Situation 1
    if not overwrite:

        # If not Situation 1
        if mz_col:
            # If an mz_col existss see if the distance between the maximum and minimum mass is low
            # If yes, Situation 2: bucket label, Neutral Mass and m/z peak will be the weighted averages of all the
            # possible peaks
            mz_col
            # Not possible currently to happem

        else:
            # Get the new bucket label
            avg_neutral_mass = np.average(other_metadata, weights=intensities.mean(axis=1))
            new_bucket = str(avg_neutral_mass) + ' Da'
            new_idxs[annotated_data.index[min_idx]] = new_bucket
            # Get the Neutral Mass
            temp_new_line[neutral_mass_col] = avg_neutral_mass

            # Storing info
            merge_description[new_bucket] = {'Nº merged peaks': len(subset_df),
                                             'Merged peaks': list(subset_df.index),
                        'Situation': 'Individual - Merging same adducts (no m/z col)'}

        # Putting the merged line in the DataFrame
        annotated_data.loc[annotated_data.index[min_idx]] = temp_new_line.copy()

    # Removing lost idxs peaks
    named_idxs = []
    for idx in lost_idxs:
        named_idxs.append(annotated_data.index[idx])
    for idx in named_idxs:
        annotated_data = annotated_data.drop(index=idx)

    # Assigning the new bucket labels
    for old_idx, new_idx in new_idxs.items():
        annotated_data = annotated_data.rename(index= {old_idx : new_idx})

    return annotated_data, merge_description


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
    tf_translation = {'Generalized Logarithmic Transformation (glog)':'glog', 'None': None}
    tf = tf_translation[PreTreatment_Method.tf_method]

    tf_kw = PreTreatment_Method.tf_kw

    # Scaling
    scaling_translation = {'Pareto Scaling':'pareto', 'Auto Scaling':'auto', 'Mean Centering':'mean_center',
                      'Range Scaling':'range', 'Vast Scaling': 'vast', 'Level Scaling': 'level', 'None': None}
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
    classes = pd.unique(np.array(target))
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
    com_exc_compounds.groups = {}
    for cl in target_list.color_classes.keys(): # Setting up the keys (classes)
        com_exc_compounds.groups[cl] = []

    for c, t in zip(target_list.sample_cols, target_list.target): # Setting up the values
        for g in com_exc_compounds.groups:
            if g == t:
                com_exc_compounds.groups[g].append(c)

    # Creating the DataFrame for each class by dropping features of the processed complete DataFrame if they do not appear in any sample of the class
    # One DataFrame is made considering every feature in the dataset, and another considering only annotated features
    com_exc_compounds.group_dfs = {}
    com_exc_compounds.group_dfs_ids = {}
    for g in com_exc_compounds.groups:
        for c, t in zip(target_list.sample_cols, target_list.target):
            if g == t:
                com_exc_compounds.group_dfs[g] = DataFrame_Store.processed_df.dropna(
                    subset=com_exc_compounds.groups[g], thresh=1)
                com_exc_compounds.group_dfs_ids[g] = com_exc_compounds.group_dfs[g].iloc[[
                    i for i in range(len(com_exc_compounds.group_dfs[
                    g]['Has Match?'])) if com_exc_compounds.group_dfs[g]['Has Match?'].iloc[i]]]

    # Description of the number of metabolites (and annotated metabolites) that appear in at least one sample of each class
    desc_string = ['**Nº of metabolic features per class:**', '', ]
    for g in com_exc_compounds.groups:
        group_string = f'**{g}**: **{len(com_exc_compounds.group_dfs[g])}** metabolites, from which **{len(com_exc_compounds.group_dfs_ids[g])}** have matches (annotations).'
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
        group_string = f'**{g}**: **{len(com_exc_compounds.exclusives[g].index)}** exclusive metabolites (**{len(com_exc_compounds.exclusives_id[g].index)}** annotated).'
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
        x_label = 'Nº of Metabolic Features (Nº of Matched Compounds)'
    elif com_exc_compounds.type_of_venn == 'All Metabolites': # Only show the metabolite numbers
        labels_all = labels
        x_label = 'Nº of Metabolic Features '
    else: # Only show the annotated numbers
        labels_all = labels_ids
        x_label = 'Nº of Matched Compounds'

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


def _plot_intersection_plots(com_exc_compounds, groups_dict, ups):
    "Plot an Intersection plot"
    # Plotting Intersection Plot
    f,ax = plt.subplots(1,1, constrained_layout=True, dpi=400)
    if com_exc_compounds.inter_include_counts_percentages == 'Show Nº and % of metabolites':
        include_counts, include_percentages = True, True
    elif com_exc_compounds.inter_include_counts_percentages == 'Show Nº of metabolites':
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

    ellip = hv.Ellipse(pos.iloc[0], pos.iloc[1], (width,height), orientation=theta)

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

# Plots the HCA function
def _plot_HCA(HCA_params, target_list):
    "Plots HCA Dendrogram."
    f, ax = plt.subplots(1, 1, figsize=(HCA_params.fig_x, HCA_params.fig_y), constrained_layout=True)
    plot_dendogram(HCA_params.Z,
                   target_list.target, ax=ax,
                   label_colors=target_list.color_classes,
                   x_axis_len=HCA_params.fig_x,
                   color_threshold=HCA_params.col_threshold)
    return f


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
        imp_feats[imp_feat_colname].iloc[model_feat_results[n][0]] = model_feat_results[n][1]
    # Sort from highest to lowest importance
    imp_feats = imp_feats.sort_values(by=imp_feat_colname, ascending=False)
    # Make Index be the place/position of each feature in the importance list
    imp_feats.index = range(1, len(imp_feats)+1)
    imp_feats.index.name = 'Place'

    return imp_feats


def _plot_permutation_test(perm_results, df, n_fold, metric, title='Permutation Test'):
    "Plots the permutation test results with matplotlib."

    with plt.style.context('seaborn-v0_8-whitegrid'):
        fig, ax = plt.subplots(1,1, figsize=(6,6))

        n_labels = len(df.index)
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


def _plot_PLSDA_ROC_curve(PLSDA_store, treated_df, target_list):
    "Plots the ROC Curve of PLS-DA models."

    if len(target_list.classes) == 2: # When you only have 2 classes

        # Compute ROC Curve
        resROC_PLSDA = metsta.PLSDA_ROC_cv(treated_df, target_list.target, # Data and target
                        pos_label=PLSDA_store.positive_class, # Positive class chosen
                        n_comp=PLSDA_store.n_components, # Number of components
                        scale=PLSDA_store.scale, # Set scale to True only scaling was not performed
                        n_iter=PLSDA_store.roc_n_iter, # Number of iterations to repeat
                        cv=None, n_fold=PLSDA_store.n_fold) # Number of folds

        mean_fpr = [0,] + list(resROC_PLSDA['average fpr'])
        mean_tpr = [0,] + list(resROC_PLSDA['average tpr'])
        mean_auc = resROC_PLSDA['mean AUC']

        # Plotting the ROC Curve
        fig = px.area(
            x=mean_fpr, y=mean_tpr,
            title=f'PLS-DA ROC Curve, AUC = {mean_auc:.3f}',
            labels=dict(x='False Positive Rate', y='True Positive Rate'), range_x=(-0.01, 1.001),
        )

        # Dotted line representing worse-case scenario
        fig.add_shape(
            type='line', line=dict(dash='dash'), line_width=3,
            x0=0, x1=1, y0=0, y1=1
        )

        # Setting up filename
        filename = f'PLSDA_ROCcurve_{PLSDA_store.positive_class}PosClass_{PLSDA_store.n_components}components_'
        filename = filename + f'{PLSDA_store.n_fold}-foldStratCV_{PLSDA_store.roc_n_iter}iterations_scale{PLSDA_store.scale}'

    else: # When you have more than 2 classes, plot a ROC curve for each class in a "1vsAll" fitted models

        # Create storage for results
        total_ROC_res = {}
        # Then for each class, build a custom target where samples of those classes are labelled correctly
        # And the remaining samples are labelled as 'Other'
        for cl in target_list.classes:
            temp_tg = []
            for t in target_list.target:
                if t == cl:
                    temp_tg.append(cl)
                else:
                    temp_tg.append('Other')

            # Use this makeshift target to fit "1vsAll" model and obtain ROC results
            total_ROC_res[cl] =  metsta.PLSDA_ROC_cv(treated_df, temp_tg, # Data and target
                            pos_label=cl, # Positive class is the current class
                            n_comp=PLSDA_store.n_components, # Number of components
                            scale=PLSDA_store.scale, # Set scale to True only scaling was not performed
                            n_iter=PLSDA_store.roc_n_iter, # Number of iterations to repeat
                            cv=None, n_fold=PLSDA_store.n_fold) # Number of folds

        # After having the results plot them in a Figure
        # Initialize the figure
        fig = go.Figure()
        for cl in total_ROC_res: # For each class results

            # Extract the relevant data
            mean_fpr = [0,] + list(total_ROC_res[cl]['average fpr'])
            mean_tpr = [0,] + list(total_ROC_res[cl]['average tpr'])
            mean_auc = total_ROC_res[cl]['mean AUC']

            # And plot the correspondign ROC Curve with the correct label
            fig.add_trace(go.Scatter(name=f'{cl} (AUC={mean_auc:.3f})',
                         x=mean_fpr, y=mean_tpr))

        # After plotting all curves, plot the dotted line representing worse-case scenario
        fig.add_shape(
        type='line', line=dict(dash='dash'), line_width=3,
        x0=0, x1=1, y0=0, y1=1
        )

        # Update other characteristics of the plot
        fig.update_layout(legend_title_text='Classes',
                     title='PLS-DA ROC Curves',
                    xaxis_range=(-0.01, 1.001),
                    xaxis_title='False Positive Rate',
                    yaxis_title='True Positive Rate')

        # Setting up filename
        filename = f'PLSDA_ROCcurve_1vsAllModels_{PLSDA_store.n_components}components_{PLSDA_store.n_fold}-foldStratCV'
        filename = filename + f'_{PLSDA_store.roc_n_iter}iterations_scale{PLSDA_store.scale}'

    return fig, filename


## Functions for RF section

def _optimization_n_trees_rf(RF_store, data, target):
    "Performs optimization of number of trees for Random Forest."

    # Vector with values for the parameter n_estimators
    values = {'n_estimators': range(RF_store.n_min_max_trees[0], RF_store.n_min_max_trees[1]+1, RF_store.n_interval)}

    # Setting up the RF model and Grid Search
    rf = skensemble.RandomForestClassifier(n_estimators=200)
    clf = GridSearchCV(rf, values, cv=RF_store.n_fold, scoring='accuracy') # Change cv to change cross-validation

    # Fitting RF models
    clf.fit(data, target)

    return clf.cv_results_


def _plot_RF_ROC_curve(RF_store, treated_df, target_list):
    "Plots the ROC Curve of Random Forest models."

    if len(target_list.classes) == 2: # When you only have 2 classes

        # Compute ROC Curve
        resROC_RF = metsta.RF_ROC_cv(treated_df, target_list.target, # Data and target
                                    pos_label=RF_store.positive_class, # Positive class chosen
                                    n_trees=RF_store.n_trees, # Number of trees of RF
                                    n_iter=RF_store.roc_n_iter, # Number of iterations to repeat
                                    cv=None, n_fold=RF_store.n_fold) # Number of folds

        mean_fpr = [0,] + list(resROC_RF['average fpr'])
        mean_tpr = [0,] + list(resROC_RF['average tpr'])
        mean_auc = resROC_RF['mean AUC']

        # Plotting the ROC Curve
        fig = px.area(
            x=mean_fpr, y=mean_tpr,
            title=f'Random Forest ROC Curve, AUC = {mean_auc:.3f}',
            labels=dict(x='False Positive Rate', y='True Positive Rate'), range_x=(-0.01, 1.001),
        )

        # Dotted line representing worse-case scenario
        fig.add_shape(
            type='line', line=dict(dash='dash'), line_width=3,
            x0=0, x1=1, y0=0, y1=1
        )

        # Setting up filename
        filename = f'RF_ROCcurve_{RF_store.positive_class}PosClass_{RF_store.n_trees}trees_'
        filename = filename + f'{RF_store.n_fold}-foldStratCV_{RF_store.roc_n_iter}iterations'

    else: # When you have more than 2 classes, plot a ROC curve for each class in a "1vsAll" fitted models

        # Create storage for results
        total_ROC_res = {}
        # Then for each class, build a custom target where samples of those classes are labelled correctly
        # And the remaining samples are labelled as 'Other'
        for cl in target_list.classes:
            temp_tg = []
            for t in target_list.target:
                if t == cl:
                    temp_tg.append(cl)
                else:
                    temp_tg.append('Other')

            # Use this makeshift target to fit "1vsAll" model and obtain ROC results
            total_ROC_res[cl] =  metsta.RF_ROC_cv(treated_df, temp_tg, # Data and target
                                    pos_label=cl, # Positive class chosen
                                    n_trees=RF_store.n_trees, # Number of trees of RF
                                    n_iter=RF_store.roc_n_iter, # Number of iterations to repeat
                                    cv=None, n_fold=RF_store.n_fold) # Number of folds

        # After having the results plot them in a Figure
        # Initialize the figure
        fig = go.Figure()
        for cl in total_ROC_res: # For each class results

            # Extract the relevant data
            mean_fpr = [0,] + list(total_ROC_res[cl]['average fpr'])
            mean_tpr = [0,] + list(total_ROC_res[cl]['average tpr'])
            mean_auc = total_ROC_res[cl]['mean AUC']

            # And plot the correspondign ROC Curve with the correct label
            fig.add_trace(go.Scatter(name=f'{cl} (AUC={mean_auc:.3f})',
                         x=mean_fpr, y=mean_tpr))

        # After plotting all curves, plot the dotted line representing worse-case scenario
        fig.add_shape(
        type='line', line=dict(dash='dash'), line_width=3,
        x0=0, x1=1, y0=0, y1=1
        )

        # Update other characteristics of the plot
        fig.update_layout(legend_title_text='Classes',
                     title='Random Forest ROC Curves',
                    xaxis_range=(-0.01, 1.001),
                    xaxis_title='False Positive Rate',
                    yaxis_title='True Positive Rate')

        # Setting up filename
        filename = f'RF_ROCcurve_1vsAllModels_{RF_store.n_trees}trees_{RF_store.n_fold}-foldStratCV'
        filename = filename + f'_{RF_store.roc_n_iter}iterations'

    return fig, filename



### Functions related to the Univariate analysis page of the graphical interface

def _perform_univariate_analysis(UnivarA_Store, DataFrame_Store, target_list, filt_method, filt_kw):
    "Performs Univariate Analysis."

    # See if a T-Test or Mann-Whitney test will be made
    if UnivarA_Store.univariate_test == 'T-Test (Parametric)':
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
        return univariate_df, {UnivarA_Store.test_class:univariate_df}, filt_uni_results, univariate_results, {UnivarA_Store.test_class:filt_uni_results}


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
            if UnivarA_Store.filt_method == 'Total Samples':
                f_meth = 'total_samples'
                # Adapting the filt_kw to a smaller subset of samples
                # Use percentage of the original filtering used to calculate the equivalent number of samples in subset and round UP
                # Possible Issue - since we already used the filtered dataset (because it has annotations and de-duplications),
                # the data filtering with 'total_samples' is not perfect - since a feature must pass this data filtering but also
                # the original data filtering made
                f_kw = math.ceil(UnivarA_Store.filt_kw/len(target_list.sample_cols)*sum(selection))
            elif UnivarA_Store.filt_method == 'Class Samples':
                f_meth = 'class_samples'
                f_kw = UnivarA_Store.filt_kw
            else:
                f_meth = None

            # Perform the Data Filtering and the Pre-Treament
            filt_df, _ = initial_filtering(file_temp, file_temp.columns, target=target_temp,
                                       filt_method=f_meth, filt_kw=f_kw)
            t_data,_,filt_data,_,_  = performing_pretreatment(UnivarA_Store, filt_df,
                                                                 target_temp, filt_df.columns)

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
          hover_data=["p-value", "FDR adjusted p-value", f'FC ({UnivarA_Store.test_class} / {UnivarA_Store.control_class})'],
          labels={'log2FC': f'log2 (Fold Change)',
               '-log10(Adj. p-value)': '- log10 (Adjusted (Benjamini-Hochberg) p-value)'})

    # Set Threshold lines
    fig.add_vline(x=np.log2(UnivarA_Store.current_univ_params["Fold Change Threshold"]), line_width=2, line_dash="dash", line_color="black")
    fig.add_vline(x=-np.log2(UnivarA_Store.current_univ_params["Fold Change Threshold"]), line_width=2, line_dash="dash", line_color="black")
    fig.add_hline(y=-np.log10(UnivarA_Store.current_univ_params["p-value"]), line_width=2, line_dash="dash", line_color="black")

    fig.update_layout(
        title_text=f'Volcano Plot - {UnivarA_Store.current_univ_params["Test Class"]}/{UnivarA_Store.current_univ_params["Control Class"]}',
        title_x=0.45)

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


### Functions related to the Data Diversity Visualization Plots page of the graphical interface

# Slightly altered version from what is in metanalysis_standard.py
def create_element_counts(data, formula_subset=['Formula',], compute_ratios=True,
                          series=('CHO', 'CHOS', 'CHON', 'CHNS', 'CHONS', 'CHOP', 'CHONP','CHONSP')):
    """Create DataFrame from element counts and concat to original DataFrame.

       Optionally, the ratios of H/C and O/C and element composition series are also computed"""

    # safe guard: remove empty formulae
    formulae = data[formula_subset]

    # count elements
    forms_list = []
    idxs_list = []
    for col in formulae.columns:
        formulae_col = formulae[[col]].dropna()
        for l in formulae_col.index:
            fs = formulae_col.loc[l].iloc[0]
            if type(fs) == str: # For Smart Formula and str based annotations
                forms_list.append(fs)
                idxs_list.append(l)

            else: # For meta_cols and list based annotations
                l_unique = set(fs)
                for f in l_unique:
                    #print(f)
                    forms_list.append(f)
                    idxs_list.append(l)
    ecounts_list = []
    for f in forms_list:
        ecounts_list.append(metsta.element_composition(f))

    result = pd.DataFrame(ecounts_list).fillna(0).astype(int)
    result['idxs'] = idxs_list

    # compute ratios for VK plots
    if compute_ratios:
        result['H/C'] = result['H'] / result['C']
        result['O/C'] = result['O'] / result['C']

    # compute series from compositions

    sorted_series = [''.join(sorted(list(s))) for s in series]
    result_series = []

    for composition in ecounts_list:
        nonzero = ''.join(sorted([k for k, c in composition.items() if c > 0]))

        if nonzero in sorted_series:
            result_series.append(series[sorted_series.index(nonzero)])
        else:
            result_series.append('other')

    result['Series'] = pd.Series(result_series, index=result.index)
    if type(formula_subset) != str:
        result = result.set_index('idxs')
        result = result.drop_duplicates()

    return result


def _plot_VK_diagrams_individual(filt_df, dataviz_store, norm_full_df, target, group):
    "Plot Van Krevelen Diagrams based on parameters passed."

    # Calculating H/C and O/C ratios and series classes for data diversity plots.
    forms = filt_df.dropna(subset=dataviz_store.vk_formula_to_consider, how='all')
    elems = create_element_counts(forms, formula_subset=dataviz_store.vk_formula_to_consider)

    # Set up the list of ranks and the two sloped before and after the midpoint
    rank_values = norm_full_df[elems.index].loc[np.array(target) == group].mean().rank(ascending=False)
    logInt_values = np.log(norm_full_df[elems.index].loc[np.array(target) == group].mean())
    elems['Rank'] = rank_values
    elems['logInt'] = logInt_values

    # If rank was selected
    if dataviz_store.vk_highlight_by == 'Rank':
        # Creates the 2 slopes and saves appropriate information
        tsn = TwoSlopeNorm(vmin=rank_values.min(),
                           vcenter=(1-dataviz_store.vk_midpoint)*(len(elems.index)+1),
                           vmax=rank_values.max())
        sorted_tsn = sorted(tsn(rank_values), reverse=False)
        elems[dataviz_store.vk_highlight_by + 's'] = 1-tsn(rank_values)
        sorted_vals = sorted(rank_values, reverse=True)

    # If log Int was selected
    elif dataviz_store.vk_highlight_by == 'logInt':
        # Creates the 2 slopes and saves appropriate information
        tsn = TwoSlopeNorm(vmin=logInt_values.min(),
                           vcenter=logInt_values.sort_values().iloc[int(dataviz_store.vk_midpoint*len(elems.index))],
                           vmax=logInt_values.max())
        sorted_tsn = sorted(tsn(logInt_values))
        elems[dataviz_store.vk_highlight_by + 's'] = tsn(logInt_values)
        sorted_vals = sorted(logInt_values, reverse=True)

    # In case none of them was selected:
    else:
        # Plot a simple Van Krevelen and end the function
        fig = px.scatter(elems, x="O/C", y="H/C", size=[1,]*len(elems),
                         size_max=dataviz_store.vk_max_dot_size, opacity=0.7,
                         range_x=[-0.05,1.5], range_y=[0.1,2.2],
                    hover_data={'Rank':True, 'logInt':True}, title='Van Krevelen - ' + group)
        filename = f'VK_plot_{group}_formulacolumns'
        # Create appropriate filename
        for cl in dataviz_store.vk_formula_to_consider:
            filename = filename + f'_{cl}'

        return fig, filename

    # See if dots should be colored based on intensity
    if dataviz_store.vk_colour:
        color = dataviz_store.vk_highlight_by + 's'
    else:
        color = None

    # See if dots should have different sizes based on intensity
    if dataviz_store.vk_size:
        s = dataviz_store.vk_highlight_by + 's'
    else:
        s = [1,]*len(elems)

    # Plot Van Krevelen with all the specified parameters
    fig = px.scatter(elems, x="O/C", y="H/C", color=color,
                 color_continuous_scale='Portland', size=s,
                 size_max=dataviz_store.vk_max_dot_size, opacity=0.7,
                 range_x=[-0.05,1.5], range_y=[0.1,2.2],
                    hover_data={dataviz_store.vk_highlight_by + 's':False, # remove
                             'Rank':True, 'logInt':True,
                            }, title='Van Krevelen - ' + group)

    # Set up the colorbar translating the two slopes values into the actual ranks or intensities
    if dataviz_store.vk_show_colorbar:
        a = 0
        tickvals = [] # Values to substitute
        ticktext = [] # Text to replace them with
        extra = [] # For log Int
        # For each tickval, get the ticktext
        for v in np.arange(0,1.01,0.1):
            while v > sorted_tsn[a]:
                a += 1
            tickvals.append(sorted_tsn[a])
            if dataviz_store.vk_highlight_by == 'Rank':
                ticktext.append(len(sorted_vals)-sorted_vals[a]+1)
            else:
                ticktext.append(len(sorted_vals)-a)
                extra.append(sorted_vals[a])

        # Update the colorbar with the correct ticks and text based on the methodology chosen
        if dataviz_store.vk_highlight_by == 'Rank': # Rank
            ticktext = [ticktext[-i] for i in range(1, len(ticktext)+1)]
            fig.update_coloraxes(colorbar_title=dataviz_store.vk_highlight_by, colorbar=dict(
                                      tickvals=tickvals,
                                      ticktext=ticktext))
        else: # Log Int
            true_ticktext = [] # Create a mix of intensity values and their rank
            for i in range(len(ticktext)):
                true_ticktext.append(f'{extra[-1-i]:.3f} ({ticktext[i]})')
            fig.update_coloraxes(colorbar_title=dataviz_store.vk_highlight_by +' (Rank)', colorbar=dict(
                                      tickvals=tickvals,
                                      ticktext=true_ticktext))

    # If no colorbar is to be shown
    else:
        fig.update(layout_coloraxis_showscale=False)

    # Create appropriate filename
    filename = f'VK_plot_{group}_{dataviz_store.vk_highlight_by}order_(midpoint{dataviz_store.vk_midpoint})_formulacolumns'
    for cl in dataviz_store.vk_formula_to_consider:
        filename = filename + f'_{cl}'

    return fig, filename


def vk_add_class_rectangles(fig):
    "Add rectangles delimiting zones for Peptide, Lignins, Tannins, Nucleotides and Phytochemical areas and corresponding legend."

    # Add a rectangle for each 'class'
    # Lipids
    fig.add_shape(type="rect", line_width=2,
        x0=0, y0=1.32, x1=0.6, y1=1.32+1.78,
        line=dict(color="olive"), line_dash='dash',
        showlegend=True, name="Lipids",
    )

    # Amino-Sugars
    fig.add_shape(type="rect", line_width=2,
        x0=0.56, y0=1.55, x1=0.56+0.17, y1=1.55+0.41,
        line=dict(color="salmon"), line_dash='dash',
        showlegend=True, name="Amino-Sugars",
    )

    # Carbohydrates
    fig.add_shape(type="rect", line_width=2,
        x0=0.72, y0=1.53, x1=0.72+0.35, y1=1.53+0.74,
        line=dict(color="Orange"), line_dash='dash',
        showlegend=True, name="Carbohydrates",
    )

    # Peptides
    fig.add_shape(type="rect", line_width=2,
        x0=0.23, y0=1.48, x1=0.23+0.31, y1=1.48+0.51,
        line=dict(color="Green"), line_dash='dash',
        showlegend=True, name="Peptide",
    )

    # Lignins
    fig.add_shape(type="rect", line_width=2,
        x0=0.26, y0=0.77, x1=0.26+0.36, y1=0.77+0.69,
        line=dict(color="Blue"), line_dash='dash',
        showlegend=True, name="Lignins",
    )

    # Tannins
    fig.add_shape(type="rect", line_width=2,
        x0=0.62, y0=0.55, x1=0.62+0.28, y1=0.55+0.7,
        line=dict(color="Purple"), line_dash='dash',
        showlegend=True, name="Tannins",
    )

    # Nucleotides
    fig.add_shape(type="rect", line_width=2,
        x0=0.5, y0=1, x1=0.5+0.7, y1=1+0.8,
        line=dict(color="Red"), line_dash='dash',
        showlegend=True, name="Nucleotides",
    )

    # Phytochemical
    fig.add_shape(type="rect", line_width=2,
        x0=0, y0=0.2, x1=1.15, y1=0.2+1.32,
        line=dict(color="Black"), line_dash='dash',
        showlegend=True, name="Phythochemical",
    )

    # Update and put legend on the plot
    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="right",
        x=0.99
    ))

    return fig


def _plot_KMD_plot_individual(filt_df, dataviz_store, group, neutral_mass_col):
    "Plots Kendrick Mass Defect Plots based on parameters passed."

    # Calculate points for the scatter plot, choose if rounding should happen by rounding up or to the nearest integer
    if dataviz_store.kmd_mass_rounding == 'Up':
        rounding = 'up'
    else:
        rounding = 'nearest'
    nominal, fraction = metsta.calc_kmd(filt_df, rounding=rounding, neutral_mass_col=neutral_mass_col)

    # Create DF that will be used to plot
    new_df = pd.DataFrame()
    new_df['Kendrick Nominal Mass'] = nominal
    new_df['Kendrick Mass Defect'] = fraction

    # See if dots should be coloured by the chemical composition series
    if len(dataviz_store.kmd_formula_to_consider) == 0: # If not
        c = [None, ]*len(new_df)
        new_df['Comp. Class'] = c

    else: # If Yes
        # Calculating H/C and O/C ratios and series classes for data diversity plots.
        forms = filt_df.dropna(subset=dataviz_store.kmd_formula_to_consider, how='all')
        elems = create_element_counts(forms, formula_subset=dataviz_store.kmd_formula_to_consider)

        class_labels = []
        n_form_per_peak = pd.Series(elems.index).value_counts()

        for i in filt_df.index: # For each peak to consider
            # If it has formula assigned by the formula columns selected
            if i in elems.index:
                # See if it has more than one formula
                if n_form_per_peak.loc[i] > 1:
                    # If it has, see if they belong to the same class series or not
                    cl = set(elems.loc[i, 'Series'])
                    if len(cl) == 1: # If yes assign it
                        class_labels.append(elems.loc[i, 'Series'].iloc[0])
                    else: # If not, say it is Ambiguous
                        class_labels.append('Ambiguous')

                # In case it only has one formula
                else:
                    class_labels.append(elems['Series'].loc[i])

            # If it has not formula assigned by the formula columns selected
            else:
                class_labels.append('No Formula')

        new_df['Comp. Class'] = class_labels
        c = class_labels

    # See if dots should have different sizes based on intensity
    s = [1,]*len(new_df)

    # Plot Kendrick Mass Defect Plot with all the specified parameters
    fig = px.scatter(new_df, x='Kendrick Nominal Mass', y='Kendrick Mass Defect', color='Comp. Class', size=s,
                 size_max=dataviz_store.kmd_max_dot_size, opacity=0.5,
                 range_x=[0,1250],
                 category_orders={"Comp. Class": ['CHO', 'CHOS', 'CHON', 'CHNS', 'CHONS', 'CHOP', 'CHONP','CHONSP',
                                                      'other', 'No Formula', 'Ambiguous']},
                 title = 'Kendrick Mass Defect - ' + group
                    )

    if len(dataviz_store.kmd_formula_to_consider) == 0:
        fig.update_layout(showlegend=False)
        filename = f'KMD_plot_{group}_rounded{dataviz_store.kmd_mass_rounding}'
        return fig, filename

    filename = f'KMD_plot_{group}_rounded{dataviz_store.kmd_mass_rounding}_formulacolumns'
    for cl in dataviz_store.vk_formula_to_consider:
        filename = filename + f'_{cl}'

    return fig, filename



### Functions related to the Compound Finder search tool page of the graphical interface

def build_annotation_to_idx_dict(metadata_df, col_list):
    "Builds a dictionary with all the annotations found in the provided columns as keys and the idxs where they are as values."

    # Filter the DataFrame given to only the columns provided
    filt_df = metadata_df[col_list].dropna(how='all')

    ann_to_idxs_dict = {}

    # For each bucket label with the type of compound chosen
    for idx in filt_df.index:
        for col in col_list: # For each column corresponding to the type of column
            # Get the annotations
            annots = filt_df.loc[idx, col]

            # If it is a single string annotation
            if type(annots) == str:
                if annots not in ann_to_idxs_dict.keys(): # And not added to the list
                    # Add to the dict
                    ann_to_idxs_dict[annots] = [idx, ]
                else: # If already added to the list
                    # Add the idx to the dict if it is different from what is there already
                    if idx not in ann_to_idxs_dict[annots]:
                        ann_to_idxs_dict[annots].append(idx)

            # If it is a list-like annotation
            elif type(annots) != float:
                for c in annots: # Go through every annotated compound
                    # Add to the dict when it is not already there
                    if c not in ann_to_idxs_dict.keys():
                        ann_to_idxs_dict[c] = [idx, ]
                    else: # If already added to the list
                        # Add the idx to the dict if it is different from what is there already
                        if idx not in ann_to_idxs_dict[c]:
                            ann_to_idxs_dict[c].append(idx)

    return ann_to_idxs_dict


def plot_sample_bar_plot(comp_finder, target_list, com_exc_compounds):
    "Plot the normalized intensity by sample bar plot."

    # Defining x and y labels for horizontal and vertical barplots
    if comp_finder.sample_bar_plot_type == 'Horizontal': # Horizontal bar plot
        xlabel = 'Normalized Intensity'
        ylabel = 'Samples'
    else: # Vertical bar plot
        xlabel = 'Samples'
        ylabel = 'Normalized Intensity'

    # Different plot behaviour if there is only one index to plot or more than 1
    # In case there is one
    if len(comp_finder.id_df) == 1:
        # Initialize figure
        fig = go.Figure()

        if comp_finder.sample_bar_plot_type == 'Horizontal': # Horizontal bar plot
            for g in target_list.color_classes: # Add a trace per group
                fig.add_trace(go.Bar(name=g,
                                     x=comp_finder.id_df[com_exc_compounds.groups[g]].values[0], y=com_exc_compounds.groups[g],
                                     orientation='h', marker_color=target_list.color_classes[g]))

        else: # Vertical bar plot
            for g in target_list.color_classes: # Add a trace per group
                fig.add_trace(go.Bar(name=g,
                                     x=com_exc_compounds.groups[g], y=comp_finder.id_df[com_exc_compounds.groups[g]].values[0],
                                     orientation='v', marker_color=target_list.color_classes[g]))

        # Update other characteristics of the plot
        fig.update_layout(legend_title_text='Classes',
                     title=comp_finder.id_type + ' - ' + comp_finder.id_comp,
                    xaxis_title=xlabel,
                    yaxis_title=ylabel)


    # In case there is more than one index
    else:
        # Initialize a subplot with n cols equal to the number of idxs
        if comp_finder.sample_bar_plot_type == 'Horizontal': # Horizontal bar plot
            fig = make_subplots(rows=1, cols=len(comp_finder.id_df), subplot_titles=comp_finder.id_df.index)

            for a in range(len(comp_finder.id_df)): # For each index
                for g in target_list.color_classes: # Plot the bar in the corresponding subplot
                    fig.add_trace(go.Bar(name=g + ' - col ' + str(a+1),
                                         x=comp_finder.id_df[com_exc_compounds.groups[g]].values[a], y=com_exc_compounds.groups[g],
                                         orientation='h', marker_color=target_list.color_classes[g]),
                                row=1, col=a+1)

                # Update the axis titles
                fig.update_xaxes(title_text=xlabel, row=1, col=a+1,
                                 range=(0, comp_finder.id_df[target_list.sample_cols].max().max()*1.1))
            # Update the axis titles
            fig.update_yaxes(title_text=ylabel, row=1, col=1)


        else: # Vertical bar plot
            fig = make_subplots(rows=1, cols=len(comp_finder.id_df), subplot_titles=comp_finder.id_df.index, shared_yaxes=True)

            for a in range(len(comp_finder.id_df)): # For each index
                for g in target_list.color_classes: # Plot the bar in the corresponding subplot
                    fig.add_trace(go.Bar(name=g + ' - col ' + str(a+1),
                                         x=com_exc_compounds.groups[g], y=comp_finder.id_df[com_exc_compounds.groups[g]].values[a],
                                         orientation='v', marker_color=target_list.color_classes[g]),
                                row=1, col=a+1)

                # Update the axis titles
                fig.update_xaxes(title_text=xlabel, row=1, col=a+1)
            # Update the axis titles
            fig.update_yaxes(title_text=ylabel, row=1, col=1)

        # Update other characteristics of the plot
        fig.update_layout(legend_title_text='Classes',
                         title=comp_finder.id_type + ' - ' + comp_finder.id_comp)

    return fig


def plot_class_bar_plot(comp_finder, target_list, com_exc_compounds):
    "Plot the normalized intensity by class bar plot."

    # Defining x and y labels for horizontal and vertical barplots
    if comp_finder.class_bar_plot_type == 'Horizontal': # Horizontal bar plot
        xlabel = 'Avg. Normalized Intensity (± St. Deviation)'
        ylabel = 'Classes'
    else: # Vertical bar plot
        xlabel = 'Classes'
        ylabel = 'Avg. Normalized Intensity (± St. Deviation)'

    avg_cols = [col for col in comp_finder.id_df.columns if 'Average' in col]

    # Ignore missing values or treat them as 0
    if comp_finder.ignore_missing_values:
        finder = comp_finder.id_df.copy()
    else:
        finder = comp_finder.id_df.copy().replace({np.nan:0})
        for g in com_exc_compounds.groups:
            finder[g+' Average'] = finder[finder.columns.intersection(com_exc_compounds.groups[g])].mean(axis=1)
            finder[g+' std'] = finder[finder.columns.intersection(com_exc_compounds.groups[g])].std(axis=1)


    # Different plot behaviour if there is only one index to plot or more than 1
    # In case there is one
    if len(finder) == 1:
        # Initialize figure
        fig = go.Figure()

        if comp_finder.class_bar_plot_type == 'Horizontal': # Horizontal bar plot
            for g_avg in avg_cols: # Add a trace per group
                g = g_avg.split(' Average')[0]
                fig.add_trace(go.Bar(name=g,
                            x=finder[g_avg].values, y=[g], error_x=dict(array=finder[g+' std']),
                            orientation='h', marker_color=target_list.color_classes[g]))

        else: # Vertical bar plot
            for g_avg in avg_cols: # Add a trace per group
                g = g_avg.split(' Average')[0]
                fig.add_trace(go.Bar(name=g,
                            x=[g], y=finder[g_avg].values, error_y=dict(array=finder[g+' std']),
                            orientation='v', marker_color=target_list.color_classes[g]))

        # Update other characteristics of the plot
        fig.update_layout(legend_title_text='Classes',
                     title=comp_finder.id_type + ' - ' + comp_finder.id_comp,
                    xaxis_title=xlabel,
                    yaxis_title=ylabel)


    # In case there is more than one index
    else:
        # Initialize a subplot with n cols equal to the number of idxs
        if comp_finder.class_bar_plot_type == 'Horizontal': # Horizontal bar plot
            fig = make_subplots(rows=1, cols=len(finder), subplot_titles=finder.index)

            for a in range(len(finder)): # For each index
                for g_avg in avg_cols: # Plot the bar in the corresponding subplot
                    g = g_avg.split(' Average')[0]
                    fig.add_trace(go.Bar(name=g + ' - col ' + str(a+1),
                                         x=[finder[g_avg].iloc[a]], y=[g],
                                         error_x=dict(array=[finder[g+' std'].iloc[a]]),
                                         orientation='h', marker_color=target_list.color_classes[g]),
                                row=1, col=a+1)

                # Update the axis titles
                fig.update_xaxes(title_text=xlabel, row=1, col=a+1,
                                 range=(0, finder[target_list.sample_cols].max().max()*1.1))
            # Update the axis titles
            fig.update_yaxes(title_text=ylabel, row=1, col=1)


        else: # Vertical bar plot
            fig = make_subplots(rows=1, cols=len(finder), subplot_titles=finder.index, shared_yaxes=True)

            for a in range(len(finder)): # For each index
                for g_avg in avg_cols: # Plot the bar in the corresponding subplot
                    g = g_avg.split(' Average')[0]
                    fig.add_trace(go.Bar(name=g + ' - col ' + str(a+1),
                                         x=[g], y=[finder[g_avg].iloc[a]],
                                         error_y=dict(array=[finder[g+' std'].iloc[a]]),
                                         orientation='v', marker_color=target_list.color_classes[g]),
                                row=1, col=a+1)

                # Update the axis titles
                fig.update_xaxes(title_text=xlabel, row=1, col=a+1)
            # Update the axis titles
            fig.update_yaxes(title_text=ylabel, row=1, col=1)

        # Update other characteristics of the plot
        fig.update_layout(legend_title_text='Classes',
                         title=comp_finder.id_type + ' - ' + comp_finder.id_comp)

    return fig


def plot_class_boxplot(comp_finder, target_list, com_exc_compounds):
    "Plot the normalized intensity by sample bar plot."

    # Defining x and y labels for horizontal and vertical barplots
    if comp_finder.sample_bar_plot_type == 'Horizontal': # Horizontal bar plot
        xlabel = 'Avg. Normalized Intensity'
        ylabel = 'Classes'
    else: # Vertical bar plot
        xlabel = 'Classes'
        ylabel = 'Avg. Normalized Intensity'

    if comp_finder.class_boxplot_points == 'Only outliers':
        boxpoints = 'outliers'
    else:
        boxpoints = 'all'

    # Ignore missing values or treat them as 0
    if comp_finder.ignore_missing_values:
        finder = comp_finder.id_df.copy()
    else:
        finder = comp_finder.id_df.copy().replace({np.nan:0})
        for g in com_exc_compounds.groups:
            finder[g+' Average'] = finder[finder.columns.intersection(com_exc_compounds.groups[g])].mean(axis=1)
            finder[g+' std'] = finder[finder.columns.intersection(com_exc_compounds.groups[g])].std(axis=1)

    # Different plot behaviour if there is only one index to plot or more than 1
    # In case there is one
    if len(finder) == 1:
        # Initialize figure
        fig = go.Figure()

        for g in target_list.color_classes: # Add a trace per group with correct colour
            fig.add_trace(go.Box(name=g,
                                 y=finder[com_exc_compounds.groups[g]].values[0],
                                 marker_color=target_list.color_classes[g], boxpoints=boxpoints))

        # Update other characteristics of the plot
        fig.update_layout(legend_title_text='Classes',
                     title=comp_finder.id_type + ' - ' + comp_finder.id_comp,
                    xaxis_title=xlabel,
                    yaxis_title=ylabel)


    # In case there is more than one index
    else:
        # Initialize a subplot with n cols equal to the number of idxs
        fig = make_subplots(rows=1, cols=len(finder), subplot_titles=finder.index, shared_yaxes=True)

        for a in range(len(finder)): # For each index
            for g in target_list.color_classes: # Plot the bar in the corresponding subplot
                fig.add_trace(go.Box(name=g + ' - col ' + str(a+1),
                                     y=finder[com_exc_compounds.groups[g]].values[a],
                                     marker_color=target_list.color_classes[g], boxpoints=boxpoints),
                            row=1, col=a+1)

            # Update the axis titles
            fig.update_xaxes(title_text=xlabel, row=1, col=a+1)
        # Update the axis titles
        fig.update_yaxes(title_text=ylabel, row=1, col=1)

        # Update other characteristics of the plot
        fig.update_layout(legend_title_text='Classes',
                         title=comp_finder.id_type + ' - ' + comp_finder.id_comp)

    return fig