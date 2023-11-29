
## Needed imports
import pandas as pd
import numpy as np
import panel as pn
import param
import seaborn as sns
import holoviews as hv
import plotly.express as px
import plotly.graph_objects as go
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hier
import matplotlib.pyplot as plt
import matplotlib as mpl
from upsetplot import from_contents

# File with functions to auxiliate the graphical interface
import interface_aux_functions as iaf

# metanalysis_standard.py file
import metanalysis_standard as metsta
import multianalysis as ma

# The initial pages, especially the read file one does not have the nomenclature that I started using later on
# for the different widgets as well as organization
pn.extension('plotly', 'floatpanel', 'katex', 'gridstack', notifications=True)
hv.extension('plotly')
pn.config.sizing_mode="stretch_width"

# To allow the interaction of different packages
pd.DataFrame.iteritems = pd.DataFrame.items

# TODO: Make a way to choose folder where all figures and tables downloaded go to
# TODO: Change UpSetPlot name to IntersectionPlot everywhere
# TODO: Updating packages made a series of future deprecation warnings appear - adapt code to them

# Define pages as classes
# Initial Pages class building with barebones for each class
class OpeningPage:
    def __init__(self):
        self.content = pn.Column("# Welcome to MetsTA!", acronym_section,
                                 "# The Go-To place for your extreme-resolution metabolomics data analysis need.",
                                 pn.Row(homepage_page, pn.pane.Image('Picture_Test.png', height=400)))
    def view(self):
        return self.content

class DataReading:
    def __init__(self):
        self.content = pn.Column("# Section 1: Data Input", """Inputting your Excel or csv MetaboScape file.
                                 If your file is not from MetaboScape, it should have the _m/z_ peak column be called 'Bucket Label'.
                                 After confirming your dataset and moving to the next page, you will not be able to change it unless you restart the program.""",
                                section1page)

    def view(self):
        return self.content

class DataMetadata:
    def __init__(self):
        self.content = pn.Column("# Section 1.1: Selecting Metadata columns and defining your target",
                                page1_1)

    def view(self):
        return self.content
    
class DataFiltering:
    def __init__(self):
        self.content = pn.Column("# Section 1.2: Selecting Data Filtering Method",
                                page1_2)

    def view(self):
        return self.content


class DataAnnotation:
    def __init__(self):
        self.content = pn.Column("# Section 2: Data Annotation", 
    """Perform Annotations based on available databases. **Cannot perform annotation without a selected Neutral Mass column.**
    You can annotate with multiple databases. However, each database is annotated individually.
    Annotation works by assigning to a m/z peak / feature all metabolites of a database that are within the provided error margin.
    Annotation from two different databases might annotate different metabolites for the same m/z peak / feature.
    Thus, each database annotation will generate 3 columns added to the metadata: one with the IDs of the metabolites annotated, another with their formula and another with their name.""",
                                page2)

    def view(self):
        return self.content


class DataPreTreatment:
    def __init__(self):
        
        self.content = pn.Column("# Section 3: Data Pre-Treatment", "Choose the pre-treatment to apply to the data.",
                                page3)

    def view(self):
        return self.content


class ClassColours:
    def __init__(self):

        self.content = pn.Column("# Section 3.1: Select Colours for each Class", "These colours will be used in the different figures made hereafter.",
                                 "### Choose the colours for each class", page4)

    def view(self):
        return self.content


class TransitionalPage:
    def __init__(self):

        self.content = pn.Column("# Select which Analysis you would like to do:", transitional_page)

    def view(self):
        return self.content

class CommonExclusivePage:
    def __init__(self):

        self.content = pn.Column("# Seeing Common and Exclusive Compounds Between Biological Classes",
                                 "This includes an overview analysis as well as Venn Diagrams and UpSetPlots",
                                 "##### Known problem: Repeating IDs leading to problems in Venn Diagrams and UpSetPlots appearing at the same time.",
                                 comexc_page)

    def view(self):
        return self.content

class UnsupervisedAnalysisPage:
    def __init__(self):

        self.content = pn.Column("# Performing Unsupervised Analysis",
                                 """In Multivariate Unsupervised Analysis, class memberships are not provided to the methodologies used. Thus, they are used to observe the intrinsic patterns and structure in the data.
                                  **Principal Component Analysis (PCA)** and **Hierarchical Clustering Analysis (HCA)**.""",
                                 "##### Known problem: Changing HCA plot characteristics many many times can lead to memory issues.",
                                 unsup_analysis_page)

    def view(self):
        return self.content


class SupervisedAnalysisPage:
    def __init__(self):

        self.content = pn.Column("# Performing Supervised Analysis",
                                 """In Multivariate Supervised Analysis, class memberships are provided to the methodologies used.
                                  Thus, they are used to build classifiers and assess how well models classify experimental samples and discriminate biological classes as well as to extract the important metabolites for said classification.
                                  **Random Forest (RF)** and **Partial Least Squares - Discriminant Analysis (PLS-DA)** classifiers and important feature extraction.""",
                                 sup_analysis_page)

    def view(self):
        return self.content


class UnivariateAnalysisPage:
    def __init__(self):

        self.content = pn.Column("# Performing Univariate Analysis",
                                 """In Univariate Analysis, each metabolite (variable) in the experimental dataset is tested individually to observe if there is a significant difference between a test and a control class.
                                  Thus, this does not take into account any interaction between metabolites as multivariate analysis does (and is expected in metabolties within a biological system).
                                  However, it provides the metabolites which are differentially expressed between 2 classes.
                                  It is more suited to when there are a low number of missing values in your data since calculation of the fold change between the classes is severely affected with missing values since it is calculated after missing values imputation and normalization.
                                  Thus, in these cases, fold change values should be taken with a grain (or multiple grains that are actually more like rocks than grains) of salt.
                                  Including Univariate and Fold-Change Analysis""",
                                 univar_analysis_page)

    def view(self):
        return self.content


class DataVisualizationPage:
    def __init__(self):

        self.content = pn.Column("# Data Diversity Visualization Plots",
                                 """This page allows you to plot Van Krevelen, Kendrick Mass Defect  and Chemical Composition Series Plots to visualize the data diversity in your dataset.
                                 **Note**: Press the button to generate the data and update the different sections with relevant information from your workflow.
                                 """, '', '#### Known Issue: Legend in Van Krevelen Plots does not appear. Randomly some Kendrick Mass Defect Plots have bigger points than others.',
                                 data_viz_page)

    def view(self):
        return self.content


class BinSimPage:
    def __init__(self):

        self.content = pn.Column("# Performing BinSim Analysis",
                                 "Include all previous analysis specifically for BinSim treated data",
                                 binsim_analysis_page)

    def view(self):
        return self.content


class CompoundFinderPage:
    def __init__(self):

        self.content = pn.Column("# Find a specific compound", "Observe barplots and boxplots of that specific feature",
                                 comp_finder_page)

    def view(self):
        return self.content


# Homepage
# This section only includes a short HTML description of the software

acronym_string = '''<strong>Met</strong>abolomic<strong>s</strong> <strong>T</strong>reatment and
<strong>A</strong>nalysis Software'''

homepage_string = '''This software is geared towards FT-ICR-MS extreme-resolution <strong>mass spectrometry data</strong>.
<br>
<br>
It is still in the <strong>early stages of development</strong>. As such, much progress in the different areas
of the software from presentation to available tools and features is expected to happen regularly.
<br>
<br>
- It is suited to read the output generated by the MetaboScape 5.0 software (Bruker Daltonics). It can however read any
other type of data if it follows the formatting rules as declared in the "Data Reading" page.
<br>
- Thus, it can handle not only Direct Infusion FT-ICR-MS data but also LC-MS and GC-MS data. It will not use the time
retention parameters of the latter 2 and as such may not be the most suitable tool to analyse said data.
<br>
<br>
- This software includes reading the data, data filtering, data annotation based on mass (or <em>m/z</em> peak)
deviation, data pre-treatment, univariate, multivariate unsupervised and supervised analysis, common and exclusive
compound analysis and data diversity visualization tools (Van Krevelen and Chemical Composition Series).
<br>
- It also includes statistical analysis performe on data treated with Binary Simplification (BinSim),
that is, feature occurrence data that can be useful for extreme-resolution data with plenty of missing values.
<br>
- Data Annotation is made based on the metabolite database given by the user. This database should be formatted according
to the indications in the "Data Annotation" page.<br>
- Every Figure generated by the software can be downloaded and a final report of the analysis can also be generated.
<strong>(TODO)</strong>
<br>
<br>
<hr />
<br>
- Software made by: <strong>FT-ICR and Structural Mass Spectrometry Laboratory, BioISI, Faculdade de Ciências da
Universidade de Lisboa</strong> (<a href="http://ft-icr.rd.ciencias.ulisboa.pt/" target="_blank" rel="nofollow">
ft-icr.rd.ciencias.ulisboa.pt/
</a>)
<br>
<br>
If you use this software in your untargeted metabolomcis data analysis, we would be grateful if you cite our paper
introducing it:
<br>
<br>
- Paper Link: <strong>(LINK)</strong><br>
- Citation: <strong>(CITATION)</strong>
'''
acronym_section = pn.pane.HTML(acronym_string, styles={'font-size': 'large'})
homepage_page = pn.pane.HTML(homepage_string)



# Create the main area and display the first page
main_area = OpeningPage().content



# Param Class to store all DataFrames
# TODO: Put previous dataframes - filtered_df and annotated_df here as well

# Contains before treatment data, treated_data, processed_data, univariate_data, meta_data, bin_data
class DataFrame_Storage(param.Parameterized):
    "Class to contain all the more relevant DataFrames for statistical analysis."

    # Read DataFrame
    read_df = param.DataFrame()

    # Starting DataFrame (after annotation)
    original_df = param.DataFrame()

    # DataFrame treated
    treated_df = param.DataFrame()

    # Metadata
    metadata_df = param.DataFrame()

    # DataFrame for exclusive/common compounds
    processed_df = param.DataFrame()

    # DataFrame for univariate analysis
    univariate_df = param.DataFrame()

    # BinSim treated DataFrame
    binsim_df = param.DataFrame()

    def __init__(self, **params):

        super().__init__(**params)

        self.controls = pn.Param(self, parameters=['treated_df'], name='Pre-Treatment Selection')

    def concat_annots(_, MS_df, annot_df):
        "Joins m/z peak data with annotation data in a single DataFrame."
        return pd.concat((MS_df, annot_df), axis=1)

# Initializing the Store
DataFrame_Store = DataFrame_Storage()



# TODO: Make Reset button to read other datasets.


# Page 1 - Reading File

class FileReading(param.Parameterized):
    """Class to store as attributes file read."""

    temp_target = param.Dict()
    read_df = param.DataFrame(pd.DataFrame())

    def __init__(self, **params):

        super().__init__(**params)

# Initializing store for File Reading
file = FileReading()

# Widgets and reacting functions of page 1
filename = pn.widgets.FileInput(name='Choose file', accept='.csv,.xlsx,.xls')
target_included_in_file = pn.widgets.Checkbox(name='The first row of the file corresponds to the target (sample class labels).', value=False)
temp_target = param.Parameter(default={})

confirm_button_filename = pn.widgets.Button(name='Read File', button_type='primary', disabled=True)
tooltip_file = pn.widgets.TooltipIcon(
    value="""Provided file must come from MetaboScape. Alternatively, the column with the _m/z_ peaks should be labelled 'Bucket label'.""")
confirm_button_step1 = pn.widgets.Button(icon=iaf.img_confirm_button, name='Confirm - Next Step', button_type='success',
                                         disabled=True)

# Update button so it can be pressed after you put something in the filename
@pn.depends(filename.param.filename, watch=True)
def _update_confirm_button_filename(filename):
    if filename != '':
        confirm_button_filename.disabled = False
    else:
        confirm_button_filename.disabled = True

def _confirm_button_filename(event):
    "Reads the file given."

    # Read the file, updating widgets and parameters
    file.read_df, file.temp_target = iaf.read_file(filename.filename, target_included_in_file.value)

    # Enabling button for next step
    section1page[2] = pn.widgets.DataFrame(file.read_df)
    confirm_button_step1.disabled = False
# Function happens when you press the button        
confirm_button_filename.on_click(_confirm_button_filename)

# Alternatively, provide option to read an example data
load_example_df_button = pn.widgets.Button(name='Load Example Dataset', button_type='warning', disabled=False, height=50)
tooltip_example_df = pn.widgets.TooltipIcon(
    value="""Example Data consists of 15 FT-ICR-MS samples of 5 strains of the Yeast Saccharomyces cerevisiae with previsouly assigned annotations.""")

def _load_example_df_button(event):
    "Reads the example file ofthe software."

    # Read the file, updating widgets and parameters
    file.read_df, file.temp_target = iaf.read_file('5yeasts_notnorm.csv', False)

    # Enabling button for next step
    section1page[2] = pn.widgets.DataFrame(file.read_df)
    confirm_button_step1.disabled = False
# Function happens when you press the button
load_example_df_button.on_click(_load_example_df_button)


# Confirm file, show next page, disable reading files, update columns of the dataset read
def _confirm_step1(event):
    # Enabling/Disabling appropriate Widgets
    page1_1_button.disabled = False
    confirm_button_filename.disabled = True
    load_example_df_button.disabled = True
    target_included_in_file.disabled = True
    filename.disabled = True

    DataFrame_Store.read_df = file.read_df # Update DataFrame store

    # Update all the options for the Data Metadata Step - CheckBox and RadioBox Widgets
    checkbox_formula.options = list(DataFrame_Store.read_df.columns)
    checkbox_annotation.options = list(DataFrame_Store.read_df.columns)
    radiobox_neutral_mass.options = ['None'] + list(DataFrame_Store.read_df.columns)
    checkbox_others.options = list(DataFrame_Store.read_df.columns)
    checkbox_samples.options = list(DataFrame_Store.read_df.columns)

    # Update the Main page
    main_area.clear()
    show_page(pages["Data Metadata"])


# Call the appropriate functions when the buttons are pressed
confirm_button_step1.on_click(_confirm_step1)

# Setting up the page layout
section1page = pn.Column(pn.Row(pn.Column(filename, target_included_in_file), pn.Row(load_example_df_button, tooltip_example_df)),
                         pn.Row(confirm_button_filename, tooltip_file),
                         file.read_df, confirm_button_step1)




# Page 1-1 - Selecting sample columns and Target

# Making checkbox widgets for each category
checkbox_formula = pn.widgets.CheckBoxGroup(
    name='Formula', value=['Formula'], options=list(file.read_df.columns),
    inline=False, disabled=False)

checkbox_annotation = pn.widgets.CheckBoxGroup(
    name='Annotation', value=['Name'], options=list(file.read_df.columns),
    inline=False, disabled=False)

# Select only one instead of multiple - RadioBox Widget
radiobox_neutral_mass = pn.widgets.RadioBoxGroup(
    name='Neutral Mass', value='Neutral Mass', options=['None'] + list(file.read_df.columns),
    inline=False, disabled=False)

checkbox_others = pn.widgets.CheckBoxGroup(
    name='Others', options=list(file.read_df.columns),
    inline=False, disabled=False)

checkbox_samples = pn.widgets.CheckBoxGroup(
    name='Samples', options=list(file.read_df.columns),
    inline=False, disabled=True)

# Arranging the checkboxes
checkbox_arrangement = pn.Column(
    pn.Row('#### Select Formula columns:                ', '#### Select Annotations columns:            ', 
           '#### Select Neutral Mass column:            ', '#### Select Other NON-SAMPLE columns:       ',
           '#### After confirming, check sample columns:'),
    pn.Row(pn.Column(checkbox_formula, scroll=True, height=400), pn.Column(checkbox_annotation, scroll=True, height=400),
          pn.Column(radiobox_neutral_mass, scroll=True, height=400), pn.Column(checkbox_others, scroll=True, height=400),
          pn.Column(checkbox_samples, scroll=True, height=400)))

# Button to confirm the selection in checkboxes and function detailing what happens when you press it
confirm_button_column_selection = pn.widgets.Button(icon=iaf.img_confirm_button, name='Confirm Columns', button_type='success')

# Confirm column selection, update sample columns, make target editable while providing a possible target
def _update_confirm_column_selection(event):
    # Deduce sample columns
    sample_cols = []
    cols = list(DataFrame_Store.read_df.columns)
    for col in cols:
        if col not in checkbox_formula.value:
            if col not in checkbox_annotation.value:
                if col not in radiobox_neutral_mass.value:
                    if col not in checkbox_others.value:
                        sample_cols.append(col)
    target_widget.disabled = False # Make you able to type in the target
    
    checkbox_samples.value = sample_cols # Update the samples checkbox
    target_list.sample_cols = sample_cols # Save the sample cols

    # Target widget box
    if checkbox_samples.value != '':
        # If the target was in the first row of the file read, pass it to here
        if file.temp_target != {}:
            tg = [file.temp_target[s] for s in sample_cols] # Update to target present in file provided
            target_widget.placeholder = ','.join(tg)
            target_widget.value = ','.join(tg)

        # If not, attempt to make 'an automatic target' as a suggestion
        else:
            tg = [i.split('_')[0] for i in checkbox_samples.value]
            target_widget.placeholder = ','.join(tg)
            target_widget.value = ','.join(tg)

    else:
        target_widget.placeholder = target_placeholder
    
confirm_button_column_selection.on_click(_update_confirm_column_selection)

# Make the target widget
target_placeholder = "A,A,A,A,A,B,B,B,B,B"
target_widget = pn.widgets.TextAreaInput(name='Target - Class Labels', placeholder=target_placeholder, 
                                         max_length=5000, height=100, disabled=True)
target_tooltip = pn.widgets.TooltipIcon(value="Provide the class labels of your samples. No spaces between labels.")
# Create button widgets
confirm_button_target = pn.widgets.Button(icon=iaf.img_confirm_button, name='Confirm Target', button_type='success',
                                         disabled=True)
confirm_button_next_step_1_1 = pn.widgets.Button(icon=iaf.img_confirm_button,
                                                 name='Next Step - Data Filtering and Characteristics',
                                                 button_type='success', disabled=True)


# Make button be pressable when you have a target
@pn.depends(target_widget.param.value, watch=True)
def _update_read_target_button(target_widget):
    if target_widget != '':
        confirm_button_target.disabled = False
    else:
        confirm_button_target.disabled = True


# Make sure that target makes sense in comparison to the number of sample columns
def _update_confirm_target(event):
    target = target_widget.value.split(',')
    sample_cols = target_list.sample_cols
    if len(sample_cols) != len(target):
        pn.state.notifications.error(
            f'Number of class labels ({len(target)}) is different than the number of sample columns ({len(sample_cols)}).')
    else:
        confirm_button_next_step_1_1.disabled = False

# Call the function
confirm_button_target.on_click(_update_confirm_target)

# Going to the next step function
def _confirm_button_next_step_1_1(event):
    page1_2_button.disabled = False
    #page2_button.disabled = False
    main_area.clear()
    show_page(pages["Data Filtering"])

# Call the function    
confirm_button_next_step_1_1.on_click(_confirm_button_next_step_1_1)

# Organizing the page layout
page1_1 = pn.Column()
page1_1.append(checkbox_arrangement)
page1_1.append(confirm_button_column_selection)
page1_1.append(pn.Row(target_widget, target_tooltip))
page1_1.append(confirm_button_target)
page1_1.append(confirm_button_next_step_1_1)




# Page 1-2 - Feature Filtering and Data Characteristics

# Widgets for feature filtering
filt_method = pn.widgets.Select(name="Feature Filter Method", value="Total Samples",
                                       options=['Total Samples', 'Class Samples', None])
filt_kw = pn.widgets.IntSlider(name="Feature Filter Keyword", value=2, start=0, end=len(target_widget.value.split(',')))
limits_filt = {"Total Samples": len(target_widget.value.split(',')),
               "Class Samples": pd.Series(target_widget.value.split(',')).value_counts().min(), None: 1}

# Change the IntSlider limits based on the number of class labels
@pn.depends(target = target_widget.param.value,
            method = filt_method.param.value, watch=True)
def _update_filt_kw_limits(target, method):
    if method == 'Total Samples':
        filt_kw.end = len(target.split(','))
    elif method == 'Class Samples':
        filt_kw.end = pd.Series(target.split(',')).value_counts().min()
    else:
        filt_kw.end = 1

filt_method_tooltip = pn.widgets.TooltipIcon(value=
    """'Total Samples' requires a feature to appear in at least x samples in the whole dataset to be retained.
    'Class Samples' requires a feature to appear in at least x samples of at least 1 class to be retained.""")
filt_kw_tooltip = pn.widgets.TooltipIcon(value=
    """How many samples a feature has to appear to be retained based on the method chosen before.""")

# Preparing DataFrames
filtered_df = pn.widgets.DataFrame(pd.DataFrame(), name='Filtered DataFrame')
characteristics_df = pn.widgets.DataFrame(pd.DataFrame(), name='Characteristics DataFrame')

# Button to perform filtering
confirm_button_initial_filtering = pn.widgets.Button(icon=iaf.img_confirm_button, name='Perform Filtering',
                                                     button_type='success', disabled=False)

# Button to next step
confirm_button_next_step_2 = pn.widgets.Button(icon=iaf.img_confirm_button, name='Next Step - Annotation',
                                                     button_type='success', disabled=False)


# Call filtering function and extend the page with data characteristics and results of filtering
def _confirm_button_initial_filtering(event):
    "Perform feature filtering."

    sample_cols = target_list.sample_cols # Select sample columns
    target = target_widget.value.split(',') # See target
    # Perform filtering
    if filt_method.value == 'Total Samples':
        f_meth = 'total_samples'
    elif filt_method.value == 'Class Samples':
        f_meth = 'class_samples'
    filtered_df.value, characteristics_df.value = iaf.initial_filtering(DataFrame_Store.read_df,
                                    sample_cols, target=target, filt_method=f_meth, filt_kw=filt_kw.value)

    annotated_df.value = pd.DataFrame(index=filtered_df.value.index)

    # Locking in the parameters used for feature filtering
    UnivarA_Store.locking_filtering_params(filt_method, filt_kw)

    # Setup the page if not setup yet
    if len(page1_2) == 3:
        page1_2.extend(['#### Characteristics of the Dataset',characteristics_df,'#### Filtered Dataset',filtered_df,
                       confirm_button_next_step_2])

# Call the function
confirm_button_initial_filtering.on_click(_confirm_button_initial_filtering)

# Go to next step function and calling it
def _confirm_button_next_step_1_2(event):
    "Ends step 1-2 and goes to Data Annotation page."
    page2_button.disabled = False
    main_area.clear()
    show_page(pages["Data Annotation"])
confirm_button_next_step_2.on_click(_confirm_button_next_step_1_2)

# Initial page layout
page1_2 = pn.Column(pn.Row(filt_method, filt_method_tooltip), pn.Row(filt_kw, filt_kw_tooltip),
                    confirm_button_initial_filtering)




# Page 2 - Annotation of Metabolites

# Widgets for selecting number of databases
n_databases_show = pn.widgets.IntInput(name='Nº of Databases to annotate', value=1, step=1, start=0, end=5)
n_databases = pn.widgets.IntInput(name='Nº of Databases to annotate', value=1, step=1, start=0, end=5)
tooltip_n_databases = pn.widgets.TooltipIcon(value="Select how many (0-5) databases you want to use for annotation.")
# Button to perform filtering
confirm_button_n_databases = pn.widgets.Button(icon=iaf.img_confirm_button, name='Select Databases',
                                                     button_type='success', disabled=False)

# Class to organize how a single Database section will be shown and that can be repeated
class DatabaseSection():
    "Set up the parameters needed to read one metabolite database."

    def __init__(self):
        # Setting up the relevant widgets
        static_db_file = pn.widgets.StaticText(name='Database File', value='', styles={'align':'center'})
        db_file_input = pn.widgets.TextInput(placeholder='hmdb.csv')
        
        static_db_abv = pn.widgets.StaticText(name='Abbreviation', value='')
        db_abv_input = pn.widgets.TextInput(placeholder='HMDB')

        static_db_IDcol = pn.widgets.StaticText(name='DB ID - Index Column', value='')
        db_IDcol_input = pn.widgets.TextInput(placeholder='accession')

        static_db_annotation = pn.widgets.StaticText(name='DB Annotation Column', value='')
        db_annotation_input = pn.widgets.TextInput(placeholder='Name')

        static_db_formula = pn.widgets.StaticText(name='DB Formula Column', value='')
        db_formula_input = pn.widgets.TextInput(placeholder='Formula')

        #static_db_mass = pn.widgets.StaticText(name='DB Mass Column', value='')
        #db_mass_input = pn.widgets.TextInput(placeholder='monoisotopic_molecular_weight')
        
        confirm_button_db = pn.widgets.Button(name='Read Database', button_type='success', disabled=True)
        
        db = pn.widgets.DataFrame(pd.DataFrame(), name='Database') # Store the database
        
        # Make the button pressable when every field has something
        @pn.depends(a=db_file_input.param.value,
                    b=db_abv_input.param.value,
                    c=db_IDcol_input.param.value,
                    d=db_annotation_input.param.value,
                    e=db_formula_input.param.value,
                    watch=True)
        def _update_confirm_button_db(a,b,c,d,e):
            if a != '':
                if b!= '':
                    if c!= '':
                        if d != '':
                            if e != '':
                                confirm_button_db.disabled = False

        # React to press of the buttons by reading the database and updating all status and values of the class with the read data
        # Also make sure database layout section continues organized
        @pn.depends(a=db_file_input.param.value,
                    b=db_abv_input.param.value,
                    c=db_IDcol_input.param.value,
                    d=db_annotation_input.param.value,
                    e=db_formula_input.param.value,
                    button=confirm_button_db.param.clicks,
                    watch=True)
        def _press_confirm_button_db(a,b,c,d,e,button):
            if button != 0:
                # Adjust the names
                filename=a
                abv=b
                ID_col=c
                name_col=d
                formula_col=e
                button=button
                self.file = db_file_input.value
                self.abv = db_abv_input.value
                self.IDcol = db_IDcol_input.value
                self.annotation = db_annotation_input.value
                self.formula = db_formula_input.value

                # Read the database
                db.value = iaf.read_database(filename, abv, ID_col, name_col, formula_col)
                confirm_button_db.param.clicks = 0
                self.db = db

                # Description of the database and indication of successful database reading
                if len(self.content) == 6:
                    self.content.append(f'Database {self.abv} has {len(self.db.value)} metabolites.')
                else:
                    self.content[6] = f'Database {self.abv} has {len(self.db.value)} metabolites.'
                self.read.value = True

                # If re-reading the databases, eliminate elements after database reading from the page to re-run them after.
                while len(page2)>4:
                    page2.pop(-1)
                confirm_button_databases_read.disabled = False
                confirm_button_annotation_perform.disabled = True

        # Initial parameters of the Database section so we can grab values easier down the line
        self.content = pn.Column(pn.Row(static_db_file, db_file_input),
                          pn.Row(static_db_abv, db_abv_input),
                          pn.Row(static_db_IDcol, db_IDcol_input),
                          pn.Row(static_db_annotation, db_annotation_input),
                          pn.Row(static_db_formula, db_formula_input),
                                confirm_button_db)
        self.file = db_file_input.value
        self.abv = db_abv_input.value
        self.IDcol = db_IDcol_input.value
        self.annotation = db_annotation_input.value
        self.formula = db_formula_input.value
        #self.mass = db_mass_input.value
        self.db = db
        self.read = pn.widgets.Switch(name='Switch', value=False) # Important for verification later
                
# Have the maximum 5 database sections ready and organized
DB_dict = {'1':DatabaseSection(), '2':DatabaseSection(), '3':DatabaseSection(), '4':DatabaseSection(),
           '5':DatabaseSection()}
        
dbs_arrangement_all = pn.Row(DB_dict['1'].content, DB_dict['2'].content, DB_dict['3'].content,
                   DB_dict['4'].content, DB_dict['5'].content)

# Make the designated number of database sections appear
def _confirm_button_n_databases(event):
    n_databases.value = n_databases_show.value
    titles = pn.Row()
    dbs_arrangement = pn.Row()
    for i in range(n_databases.value):
        dbs_arrangement.append(dbs_arrangement_all[i])
        titles.append(f'#### Database {i+1}')
    
    # Keep the page layout organized
    if len(page2) == 1:
        page2.append(titles)
        page2.append(dbs_arrangement)
    else:
        while len(page2)>4:
            page2.pop(-1)
        page2[1] = pn.Row(titles)
        page2[2] = dbs_arrangement

    if n_databases.value == 0: # Case where no database is going to be used for annotation
        page2.append(confirm_button_databases_read)
        confirm_button_databases_read.disabled = False

confirm_button_n_databases.on_click(_confirm_button_n_databases)

# Initial page layout
page2 = pn.Column(pn.Row(n_databases_show, tooltip_n_databases, confirm_button_n_databases))
confirm_button_databases_read = pn.widgets.Button(icon=iaf.img_confirm_button, name='Confirm Databases',
                                                 button_type='success', disabled=False)

# Make button to confirm databases appear after all databases are read
@pn.depends(db1=DB_dict['1'].read.param.value,
            db2=DB_dict['2'].read.param.value,
            db3=DB_dict['3'].read.param.value,
            db4=DB_dict['4'].read.param.value,
            db5=DB_dict['5'].read.param.value, watch=True)
def _press_confirm_button_db(db1, db2, db3, db4, db5):
    all_read = False
    for i in range(n_databases.value):
        n_data = str(i+1)
        if DB_dict[n_data].read.value == True:
            if i == n_databases.value - 1:
                all_read = True
            continue
        else:
            break
    if all_read:
        page2.append(confirm_button_databases_read)

# Annotation parameters widgets
annotation_margin_method_radio = pn.widgets.RadioBoxGroup(name='Annotation Margin Method', value='PPM Deviation',
                                                options=['PPM Deviation', 'Absolute Dalton Deviation'], inline=True)
annotation_ppm_deviation = pn.widgets.IntInput(name='Maximum PPM Deviation', value=1, step=1, start=1)
annotation_Da_deviation = pn.widgets.FloatInput(name='Maximum Absolute Dalton Deviation', value=0.001, step=0.001, 
                                                page_step_multiplier=10)
tooltip_annotation = pn.widgets.TooltipIcon(
    value="Choose the maximum allowed deviation (in PPM or Da) for annotating a metabolite.")
confirm_button_annotation_perform = pn.widgets.Button(name='Perform Annotation',
                                                     button_type='success', disabled=False)

# Update the parameter input given based on the methodology chosen
@pn.depends(annotation_margin_method_radio.param.value, watch=True)
def _annotation_margin_method(method):
    if method == 'PPM Deviation':
        annotation_param_selection[1][0] = annotation_ppm_deviation
    else:
        annotation_param_selection[1][0] = annotation_Da_deviation

# Organizing section of parameters for annotation
annotation_param_selection = pn.Column(annotation_margin_method_radio, 
                                       pn.Row(annotation_ppm_deviation, tooltip_annotation),
                                       confirm_button_annotation_perform)

# Make the annotation part of the layout appear and disable button to confirm databases
def _confirm_button_databases_read(event):
    confirm_button_databases_read.disabled = True
    confirm_button_annotation_perform.disabled = False
    annotated_df.value = pd.DataFrame(index=filtered_df.value.index) # Setup the annotation df

    page2.append(annotation_param_selection)
confirm_button_databases_read.on_click(_confirm_button_databases_read)

# Widgets for annotation part
performing_annotation_arrangement = pn.Column() # Start with empty page for the annotation widgets
tqdm_database = {i+1:pn.widgets.Tqdm() for i in range(n_databases.end)}
verbose_annotated_compounds = {i+1:pn.widgets.StaticText(name='', value=f'') for i in range(n_databases.end)}
annotated_df = pn.widgets.DataFrame(pd.DataFrame(index=filtered_df.value.index))

# Function to perform metabolite annotation (also contributes to updating the page)
def metabolite_annotation():
    "Perform metabolite annotation for every database loaded."

    # For each database, perform annotation adding a section with 4 columns (ID, metabolite name, formula and number of matches)
    for i in range(n_databases.value):
        # Get the correct database
        db_to_use = str(i+1)

        # Prepare columns
        matched_ids_col = 'Matched '+DB_dict[db_to_use].abv+' IDs'
        matched_names_col = 'Matched '+DB_dict[db_to_use].abv+' names'
        matched_formulas_col = 'Matched '+DB_dict[db_to_use].abv+' formulas'
        match_count_col = DB_dict[db_to_use].abv+' match count'
        annotated_df.value[matched_ids_col] = ''
        annotated_df.value[matched_names_col] = ""
        annotated_df.value[matched_formulas_col] = ""
        annotated_df.value[match_count_col] = ""

        # Update the page layout correctly, whether it is a repeat annotation or new one
        if len(performing_annotation_arrangement) == i:
            performing_annotation_arrangement.append(pn.Row(f'Annotating {DB_dict[db_to_use].abv} Database:',
                                                            tqdm_database[i+1], verbose_annotated_compounds[i+1]))
        else:
            performing_annotation_arrangement[i] = pn.Row(f'Annotating {DB_dict[db_to_use].abv} Database:',
                                                            tqdm_database[i+1], verbose_annotated_compounds[i+1])

        # Annotation for each metabolite
        for a in tqdm_database[i+1](range(filtered_df.value.shape[0])):
            matched_ids = []
            matched_names = []
            matched_formulas = []

            # See candidates for annotation to add
            # Option 1 - Based on maximum PPM Deviation
            if annotation_margin_method_radio.value == 'PPM Deviation':
                ppm_margin = annotation_ppm_deviation.value

                candidates_for_annotation = abs((DB_dict[db_to_use].db.value['Mass']-filtered_df.value[
                    radiobox_neutral_mass.value].iloc[a])/filtered_df.value[radiobox_neutral_mass.value].iloc[a])*10**6
                candidates_for_annotation = candidates_for_annotation[candidates_for_annotation<ppm_margin]

            # Option 2 - Based on maximum Dalton Deviation
            elif annotation_margin_method_radio.value == 'Absolute Dalton Deviation':
                Da_margin = annotation_Da_deviation.value

                candidates_for_annotation = abs(
                    DB_dict[db_to_use].db.value['Mass']-filtered_df.value[radiobox_neutral_mass.value][a])
                candidates_for_annotation = candidates_for_annotation[candidates_for_annotation<Da_margin]

            # Store candidates
            for m in candidates_for_annotation.index:
                matched_ids.append(m)
                matched_names.append(DB_dict[db_to_use].db.value[DB_dict[db_to_use].annotation][m])
                matched_formulas.append(DB_dict[db_to_use].db.value[DB_dict[db_to_use].formula][m])

            # Add the annotation candidates
            if len(matched_ids) > 0:
                annotated_df.value.at[annotated_df.value.index[a], matched_ids_col] = matched_ids
                annotated_df.value.at[annotated_df.value.index[a], matched_names_col] = matched_names
                annotated_df.value.at[annotated_df.value.index[a], matched_formulas_col] = matched_formulas
                annotated_df.value.at[annotated_df.value.index[a], match_count_col] = len(matched_ids)
            else:
                annotated_df.value.at[annotated_df.value.index[a], matched_ids_col] = np.nan
                annotated_df.value.at[annotated_df.value.index[a], matched_names_col] = np.nan
                annotated_df.value.at[annotated_df.value.index[a], matched_formulas_col] = np.nan
                annotated_df.value.at[annotated_df.value.index[a], match_count_col] = np.nan

        verbose_annotated_compounds[
            i+1].value = f'Annotated {annotated_df.value[matched_ids_col].notnull().sum()} compounds.'

# Perform annotation, update page layout
def _press_confirm_annotation_perform(event):
    "Perform metabolite annotation and update the page layout accordingly."
    while len(performing_annotation_arrangement)>0:
        performing_annotation_arrangement.pop(-1)
    page2.append(performing_annotation_arrangement)
    confirm_button_annotation_perform.disabled=True
    # Perform metabolite annotation
    metabolite_annotation()
    page2.append(confirm_button_next_step_3)

confirm_button_annotation_perform.on_click(_press_confirm_annotation_perform)

# Button to next step
confirm_button_next_step_3 = pn.widgets.Button(icon=iaf.img_confirm_button, name='Next Step - Data Pre-Treatment',
                                                     button_type='success', disabled=False)

# Go to next step function and calling it
def _confirm_button_next_step_3(event):
    page3_button.disabled = False
    # Join Filtered and Annotated dfs to make the original DataFrame for pre-treatment
    DataFrame_Store.original_df = DataFrame_Store.concat_annots(filtered_df.value, annotated_df.value)
    # Creates in the DataFrame a column ('Has Match?') that indicates if a feature was annotated either previously to being
    # inputted in this software or using the data annotation of this software
    iaf.creating_has_match_column(DataFrame_Store, n_databases, checkbox_annotation)

    # Update to show the Data Pre-Treatment page
    main_area.clear()
    show_page(pages["Data Pre-Treatment"])
confirm_button_next_step_3.on_click(_confirm_button_next_step_3)




#### This marks the separation to use mainly param instead of only panel

# Page 3 - Data Pre-Treatment

# Param to encompass the choice of Pre-Treatment methods
# TODO: Introduce way to read reference sample for PQN normalization
class PreTreatment(param.Parameterized):
    """Class to store as attributes the parameters chosen to perform pre-treatments.

       Each type of pre-treatment has 2 main parameters:
        '_method': indicates the methodology to be used in each pre-treatment step.
        '_kw': indicates a keyword parameter to use relative to the methodology chosen in each pre-treatment step."""

    # Update the keyword slider based on methodology chosen

    # Missing Value Imputation
    mvi_method = param.Selector(default="Minimum of Sample")
    mvi_kw = param.Number(default=0.2, bounds=(0,1))

    # Normalization
    norm_method = param.Selector(default="Reference Feature")
    norm_kw = param.Selector(default=None)

    # Transformation
    tf_method = param.Selector(default="Generalized Logarithmic Transformation (glog)")
    tf_kw = param.Number(default=None)

    # Scaling
    scaling_method = param.Selector(default="Pareto Scaling")
    scaling_kw = param.String(default="")

    # Confirm Pre-Treatment Selection
    confirm_button = param.Boolean(default=False)

    # Function to confirm Pre-Treatment Selection and Updating DataFrames
    # TODO: Make Metadata appear in a cleaner way
    def _confirm_button_press(self, event):
        "Perform pre-treatment and update page layout."
        treat, proc, uni, meta, bin = iaf.performing_pretreatment(self, DataFrame_Store.original_df,
                                                                   target_list.target, target_list.sample_cols)
        DataFrame_Store.treated_df, DataFrame_Store.processed_df, DataFrame_Store.univariate_df = treat, proc, uni
        DataFrame_Store.metadata_df, DataFrame_Store.binsim_df = meta, bin

        # Locking in pre-treatment parameters chosen
        UnivarA_Store.locking_pretreatment_params(self)

        page3[:2,2:5] = pn.Tabs(('Treated Data', DataFrame_Store.treated_df),
            ('Metadata', DataFrame_Store.metadata_df.T),
            ('BinSim Treated Data', DataFrame_Store.binsim_df), height=600, dynamic=True)
        confirm_button_next_step_4.disabled = False
        save_data_dataframes_button.disabled = False


    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets = {
            'mvi_method': pn.widgets.Select(name="Missing Value Imputation Method", value = 'Minimum of sample',
                            options=['Minimum of Sample', 'Minimum of Feature', 'Minimum of Data', 'Zero'],
                description='What Minimum to use as reference and what fraction of it for constant value imputation.'),
            'mvi_kw': pn.widgets.FloatSlider(name="Missing Value Imputation Fraction of Minimum Chosen",
                                             value=0.2, step=0.01),

            'norm_method': pn.widgets.Select(name="Normalization Method", value = 'Write reference feature here',
                            options=['Reference Feature', 'Total Intensity Sum', 'PQN', 'Quantile', 'None']),
            'norm_kw': pn.widgets.MultiChoice(name="Normalization Reference Sample or Feature", max_items=1,
                                              option_limit=8, search_option_limit=8, disabled=True,
                        description="""If you do not find your reference feature (Norm. by Ref. Feat.) or your sample (if you are choosing a sample for PQN), start writing it in the box below until it appears.
                        Mean or Median use the mean or median of every sample as reference, respectively."""),

            'tf_method': pn.widgets.Select(name="Transformation Method",
                                           value = 'Generalized Logarithmic Transformation (glog)',
                            options=['Generalized Logarithmic Transformation (glog)', None]),
            'tf_kw': pn.widgets.FloatInput(name="Transformation Lambda Parameter", value=None,
                    description="""Lambda value in glog transformation. Leave empty or 0 for usual log transformation.
                    Ignored if Transformation Method is None.
            glog equation is: log\N{SUBSCRIPT TWO}((y + \u221A(y\N{SUPERSCRIPT TWO} + lambda\N{SUPERSCRIPT TWO}) ) /2)
                    """),

            'scaling_method': pn.widgets.Select(name="Scaling Method", value = 'Pareto Scaling',
                            options=['Pareto Scaling', 'Auto Scaling', 'Mean Centering', 'Range Scaling',
                                     'Vast Scaling', 'Level Scaling', None]),
            'scaling_kw': pn.widgets.Select(name="Scaling Factor (for Level Scaling only)", value='', disabled=True,
                                           options=['Average', 'Median'], description='After mean-centering, divide each value by its respective feature average or median.'),

            'confirm_button': pn.widgets.Button(name="Confirm Pre-Treatment", button_type='primary'),
        }
        self.controls = pn.Param(self, widgets=widgets, name='Pre-Treatment Selection')

# Initializing Pre-Treatment Param
PreTreatment_Method = PreTreatment()

# Update options for missing value imputation keyword based on missing value imputation
limits_mvi = {"Minimum of Sample": 1, "Minimum of Feature": 1, "Minimum of Data": 1, "Zero": 0}
@pn.depends(PreTreatment_Method.controls.widgets['mvi_method'].param.value, watch=True)
def _update_limit_mvi(mvi_method):
    if mvi_method == None:
        PreTreatment_Method.controls.widgets['mvi_kw'].value = 0
    PreTreatment_Method.controls.widgets['mvi_kw'].end = limits_mvi[
        PreTreatment_Method.controls.widgets['mvi_method'].value]

# Update options for normalization keyword based on normalization
options_norm = {"Reference Feature": None, "Total Intensity Sum": '', "PQN": ["mean", "median"],
                "Quantile": ["mean", "median"], 'None': ''}
@pn.depends(PreTreatment_Method.controls.widgets['norm_method'].param.value, watch=True)
def _update_options_norm(norm_method):
    if norm_method in ['Total Intensity Sum', 'None']:
        PreTreatment_Method.controls.widgets['norm_kw'].disabled = True
        PreTreatment_Method.controls.widgets['norm_kw'].placeholder = ''
        PreTreatment_Method.controls.widgets['norm_kw'].value = []

    elif norm_method == 'Reference Feature':
        PreTreatment_Method.controls.widgets['norm_kw'].value = []
        PreTreatment_Method.controls.widgets['norm_kw'].options = list(filtered_df.value.index)
        PreTreatment_Method.controls.widgets['norm_kw'].placeholder = 'Reference Feature'
        PreTreatment_Method.controls.widgets['norm_kw'].disabled = False

    elif norm_method == 'PQN':
        PreTreatment_Method.controls.widgets['norm_kw'].value = []
        PreTreatment_Method.controls.widgets['norm_kw'].options = options_norm[norm_method] + list(target_list.sample_cols)
        PreTreatment_Method.controls.widgets['norm_kw'].placeholder = ''
        PreTreatment_Method.controls.widgets['norm_kw'].disabled = False

    else:
        PreTreatment_Method.controls.widgets['norm_kw'].value = []
        PreTreatment_Method.controls.widgets['norm_kw'].options = options_norm[norm_method]
        PreTreatment_Method.controls.widgets['norm_kw'].placeholder = ''
        PreTreatment_Method.controls.widgets['norm_kw'].disabled = False

# Update options for scaling keyword based on scaling
@pn.depends(PreTreatment_Method.controls.widgets['scaling_method'].param.value, watch=True)
def _update_options_scaling(scaling_method):
    if scaling_method == 'Level Scaling':
        PreTreatment_Method.controls.widgets['scaling_kw'].disabled = False
        PreTreatment_Method.controls.widgets['scaling_kw'].value = 'Average'

    else:
        PreTreatment_Method.controls.widgets['scaling_kw'].disabled = True
        PreTreatment_Method.controls.widgets['scaling_kw'].value = ''

# Click button to confirm Pre-Treatment
PreTreatment_Method.controls.widgets['confirm_button'].on_click(PreTreatment_Method._confirm_button_press)


# Button to next step
confirm_button_next_step_4 = pn.widgets.Button(icon=iaf.img_confirm_button, name='Next Step - Class Colours',
                                                     button_type='success', disabled=True)

# Go to next step function and calling it
def _confirm_button_next_step_4(event):
    page4_button.disabled = False

    # Filling the target storage with the correct target and default colours
    target_list.target = target_widget.value.split(',')
    target_list.classes = list(pd.unique(np.array(target_list.target)))
    target_list(target_widget.value.split(','), colours)

    # Setting up the layout of page 4
    n_classes = len(target_list.color_classes)
    for row in range(0, n_classes, 5):
        new_row = pn.Row()
        for col in range(5):
            key = list(target_list.color_classes.keys())[row+col]
            new_row.append(pn.widgets.ColorPicker(value=target_list.color_classes[key], name=str(key)))

        page4.append(new_row)
    page4.append(confirm_button_next_step_transitionalpage)

    # Update the main area
    main_area.clear()
    show_page(pages["Class Colours"])
confirm_button_next_step_4.on_click(_confirm_button_next_step_4)

# Widget to save dataframes in .csv format (specifically annotated data, treated data and metadata)
save_data_dataframes_button = pn.widgets.Button(name='Save representative Dataframes as .csv files (in current folder)',
                                                button_type='warning', icon=iaf.download_icon, disabled=True)

# When pressing the button, downloads the dataframes
def _save_data_dataframes_button(event):
    DataFrame_Store.original_df.to_csv('annotated_df.csv')
    DataFrame_Store.treated_df.to_csv('treated_df.csv')
    pd.concat((DataFrame_Store.metadata_df, DataFrame_Store.treated_df.T), axis=1).to_csv('complete_treated_df.csv')
    pn.state.notifications.success(f'DataFrames successfully saved.')

save_data_dataframes_button.on_click(_save_data_dataframes_button)

# Organize Page Layout
page3 = pn.GridSpec(mode='override')
page3[:2,0:2] = PreTreatment_Method.controls
page3[:2,2:5] = pn.Tabs(('Treated Data', DataFrame_Store.treated_df),
        ('Metadata', DataFrame_Store.metadata_df),
       ('BinSim Treated Data', DataFrame_Store.binsim_df), height=600, dynamic=True)
page3[2, :] = pn.Column(confirm_button_next_step_4,
                        '## Optionally save data as .csv files',
                        '''The button saves data after annotation (without pre-treatment), the treated intensity data without metadata and with the metadata.
                        These files are saved in the current folder as annotated_df.csv, treated_df.csv, complete_treated_df.csv.''',
                        save_data_dataframes_button)




# Page 4 - Choosing Class Colours

# Default colours - 10 possible colours
colours = sns.color_palette('tab10', 10)


# Class to store target and to store target colours
class TargetStorage(param.Parameterized):
    "Class to store information regarding target, classes and colours associated with the classes."

    # Function to do when calling the function
    updater = param.Callable(iaf.TargetStorage_filling)

    # Setting up parameters/attributes
    target = param.List(default=target_widget.value.split(',')) # Target List
    classes = param.List(default=[]) # Classes List
    color_classes = param.Dict(default={}) # Colours assigned to classes
    sample_cols = param.List(default=checkbox_samples.value) # Columns related to samples

    def __init__(self, **params):
        super().__init__(**params)

    def __call__(self, target, colours):
        self.color_classes = self.updater(target, colours)

        return self

# Initialize the target store
target_list = TargetStorage()

# Button to next step and to confirm colours
confirm_button_next_step_transitionalpage = pn.widgets.Button(icon=iaf.img_confirm_button, name='Next Step - Analysis',
                                                     button_type='success', disabled=False)
# Confirm colours, go to next step function and calling it, enabling all buttons for analysis and performing initial computations for each analysis
def _confirm_button_next_step_5(event):

    # Pass the colours picked to the target storage
    n_classes = len(target_list.color_classes)
    for row in range(0, n_classes, 5):
        for col in range(5):
            key = list(target_list.color_classes.keys())[row+col]
            target_list.color_classes[key] = page4[row//5][col].value

    # Enable all statistical analysis related buttons
    page5_button.disabled = False
    page6_button.disabled = False
    page7_button.disabled = False
    page8_button.disabled = False
    page9_button.disabled = False
    page10_button.disabled = False
    page11_button.disabled = False

    # Initial Calculations for PCA and storing initial plots
    principaldf, var, loadings = metsta.compute_df_with_PCs_VE_loadings(DataFrame_Store.treated_df,
                                       n_components=n_components_compute.value,
                                       whiten=True, labels=target_list.target, return_var_ratios_and_loadings=True)
    PCA_params.pca_scores = principaldf
    PCA_params.explained_variance = var
    PCA_params.pca_loadings = pd.DataFrame(loadings)
    PCA_params.PCA_plot[0] = iaf._plot_PCA(PCA_params, target_list)
    middle_page_PCA[0,1:3] = pn.pane.Plotly(PCA_params.PCA_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_plot', 'scale':4,}})
    PCA_params.exp_var_fig_plot[0] = iaf._plot_PCA_explained_variance(PCA_params)
    end_page_PCA[0] = pn.pane.Plotly(PCA_params.exp_var_fig_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_exp_var_plot', 'scale':4,}})
    PCA_params.scatter_PCA_plot[0] = iaf._scatter_PCA_plot(PCA_params, target_list)
    end_page_PCA[1] = pn.pane.Plotly(PCA_params.scatter_PCA_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_scatter_plot', 'scale':4,}})

    # Initial calculations for HCA and storing initial plots
    HCA_params.dists = dist.pdist(DataFrame_Store.treated_df, metric=HCA_params.dist_metric)
    HCA_params.Z = hier.linkage(HCA_params.dists, method=HCA_params.link_metric)
    HCA_params.HCA_plot[0] = _plot_HCA()
    page_HCA[0:6,1:4] = HCA_params.HCA_plot[0]

    # Updating Widgets for Unsupervised Analysis
    UnivarA_Store._update_widgets()

    # Updating Widgets for Supervised Analysis
    PLSDA_store.n_folds_limits(target_list)
    RF_store.n_folds_limits(target_list)

    # Updating the layout for the transitional page
    main_area.clear()
    show_page(pages["Transitional Page"])
confirm_button_next_step_transitionalpage.on_click(_confirm_button_next_step_5)

# Setting up empty page
page4 = pn.Column()




# Transitional Page Layout

# Buttons for the main analyses steps
ComExc_A = pn.widgets.Button(name='Common/Exclusive Comp.', button_type='default')
Unsup_A = pn.widgets.Button(name='Unsupervised Analysis', button_type='default')
Sup_A = pn.widgets.Button(name='Supervised Analysis', button_type='default')
Univariate_A = pn.widgets.Button(name='Univariate Analysis', button_type='default')
DataViz_A = pn.widgets.Button(name='Data Visualization', button_type='default')

# Buttons for the secondary analyses steps
BinSim_A = pn.widgets.Button(name='BinSim Specific Analysis', button_type='default')
CompFinder_A = pn.widgets.Button(name='Compound Finder', button_type='default')
ToBeAdded_A = pn.widgets.Button(name='More to be added', button_type='danger', disabled=True)

# Functions for pressing each button
def _confirm_button_ComExc_A(event):
    main_area.clear()
    show_page(pages["Common and Exclusive Compounds"])
ComExc_A.on_click(_confirm_button_ComExc_A)

def _confirm_button_Unsup_A(event):
    main_area.clear()
    show_page(pages["Unsupervised Analysis"])
Unsup_A.on_click(_confirm_button_Unsup_A)

def _confirm_button_Sup_A(event):
    main_area.clear()
    show_page(pages["Supervised Analysis"])
Sup_A.on_click(_confirm_button_Sup_A)

def _confirm_button_Univariate_A(event):
    main_area.clear()
    show_page(pages["Univariate Analysis"])
Univariate_A.on_click(_confirm_button_Univariate_A)

def _confirm_button_DataViz_A(event):
    main_area.clear()
    show_page(pages["Data Visualization"])
DataViz_A.on_click(_confirm_button_DataViz_A)

def _confirm_button_BinSim_A(event):
    main_area.clear()
    show_page(pages["BinSim Analysis"])
BinSim_A.on_click(_confirm_button_BinSim_A)

def _confirm_button_CompFinder_A(event):
    main_area.clear()
    show_page(pages["Compound Finder"])
CompFinder_A.on_click(_confirm_button_CompFinder_A)

# Transitional page Layout
transitional_page = pn.Column(pn.Row(ComExc_A, Unsup_A, Sup_A, Univariate_A, DataViz_A),
                             '## Other Options:',
                             pn.Row(BinSim_A, CompFinder_A, ToBeAdded_A))




# Page for Common and Exclusive Compounds
# TODO: Something is funky with the updates of the different pages, maybe due to overlapping pythons IDs. The page doesn't update correctly sometimes.
# Root cause is unknown

# Param Class to store parameters and data regarding Common and Exclusive Compounds
class ComExc_Storage(param.Parameterized):
    "Class to store all information on common and exclusive compounds and to plot Venn diagrams and UpSetPlots."

    # Dictionaries to group information
    groups = param.Dict(default={})
    group_dfs = param.Dict(default={})
    group_dfs_ids = param.Dict(default={})
    groups_description = param.String(default='')

    # DataFrame with the common metabolites
    common_all = param.DataFrame()
    common_all_id = param.DataFrame()

    # Dictionary with DataFrames for exclusive compounds
    exclusives = param.Dict(default={})
    exclusives_id = param.Dict(default={})

    com_exc_desc = param.String(default='')

    # Parameters and DataFrame of chosen sub-section of classes
    class_subset = param.List(default=list())
    df_type = param.String(default='Metabolites in chosen classes')
    annot = param.String(default='See annotated and non-annotated compounds')
    specific_cl_df = param.DataFrame()
    class_specific_cl_desc = param.Number(default=0)

    # Venn Diagram parameters
    venn_class_subset = param.List(default=list())
    venn_alpha = param.Number(default=0.3, bounds=(0,1.0))
    type_of_venn = param.String(default='All Metabolites (Annotated)')
    dpi_venn = param.Number(default=200)

    # UpsetPlot parameters
    upset_class_subset = param.List(default=list())
    upset_include_counts_percentages = param.String(default='Show Nº of metabolites')
    dpi_upset = param.Number(default=200)

    # Storing figure
    Venn_plot = param.List(default=['Pane for Venn Diagram'])
    UpSetPlot = param.List(default=['Pane for UpSetPlot1', 'Pane for UpSetPlot2'])


    # Update the subset df plot based on specifications chosen
    @param.depends('class_subset', 'df_type', 'annot', watch=True)
    def _update_specific_cl_df(self):
        "Update the subset DataFrame based on inputs given."

        df_list = []

        # All features of the dataset
        if self.annot == 'See annotated and non-annotated compounds':
            for cl in self.class_subset: # See the DataFrames of each class considered
                df_list.append(self.group_dfs[cl])
            # Calculate common features between them
            if len(df_list) != 0:
                df_common = metsta.common(df_list)

                if self.df_type == 'Metabolites in chosen classes and no other class':
                    # For this specific case
                    # Remove common compounds that also appear in other not selected classes
                    for s in self.group_dfs:
                        if s not in self.class_subset:
                            exclude = self.group_dfs[s]
                            df_common = df_common.loc[~(df_common.index.isin(exclude.index))]

            else:
                df_common = pd.DataFrame()

        # Only annotated features of the dataset
        else:
            for cl in self.class_subset: # See the DataFrames of each class considered
                df_list.append(self.group_dfs_ids[cl])
            # Calculate common features between them
            if len(df_list) != 0:
                df_common = metsta.common(df_list)

                if self.df_type == 'Metabolites in chosen classes and no other class':
                    # For this specific case
                    # Remove common compounds that also appear in other not selected classes
                    for s in self.group_dfs_ids:
                        if s not in self.class_subset:
                            exclude = self.group_dfs_ids[s]
                            df_common = df_common.loc[~(df_common.index.isin(exclude.index))]

            else:
                df_common = pd.DataFrame()

        # Updating information on class and panel sections
        self.specific_cl_df = df_common
        self.specific_cl_desc = len(df_common.index)
        subsetdf_comexc_section_page[0:4,1:3] = pn.widgets.DataFrame(self.specific_cl_df)
        subsetdf_comexc_section_page[2,0].value = self.specific_cl_desc


    # Update the Venn Diagram based on specifications chosen
    @param.depends('venn_class_subset', 'venn_alpha', 'type_of_venn', 'dpi_venn', watch=True)
    def _update_Venn_diagram(self):
        "Update the Venn Diagram based on parameters given."
        if 1 < len(self.venn_class_subset) < 7:
            self.Venn_plot = []
            self.Venn_plot.append(iaf._plot_Venn_diagram(self, target_list))
            venn_page[0,1:3] = pn.pane.Matplotlib(self.Venn_plot[0], height=800)
            plt.close()
        else:
            pn.state.notifications.info(f'Venn Diagram can only be made with 2 to 6 different classes. You currently have {len(self.venn_class_subset)} classes.')
        #if len(end_page_comexc) >= 2:
        #    end_page_comexc[1] = venn_page


    # Update the UpSetPlot based on specifications chosen
    @param.depends('upset_class_subset', 'upset_include_counts_percentages', 'dpi_upset', watch=True)
    def _update_Upsetplot(self):
        "Update the UpSetPlots based on parameters given."
        if 1 < len(self.upset_class_subset):
            self.UpSetPlot = []

            # Select the relevant data - all metabolites and annotated metabolites
            groups_dict = {}
            groups_dict_ids = {}
            for df in self.upset_class_subset:
                groups_dict[df] = self.group_dfs[df].index
                groups_dict_ids[df] = self.group_dfs_ids[df].index

            ups = from_contents(groups_dict)
            self.UpSetPlot.append(iaf._plot_upsetplots(self, groups_dict, ups))
            upset_page[2] = pn.pane.Matplotlib(self.UpSetPlot[0], height=300)
            plt.close()

            ups_ids = from_contents(groups_dict_ids)

            self.UpSetPlot.append(iaf._plot_upsetplots(self, groups_dict_ids, ups_ids))

            upset_page[5] = pn.pane.Matplotlib(self.UpSetPlot[1], height=300)
            plt.close()
        else:
            pn.state.notifications.info(f'UpSetPlot can only be made with 2 or more classes. You currently have {len(self.upset_class_subset)} classes.')


    def _update_widgets(self):
        "Update widgets to consider the update list of classes in the target."
        widgets = {
            'class_subset': pn.widgets.CheckBoxGroup(name='Classes', value=[], options=target_list.classes,
                                inline=False, disabled=False),
            'df_type': pn.widgets.RadioBoxGroup(name='Type of DataFrame', value='Metabolites in chosen classes',
                                options=['Metabolites in chosen classes',
                                         'Metabolites in chosen classes and no other class'],
                                inline=False, disabled=False),
            'annot': pn.widgets.RadioBoxGroup(name='Annotated', value='See annotated and non-annotated compounds',
                                options=['See annotated and non-annotated compounds',
                                         'Only see annotated compounds'],
                                inline=False, disabled=False),
            'venn_class_subset': pn.widgets.CheckBoxGroup(name='Classes', value=target_list.classes,
                                options=target_list.classes, inline=False, disabled=False),
            'venn_alpha': pn.widgets.FloatInput(name='Transparency of Venn Diagram Circles', value=0.3, start=0.0, end=1.0,
                                step = 0.01),
            'type_of_venn': pn.widgets.Select(name='Consider what type of metabolites',
                                value='All Metabolites (Annotated)',
                                options=['All Metabolites (Annotated)', 'All Metabolites','Annotated Metabolites']),
            'dpi_venn': pn.widgets.IntInput(name="DPI (Resolution)",
                                    value=200, step=10, start=100, disabled=False,
                                    description='Set the resolution of diagram'),
            'upset_class_subset': pn.widgets.CheckBoxGroup(name='Classes', value=target_list.classes,
                                options=target_list.classes, inline=False, disabled=False),
            'upset_include_counts_percentages': pn.widgets.RadioBoxGroup(name='Show:',
                                value='Show Nº of metabolites', options=['Show Nº and % of metabolites',
                                'Show Nº of metabolites', 'Show neither'], inline=False, disabled=False),
            'dpi_upset': pn.widgets.IntInput(name="DPI (Resolution)", value=200, step=10, start=100, disabled=False,
                                    description='Set the resolution of UpSetPlots'),
        }
        self.class_subset = []
        self.venn_class_subset = target_list.classes
        self.upset_class_subset = target_list.classes

        # Control panel for the overview section, for the Venn diagram section and for the Upsetplot section
        return (pn.Param(self, parameters=['class_subset', 'df_type', 'annot'], widgets=widgets, name='Subset of Data to See'),
                pn.Param(self, parameters=['venn_class_subset', 'venn_alpha', 'type_of_venn', 'dpi_venn'], widgets=widgets,
                         name='Parameters to draw Venn Diagram'),
                pn.Param(self, parameters=['upset_class_subset', 'upset_include_counts_percentages','dpi_upset'],
                         widgets=widgets, name='Parameters to draw UpSetPlots', default_layout=pn.Row))


    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets = {
            'class_subset': pn.widgets.CheckBoxGroup(name='Classes', value=[], options=target_list.classes,
                                inline=False, disabled=False),
            'df_type': pn.widgets.RadioBoxGroup(name='Type of DataFrame', value='Metabolites in chosen classes',
                                options=['Metabolites in chosen classes',
                                         'Metabolites in chosen classes and no other class'],
                                inline=False, disabled=False),
            'annot': pn.widgets.RadioBoxGroup(name='Annotated', value='See annotated and non-annotated compounds',
                                options=['See annotated and non-annotated compounds',
                                         'Only see annotated compounds'],
                                inline=False, disabled=False),
            'venn_class_subset': pn.widgets.CheckBoxGroup(name='Classes', value=target_list.classes,
                                options=target_list.classes, inline=False, disabled=False),
            'venn_alpha': pn.widgets.FloatInput(name='Transparency of Venn Diagram Circles', value=0.3, start=0.0, end=1.0,
                                step = 0.01),
            'type_of_venn': pn.widgets.Select(name='Consider what type of metabolites',
                                value='All Metabolites (Annotated)',
                                options=['All Metabolites (Annotated)', 'All Metabolites','Annotated Metabolites']),
            'dpi_venn': pn.widgets.IntInput(name="DPI (Resolution)",
                                    value=200, step=10, start=100, disabled=False,
                                    description='Set the resolution of diagram'),
            'upset_class_subset': pn.widgets.CheckBoxGroup(name='Classes', value=target_list.classes,
                                options=target_list.classes, inline=False, disabled=False),
            'upset_include_counts_percentages': pn.widgets.RadioBoxGroup(name='Show:',
                                value='Show Nº of metabolites', options=['Show Nº and % of metabolites',
                                'Show Nº of metabolites', 'Show neither'], inline=False, disabled=False),
            'dpi_upset': pn.widgets.IntInput(name="DPI (Resolution)", value=200, step=10, start=100, disabled=False,
                                    description='Set the resolution of UpSetPlots'),
        }

        # Control panel for the overview section
        self.controls = pn.Param(self,
                                 parameters=['class_subset', 'df_type', 'annot'],
                                 widgets=widgets, name='Subset of Data to See')
        # Control panel for the Venn diagram section
        self.venn_controls = pn.Param(self,
                                 parameters=['venn_class_subset', 'venn_alpha', 'type_of_venn', 'dpi_venn'],
                                 widgets=widgets, name='Parameters to draw Venn Diagram')
        # Control panel for the UpSetPlots section
        self.upsetplot_controls = pn.Param(self, parameters=['upset_class_subset', 'upset_include_counts_percentages',
                                                             'dpi_upset'],
                                 widgets=widgets, name='Parameters to draw UpSetPlots', default_layout=pn.Row)


# Initialize common and exclusive compound storage
com_exc_compounds = ComExc_Storage()

# Widgets for the page
compute_ComExc_button = pn.widgets.Button(name='Compute', button_type='success', height=50)
checkbox_com_exc = pn.widgets.CheckBoxGroup(name='Include:', value=['Venn Diagram', 'UpSetPlot'],
                                            options=['Venn Diagram', 'UpSetPlot'], inline=True)

# When pressing the button, performs common and exclusive compound calculations, and sets up the different pages in the tabs section
def _compute_ComExc_button(event):
    "Actions after pressing the button"

    iaf._group_compounds_per_class(com_exc_compounds, target_list, DataFrame_Store) # Add compounds per class dfs
    iaf._compute_com_exc_compounds(com_exc_compounds) # Add common to all and exclusive to each class compounds dfs

    new_controls = com_exc_compounds._update_widgets() # Update Widgets
    com_exc_compounds.controls = new_controls[0]
    com_exc_compounds.venn_controls = new_controls[1]
    com_exc_compounds.upsetplot_controls = new_controls[2]

    subsetdf_comexc_section_page[0:2,0] = com_exc_compounds.controls
    venn_page[0,0][0] = com_exc_compounds.venn_controls
    upset_page[0] = com_exc_compounds.upsetplot_controls

    # Setting up the overview page
    if len(overview_page) == 0:
        overview_page.append(pn.pane.Markdown(com_exc_compounds.groups_description))
        overview_page.append(pn.pane.Markdown(com_exc_compounds.com_exc_desc))
        overview_page.append(pn.Row(save_comexc_dfs_button, tooltip_comexc_dfs))
        overview_page.append(subsetdf_comexc_section_page)
    else:
        overview_page[0] = pn.pane.Markdown(com_exc_compounds.groups_description)
        overview_page[1] = pn.pane.Markdown(com_exc_compounds.com_exc_desc)
        overview_page[2] = pn.Row(save_comexc_dfs_button, tooltip_comexc_dfs)
        overview_page[3] = subsetdf_comexc_section_page

    # Organizing the Accordion
    if len(end_page_comexc) == 0:
        end_page_comexc.append(('Overview', overview_page))

    while len(end_page_comexc) > 1:
        end_page_comexc.pop(1)

    if 'Venn Diagram' in checkbox_com_exc.value:
        end_page_comexc.append(('Venn Diagram', venn_page))
    if 'UpSetPlot' in checkbox_com_exc.value:
        end_page_comexc.append(('UpSetPlot', upset_page))

    com_exc_compounds._update_specific_cl_df() # Update the DataFrame shown in the overview tab
    com_exc_compounds.venn_class_subset = com_exc_compounds.venn_controls[1].value # Updating starting subset value for Venn diagram
    com_exc_compounds.upset_class_subset = com_exc_compounds.upsetplot_controls[0][1].value # Updating starting subset value for UpSetPlots


# Action when pressing the button
compute_ComExc_button.on_click(_compute_ComExc_button)

# First section of the page
initial_page_comexc = pn.GridSpec(mode='override')
initial_page_comexc[0,0] = '### Compute Common and Exclusive Compound Calculations'
initial_page_comexc[1,0] = pn.Row(pn.widgets.StaticText(name='Include', value='', styles={'font-size': 'large'}),
                                  checkbox_com_exc)
initial_page_comexc[:2,1] = compute_ComExc_button

# 1st Tab
overview_page = pn.Column()

# Widget related to creating and saving the common and exclusive compound Excel
save_comexc_dfs_button = pn.widgets.Button(name='Save Excel with annotated comp. common to all or exclusive to each class',
                                           button_type='warning', icon=iaf.download_icon)
tooltip_comexc_dfs = pn.widgets.TooltipIcon(value="""File is created only considering **Annotated** compounds by this software or previously.
    Excel created has n+1 sheets where n is the number of classes.
    One sheet is for the compounds common to all classes. The rest is one per class showing their exclusive compounds.""")

# When pressing the button, creates and saves the common and exclusive compound Excel
def _save_comexc_dfs_button(event):
    "Creates and saves to Excelthe common and exclusive annotated compound dataframes."
    # Building the common and exclusive dataframes
    common_df, exclusive_dfs = iaf.build_common_exclusive_dfs_to_save(com_exc_compounds, target_list, checkbox_annotation,
                                                              checkbox_formula, radiobox_neutral_mass, checkbox_others)
    iaf.common_exclusive_compound_excel_writer(common_df, exclusive_dfs)

    pn.state.notifications.success(f'Common and Exclusive Ann. Compound Excel successfully saved.')

save_comexc_dfs_button.on_click(_save_comexc_dfs_button)


# Widget to save dataframe (with characteristics chosen) in .csv format
save_comexc_df_button = pn.widgets.Button(name='Save shown Dataframe as .csv (in current folder)',
                                                button_type='warning', icon=iaf.download_icon)

# When pressing the button, downloads the dataframe
def _save_comexc_df_button(event):
    # Building the datafile name
    # Based on what type of metabolites were chosen
    if com_exc_compounds.annot == 'See annotated and non-annotated compounds':
        annot_string = 'All_metabolites'
    else:
        annot_string = 'Annotated_metabolites'

    # Based on type of dataframe chosen
    if com_exc_compounds.df_type == 'Metabolites in chosen classes':
        df_type_string = '_in_all_chosen'
    else:
        df_type_string = '_only_in_all_chosen'

    # Based on the biological classes chosen
    class_subset_string = ''
    for cl in com_exc_compounds.class_subset:
        class_subset_string = class_subset_string + '_'+cl

    # Final name
    filename_string = annot_string + df_type_string + class_subset_string + '_classes.csv'
    # Saving the file
    com_exc_compounds.specific_cl_df.to_csv(filename_string)
    pn.state.notifications.success(f'{filename_string} successfully saved.')

save_comexc_df_button.on_click(_save_comexc_df_button)

# Create specific subset df section of the page  for the overview Tab
subsetdf_comexc_section_page = pn.GridSpec(mode='override')
subsetdf_comexc_section_page[0:2,0] = com_exc_compounds.controls
subsetdf_comexc_section_page[2,0] = pn.indicators.Number(name='Nº of Metabolites', font_size='14pt', title_size='14pt',
                                                         value=com_exc_compounds.class_specific_cl_desc)
subsetdf_comexc_section_page[3,0] = save_comexc_df_button
subsetdf_comexc_section_page[0:4,1:3] = pn.widgets.DataFrame(com_exc_compounds.specific_cl_df)


# Widget to save Venn Diagram (needed since it is a matplotlib plot instead of a plotly plot)
save_Venn_diag_button = pn.widgets.Button(name='Save as a png (in current folder)', button_type='success',
                                         icon=iaf.download_icon)
# When pressing the button, downloads the figure
def _save_Venn_diag_button(event):
    com_exc_compounds.Venn_plot[0].savefig('Venn_diagram.png', dpi=com_exc_compounds.dpi_venn)
    pn.state.notifications.success(f'Venn Diagram successfully saved.')
save_Venn_diag_button.on_click(_save_Venn_diag_button)

# Create specific section of the page for the Venn Diagram tab
venn_page = pn.GridSpec(mode='override')
venn_page[0,0] = pn.Column(com_exc_compounds.venn_controls, save_Venn_diag_button)
venn_page[0,1:3] = com_exc_compounds.Venn_plot[0]


# Widgets to save UpSetPlots (needed since they are matplotlib plots instead of a plotly plots)
save_UpSetPlot_allmets_button = pn.widgets.Button(name='Save as a png (in current folder)', button_type='success',
                                         icon=iaf.download_icon)
save_UpSetPlot_annotatedmets_button = pn.widgets.Button(name='Save as a png (in current folder)', button_type='success',
                                         icon=iaf.download_icon)
# When pressing the button, downloads the figures
def _save_UpSetPlot_allmets_button(event):
    com_exc_compounds.UpSetPlot[0].savefig('UpSetPlot_all_metabolites.png', dpi=com_exc_compounds.dpi_upset)
    pn.state.notifications.success(f'Intersection Plot (all metabolites) successfully saved.')
save_UpSetPlot_allmets_button.on_click(_save_UpSetPlot_allmets_button)

def _save_UpSetPlot_annotatedmets_button(event):
    com_exc_compounds.UpSetPlot[1].savefig('UpSetPlot_annotated_metabolites.png', dpi=com_exc_compounds.dpi_upset)
    pn.state.notifications.success(f'Intersection Plot (annotated metabolites) successfully saved.')
save_UpSetPlot_annotatedmets_button.on_click(_save_UpSetPlot_annotatedmets_button)

# Create specific section of the page for the UpSetPlot tab
upset_page = pn.Column(com_exc_compounds.upsetplot_controls, # Control parameters for UpSetPlot
                      '## UpSetPlot with all metabolites of the dataset',
                       com_exc_compounds.UpSetPlot[0],
                       save_UpSetPlot_allmets_button,
                       '## UpSetPlot with only annotated metabolites in the dataset',
                       com_exc_compounds.UpSetPlot[1],
                       save_UpSetPlot_annotatedmets_button)


end_page_comexc = pn.Accordion(toggle=True)

comexc_page = pn.Column(initial_page_comexc, end_page_comexc)




# Page for Unsupervised Analysis

# Tab for PCA analysis
# TODO: PCA is straight away computed with 10 components. If your data has less than 10 features, this will lead to an error and everything will fail.
# Should this possibility be taken into account?

# Param Class to store parameters and data regarding PCA
class PCA_Storage(param.Parameterized):
    "Class use to store PCA and PCA related plots parameters, dataframes and figures."
    # Update the keyword slider based on methodology chosen

    # Components to compute PCA
    n_components = param.Number(default=10)

    # N dimensions to show in plot
    n_dimensions = param.String(default='2 Components')

    # PCAs to plot
    PCx = param.String(default='PC 1')
    PCy = param.String(default='PC 2')
    PCz = param.String(default='PC 3')

    # Draw Ellipses
    ellipse_draw = param.Boolean(default=True)
    confidence = param.Number(default=0.95)
    confidence_std = param.Number(default=2)

    # Other params
    dot_size = param.Number(default=5)

    # Computed PCA df
    pca_scores = param.DataFrame()
    explained_variance = param.Array()
    pca_loadings = param.DataFrame()

    # Storing figures
    PCA_plot = param.List(default=['a'])
    exp_var_fig_plot = param.List(default=['a'])
    scatter_PCA_plot = param.List(default=['a'])

    # Update the PCA plot
    @param.depends('n_dimensions', 'PCx', 'PCy', 'PCz', 'ellipse_draw', 'confidence', 'confidence_std', 'dot_size', watch=True)
    def _update_PCA_plot(self):
        self.PCA_plot[0] = iaf._plot_PCA(PCA_params, target_list)
        middle_page_PCA[0,1:3] = pn.pane.Plotly(self.PCA_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_plot', 'scale':4}})

    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets = {
            'n_components': pn.widgets.IntInput(name='Number of Components to Compute:',
                                           value=10, step=1, start=2, end=20,
                                          description='Select 2-20 components.'),
            'n_dimensions': pn.widgets.RadioBoxGroup(name="n_dimensions",
                                                     value = '2 Components',
                            options=['2 Components', '3 Components'], inline = True),
            'PCx': pn.widgets.Select(name="Principal Component in X axis",
                                    value='PC 1', options=['PC '+str(i+1) for i in range(self.n_components)]),
            'PCy': pn.widgets.Select(name="Principal Component in Y axis",
                                    value='PC 2', options=['PC '+str(i+1) for i in range(self.n_components)]),
            'PCz': pn.widgets.Select(name="Principal Component in Z axis",
                                    value='PC 3', options=['PC '+str(i+1) for i in range(self.n_components)],
                                     disabled=True),
            'ellipse_draw': pn.widgets.Checkbox(name="Draw Confidence Ellipses (Only for 2-D)",
                                    value=True),
            'confidence': pn.widgets.FloatInput(name="Ellipse Confidence Level",
                                    value=0.95, step=0.01, start=0, end=1, disabled=False,
                                    description='E.g. 0.95 confidence level means a 95% confidence level of the sigma error ellipse based on the covariance matrix'),
            'confidence_std': pn.widgets.FloatInput(name="Ellipse Confidence Level (standard deviation unit)",
                                    value=2, step=0.50, start=1, end=5,
                description='''Only used if Level of Confidence of Ellipse is 0.
                E.g. 1 stands for 68.3% and 2 for 95.4% confidence level of the sigma error ellipse based on the covariance matrix''',
                                    disabled=True),
            'dot_size': pn.widgets.IntSlider(name='Size of Points in Projection', start=1, end=20, value=5, step=1),
        }
        self.controls = pn.Param(self, parameters=['n_dimensions', 'PCx', 'PCy', 'PCz',
                                                   'ellipse_draw', 'confidence', 'confidence_std', 'dot_size'],
                                 widgets=widgets, name='Parameters for PCA Projection plot')

# Running initial param to store PCA details
PCA_params = PCA_Storage()


# Extra widgets for the page
compute_PCA_button = pn.widgets.Button(name='Compute', button_type='success')
n_components_compute = pn.widgets.IntInput(name='Number of Components to Compute:',
                                           value=10, step=1, start=2, end=20,
                                          description='Select 2-20 components.')

# When pressing the button, runs again the PCA with the designated number of components
def _compute_PCA_button(event):
    principaldf, var, loadings = metsta.compute_df_with_PCs_VE_loadings(DataFrame_Store.treated_df,
                                       n_components=n_components_compute.value,
                                       whiten=True, labels=target_list.target, return_var_ratios_and_loadings=True)
    PCA_params.pca_scores = principaldf
    PCA_params.explained_variance = var
    PCA_params.pca_loadings = pd.DataFrame(loadings)
    PCA_params.n_components = n_components_compute.value
    PCA_params.controls.widgets['n_components'].value = n_components_compute.value

compute_PCA_button.on_click(_compute_PCA_button)


# Function to see if Z-axis can be edited (3D) or not (2D)
@pn.depends(PCA_params.controls.widgets['n_dimensions'].param.value, watch=True)
def _update_PCz_disabled(dimensions):
    if dimensions == '2 Components':
        PCA_params.controls.widgets['PCz'].disabled = True
    elif dimensions == '3 Components':
        PCA_params.controls.widgets['PCz'].disabled = False


# Function updating the possible components to choose based on number of components of the PCA and updating the explained variance plot
@pn.depends(PCA_params.controls.widgets['n_components'].param.value, watch=True)
def _update_PC_options(components):
    "Controlling options and plots based on the number of components PCA was computed with."
    # Controlling widget options for components to show in the projection plot
    PCA_params.controls.widgets['PCx'].options = ['PC '+str(i+1) for i in range(components)]
    PCA_params.controls.widgets['PCy'].options = ['PC '+str(i+1) for i in range(components)]
    PCA_params.controls.widgets['PCz'].options = ['PC '+str(i+1) for i in range(components)]

    # Updating the explained variance figure with thee number of components PCA was computed
    PCA_params.exp_var_fig_plot[0] = iaf._plot_PCA_explained_variance(PCA_params)
    end_page_PCA[0] = pn.pane.Plotly(PCA_params.exp_var_fig_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_exp_var_plot', 'scale':4,}})

    # Control the many PCA scatter plots to include less than 6 components if PCA was computed with less than 6 components
    if len(PCA_params.scatter_PCA_plot[0].data[0]['dimensions']) < 6:
        PCA_params.scatter_PCA_plot[0] = iaf._scatter_PCA_plot(PCA_params, target_list)
        end_page_PCA[1] = pn.pane.Plotly(PCA_params.scatter_PCA_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_scatter_plot', 'scale':4,}})
    else:
        if components < 6:
            PCA_params.scatter_PCA_plot[0] = iaf._scatter_PCA_plot(PCA_params, target_list)
            end_page_PCA[1] = pn.pane.Plotly(PCA_params.scatter_PCA_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_scatter_plot', 'scale':4,}})


# Function enabling/disabling the confidence level parameters for ellipses
@pn.depends(PCA_params.controls.widgets['ellipse_draw'].param.value,
            watch=True)
def _update_ellipse_options(ellipse):
    if ellipse:
        PCA_params.controls.widgets['confidence'].disabled = False
        if PCA_params.controls.widgets['confidence'].value == 0:
            PCA_params.controls.widgets['confidence_std'].disabled = False
    else:
        PCA_params.controls.widgets['confidence'].disabled = True
        PCA_params.controls.widgets['confidence_std'].disabled = True


# Function enabling/disabling the confidence level based on std parameter for ellipses
@pn.depends(PCA_params.controls.widgets['confidence'].param.value,
            watch=True)
def _update_ellipse_std_options(confidence):
    if confidence == 0:
        PCA_params.controls.widgets['confidence_std'].disabled = False
    else:
        PCA_params.controls.widgets['confidence_std'].disabled = True


# First section of the page
initial_page_PCA = pn.GridSpec(mode='override')
initial_page_PCA[0,:2] = n_components_compute
initial_page_PCA[0,2] = compute_PCA_button

# Middle section of the page
middle_page_PCA = pn.GridSpec(mode='override')
middle_page_PCA[0,0] = PCA_params.controls
middle_page_PCA[0,1:3] = 'To plot a PCA'

# Final section of the page
end_page_PCA = pn.Column('To plot explained variance figure', 'To plot matrices of PCA projections')
page_PCA = pn.Column(initial_page_PCA, middle_page_PCA, end_page_PCA)



# Tab for HCA analysis
# Param Class to store parameters and data regarding HCA
# TODO: Each time you modify the HCA, a new matplotlib plot is created, thus this could lead to memory issues
class HCA_Storage(param.Parameterized):

    # Distance metric to compute HCA
    dist_metric = param.String(default='euclidean')

    # Linkage metric
    link_metric = param.String(default='ward')

    # Distance matrix and Z linkage matrix
    dists = param.Array()
    Z = param.Array()

    # Figure details
    fig_text = param.String('Figure Details')
    fig_x = param.Number(default=4)
    fig_y = param.Number(default=4)
    col_threshold = param.Number(default=0)
    dpi = param.Number(default=200)

    # Storing figure
    HCA_plot = param.List(default=['Pane to Plot a HCA Dendrogram'])

    # Update the HCA plot
    @param.depends('dist_metric', 'link_metric', watch=True)
    def _update_HCA_compute_plot(self):
        # Recompute Distance and Linkage matrices
        self.dists = dist.pdist(DataFrame_Store.treated_df, metric=self.dist_metric)
        self.Z = hier.linkage(self.dists, method=self.link_metric)
        # Plot the new HCA
        self.HCA_plot[0] = _plot_HCA()
        page_HCA[0:6,1:4] = pn.pane.Matplotlib(self.HCA_plot[0], dpi=self.dpi)

    @param.depends('fig_text', 'fig_x', 'fig_y', 'col_threshold', 'dpi', watch=True)
    def _update_HCA_plot(self):
        self.HCA_plot[0] = _plot_HCA()
        page_HCA[0:6,1:4] = pn.pane.Matplotlib(self.HCA_plot[0], dpi=self.dpi)

    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets = {
            'dist_metric': pn.widgets.Select(name='Distance Metric:',
                                           value='euclidean', options=['euclidean', 'cityblock', 'minkowski', 'seuclidean',
                    'sqeuclidean', 'braycurtis', 'canberra', 'chebyshev', 'correlation', 'cosine', 'dice', 'hamming',
                    'jaccard', 'jensenshannon', 'kulczynski1', 'mahalanobis', 'matching', 'rogerstanimoto', 'russellrao',
                    'sokalmichener', 'sokalsneath', 'yule'],
                    description='Select distance metric. See details in https://docs.scipy.org/doc/scipy/reference/spatial.distance.html.'),
            'link_metric': pn.widgets.Select(name='Linkage:',
                                           value='ward', options=['ward', 'average', 'centroid', 'single', 'complete',
                                                                  'weighted', 'median'],
                    description='Select Linkage. See details in https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html.'),
            'fig_text': pn.widgets.StaticText(name="", value='Figure Details',
                                             styles={'font-size': 'large'}),
            'fig_x': pn.widgets.IntInput(name="X-axis", value=4, step=1, start=1, end=20, disabled=False,
                                    description='X axis length for HCA plot.'),
            'fig_y': pn.widgets.IntInput(name="Y-axis", value=4, step=1, start=1, end=20, disabled=False,
                                    description='Y axis length for HCA plot.'),
            'col_threshold': pn.widgets.IntInput(name="Color Threshold",
                                    value=0, step=1, start=0, disabled=False,
                                    description='''Select a distance threshold from where clusters are coloured.
                                    E.g. 70 would colour the branches that split under a distance of 70 in the same colour.'''),
            'dpi': pn.widgets.IntInput(name="DPI (Resolution)",
                                    value=200, step=10, start=100, disabled=False,
                                    description='Set the resolution of figure'),
        }
        self.controls = pn.Param(self,
                                 parameters=['dist_metric', 'link_metric', 'fig_text', 'fig_x', 'fig_y', 'col_threshold','dpi'],
                                 widgets=widgets, name='Parameters for HCA')

# Initializes HCA Store
HCA_params = HCA_Storage()

# Plots the HCA function
def _plot_HCA():
    f, ax = plt.subplots(1, 1, figsize=(HCA_params.fig_x, HCA_params.fig_y), constrained_layout=True)
    iaf.plot_dendogram(HCA_params.Z,
                   target_list.target, ax=ax,
                   label_colors=target_list.color_classes,
                   x_axis_len=HCA_params.fig_x,
                   color_threshold=HCA_params.col_threshold)
    return f


# Widget to save HCA plot (needed since it is a matplotlib plot instead of a plotly plot)
save_HCA_plot_button = pn.widgets.Button(name='Save as a png (in current folder)', button_type='success',
                                         icon=iaf.download_icon)
# When pressing the button, downloads the figure
def _save_HCA_plot_button(event):
    HCA_params.HCA_plot[0].savefig('HCA_plot.png', dpi=HCA_params.dpi)
    pn.state.notifications.success(f'Dendrogram successfully saved.')
save_HCA_plot_button.on_click(_save_HCA_plot_button)

# Organization for the HCA page
page_HCA = pn.GridSpec(mode='override')
page_HCA[0:5,0] = HCA_params.controls
page_HCA[0:6,1:4] = HCA_params.HCA_plot[0]
page_HCA[5,0] = save_HCA_plot_button

# Page with PCA and HCA analysis
unsup_analysis_page = pn.Tabs(('PCA', page_PCA), ('HCA', page_HCA))




# Page for Supervised Analysis

# TODO: Add ROC curves to both analyses

# Section for the PLS-DA analysis
plsda_opening_string = """Here, <strong>Partial Least Squares - Discriminant Analysis</strong> models can be built.
<br>
<br>
As the number of components (Latent Variables) to use in a model such as this one is of essential importance, first we will
start with an optimization of this number. This is performed by fitting a PLS model in your data from a chosen initial
number of starting components to a chosen final number of components. The <strong>Q<sup>2</sup> and R<sup>2</sup></strong>
scores are then estimated by stratified cross-validation (number of folds <strong>K</strong> is also chosen by the user).
The chosen number of components should be the number <strong>where Q<sup>2</sup> specifically</strong> stops increasing. A
figure is provided to evaluate this.
<br>
<br>
- <strong>Q<sup>2</sup></strong> - Mean squared error of PLS-DA predictions based on the test samples, thus it is ideal
to test if the model will overfit. This will increase until a certain number of components that should be chosen. Then
it usually stabilizes but from a certain point it might start to decrease which would mean the model is overfitting.
<br>
- <strong>R<sup>2</sup></strong> - Mean squared error of PLS-DA predictions based on the training samples used to train
the model (it will be higher than <strong>Q<sup>2</sup></strong> but it should not be used to choose the number of
components. This metric always increases with more components (estimation errorss might lead to temporary decreases)
which means it will overfit the model eventually.
<br>
<br>
After this optimization, a PLS model can be fitted with the chosen parameters and number of components. Model performance
can be estimated by accuracy, recall, precision and F1-score (the latter 3 weighted by class) using stratified fold
cross-validation. This cross-validation will also be used to estimate feature importance - variable importance in projection
(<strong>VIP</strong> scores). These are calculated based on Keiron Teilo O'Shea provided code in <a
href="https://www.researchgate.net/post/How-can-I-compute-Variable-Importance-in-Projection-VIP-in-Partial-Least-Squares-PLS"
target="_blank" rel="nofollow">
https://www.researchgate.net/post/How-can-I-compute-Variable-Importance-in-Projection-VIP-in-Partial-Least-Squares-PLS
</a>. This calculation can be slow when you have a high number of features. Importances estimated are then averaged across
the different folds. This can be iterated multiple times using randomized cross-validation folds. Model performance can
also be estimated with a Receiver Operating Characteristic (ROC) curve when only 2 classes are present.
<br>
<br>
Furthermore, a projection of the samples based on the PLS regression model using all available experimental samples on the
chosen Latent Variables akin to the PCA projections in the unsupervised analysis page is provided. Finally, Permutation
Testing (<strong>slow</strong>) can be performed to observe if the model and results obtained are significative.
"""


# Param Class to store parameters and data regarding PLS-DA
class PLSDA_Storage(param.Parameterized):
    "Class to store PLS-DA models, parameters and results."

    # PLS-DA Optimization
    n_min_max_components = param.Range(default=(2,20))
    n_fold = param.Number(default=5)
    scale = param.Boolean(default=False)
    confirm_optim_button = param.Boolean(default=False)
    optim_q2 = param.List(default=[''])
    optim_r2 = param.List(default=[''])
    rec_components = param.Number(default=np.nan)

    # PLS-DA model
    n_components = param.Number(default=5)
    n_iterations = param.Number(default=10)
    static_text_metrics = param.String(default='Choose what model performance metrics to use')
    metrics_to_use = param.List(default=['Accuracy', 'F1-Score (weighted)', 'Precision (weighted)', 'Recall (weighted)'])
    imp_feature_metric = param.String(default='VIP')
    confirm_plsda_button = param.Boolean(default=False)
    n_results = param.DataFrame()
    feat_impor = param.DataFrame()

    # Params for PLS Projection
    models = param.List(default=[''])
    x_scores = param.DataFrame()

    # Params Used Store
    current_plsda_params = param.Dict()

    # PLS Projection plot parameters
    n_dimensions = param.String(default='2 Components') # N dimensions to show in plot
    # Latent Variables to plot
    LVx = param.String(default='LV 1')
    LVy = param.String(default='LV 2')
    LVz = param.String(default='LV 3')
    # Draw Ellipses
    ellipse_draw = param.Boolean(default=True)
    confidence = param.Number(default=0.95)
    confidence_std = param.Number(default=2)
    # Extra params
    dot_size = param.Number(default=5)

    # Permutation Test
    n_perm = param.Number(default=500)
    perm_metric = param.String(default='Accuracy')
    dpi = param.Number(default=200)
    confirm_button_permutation = param.Boolean(default=False)
    save_figure_button_permutation = param.Boolean(default=False)
    current_plsda_params_permutation = param.Dict()

    # ROC Curve

    # Storing figures
    optim_figure = param.List(default=['To Plot the Optimization PLS Figure'])
    PLS_plot = param.List(default=['To Plot the PLS Projection Figure'])
    ROC_figure = param.List(default=['a'])
    perm_figure = param.List(default=['To Plot the Permutation Test Figure'])

    # Update number of folds limits
    def n_folds_limits(self, target_list):
        "Updates the limits of the cross-validation fold number based on the number of samples per class."
        # Updating the end value limits
        min_samples_in_class = pd.Series(target_list.target).value_counts().min()
        self.controls_optim.widgets['n_fold'].end = min_samples_in_class
        self.controls.widgets['n_fold'].end = min_samples_in_class

        # Updating the proper value if it is above the highest possible number
        if min_samples_in_class < self.controls_optim.widgets['n_fold'].value:
            self.controls_optim.widgets['n_fold'].value = min_samples_in_class
            self.controls.widgets['n_fold'].value = min_samples_in_class
            self.n_fold = min_samples_in_class


    # Function to confirm the optimization
    def _confirm_optim_button(self, event):
        "Performs optimization, plots the optimization figure and updated corresponding layout."

        # Loading Widget while the PLS-DA optimization is being performed
        pls_optim_section[1] = pn.indicators.LoadingSpinner(value=True, size=90,
                                                            name='Performing Optimization of PLS-DA Components...')

        # Round to integer the values in the range
        self.n_min_max_components = int(np.round(self.n_min_max_components[0])), int(np.round(self.n_min_max_components[1]))

        # Performs optimization
        PLS_optim = metsta.optim_PLSDA_n_components(DataFrame_Store.treated_df, target_list.target, # Data and target
                                encode2as1vector=True,
                                max_comp=self.n_min_max_components[1], # Max. number of components to search
                                min_comp=self.n_min_max_components[0], # Min. number of components to search
                                kf=None, n_fold=self.n_fold, # Nº of folds
                                scale=self.scale)

        # Transforms the results to be suited to be used in plotly
        # DataFrame with Q2 and R2 values in the same column, a column to identify which value is what and a
        # column with the number of components that led to that result
        pls_optim_values = pd.DataFrame(PLS_optim,
                                index=range(self.n_min_max_components[0], self.n_min_max_components[1]+1))
        pls_joined = pd.DataFrame(pd.concat((pls_optim_values['CVscores'], pls_optim_values['CVR2scores'])),
                                  columns=['Scores']) # Q2 and R2 values in the same column
        dif_max_min = self.n_min_max_components[1] - self.n_min_max_components[0] + 1
        pls_joined['Class'] = ['Q2',]*dif_max_min + ['R2',]*dif_max_min # Identification of Q2 or R2 value
        pls_joined['Nº of Components'] = pls_joined.index # Nº of Components column

        # Plot the figure
        rec_comp, fig = iaf._plot_PLS_optimization_components_fig(pls_joined)
        self.optim_figure[0] = fig
        self.rec_components = rec_comp
        # Update the layout
        pls_optim_section[1] = pn.pane.Plotly(self.optim_figure[0],
                                              config = {'toImageButtonOptions': {'filename': 'PLS_optim_plot', 'scale':4,}})
        pls_optim_section[0][1].value = self.rec_components


    # Function to fit the PLS-DA model and obtain results
    def _confirm_plsda_button(self, event):
        "Fits a PLS-DA model, gives model performance estimations and feature importance lists, updating the layout."

        # Loading Widget while the PLS-DA model is being fitted and assessed
        pls_results_section[0][1][1] = pn.indicators.LoadingSpinner(value=True, size=90,
                                                            name='Fitting and Assessing PLS-DA model...')

        # Sees what type of feature importance metric will be used
        if self.imp_feature_metric == 'VIP':
            feat_type = 'VIP'
        elif self.imp_feature_metric == 'Coefficients':
            feat_type = 'Coef'
        else:
            feat_type = 'Weights'

        # Fits a PLS-DA model under a stratified cross validation scheme
        PLSDA_results = metsta.PLSDA_model_CV(DataFrame_Store.treated_df, target_list.target, # Data and target
                       n_comp=self.n_components, # Number of components of PLS-DA model - very important
                       kf = None, n_fold=self.n_fold, # Nº of folds
                       iter_num=self.n_iterations, # Number of iterations of cross-validation to do
                       encode2as1vector=True,
                       scale=self.scale, # Set scale to True only if you did not do scaling in pre-treatments
                       feat_type=feat_type) # Feature Importance Metric to use, default is VIP scores

        # Exclude metrics that were not chosen to be shown
        if 'Accuracy' not in self.metrics_to_use:
            PLSDA_results.pop('accuracy')
        if 'F1-Score (weighted)' not in self.metrics_to_use:
            PLSDA_results.pop('F1-scores')
        if 'Precision (weighted)' not in self.metrics_to_use:
            PLSDA_results.pop('precision')
        if 'Recall (weighted)' not in self.metrics_to_use:
            PLSDA_results.pop('recall')

        # Results Table
        results_summary = pd.DataFrame(columns=['Value', 'Standard Deviation'])
        for k,v in PLSDA_results.items():
            if k != 'Q2' and k != 'imp_feat':
                results_summary.loc[k] = np.mean(v), np.std(v)
        self.n_results = results_summary

        # Important Feature Table
        self.feat_impor = iaf.creating_importance_feat_table(self.imp_feature_metric, DataFrame_Store, PLSDA_results['imp_feat'])

        # Store Parameters used for the model
        self.current_plsda_params = {'n_components': self.n_components, 'n_folds': self.n_fold,
                                     'n_iterations': self.n_iterations, 'scale': self.scale,
                                     'feat_imp': self.imp_feature_metric}

        # Update the layout
        pls_results_section[0][1][1] = pn.pane.DataFrame(self.n_results)
        pls_results_section[3] = pn.pane.DataFrame(self.feat_impor, height=600)
        plsda_feat_imp_show_annots_only.value = False
        save_plsda_feat_imp_button.disabled = False

        # Fit a PLS-DA model with all samples and store model and x_scores
        self.models[0], self.x_scores = ma.fit_PLSDA_model(
            DataFrame_Store.treated_df, target_list.target, # Data and target
            n_comp=self.n_components, return_scores=True,
            scale=self.scale,
            encode2as1vector=True, lv_prefix='LV ', label_name='Label')

        # Name of the file
        filename_string = f'PLS_plot_({self.n_components}comp)'
        if self.n_dimensions == '2 Components':
            if self.ellipse_draw:
                if self.confidence != 0:
                    filename_string = filename_string + f'_ellipse({self.confidence*100}%confidence)'
                else:
                    filename_string = filename_string + f'_ellipse({self.confidence_std}std)'

        self.PLS_plot[0] = iaf._plot_PLS(self, target_list)
        pls_proj_page[0,1:3] = pn.pane.Plotly(self.PLS_plot[0], config={'toImageButtonOptions': {'filename': filename_string, 'scale':4}})

        if len(pls_results_section) == 5:
            pls_results_section.append('### PLS Projection Section')
            pls_results_section.append(pls_proj_page)
            pls_results_section.append('### PLS-DA Permutation Test Section')
            pls_results_section.append(pn.Column(pn.pane.HTML(permutation_test_description),
                                                pn.pane.LaTeX(pvalue_equation_string, styles={'font-size': '14pt'})))
            pls_results_section.append(pn.Row(self.controls_permutation,
                                        self.perm_figure[0]))


    # Update the PLS Projection plot
    @param.depends('n_dimensions', 'LVx', 'LVy', 'LVz', 'ellipse_draw', 'confidence', 'confidence_std', 'dot_size', watch=True)
    def _update_PLS_plot(self):
        # Name of the file
        filename_string = f'PLS_plot_({self.n_components}comp)'
        if self.n_dimensions == '2 Components':
            if self.ellipse_draw:
                if self.confidence != 0:
                    filename_string = filename_string + f'_ellipse({self.confidence*100}%confidence)'
                else:
                    filename_string = filename_string + f'_ellipse({self.confidence_std}std)'

        self.PLS_plot[0] = iaf._plot_PLS(self, target_list)
        pls_proj_page[0,1:3] = pn.pane.Plotly(self.PLS_plot[0], config={'toImageButtonOptions': {'filename': filename_string, 'scale':4}})


    # Function to confirm the permutation test to perform
    def _confirm_button_permutation(self, event):
        "Performs a permutation test for the PLS-DA model and updates the layout."

        # Loading Widget while the Permutation Test is being performed
        pls_results_section[9][1] = pn.indicators.LoadingSpinner(value=True, size=90,
                                                            name='Performing Permutation Test... (It can take a while)')

        # Sees what type of feature importance metric will be used
        if self.perm_metric == 'Accuracy':
            perm_metric = 'accuracy'
        elif self.perm_metric == 'F1-score (weighted)':
            perm_metric = 'f1_weighted'
        elif self.perm_metric == 'Precision (weighted)':
            perm_metric = 'precision_weighted'
        else:
            perm_metric = 'recall_weighted'

        # Performs the Permutation Test on the PLS-DA model under a stratified cross validation scheme
        perm_results_PLSDA = metsta.permutation_PLSDA(
            DataFrame_Store.treated_df, target_list.target, # Data and target
            n_comp=self.n_components, # Number of components
            iter_num=self.n_perm, # Nº of permutations to do in your test - around 500 should be enough
            cv=None, n_fold=self.n_fold, # Choose the number of folds
            random_state=None, # Random seed given to make the permutations rng class labels
            encode2as1vector=True, scale=self.scale, # Set scale to True only if you did not do scaling in pre-treatments
            metric=perm_metric) # Choose a metric to use to evaluate if the model is significant

        # Store Parameters used for the model
        self.current_plsda_params_permutation = {'n_components': self.n_components, 'n_folds': self.n_fold,
                                     'n_permutations': self.n_perm, 'scale': self.scale,
                                     'perm_metric': self.perm_metric}

        # Plot the permutation test
        self.perm_figure[0] = iaf._plot_permutation_test(perm_results_PLSDA, DataFrame_Store, self.n_fold,
                                                     self.perm_metric, 'PLS-DA Permutation Test')

        # Update the layout
        pls_results_section[9][1] = pn.pane.Matplotlib(self.perm_figure[0], height=600)
        self.controls_permutation.widgets['save_figure_button_permutation'].disabled = False


    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets_optim = {
            'n_min_max_components': pn.widgets.EditableRangeSlider(
                name='Range of Components to test for optimization',
                start=1, end=30, fixed_start=1, fixed_end=40, value=(2, 20), step=1),
            'n_fold': pn.widgets.IntInput(name="Number of folds for stratified cross-validation",
                value=5, start=2, end=20,
                description='''Value cannot be higher than the number of samples of your lowest sample number class.
                Usual values are 5, 7 or 10. A good value would be a divisor of the number of samples in each of the classes.'''),
            'scale': pn.widgets.Checkbox(name='Check if you did not perform scaling in the pre-treatment section. Performs auto scale.',
                                        value=False),
            'confirm_optim_button': pn.widgets.Button(name="Confirm Optimization Parameters", button_type='primary')
        }

        widgets = {
            'n_components': pn.widgets.IntSlider(name='Number of Components for PLS-DA model',
                start=1, end=40, value=5, step=1, styles={'font-weight': 'bold'}),
            'n_iterations': pn.widgets.IntInput(name='Number of Times to repeat analysis',
                start=1, value=10, step=1),
            'n_fold': pn.widgets.IntInput(name="Number of folds for stratified cross-validation",
                value=5, start=2, end=20,
                description='''Value cannot be higher than the number of samples of your lowest sample number class.
                Usual values are 5, 7 or 10. A good value would be a divisor of the number of samples in each of the classes.'''),
            'scale': pn.widgets.Checkbox(name='Check if you did not perform scaling in the pre-treatment section. Performs auto scale.',
                                        value=False),
            'static_text_metrics': pn.widgets.StaticText(name='', value='Choose what model performance metrics to use'),
            'metrics_to_use': pn.widgets.CheckBoxGroup(name='Choose what metrics to use',
                value=['Accuracy', 'F1-Score (weighted)', 'Precision (weighted)', 'Recall (weighted)'],
                options=['Accuracy', 'F1-Score (weighted)', 'Precision (weighted)', 'Recall (weighted)'],
                inline=False),
            #'static_text_imp_feat_metric': pn.widgets.StaticText(name='', value='Choose what feature importance metric to use'),
            'imp_feature_metric': pn.widgets.Select(name='Choose what feature importance metric to use', value='VIP',
                options=['VIP', 'Coefficients', 'X-Weights'], disabled_options=['Coefficients'],
                description='VIP is the most used metric but it also is the slowest to compute (Coefficients option currently not working).'),
            'confirm_plsda_button': pn.widgets.Button(name="Fit PLS-DA model", button_type='primary')
        }

        widgets_PLS_proj = {'n_components': pn.widgets.IntInput(name='Number of Components for PLS-DA model',
                start=1, end=40, value=5, step=1, disabled=True,
                description='This parameter was chosen when fitting the PLS-DA model and assessing its performance'),
            'n_dimensions': pn.widgets.RadioBoxGroup(name="n_dimensions",
                value = '2 Components', options=['2 Components', '3 Components'], inline = True),
            'LVx': pn.widgets.Select(name="Latent Variable in X axis",
                value='LV 1', options=['LV '+str(i+1) for i in range(self.n_components)]),
            'LVy': pn.widgets.Select(name="Latent Variable in Y axis",
                value='LV 2', options=['LV '+str(i+1) for i in range(self.n_components)]),
            'LVz': pn.widgets.Select(name="Latent Variable in Z axis",
                value='LV 3', options=['LV '+str(i+1) for i in range(self.n_components)], disabled=True),
            'ellipse_draw': pn.widgets.Checkbox(name="Draw Confidence Ellipses (Only for 2-D)",
                                    value=True),
            'confidence': pn.widgets.FloatInput(name="Ellipse Confidence Level",
                                    value=0.95, step=0.01, start=0, end=1, disabled=False,
                                    description='E.g. 0.95 confidence level means a 95% confidence level of the sigma error ellipse based on the covariance matrix'),
            'confidence_std': pn.widgets.FloatInput(name="Ellipse Confidence Level (standard deviation unit)",
                                    value=2, step=0.50, start=1, end=5,
                description='''Only used if Level of Confidence of Ellipse is 0.
                E.g. 1 stands for 68.3% and 2 for 95.4% confidence level of the sigma error ellipse based on the covariance matrix''',
                                    disabled=True),
            'dot_size': pn.widgets.IntSlider(name='Size of Points in Projection', start=1, end=20, value=5, step=1),
        }

        widgets_PLS_perm = {'n_components': pn.widgets.IntInput(name='Number of Components for PLS-DA model',
                start=1, end=40, value=5, step=1, disabled=True,
                description='This parameter was chosen when fitting the PLS-DA model and assessing its performance'),
            'n_fold': pn.widgets.IntInput(name="Number of folds for stratified cross-validation",
                value=5, start=2, end=20, disabled=True,
                description='This parameter was chosen when fitting the PLS-DA model and assessing its performance'),
            'scale': pn.widgets.Checkbox(name='Check if you did not perform scaling in the pre-treatment section. Performs auto scale.',
                value=False, disabled=True),
            'n_perm': pn.widgets.EditableIntSlider(name="Nº of Permutations to perform", start=100, value=500, end=2000,
                step=100),
            'perm_metric': pn.widgets.Select(name="Model Performance Metric", value='Accuracy',
                options=['Accuracy', 'F1-score (weighted)', 'Precision (weighted)', 'Recall (weighted)'],
                description='''Model performanced is estimated by default with **accuracy**.
                However, if your data is imbalanced (that is, at least one class has much fewer samples than another), accuracy is not the best metric.
                Consider using another metric such as the **F1 score**.'''),
            'dpi': pn.widgets.IntInput(name="DPI (Resolution)", value=200, step=10, start=100, disabled=False,
                                    description='Set the resolution of Permutation Test Figure'),
            'confirm_button_permutation': pn.widgets.Button(name="Perform Permutation Test (Slow)", button_type='primary', icon=iaf.hourglass_icon),
            'save_figure_button_permutation': pn.widgets.Button(name="Save as a png (in current folder)",
                button_type='success', icon=iaf.download_icon, disabled=True),
        }

        self.controls = pn.Param(self, parameters=['n_components', 'n_iterations', 'n_fold', 'scale', 'static_text_metrics',
                                                   'metrics_to_use', 'imp_feature_metric', 'confirm_plsda_button'],
                                 widgets=widgets, name='Parameters for PLS-DA model fitting')

        self.controls_optim = pn.Param(self, parameters=['n_min_max_components', 'n_fold', 'scale', 'confirm_optim_button'],
                                 widgets=widgets_optim, name='Parameters for optimization of PLS-DA nº of components')

        self.controls_projection = pn.Param(self, parameters=['n_components', 'n_dimensions', 'LVx', 'LVy', 'LVz',
                                                             'ellipse_draw', 'confidence', 'confidence_std', 'dot_size'],
                                 widgets=widgets_PLS_proj, name='Parameters for PLS Projection plot')

        self.controls_permutation = pn.Param(self, parameters=['n_components', 'n_fold', 'scale', 'n_perm', 'perm_metric',
                                                               'dpi', 'confirm_button_permutation',
                                                               'save_figure_button_permutation'],
                                 widgets=widgets_PLS_perm, name='Parameters for PLS-DA Permutation Test')

# Running initial param to store PLSDA details
PLSDA_store = PLSDA_Storage()

# Click button to confirm PLS Optimization
PLSDA_store.controls_optim.widgets['confirm_optim_button'].on_click(PLSDA_store._confirm_optim_button)

# Click button to fit the PLS-DA model and obtain model performance metrics
PLSDA_store.controls.widgets['confirm_plsda_button'].on_click(PLSDA_store._confirm_plsda_button)

# Click button to perform PLS-DA Permutation Test
PLSDA_store.controls_permutation.widgets['confirm_button_permutation'].on_click(PLSDA_store._confirm_button_permutation)
# Function to save the Permutation Test figure as png
def _save_figure_button_permutation_PLSDA(event):
    filename_string = f'PLS-DA_permutation_test_{PLSDA_store.current_plsda_params_permutation["n_permutations"]}perm_'
    filename_string = filename_string + f'{PLSDA_store.current_plsda_params_permutation["n_components"]}comp_'
    filename_string = filename_string + f'{PLSDA_store.current_plsda_params_permutation["n_folds"]}-foldstratCV_scale'
    filename_string = filename_string + f'{PLSDA_store.current_plsda_params_permutation["scale"]}_metric'
    filename_string = filename_string + f'{PLSDA_store.current_plsda_params_permutation["perm_metric"]}'
    PLSDA_store.perm_figure[0].savefig(filename_string+'.png', dpi=PLSDA_store.dpi)
    pn.state.notifications.success(f'Figure {filename_string} successfully saved.')
# Click button to save the aforementioned figure
PLSDA_store.controls_permutation.widgets['save_figure_button_permutation'].on_click(_save_figure_button_permutation_PLSDA)

# Widget to add recommended number of components to page
rec_comp_indicator_widget = pn.indicators.Number(name='Recommended Components (based on max. Q2)', font_size='14pt', title_size='14pt',
                                    value=PLSDA_store.rec_components)

# Organizing the optimization section of the page
pls_optim_section = pn.Row(pn.Column(PLSDA_store.controls_optim, rec_comp_indicator_widget,
                                    'A much lower number of components with a similar Q2 may be preferable than the number shown'),
                           PLSDA_store.optim_figure[0])


# Results Section of the PLS-DA page
# Specific Widget for PLS results section of the page, shows DataFrame with only annotated metabolites or all metabolites
plsda_feat_imp_show_annots_only = pn.widgets.Checkbox(name='Only show annotated metabolites in feature importance table',
                                                       value=False)

# Change the DataFrame shown based on checkbox
@pn.depends(plsda_feat_imp_show_annots_only.param.value, watch=True)
def _layout_plsda_feat_import_dataframe(plsda_feat_imp_show_annots_only):
    "Update the layout based on if we are showing all metabolites or only annotated ones."
    # Select DataFrame
    if plsda_feat_imp_show_annots_only:
        df_to_show = PLSDA_store.feat_impor[PLSDA_store.feat_impor['Has Match?']]
    else:
        df_to_show = PLSDA_store.feat_impor

    # Update the layout
    pls_results_section[3] = pn.pane.DataFrame(df_to_show, height=600)

# Widget to save dataframe with features ordered by importance
save_plsda_feat_imp_button = pn.widgets.Button(name='Save PLS-DA Feature Importance table obtained as .xlsx (in current folder)',
                                                button_type='warning', icon=iaf.download_icon, disabled=True)

# When pressing the button, downloads the dataframe (builds the appropriate filename)
def _save_plsda_feat_imp_button(event):
    "Save PLS-DA Feature Importance results as an Excel."
    try:
        # Building the datafile name
        filename_string = f'PLS-DA_FeatImp_{PLSDA_store.imp_feature_metric}_model_params_components{PLSDA_store.n_components}'
        filename_string = filename_string + f'_{PLSDA_store.n_fold}-foldstratCV_iterations{PLSDA_store.n_iterations}'
        filename_string = filename_string + f'_scale{PLSDA_store.scale}.xlsx'

        # Saving the file
        PLSDA_store.feat_impor.to_excel(filename_string)
        pn.state.notifications.success(f'{filename_string} successfully saved.')
    except:
        pn.state.notifications.error(f'File could not be saved.')

save_plsda_feat_imp_button.on_click(_save_plsda_feat_imp_button)


# Set of Functions controlling parameters for PLS Projection
# Function to see if Z-axis can be edited (3D) or not (2D)
@pn.depends(PLSDA_store.controls_projection.widgets['n_dimensions'].param.value, watch=True)
def _update_LVz_disabled(dimensions):
    if dimensions == '2 Components':
        PLSDA_store.controls_projection.widgets['LVz'].disabled = True
    elif dimensions == '3 Components':
        PLSDA_store.controls_projection.widgets['LVz'].disabled = False


# Function updating the possible Latent variables to choose based on number of LVs of the PLS-DA
@pn.depends(PLSDA_store.controls_projection.widgets['n_components'].param.value, watch=True)
def _update_LV_options(components):
    PLSDA_store.controls_projection.widgets['LVx'].options = ['LV '+str(i+1) for i in range(components)]
    PLSDA_store.controls_projection.widgets['LVy'].options = ['LV '+str(i+1) for i in range(components)]
    PLSDA_store.controls_projection.widgets['LVz'].options = ['LV '+str(i+1) for i in range(components)]


# Function enabling/disabling the confidence level parameters for ellipses
@pn.depends(PLSDA_store.controls_projection.widgets['ellipse_draw'].param.value,
            watch=True)
def _update_plsda_ellipse_options(ellipse):
    if ellipse:
        PLSDA_store.controls_projection.widgets['confidence'].disabled = False
        if PLSDA_store.controls_projection.widgets['confidence'].value == 0:
            PLSDA_store.controls_projection.widgets['confidence_std'].disabled = False
    else:
        PLSDA_store.controls_projection.widgets['confidence'].disabled = True
        PLSDA_store.controls_projection.widgets['confidence_std'].disabled = True


# Function enabling/disabling the confidence level based on std parameter for ellipses
@pn.depends(PLSDA_store.controls_projection.widgets['confidence'].param.value,
            watch=True)
def _update_plsda_ellipse_std_options(confidence):
    if confidence == 0:
        PLSDA_store.controls_projection.widgets['confidence_std'].disabled = False
    else:
        PLSDA_store.controls_projection.widgets['confidence_std'].disabled = True


# Permutation Test associated widgets
# This description and equation will also be used for the Random Forest section
permutation_test_description = '''Permutation tests permutate  the class labels of the samples, that is, all classes will
be randomized while maintaining the same number of samples per class and classes. For each permutation, model performance is
assessed by stratified cross-validation exactly as the previous sections. These performances are then compared with the
performance of a single iteration of a model built based on the non-permutated model that appears with a red line. Thus,
this red line might not be exactly the same value as the average model performance obtained in the supervised model (PLS-DA
or Random Forest) fitting section.
<br>
<br>
Thus, this is a test to observe model performance significancy, that is, if it is better than a random model. If it is,
then the remaining results from the important features give meaningful information.
<br>
<br>
<strong>Note: Permutation tests are slow to be computed.</strong>
<br>
<br>
P-value is calculated using the following equation:'''
pvalue_equation_string = r'p-value = \(\frac{1 + \text{times permutated model has better performance than non-permutated model}}{\text{number of permutations}}\)'


# PLS projection section of the page
pls_proj_page = pn.GridSpec(mode='override')
pls_proj_page[0,0] = PLSDA_store.controls_projection
pls_proj_page[0,1:3] = 'To plot a PLS'

# Layout of the full results section (Partial, more is added when fitting PLS-DA model)
pls_results_section = pn.Column(pn.Row(PLSDA_store.controls,
                                       pn.Column('### Model Performance Metrics', PLSDA_store.n_results)),
                                '### Feature Importance Table',
                                plsda_feat_imp_show_annots_only,
                               pn.pane.DataFrame(PLSDA_store.feat_impor),
                               save_plsda_feat_imp_button,)

# Layout of the PLS-DA page
page_PLSDA = pn.Column(pn.pane.HTML(plsda_opening_string), pls_optim_section, pls_results_section)



# Section for the Random Forest analysis
rf_opening_string = """Here, <strong>Random Forest (RF)</strong> models can be built.
<br>
<br>
First, a focus on the optimization of the <strong>number of trees</strong> to use to build the Random Forest models is
available. This will allow the user to choose an optimal number of trees by observing the estimated accuracy of Random
Forest models (by stratified cross validation - number of folds <strong>K</strong> chosen by the user). The usual pattern
is that the accuracy will increase until a certain number of trees and will then fluctuate around that accuracy. Since RF
are somewhat resistant to overfitting, the value chosen will be one where the accuracy has stop increasing but not too high
as to make the model training too slow. A usual number used in our lab is 200.
<br>
<br>
After this optimization, a RF model can be fitted with the chosen parameters and number of trees. Model performance
can be estimated by accuracy, recall, precision and F1-score (the latter 3 weighted by class) using stratified fold
cross-validation. This cross-validation will also be used to estimate feature importance - <strong>Gini Importance</strong>.
Importances estimated are then averaged across the different folds. This can be iterated multiple times using randomized
cross-validation folds. Model performance can also be estimated with a Receiver Operating Characteristic (ROC) curve when
only 2 classes are present.
<br>
<br>
Finally, Permutation Testing (<strong>slow</strong>) can be performed to observe if the model and results obtained are
significative.
"""


# Param Class to store parameters and data regarding Random Forests
class RF_Storage(param.Parameterized):
    "Class to store Random Forest models, parameters and results."

    # RF Optimization
    n_min_max_trees = param.Range(default=(20,300))
    n_interval = param.Number(default=5, doc='Test')
    n_fold = param.Number(default=5)
    confirm_optim_button = param.Boolean(default=False)
    optim_scores = param.List(default=[''])
    optim_ntrees = param.List(default=[''])

    # RF model
    n_trees = param.Number(default=200)
    n_iterations = param.Number(default=10)
    static_text_metrics = param.String(default='Choose what model performance metrics to use')
    metrics_to_use = param.List(default=['Accuracy', 'F1-Score (weighted)', 'Precision (weighted)', 'Recall (weighted)'])
    confirm_rf_button = param.Boolean(default=False)
    n_results = param.DataFrame()
    feat_impor = param.DataFrame()
    models = param.List(default=[''])

    # Params Used Store
    current_rf_params = param.Dict()

    # Permutation Test
    n_perm = param.Number(default=500)
    perm_metric = param.String(default='Accuracy')
    dpi = param.Number(default=200)
    confirm_button_permutation = param.Boolean(default=False)
    save_figure_button_permutation = param.Boolean(default=False)
    current_rf_params_permutation = param.Dict()

    # ROC Curve

    # Storing figures
    optim_figure = param.List(default=['To Plot the Optimization RF Figure'])
    ROC_figure = param.List(default=['a'])
    perm_figure = param.List(default=['To Plot the Permutation Test Figure'])

    # Update number of folds limits
    def n_folds_limits(self, target_list):
        "Updates the limits of the cross-validation fold number based on the number of samples per class."
        # Updating the end value limits
        min_samples_in_class = pd.Series(target_list.target).value_counts().min()
        self.controls_optim.widgets['n_fold'].end = min_samples_in_class
        self.controls.widgets['n_fold'].end = min_samples_in_class

        # Updating the proper value if it is above the highest possible number
        if min_samples_in_class < self.controls_optim.widgets['n_fold'].value:
            self.controls_optim.widgets['n_fold'].value = min_samples_in_class
            self.controls.widgets['n_fold'].value = min_samples_in_class
            self.n_fold = min_samples_in_class


    # Function to confirm the optimization
    def _confirm_optim_button(self, event):
        "Performs optimization, plots the optimization figure and updated corresponding layout."

        # Loading Widget while the Random Forst optimization is being performed
        rf_optim_section[1] = pn.indicators.LoadingSpinner(value=True, size=90,
                                                            name='Performing Optimization of RF nº of trees...')

        # Round to integer the values in the range
        self.n_min_max_trees = int(np.round(self.n_min_max_trees[0])), int(np.round(self.n_min_max_trees[1]))

        # Perform Optimization
        rf_optim = iaf._optimization_n_trees_rf(self, DataFrame_Store.treated_df, target_list.target)
        self.optim_scores = list(rf_optim['mean_test_score'])
        self.optim_ntrees = list(rf_optim['param_n_estimators'])

        # RF optimization Figure
        rf_optim_results = pd.DataFrame([self.optim_ntrees, self.optim_scores],
                                       index=['Number of Trees', 'Model Accuracy (estimated by CV)']).T
        self.optim_figure[0] = px.line(rf_optim_results, x='Number of Trees', y='Model Accuracy (estimated by CV)',
                title='Random Forest Optimization Plot', range_y=(0, 1.05),
                range_x=(self.n_min_max_trees[0] - 10, self.n_min_max_trees[1] + 10)
                 )

        # Update the layout
        filename_string = f'RF_optim_plot_{self.n_fold}-stratCV'
        rf_optim_section[1] = pn.pane.Plotly(self.optim_figure[0],
                                              config = {'toImageButtonOptions': {'filename': filename_string, 'scale':4,}})


    # Function to fit the RF model and obtain results
    def _confirm_rf_button(self, event):
        "Fits a Random Forest model, gives model performance estimations and feature importance lists, updating the layout."

        # Loading Widget while the Random Forest model is being fitted and assessed
        rf_results_section[0][1][1] = pn.indicators.LoadingSpinner(value=True, size=90,
                                                            name='Fitting and Assessing RF model...')

        metrics_to_use = []

        # Include the metrics chosen to be shown
        if 'Accuracy' in self.metrics_to_use:
            metrics_to_use.append('accuracy')
        if 'F1-Score (weighted)' in self.metrics_to_use:
            metrics_to_use.append('f1_weighted')
        if 'Precision (weighted)' in self.metrics_to_use:
            metrics_to_use.append('precision_weighted')
        if 'Recall (weighted)' in self.metrics_to_use:
            metrics_to_use.append('recall_weighted')
        if len(self.metrics_to_use) == 0:
            metrics_to_use.append('accuracy')

        # Fits a RF model under a stratified cross validation scheme
        RF_results = metsta.RF_model(DataFrame_Store.treated_df, target_list.target, # Data and labels
                         return_cv=True, iter_num=self.n_iterations, # Number of iterations for it
                         n_trees=self.n_trees, # Number of trees in the model
                         cv=None, n_fold=self.n_fold, # Number of folds
                         metrics = metrics_to_use) # Choose the performance metrics
        if len(self.metrics_to_use) == 0:
            RF_results.pop('accuracy')

        # Results Table
        results_summary = pd.DataFrame(columns=['Value', 'Standard Deviation'])
        for k,v in RF_results.items():
            if k != 'model' and k != 'imp_feat':
                results_summary.loc[k] = np.mean(v), np.std(v)
        self.n_results = results_summary

        # Important Feature Table
        self.feat_impor = iaf.creating_importance_feat_table('Gini Importance', DataFrame_Store, RF_results['imp_feat'])

        # Store model
        self.models[0] = RF_results['model']

        # Store Parameters used for the model
        self.current_rf_params = {'n_trees': self.n_trees, 'n_folds': self.n_fold, 'n_iterations': self.n_iterations}

        # Update the layout
        rf_results_section[0][1][1] = pn.pane.DataFrame(self.n_results)
        rf_results_section[3] = pn.pane.DataFrame(self.feat_impor, height=600)
        rf_feat_imp_show_annots_only.value = False
        save_rf_feat_imp_button.disabled = False

        if len(rf_results_section) == 5:
            rf_results_section.append('### Random Forest Permutation Test Section')
            rf_results_section.append(pn.Column(pn.pane.HTML(permutation_test_description),
                                                 pn.pane.LaTeX(pvalue_equation_string, styles={'font-size': '14pt'})))
            rf_results_section.append(pn.Row(self.controls_permutation,
                                          self.perm_figure[0]))


    # Function to confirm the permutation test to perform
    def _confirm_button_permutation(self, event):
        "Performs a permutation test for the Random Forest model and updates the layout."

        # Loading Widget while the Permutation Test is being performed
        rf_results_section[7][1] = pn.indicators.LoadingSpinner(value=True, size=90,
                                                            name='Performing Permutation Test... (It can take a while)')

        # Sees what type of feature importance metric will be used
        if self.perm_metric == 'Accuracy':
            perm_metric = 'accuracy'
        elif self.perm_metric == 'F1-score (weighted)':
            perm_metric = 'f1_weighted'
        elif self.perm_metric == 'Precision (weighted)':
            perm_metric = 'precision_weighted'
        elif self.perm_metric == 'Recall (weighted)':
            perm_metric = 'recall_weighted'

        # Performs the Permutation Test on the Random Forest model under a stratified cross validation scheme
        perm_results_RF = metsta.permutation_RF(
            DataFrame_Store.treated_df, target_list.target, # Data and target
            iter_num=self.n_perm, # Nº of permutations to do in your test - around 500 should be enough
            n_trees=self.n_trees, # Number of trees in the model
            cv=None, n_fold=self.n_fold, # Choose the number of folds
            random_state=None, # Random seed given to make the permutations rng class labels
            metric=perm_metric) # Choose a metric to use to evaluate if the model is significant

        # Store Parameters used for the model
        self.current_rf_params_permutation = {'n_trees': self.n_trees, 'n_folds': self.n_fold,
                                     'n_permutations': self.n_perm, 'perm_metric': self.perm_metric}

        # Plot the permutation test
        self.perm_figure[0] = iaf._plot_permutation_test(perm_results_RF, DataFrame_Store, self.n_fold,
                                                     self.perm_metric, 'Random Forest Permutation Test')

        # Update the layout
        rf_results_section[7][1] = pn.pane.Matplotlib(self.perm_figure[0], height=600)
        self.controls_permutation.widgets['save_figure_button_permutation'].disabled = False


    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets_optim = {
            'n_min_max_trees': pn.widgets.EditableRangeSlider(
                name='Range of Trees to test for optimization', format='0',
                start=1, end=600, fixed_start=1, value=(20, 300), step=1),
            'n_interval': pn.widgets.IntInput(name='Step of number of trees for Optimization', start=1, end=20, value=5, step=1,
                description='''E.g. If 5, then the optimization is performed with the minimum number of trees, then that number + 5 and again until reaching the maximum number of trees.
                It should be a divisor of the subtraction between the maximum and minimum number of trees allowed for the optimization.'''),
            'n_fold': pn.widgets.IntInput(name="Number of folds for stratified cross-validation",
                value=5, start=2, end=20,
                description='''Value cannot be higher than the number of samples of your lowest sample number class.
                Usual values are 5, 7 or 10. A good value would be a divisor of the number of samples in each of the classes.'''),
            'confirm_optim_button': pn.widgets.Button(name="Confirm Optimization Parameters", button_type='primary')
        }

        widgets = {
            'n_trees': pn.widgets.IntSlider(name='Number of Trees for Random Forest model',
                start=1, end=600, value=200, step=1, styles={'font-weight': 'bold'}),
            'n_iterations': pn.widgets.IntInput(name='Number of Times to repeat analysis',
                start=1, value=10, step=1),
            'n_fold': pn.widgets.IntInput(name="Number of folds for stratified cross-validation",
                value=5, start=2, end=20,
                description='''Value cannot be higher than the number of samples of your lowest sample number class.
                Usual values are 5, 7 or 10. A good value would be a divisor of the number of samples in each of the classes.'''),
            'static_text_metrics': pn.widgets.StaticText(name='', value='Choose what model performance metrics to use'),
            'metrics_to_use': pn.widgets.CheckBoxGroup(name='Choose what metrics to use',
                value=['Accuracy', 'F1-Score (weighted)', 'Precision (weighted)', 'Recall (weighted)'],
                options=['Accuracy', 'F1-Score (weighted)', 'Precision (weighted)', 'Recall (weighted)'],
                inline=False),
            'confirm_rf_button': pn.widgets.Button(name="Fit Random Forest model", button_type='primary')
        }

        widgets_perm = {'n_trees': pn.widgets.IntInput(name='Number of Trees for Random Forest model',
                start=1, end=600, value=200, step=1, disabled=True,
                description='This parameter was chosen when fitting the Random Forest model and assessing its performance'),
            'n_fold': pn.widgets.IntInput(name="Number of folds for stratified cross-validation",
                value=5, start=2, end=20, disabled=True,
                description='This parameter was chosen when fitting the Random Forest model and assessing its performance'),
            'n_perm': pn.widgets.EditableIntSlider(name="Nº of Permutations to perform", start=100, value=500, end=2000,
                step=100),
            'perm_metric': pn.widgets.Select(name="Model Performance Metric", value='Accuracy',
                options=['Accuracy', 'F1-score (weighted)', 'Precision (weighted)', 'Recall (weighted)'],
                description='''Model performanced is estimated by default with **accuracy**.
                However, if your data is imbalanced (that is, at least one class has much fewer samples than another), accuracy is not the best metric.
                Consider using another metric such as the **F1 score**.'''),
            'dpi': pn.widgets.IntInput(name="DPI (Resolution)", value=200, step=10, start=100, disabled=False,
                                    description='Set the resolution of Permutation Test Figure'),
            'confirm_button_permutation': pn.widgets.Button(name="Perform Permutation Test (Slow)", button_type='primary', icon=iaf.hourglass_icon),
            'save_figure_button_permutation': pn.widgets.Button(name="Save as a png (in current folder)",
                button_type='success', icon=iaf.download_icon, disabled=True),
        }

        self.controls = pn.Param(self, parameters=['n_trees', 'n_iterations', 'n_fold', 'static_text_metrics',
                                                   'metrics_to_use', 'confirm_rf_button'],
                                 widgets=widgets, name='Parameters for Random Forest model fitting')

        self.controls_optim = pn.Param(self, parameters=['n_min_max_trees', 'n_interval', 'n_fold', 'confirm_optim_button'],
                                 widgets=widgets_optim, name='Parameters for optimization of nº of trees in Random Forest')

        self.controls_permutation = pn.Param(self, parameters=['n_trees', 'n_fold', 'n_perm', 'perm_metric',
                                                               'dpi', 'confirm_button_permutation',
                                                               'save_figure_button_permutation'],
                                 widgets=widgets_perm, name='Parameters for Random Forest Permutation Test')

# Running initial param to store RF details
RF_store = RF_Storage()

# Click button to confirm Random Forest Optimization
RF_store.controls_optim.widgets['confirm_optim_button'].on_click(RF_store._confirm_optim_button)

# Click button to fit the Random Forest model and obtain model performance metrics
RF_store.controls.widgets['confirm_rf_button'].on_click(RF_store._confirm_rf_button)

# Click button to perform Random Forest Permutation Test
RF_store.controls_permutation.widgets['confirm_button_permutation'].on_click(RF_store._confirm_button_permutation)
# Function to save the Permutation Test figure as png
def _save_figure_button_permutation_RF(event):
    filename_string = f'RF_permutation_test_{RF_store.current_rf_params_permutation["n_permutations"]}perm_'
    filename_string = filename_string + f'{RF_store.current_rf_params_permutation["n_trees"]}trees_'
    filename_string = filename_string + f'{RF_store.current_rf_params_permutation["n_folds"]}-foldstratCV_scale'
    filename_string = filename_string + f'{RF_store.current_rf_params_permutation["perm_metric"]}'
    RF_store.perm_figure[0].savefig(filename_string+'.png', dpi=RF_store.dpi)
    pn.state.notifications.success(f'Figure {filename_string} successfully saved.')
# Click button to save the aforementioned figure
RF_store.controls_permutation.widgets['save_figure_button_permutation'].on_click(_save_figure_button_permutation_RF)


# Optimization section of the page
# Organizing the optimization section of the page
rf_optim_section = pn.Row(RF_store.controls_optim, RF_store.optim_figure[0])


# Results Section of the Random Forest page
# Specific Widget for Random Forest results section, shows DataFrame with only annotated metabolites or all metabolites
rf_feat_imp_show_annots_only = pn.widgets.Checkbox(name='Only show annotated metabolites in feature importance table',
                                                       value=False)

# Change the DataFrame shown based on checkbox
@pn.depends(rf_feat_imp_show_annots_only.param.value, watch=True)
def _layout_rf_feat_import_dataframe(rf_feat_imp_show_annots_only):
    "Update the layout based on if we are showing all metabolites or only annotated ones."
    # Select DataFrame
    if rf_feat_imp_show_annots_only:
        df_to_show = RF_store.feat_impor[RF_store.feat_impor['Has Match?']]
    else:
        df_to_show = RF_store.feat_impor

    # Update the layout
    rf_results_section[3] = pn.pane.DataFrame(df_to_show, height=600)

# Widget to save dataframe with features ordered by importance
save_rf_feat_imp_button = pn.widgets.Button(name='Save Random Forest Feature Importance table obtained as .xlsx (in current folder)',
                                                button_type='warning', icon=iaf.download_icon, disabled=True)

# When pressing the button, downloads the dataframe (builds the appropriate filename)
def _save_rf_feat_imp_button(event):
    "Save Random Forest Feature Importance results as an Excel."
    try:
        # Building the datafile name
        filename_string = f'RF_FeatImp_Gini_model_params_{RF_store.n_trees}trees_{RF_store.n_fold}-foldstratCV_'
        filename_string = filename_string + f'iterations{RF_store.n_iterations}.xlsx'

        # Saving the file
        RF_store.feat_impor.to_excel(filename_string)
        pn.state.notifications.success(f'{filename_string} successfully saved.')
    except:
        pn.state.notifications.error(f'File could not be saved.')

save_rf_feat_imp_button.on_click(_save_rf_feat_imp_button)

# Layout of the full results section (Partial, more is added when fitting RF model)
rf_results_section = pn.Column(pn.Row(RF_store.controls,
                                       pn.Column('### Model Performance Metrics', RF_store.n_results)),
                                '### Feature Importance Table',
                                rf_feat_imp_show_annots_only,
                               pn.pane.DataFrame(RF_store.feat_impor),
                               save_rf_feat_imp_button,)


# Layout of the RF page
page_RF = pn.Column(pn.pane.HTML(rf_opening_string), rf_optim_section, rf_results_section)


# Layout of the full supervised analysis page
sup_analysis_page = pn.Tabs(('PLS-DA', page_PLSDA), ('RF ', page_RF))




# Page for Univariate Analysis
univ_opening_string = '''In this section, both Univariate Analysis and Fold-Change analysis are performed.
<br>
<br>
The Fold change is calculated in a dataset with missing values imputed and normalized after. <strong>This means that with
our very high number of missing values in FT-ICR-MS data, it affects the calculation of the fold change a lot. Thus, take
this fold changes values with a grain (or multiple grains that are actually more like rocks than grains) of salt.</strong>
<br>
<br>
Choose between the parametric <strong>t-test</strong> and non-parametric <strong>Mann-Whitney test</strong>.
<br>
<br>
<strong>Warning</strong>: This type of analysis is only done between 2 classes. If you have more than 2 classes, you can
also choose one as the control class and one as the test class to perform univariate analysis in the 1st tab.
<br>
<br>
Finally, when you have more than 2 classes, each unsueprvised analysis will select the samples respective to the 2 classes,
filter and pre-treatment them the same way it was done on the full dataset before performing unsupervised analysis.
'''

class UnivariateAnalysis_Store(param.Parameterized):
    """Class to store the unsupervised analysis and the locked in parameters used in data filtering and pre-treatment."""

    # Locked in attributes
    # Filtering
    filt_method = param.Selector(default="Total Samples")
    filt_kw = param.Number(default=2)

    # Pre-Treatment
    # Missing Value Imputation
    mvi_method = param.Selector(default="Minimum of Sample")
    mvi_kw = param.Number(default=0.2, bounds=(0,1))

    # Normalization
    norm_method = param.Selector(default="Reference Feature")
    norm_kw = param.Selector(default=None)

    # Transformation
    tf_method = param.Selector(default="Generalized Logarithmic Transformation (glog)")
    tf_kw = param.Number(default=None)

    # Scaling
    scaling_method = param.Selector(default="Pareto Scaling")
    scaling_kw = param.String(default="")

    # Unsupervised Analysis Main Parameters
    control_class = param.Selector(default='')
    test_class = param.Selector(default='')
    univariate_test_str = param.String('What Univariate Test to perform')
    univariate_test = param.Selector(default='T-Test (Non-Parametric)')
    expected_equal_variance = param.Boolean(default=True)
    p_value_threshold = param.Number(default=0.05)
    fold_change_threshold = param.Number(default=2)
    univariate_df = param.DataFrame()
    univariate_df_set = param.Dict()

    # Store Univariate results
    univariate_results = param.DataFrame()
    univariate_results_non_filt = param.DataFrame()
    univariate_results_set = param.Dict()

    # Confirm Univariate Test to Perform
    confirm_button_univariate = param.Boolean(default=False)

    # Volcano Plot Parameters
    color_non_sig = param.Color(default=mpl.colors.CSS4_COLORS['silver'])
    color_down_sig = param.Color(default=mpl.colors.CSS4_COLORS['deepskyblue'])
    color_up_sig = param.Color(default=mpl.colors.CSS4_COLORS['lightcoral'])
    Volcano_fig = param.List(default=['To Plot a Volcano Plot'])

    # Parameters and DataFrame of chosen sub-section of test_classes
    test_class_subset = param.List(default=list())
    show_annots_only = param.Boolean(default=False)
    specific_cl_df = param.DataFrame()
    inter_description = param.String(default='')


    def locking_filtering_params(self, filt_method_widget, filt_kw_widget):
        "Locking in parameters regarding to the data filtering made."
        self.filt_method = filt_method_widget.value
        self.filt_kw = filt_kw_widget.value


    def locking_pretreatment_params(self, PreTreatment_Method):
        "Locking in parameters regarding the data pre-treatment made."
        # Missing Value Imputation
        self.mvi_method = PreTreatment_Method.mvi_method
        self.mvi_kw = PreTreatment_Method.mvi_kw

        # Normalization
        self.norm_method = PreTreatment_Method.norm_method
        self.norm_kw = PreTreatment_Method.norm_kw

        # Transformation
        self.tf_method = PreTreatment_Method.tf_method
        self.tf_kw = PreTreatment_Method.tf_kw

        # Scaling
        self.scaling_method = PreTreatment_Method.scaling_method
        self.scaling_kw = PreTreatment_Method.scaling_kw


    def _update_widgets(self):
        "Update the needed widget values."
        self.controls.widgets['control_class'].options = target_list.classes
        self.controls.widgets['test_class'].options = target_list.classes
        if 'control' in target_list.classes:
            self.controls.widgets['control_class'].value = 'control'
            self.control_class = 'control'
        else:
            self.controls.widgets['control_class'].value = target_list.classes[0]
            self.control_class = target_list.classes[0]
        self.controls.widgets['test_class'].value = target_list.classes[1]
        self.test_class = target_list.classes[1]


    # Function to confirm and perform Univariate Analysis
    def _confirm_button_univariate(self, event):
        "Perform Univariate Analysis."

        # Initial check to see if control and test classes are different
        if self.control_class == self.test_class:
            pn.state.notifications.error('Test class is the same as Control class. Thus, analysis cannot be made.')
            return

        # Performing and storing results from Univariate Analysis
        a,b,c,d,e = iaf._perform_univariate_analysis(self, DataFrame_Store, target_list, filt_method, filt_kw)
        self.univariate_df, self.univariate_df_set = a, b
        self.univariate_results, self.univariate_results_non_filt, self.univariate_results_set = c, d, e

        # Joining the results of the univariate analysis to the metadata available
        self.univariate_results = pd.concat((self.univariate_results,
                                             DataFrame_Store.metadata_df.loc[self.univariate_results.index]),
                                            axis=1)

        # Updating the page layout
        _updating_univariate_analysis_page_layout()


    # Update the Volcano plot
    @param.depends('color_non_sig', 'color_down_sig', 'color_up_sig', watch=True)
    def _update_Volcano_plot(self):
        "Update the Volcano plot based on change in colours."
        # Repeated content from the main function to not add another attribute to the class
        results_df = UnivarA_Store.univariate_results_non_filt.copy()
        results_df['-log10(Adj. p-value)'] = -np.log10(results_df['FDR adjusted p-value'])

        expression = []
        for i in results_df.index:
            if i not in UnivarA_Store.univariate_results.index:
                expression.append('Non-Significant')
            elif results_df.loc[i, 'log2FC'] < 0:
                expression.append('Downregulated')
            else:
                expression.append('Upregulated')
        results_df['Expression'] = expression

        self.Volcano_fig[0] = iaf._plot_Volcano_plot(results_df, UnivarA_Store, )
        # Update the layout
        layout_volcano[0,1:3] = pn.pane.Plotly(UnivarA_Store.Volcano_fig[0], config={'toImageButtonOptions': {
                   'filename': f'VolcanoPlot - {UnivarA_Store.test_class}/{UnivarA_Store.control_class}', 'scale':4}})


    @param.depends('test_class_subset', 'show_annots_only', watch=True)
    def _update_intersections(self):
        "Update the intersections DataFrame and description."
        iaf._univariate_intersections(self, DataFrame_Store)
        layout_final_section[0][1] = UnivarA_Store.inter_description
        layout_final_section[1] = pn.pane.DataFrame(UnivarA_Store.specific_cl_df, height=600)


    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets = {
            'control_class': pn.widgets.Select(name='Control Class', options=target_list.classes),
            'test_class': pn.widgets.Select(name='Test Class', options=target_list.classes),
            'univariate_test_str': pn.widgets.StaticText(value='What Univariate Test to perform',
                                                         styles={'font-size': 'medium', 'font-weight': 'bold'}),
            'univariate_test': pn.widgets.RadioBoxGroup(name='What Univariate Test to perform',
                    value = 'T-Test (Non-Parametric)',
                    options=['T-Test (Non-Parametric)', 'Mann-Whitney Test (Parametric)']),
            'expected_equal_variance': pn.widgets.Checkbox(name='Consider Variance between classes as equal (only affects T-Test)', value=True),
            'p_value_threshold': pn.widgets.EditableFloatSlider(name='P-value threshold', start=0, end=1, value=0.05,
                    step=0.001),
            'fold_change_threshold': pn.widgets.FloatInput(name='Fold Change Threshold', value=2, step=0.1, start=1,
                    description='''The threshold is set so as only selecting features as significant if either the control class or the test class average is "chosen" threhold times higher than the opposing class.
                    0 will essentially skip this threshold and use only the p-value threshold. Setting 1 for p-value threshold has the same effect for the p-value step.
                    E.g: if 2, the average of a feature in control class samples has to be double or more that of the test class or vice-versa.'''),
            'confirm_button_univariate': pn.widgets.Button(name="Perform Univariate Analysis", button_type='primary'),
            'color_non_sig': pn.widgets.ColorPicker(name='Color Non Significant Metabolites',
                                                    value=mpl.colors.CSS4_COLORS['silver']),
            'color_down_sig': pn.widgets.ColorPicker(name='Color Downregulated Sig. Metabolites',
                                                    value=mpl.colors.CSS4_COLORS['deepskyblue']),
            'color_up_sig': pn.widgets.ColorPicker(name='Color Upregulated Sig. Metabolites',
                                                    value=mpl.colors.CSS4_COLORS['lightcoral']),
        }

        widgets2 = {'test_class_subset': pn.widgets.CheckBoxGroup(name='Test Classes', value=[],
                options=target_list.classes,inline=False),
                    'show_annots_only': pn.widgets.Checkbox(
                name='Only show annotated metabolites from the analysis', value=False),}

        self.controls = pn.Param(self, parameters=['control_class', 'test_class', 'univariate_test_str', 'univariate_test',
                                'expected_equal_variance', 'p_value_threshold','fold_change_threshold', 'univariate_test',
                                'confirm_button_univariate'],
                                 widgets=widgets, name='Univariate Analysis Parameters')
        self.volcano_colors = pn.Param(self, parameters=['color_non_sig', 'color_down_sig', 'color_up_sig'],
                                 widgets=widgets, name='Parameters for Volcano Plot')
        self.specific_df_controls = pn.Param(self, parameters=['test_class_subset', 'show_annots_only'],
                                 widgets=widgets2, name='Choose test classes to test against')


# Initializing Store for Univariate Analysis
UnivarA_Store = UnivariateAnalysis_Store()

# Click button to confirm Univariate Analysis
UnivarA_Store.controls.widgets['confirm_button_univariate'].on_click(UnivarA_Store._confirm_button_univariate)

# Specific Widget for middle section of the page, shows DataFrame with only annotated metabolites or all metabolites
univar_results_show_annots_only = pn.widgets.Checkbox(name='Only show annotated metabolites from the analysis', value=False)

# Change the DataFrame shown based on checkbox
@pn.depends(univar_results_show_annots_only.param.value, watch=True)
def _layout_df_dataframe(univar_results_show_annots_only):
    "Update the layout based on if we are showing all metabolites or only annotated ones."
    # Select DataFrame
    if univar_results_show_annots_only:
        df_to_show = UnivarA_Store.univariate_results[UnivarA_Store.univariate_results['Has Match?']]
    else:
        df_to_show = UnivarA_Store.univariate_results

    # Update the layout
    if len(layout_df) > 1:
        layout_df[2] = pn.pane.DataFrame(df_to_show, height=600)

# Widget to save dataframe of univariate analysis performed in .csv format
save_univariate_results_button = pn.widgets.Button(name='Save univariate analysis results as .csv (in current folder)',
                                                button_type='warning', icon=iaf.download_icon)

# When pressing the button, downloads the dataframe (builds the appropriate filename)
def _save_univariate_results_button(event):
    "Save univariate results performed"
    # Building the datafile name
    test_performed = UnivarA_Store.univariate_test.split(' ')[0]
    filename_string = f'Univar_res_{UnivarA_Store.test_class}_vs_{UnivarA_Store.control_class}_{test_performed}'
    filename_string = filename_string + f'_pvalue{UnivarA_Store.p_value_threshold}_FC{UnivarA_Store.fold_change_threshold}'
    if test_performed == 'T-Test':
        if UnivarA_Store.expected_equal_variance:
            filename_string = filename_string + f'_equalvariance.csv'
        else:
            filename_string = filename_string + f'_notequalvariance.csv'
    else:
        filename_string = filename_string + f'.csv'

    # Saving the file
    UnivarA_Store.univariate_results.to_csv(filename_string)
    pn.state.notifications.success(f'{filename_string} successfully saved.')

save_univariate_results_button.on_click(_save_univariate_results_button)

# Widget button to save dataframe shown when seeing intersections of multiple univariate analysis of
# different chosen test classes against a chosen control class
save_multiple_univariate_button = pn.widgets.Button(name='Save shown Dataframe with int. values as .csv file (in current folder)',
                                                button_type='warning', icon=iaf.download_icon)

# When pressing the button, downloads the dataframe (filename quite big)
def _save_multiple_univariate_button(event):
    # Building the datafile name
    # Based on what type of metabolites were chosen
    if UnivarA_Store.show_annots_only:
        annot_string = 'Annotated_metabolites'
    else:
        annot_string = 'All_metabolites'

    # Based on the test classes chosen
    class_subset_string = ''
    for cl in UnivarA_Store.test_class_subset:
        class_subset_string = class_subset_string + '_'+cl

    # Final name
    test_performed = UnivarA_Store.univariate_test.split(' ')[0]
    filename_string = annot_string + '_significant_in_testing_each_chosen_test_class' + class_subset_string + 'against'
    filename_string = filename_string + f'_control_{UnivarA_Store.control_class}_{test_performed}'
    filename_string = filename_string + f'_pvalue{UnivarA_Store.p_value_threshold}_FC{UnivarA_Store.fold_change_threshold}'
    if test_performed == 'T-Test':
        if UnivarA_Store.expected_equal_variance:
            filename_string = filename_string + f'_equalvariance.csv'
        else:
            filename_string = filename_string + f'_notequalvariance.csv'
    else:
        filename_string = filename_string + f'.csv'

    # Saving the file
    pd.concat((UnivarA_Store.specific_cl_df, DataFrame_Store.original_df[target_list.sample_cols]), join='inner', axis=1).to_csv(filename_string)
    pn.state.notifications.success(f'{filename_string} successfully saved.')

save_multiple_univariate_button.on_click(_save_multiple_univariate_button)


# Set up the different sections of the univariate analysis page
layout_df = pn.Column(univar_results_show_annots_only) # For DataFrame section
# For Volcano Plot section
layout_volcano = pn.GridSpec(mode='override')
layout_volcano[0,0] = UnivarA_Store.volcano_colors
layout_volcano[0,1:3] = UnivarA_Store.Volcano_fig[0]
# For intersection section
layout_final_section = pn.Row(pn.Column(UnivarA_Store.specific_df_controls, UnivarA_Store.inter_description, save_multiple_univariate_button),
       pn.pane.DataFrame(UnivarA_Store.specific_cl_df, height=600))

# Setting up the different pages
middle_page_univar = pn.Column('## DataFrame from the Univariate Analysis:',
                              layout_df,
                              '## Volcano Plot of the Univariate Analysis:',
                              layout_volcano,
                              '## Control Class versus all possible Test Classes (Intersections)',
                              layout_final_section)


def _updating_univariate_analysis_page_layout():
    "Updated the analysis of the univariate analysis page after unviariate analysis."

    # Updating the DataFrame section of the layout
    univar_results_show_annots_only.value = False

    n_sig_met = UnivarA_Store.univariate_results.shape[0]
    n_sig_annotated = UnivarA_Store.univariate_results[UnivarA_Store.univariate_results['Has Match?']].shape[0]
    if len(layout_df) == 1:
        layout_df.append(f'**{n_sig_met}** metabolites are significant, **{n_sig_annotated}** of which are annotated.')
        layout_df.append(pn.pane.DataFrame(UnivarA_Store.univariate_results, height=600))
        layout_df.append(save_univariate_results_button) # Button does not change - only needs to be added once
    else:
        layout_df[1] = f'**{n_sig_met}** metabolites are significant, **{n_sig_annotated}** of which are annotated.'
        layout_df[2] = pn.pane.DataFrame(UnivarA_Store.univariate_results, height=600)


    # Update the Volcano plot section of the layout
    # Generate and the Volcano plot data
    results_df = UnivarA_Store.univariate_results_non_filt.copy()
    results_df['-log10(Adj. p-value)'] = -np.log10(results_df['FDR adjusted p-value'])

    expression = []
    for i in results_df.index:
        if i not in UnivarA_Store.univariate_results.index:
            expression.append('Non-Significant')
        elif results_df.loc[i, 'log2FC'] < 0:
            expression.append('Downregulated')
        else:
            expression.append('Upregulated')
    results_df['Expression'] = expression

    # Plot the Volcano plot
    UnivarA_Store.Volcano_fig[0] = iaf._plot_Volcano_plot(results_df, UnivarA_Store)
    layout_volcano[0,1:3] = pn.pane.Plotly(UnivarA_Store.Volcano_fig[0], config={'toImageButtonOptions': {
                   'filename': f'VolcanoPlot - {UnivarA_Store.test_class}/{UnivarA_Store.control_class}', 'scale':4}})

    # Update the specific DF section of the layout
    UnivarA_Store.specific_df_controls.widgets['test_class_subset'].options = [
        i for i in target_list.classes if i != UnivarA_Store.control_class]
    UnivarA_Store.specific_df_controls.widgets['test_class_subset'].value = [UnivarA_Store.test_class,]
    UnivarA_Store.test_class_subset = [UnivarA_Store.test_class,]
    # See the intersections and update layout
    iaf._univariate_intersections(UnivarA_Store, DataFrame_Store)
    layout_final_section[0][1] = UnivarA_Store.inter_description
    layout_final_section[1] = pn.pane.DataFrame(UnivarA_Store.specific_cl_df, height=600)

    # Add the analysis section of the page in case it was not added yet
    if len(univar_analysis_page) == 2:
        univar_analysis_page.append(middle_page_univar)

# Setting up the initial univariate analysis page
univar_analysis_page = pn.Column(pn.pane.HTML(univ_opening_string), UnivarA_Store.controls)




# Page for Data Visualization

# Three sections: Van Krevelen Plots, Kendrick Mass Defect Plots and Chemical Composition Series
# TODO: Legend does not appear in Van Krevelen Plot - make it appear
# TODO: Make Kendrick Mass Defect Plots have different sizes based on avg. intensity like in VK plots?
# TODO: KMD Plots of different classes have different dot sizes - why?? What can even cause this?

# Param Class to store parameters and data regarding the 3 different plots
class VanKrev_KMD_CCS_Storage(param.Parameterized):
    "Class to store all information and parameters to plot Van Krevelen, Kendrick Mass Defect and Chemical Composition Series."

    # Van Krevelen Parameters
    vk_highlight_by = param.String(default='Rank')
    vk_colour = param.Boolean(True)
    vk_size = param.Boolean(True)
    vk_midpoint = param.Number(default=0.7)
    vk_max_dot_size = param.Number(default=8)
    vk_show_colorbar = param.Boolean(default=True)
    vk_draw_class_rectangle = param.Boolean(default=False)
    vk_text = param.String('Select which columns with Formulas to consider (at least 1 has to be selected):')
    vk_formula_to_consider = param.List(default=checkbox_formula.value)

    # Kendrick Mass Defect Plot Parameters
    kmd_mass_rounding = param.String(default='Up')
    kmd_max_dot_size = param.Number(default=8)
    kmd_text = param.String('Select which columns with Formulas to consider for colouring (if none is selected, dots are not coloured):')
    kmd_formula_to_consider = param.List(default=checkbox_formula.value)

    # Chemical Composition Series Plot Parameters
    ccs_bar_plot_type = param.String(default='Horizontal')
    ccs_text = param.String('Select which columns with Formulas to consider for counting (at least 1 has to be selected):')
    ccs_formula_to_consider = param.List(default=checkbox_formula.value)
    ccs_desc = param.String(default='')

    # Storing figures
    VanKrevelen_plot = param.List(default=['Pane for Van Krevelen Plot'])
    VanKrevelen_filenames = param.List(default=[])
    KendrickMD_plot = param.List(default=['Pane for Kendrick Mass Defect Plot'])
    KendrickMD_filenames = param.List(default=[])
    CCS_plot = param.List(default=['Pane for Chemical Composition Series'])

    # Update the VK Plots
    @param.depends('vk_highlight_by', 'vk_colour', 'vk_size', 'vk_midpoint', 'vk_max_dot_size', 'vk_show_colorbar',
                   'vk_draw_class_rectangle', 'vk_formula_to_consider', watch=True)
    def _compute_VK_plots(self):
        "Computes Van Krevelen Plots and updates layout and widgets."

        if len(self.ccs_formula_to_consider) == 0:
            pn.state.notifications.error('At least 1 Formula Annotation column must be provided for VK plot.')
            raise ValueError('At least 1 Formula Annotation column must be provided for VK plot.')

        # Enabling and Disabling widgets
        if self.vk_highlight_by == 'None':
            self.controls_vk.widgets['vk_colour'].disabled = True
            self.controls_vk.widgets['vk_size'].disabled = True
            self.controls_vk.widgets['vk_midpoint'].disabled = True
            self.controls_vk.widgets['vk_show_colorbar'].disabled = True
        else:
            if self.controls_vk.widgets['vk_colour'].disabled == True:
                self.controls_vk.widgets['vk_colour'].disabled = False
                self.controls_vk.widgets['vk_size'].disabled = False
                self.controls_vk.widgets['vk_midpoint'].disabled = False
                self.controls_vk.widgets['vk_show_colorbar'].disabled = False

        # Setting up stores for figures and filenames
        figures_list = []
        filenames = []

        for g in com_exc_compounds.group_dfs:
            # Getting the Van Krevelen Plots and filenames for each class and adding it to the store
            fig, f = iaf._plot_VK_diagrams_individual(com_exc_compounds.group_dfs[g], # Filtered DF for specific class
                                                  self, # Parameters to use for VK
                                                  DataFrame_Store.univariate_df, # Full normalized data
                                                  target_list.target, g) # Target and current class

            if self.vk_draw_class_rectangle:
                fig = iaf.vk_add_class_rectangles(fig)

            figures_list.append(fig)
            filenames.append(f)

        # Storing as attributes
        self.VanKrevelen_plot = figures_list
        self.VanKrevelen_filenames = filenames

        # Updating the layouts with the new VK plots
        vk_plots.clear()
        for i in range(len(self.VanKrevelen_plot)):
            # Repeating this since otherwise the legend does not appear - it still does not appear
            self.VanKrevelen_plot[i].update_layout(showlegend=True, legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99)),
            pane = pn.pane.Plotly(self.VanKrevelen_plot[i],
                           config = {'toImageButtonOptions': {'filename': self.VanKrevelen_filenames[i], 'scale':4,}})
            vk_plots.append(pane)


    # Update the KMD Plots
    @param.depends('kmd_mass_rounding', 'kmd_max_dot_size', 'kmd_formula_to_consider', watch=True)
    def _compute_KMD_plots(self):
        "Computes Kendrick Mass Defect Plots and updates layout."

        # Setting up stores for figures and filenames
        figures_list = []
        filenames = []

        # Kendrick Mass Defect Plots cannot be computed without a neutral mass column
        if radiobox_neutral_mass.value == 'None':
            return

        for g in com_exc_compounds.group_dfs:
            # Getting the Kendrick Mass Defect Plots for each class and adding it to the store
            fig, f = iaf._plot_KMD_plot_individual(com_exc_compounds.group_dfs[g], # Filtered DF for specific class
                                                   self, # Parameters to use for VK
                                                   g, # Current class
                                                   radiobox_neutral_mass.value) # Neutral Mass Column

            figures_list.append(fig)
            filenames.append(f)

        # Storing as attributes
        self.KendrickMD_plot = figures_list
        self.KendrickMD_filenames = filenames

        # Updating the layouts with the new Kendrick Mass Defect plots
        kmd_plots.clear()
        for i in range(len(self.KendrickMD_plot)):
            pane = pn.pane.Plotly(self.KendrickMD_plot[i],
                           config = {'toImageButtonOptions': {'filename': self.KendrickMD_filenames[i], 'scale':4,}})
            kmd_plots.append(pane)


    # Update the CCS Plot
    @param.depends('ccs_bar_plot_type', 'ccs_formula_to_consider', watch=True)
    def _compute_CCS_plot(self, group_dfs=com_exc_compounds.group_dfs,
                          series_order=('CHO', 'CHOS', 'CHON', 'CHNS', 'CHONS', 'CHOP', 'CHONP','CHONSP', 'other')):
        "Computes and plots the chemical compostion series plot, updating the corresponding layout."

        if len(self.ccs_formula_to_consider) == 0:
            pn.state.notifications.error('At least 1 Formula Annotation column must be provided for CCS plot.')
            raise ValueError('At least 1 Formula Annotation column must be provided for CCS plot.')

        # Initialize figure
        fig = go.Figure()

        # Initialize the DataFrame
        series_df = pd.DataFrame()

        # Description of the number of formulas assigned considered for each class and from how many peaks they came from.
        desc_string = [f'**Description with number of formulas assigned considered**', '']

        # For each class
        for g in group_dfs:
            # Calculating H/C and O/C ratios and series classes for data diversity plots.
            forms = com_exc_compounds.group_dfs[g].dropna(subset=self.ccs_formula_to_consider, how='all')
            elems = iaf.create_element_counts(forms, formula_subset=self.ccs_formula_to_consider)

            # Calculating counts of each series for each class
            counts = elems['Series'].value_counts().reindex(series_order)
            series_df[g] = counts # Store the counts in the dataframe

            # Adding the class' bars to the plot
            if self.ccs_bar_plot_type == 'Horizontal': # Horizontal bar plot
                fig.add_trace(go.Bar(name=g, x=counts, y=series_order, orientation='h'))
            else: # Vertical bar plot
                fig.add_trace(go.Bar(name=g, x=series_order, y=counts, orientation='v'))

            # Adding to the description
            n_peaks = pd.Series(elems.index).value_counts().shape[0]
            group_string = f'**{g}**: **{counts.sum()}** formulas were considered from a total of **{n_peaks}** peaks.'
            desc_string.append(group_string)

        # Defining x and y labels for horizontal and vertical barplots
        if self.ccs_bar_plot_type == 'Horizontal': # Horizontal bar plot
            xlabel = 'Nº of Formulas'
            ylabel = 'Composition Series'
        else: # Vertical bar plot
            xlabel = 'Composition Series'
            ylabel = 'Nº of Formulas'

        # Defining major figure characteristics
        fig.update_layout(
            title="Chemical Composition Series",
            xaxis_title=xlabel,
            yaxis_title=ylabel,
            legend_title="Classes")

        # Setting up the attributes
        self.CCS_plot[0] = fig
        self.ccs_desc = '<br />'.join(desc_string) # <br /> leads to line breaks in Markdown
        filename = 'CCS_Plot_formulacolumns_'
        # Create appropriate filename
        for cl in self.ccs_formula_to_consider:
            filename = filename + f'_{cl}'

        # Updating chemical composition series layout
        if len(ccs_page) == 2:
            ccs_page.append(self.ccs_desc)
            ccs_page.append(pn.pane.Plotly(self.CCS_plot[0], height=600,
                            config = {'toImageButtonOptions': {'filename': filename, 'scale':4,}}))
            ccs_page.append(pn.widgets.DataFrame(series_df.T))
        else:
            ccs_page[2] = self.ccs_desc
            ccs_page[3] = pn.pane.Plotly(self.CCS_plot[0], height=600,
                            config = {'toImageButtonOptions': {'filename': filename, 'scale':4,}})
            ccs_page[4] = pn.widgets.DataFrame(series_df.T)


    def update_widgets(self):
        "Update the needed widget values."
        # Calculate specific class DataFrames in case it has not been calculated before
        if com_exc_compounds.group_dfs == {}:
            iaf._group_compounds_per_class(com_exc_compounds, target_list, DataFrame_Store) # Add compounds per class dfs

        formula_cols = checkbox_formula.value + [i for i in DataFrame_Store.metadata_df.columns if i.startswith('Matched') and i.endswith('formulas')]

        # VK Plots
        self.controls_vk.widgets['vk_formula_to_consider'].options = formula_cols
        self.controls_vk.widgets['vk_formula_to_consider'].value = formula_cols
        self._compute_VK_plots()

        # KMD Plots
        self.controls_kmd.widgets['kmd_formula_to_consider'].options = formula_cols
        self.controls_kmd.widgets['kmd_formula_to_consider'].value = formula_cols
        self._compute_KMD_plots()

        # CCS
        self.controls_ccs.widgets['ccs_formula_to_consider'].options = formula_cols
        self.controls_ccs.widgets['ccs_formula_to_consider'].value = formula_cols
        self._compute_CCS_plot(com_exc_compounds.group_dfs)


    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets_vk = {
            'vk_highlight_by': pn.widgets.Select(name='What method to use to highlight points based on avg. intensity:',
                                value='Rank', options=['Rank', 'logInt', 'None'],
                                description='''Van Krevelen points (peaks) are coloured or sized based on their average intensity.
                                Rank - by the rank of their average intensity compared to others.
                                logInt - by the logarithm (base 10) of their averaged intensity compared to others.
                                Midpoint marks the point from which a peak is considered low or high intensity.
                                e.g. using 0.7 and "Rank" means that, if selected, 70% of points (lowest intensity) will have the low intensity colour and will be smaller and 30% the higher intensity colour and will be bigger.
                                '''),
            'vk_colour': pn.widgets.Checkbox(name='Colour points by intensity', value=True),
            'vk_size': pn.widgets.Checkbox(name='Size points by intensity', value=True),
            'vk_midpoint': pn.widgets.FloatSlider(name='Midpoint (% of points with "low" intensity)',
                                value=0.7, start=0.0, end=1.0, step=0.05),
            'vk_max_dot_size': pn.widgets.IntSlider(name='Max. dot size (or dot size in general)',
                                value=8, start=1, end=20, step=1),
            'vk_show_colorbar': pn.widgets.Checkbox(name='Show colorbar', value=True),
            'vk_draw_class_rectangle': pn.widgets.Checkbox(name='Draw Peptide, Lignins, Tannins, Nucleotides and Phytochemical area rectangles', value=False),
            'vk_text': pn.widgets.StaticText(name='',
                                value='Select which columns with Formulas to consider (at least 1 has to be selected):',
                                styles={'font-weight': 'bold'}),
            'vk_formula_to_consider': pn.widgets.CheckBoxGroup(name='Columns', value=checkbox_formula.value,
                                options=checkbox_formula.value, inline=False),
        }

        widgets_kmd = {
            'kmd_mass_rounding': pn.widgets.Select(name='Choose how to round experimental masses to nominal:',
                                value='Up', options=['Up', 'Nearest'],
                                description='''Kendrick Nominal Mass is obtained by rounding up or to the nearest integer.
                                '''),
            'kmd_max_dot_size': pn.widgets.IntSlider(name='Dot size',
                                value=8, start=1, end=20, step=1),
            'kmd_text': pn.widgets.StaticText(name='',
                                value='Select which columns with Formulas to consider for colouring (if none is selected, dots are not coloured):',
                                styles={'font-weight': 'bold'}),
            'kmd_formula_to_consider': pn.widgets.CheckBoxGroup(name='Columns', value=checkbox_formula.value,
                                options=checkbox_formula.value, inline=False),
        }

        widgets_ccs = {
            'ccs_bar_plot_type': pn.widgets.Select(name='Type of bar plot:',
                                value='Horizontal', options=['Horizontal', 'Vertical']),
            'ccs_text': pn.widgets.StaticText(name='',
                                value='Select which columns with Formulas to consider for counting (at least 1 has to be selected):',
                                styles={'font-weight': 'bold'}),
            'ccs_formula_to_consider': pn.widgets.CheckBoxGroup(name='Columns', value=checkbox_formula.value,
                                options=checkbox_formula.value, inline=False),
        }

        # Control panel for the Van Krevelen section
        self.controls_vk = pn.Param(self,
                                 parameters=['vk_highlight_by', 'vk_colour', 'vk_size', 'vk_midpoint', 'vk_max_dot_size',
                                            'vk_show_colorbar', 'vk_draw_class_rectangle', 'vk_text', 'vk_formula_to_consider'],
                                 widgets=widgets_vk, name='Van Krevelen Plot Parameters')

        # Control panel for the Kendrick Mass Defect section
        self.controls_kmd = pn.Param(self,
                                 parameters=['kmd_mass_rounding', 'kmd_max_dot_size',
                                            'kmd_text', 'kmd_formula_to_consider'],
                                 widgets=widgets_kmd, name='Kendrick Mass Defect Plot Parameters')

        # Control panel for the Chemical Composition Series section
        self.controls_ccs = pn.Param(self,
                                 parameters=['ccs_bar_plot_type', 'ccs_text', 'ccs_formula_to_consider'],
                                 widgets=widgets_ccs, name='Chemical Composition Series Parameters')


# Initializing the store
dataviz_store = VanKrev_KMD_CCS_Storage()

# Button to compute the data diversity visualization plots for the first time
compute_vk_kmd_ccs_button = pn.widgets.Button(
    name='Compute Van Krevelen, Kendrick Mass Defect and Chemical Composition Series Plots', button_type='success')
def _compute_vk_kmd_ccs_button(event):
    dataviz_store.update_widgets()
compute_vk_kmd_ccs_button.on_click(_compute_vk_kmd_ccs_button)


# Van Krevelen Plot Section
vk_opening_string = '''<strong>Van Krevelen Plot section</strong>
<br>
<br>
This section plots a Van Krevelen Plot for each of the classes under analysis. This is made by only considering metabolites
(features) that appear at least in one sample of said class.
<br>
Only metabolites with assigned formulas are considered. You can choose which formula annotation you want to use whether it
is a previous annotation performed, an annotation performed in this software, or even multiple different annotations
(<strong>Note: you HAVE to select at least one annotation</strong>). If multiple formulas can be assigned to a metabolite
whether within the same database annotation or different, they are both considered and plotted in the Van Krevelen.
<br>
<br>
You can also decide to colour the dots and make them different sized based on their average intensity in the corresponding
class using either a rank-based (more linear) approach or the logarithm of the averages intensity (less linear). You can
also choose a "midpoint", where the dots with avg. intensity below the midpoint % of the intensity of all features have a
more bluish colour, while those above have a redder colour.
'''

vk_plots = pn.Column()

vk_page = pn.Column(pn.pane.HTML(vk_opening_string), dataviz_store.controls_vk, vk_plots)


# Kendrick Mass Defect Plot Section
kmd_opening_string = '''<strong>Kendrick Mass Defect Plot section</strong>
<br>
<br>
This section plots a Kendrick Mass Defect Plot for each of the classes under analysis. This is plotted based on Neutral
Mass column selected.
<br>
Kendrick Nominal Mass and Mass Defect can be calculated in two ways:
<br>
- <strong>'Up'</strong>: by rounding the neutral mass up, thus the Mass Defects will vary between 0 and 1.
<br>
- <strong>'Nearest'</strong>: by rounding the neutral mass to the nearest integer, thus the Mass Defects will vary between -0.5 and 0.5.
<br>
<br>
You can decide to colour the dots based on the chemical composition series they belong to ('CHO', 'CHON', 'CHOP', 'CHONSP',
etc.) by selecting the formula columns of previous annotations, of the annotations performed in the software or a
combination of them. When, for one peak, there are multiple candidate formulas and they do not belong to the same chemical
composition series, they get assigned as <strong>'Ambiguous'</strong>. If <strong>None</strong> is selected, then the
points are not coloured by chemical composition series.
'''

kmd_plots = pn.Column()

kmd_page = pn.Column(pn.pane.HTML(kmd_opening_string), dataviz_store.controls_kmd, kmd_plots)


# Chemical Composition Series Section
ccs_opening_string = '''<strong>Chemical Composition Series section</strong>
<br>
<br>
This section plots a Chemical Composition Series bar plot indicating how many formulas of different chemical series ('CHO',
'CHON', 'CHOP', 'CHONSP', etc.) were assigned to samples of each class.
<br>
<br>
You can decide if you want to use formula annotations made previously, made in this software or a combination. <strong>If
multiple formulas</strong> can be assigned to the same m/z peak whether within the same formula annotation made or between
different annotations, <strong>each one will be counted in this plot</strong>. That is, if a peak in a class has 3 possible
candidate formulas, 2 belonging to the 'CHO' series and another to the 'CHOP' series; then 2 formulas will be added to the
'CHO' series and 1 to the 'CHOP' series. Thus, we are considering that the 3 elementary formulas are represented by that
peak (probably an overestimation). In all cases, it is rare to find multiple candidate formulas for the same m/z peak,
especially with extreme-resolution data. To provide an idea of how extensive this effect is, a description is also provided
detailing how many formulas are being considered for each class and from how many different features (<em>m/z</em> peaks)
they came from.
'''

ccs_page = pn.Column(pn.pane.HTML(ccs_opening_string), dataviz_store.controls_ccs)


data_viz_page = pn.Column(compute_vk_kmd_ccs_button, pn.Tabs(('Van Krevelen Plot', vk_page), ('Kendrick Mass Defect Plot', kmd_page),
                        ('Chemical Composition Series', ccs_page)))




# Page for BinSim Analysis
binsim_analysis_page = pn.Column()

# Page for Compound Finder
comp_finder_page = pn.Column()




# Overall layout of the program and initialization
# The pages we have
pages = {
    "Index": OpeningPage(),
    "Data Reading": DataReading(),
    "Data Metadata": DataMetadata(),
    "Data Filtering": DataFiltering(),
    "Data Annotation": DataAnnotation(),
    "Data Pre-Treatment": DataPreTreatment(),
    "Class Colours": ClassColours(),
    "Transitional Page": TransitionalPage(),
    "Common and Exclusive Compounds": CommonExclusivePage(),
    "Unsupervised Analysis": UnsupervisedAnalysisPage(),
    "Supervised Analysis": SupervisedAnalysisPage(),
    "Univariate Analysis": UnivariateAnalysisPage(),
    "Data Visualization": DataVisualizationPage(),
    "BinSim Analysis": BinSimPage(),
    "Compound Finder": CompoundFinderPage()
}

# Function to show the selected page - needs update (may cause bug)
def show_page(page_instance):
    main_area.clear()
    def clear_again():
        main_area.clear()
    clear_again()
    #main_area.get_param_values()
    main_area.append(page_instance.view())


# Create and organize the sidebar
# The sidebar will include buttons to every main page and will be used as a tool for navigation

# Define buttons to navigate between pages
index_button = pn.widgets.Button(name="Home", button_type="primary", icon=iaf.home_icon, disabled=False)
page1_button = pn.widgets.Button(name="Data Reading", button_type="primary")
page1_1_button = pn.widgets.Button(name="Data Metadata", button_type="primary", disabled=True)
page1_2_button = pn.widgets.Button(name="Data Filtering", button_type="primary", disabled=True)
page2_button = pn.widgets.Button(name="Data Annotation", button_type="primary", disabled=True)
page3_button = pn.widgets.Button(name="Data Pre-Treatment", button_type="primary", disabled=True)
page4_button = pn.widgets.Button(name="Class Colours", button_type="primary", disabled=True)
page5_button = pn.widgets.Button(name="Common/Exclusive Comp.", button_type="default", disabled=True)
page6_button = pn.widgets.Button(name="Unsupervised Analysis", button_type="default", disabled=True)
page7_button = pn.widgets.Button(name="Supervised Analysis", button_type="default", disabled=True)
page8_button = pn.widgets.Button(name="Univariate Analysis", button_type="default", disabled=True)
page9_button = pn.widgets.Button(name="Data Visualization", button_type="default", disabled=True)
page10_button = pn.widgets.Button(name="BinSim Analysis", button_type="default", disabled=True)
page11_button = pn.widgets.Button(name="Compound Finder", button_type="default", disabled=True)
RESET_button = pn.widgets.Button(name="RESET", button_type="danger", disabled=False)

# Set up button click callbacks
index_button.on_click(lambda event: show_page(pages["Index"]))
page1_button.on_click(lambda event: show_page(pages["Data Reading"]))
page1_1_button.on_click(lambda event: show_page(pages["Data Metadata"]))
page1_2_button.on_click(lambda event: show_page(pages["Data Filtering"]))
page2_button.on_click(lambda event: show_page(pages["Data Annotation"]))
page3_button.on_click(lambda event: show_page(pages["Data Pre-Treatment"]))
page4_button.on_click(lambda event: show_page(pages["Class Colours"]))
page5_button.on_click(lambda event: show_page(pages["Common and Exclusive Compounds"]))
page6_button.on_click(lambda event: show_page(pages["Unsupervised Analysis"]))
page7_button.on_click(lambda event: show_page(pages["Supervised Analysis"]))
page8_button.on_click(lambda event: show_page(pages["Univariate Analysis"]))
page9_button.on_click(lambda event: show_page(pages["Data Visualization"]))
page10_button.on_click(lambda event: show_page(pages["BinSim Analysis"]))
page11_button.on_click(lambda event: show_page(pages["Compound Finder"]))

# Reset panel widgets
reset_panel_float_text = pn.widgets.StaticText(name='', value='Do you want to reset all parameters of the software?')
reset_panel_float_yes_button = pn.widgets.Button(name='Yes', button_type='danger')
reset_panel_float_no_button = pn.widgets.Button(name='No', button_type='default')
reset_floatpanel = pn.layout.FloatPanel(reset_panel_float_text, reset_panel_float_yes_button, reset_panel_float_no_button,
                                      name='Reset', margin=20, position='center', config={"headerControls": {"close": "remove"}},
                                      contained=False)

def RESET(event):
    reset_panel = pn.Column(reset_floatpanel, height=200)
    #main_area.clear()
    if len(main_area) == 1:
        main_area.append(reset_panel)
    else:
        reset_floatpanel.status = 'normalized'

def No_Reset(event):
    reset_floatpanel.status = 'closed'
    main_area.pop(-1)

#def Yes_Reset(event):


RESET_button.on_click(RESET)
reset_panel_float_no_button.on_click(No_Reset)
#reset_panel_float_yes_button.on_click(Yes_Reset)

sidebar = pn.Column(index_button, '## Data Pre-Processing and Pre-Treatment', page1_button, page1_1_button, page1_2_button, page2_button, page3_button, page4_button,
                   '## Statistical Analysis', page5_button, page6_button, page7_button, page8_button, page9_button, page10_button,
                   page11_button, '## To Reset (TODO)', RESET_button)

app = pn.template.FastListTemplate(title='Testing MetsTA', sidebar=[sidebar], main=[main_area])

app.show()