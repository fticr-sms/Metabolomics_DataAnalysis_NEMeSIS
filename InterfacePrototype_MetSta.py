
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
import pickle

# File with functions to auxiliate the graphical interface
import interface_aux_functions as iaf
import description_strings as desc_str
from report_generation import ReportGenerator

# metanalysis_standard.py file
import metanalysis_standard as metsta
import multianalysis as ma

# The initial pages, especially the read file one does not have the nomenclature that I started using later on
# for the different widgets as well as organization
pn.extension('plotly', 'floatpanel', 'katex', notifications=True)
hv.extension('plotly')
pn.config.sizing_mode="stretch_width"
mpl.use('agg')

# To allow the interaction of different packages
pd.DataFrame.iteritems = pd.DataFrame.items

# TODO: Make a way to choose folder where all figures and tables downloaded go to
# TODO: Updating packages made a series of future deprecation warnings appear - adapt code to them
# TODO: Make PCA and projection plots save the current components names in the filename


# Define pages as classes
# Initial Pages class building with barebones for each class
class OpeningPage:
    def __init__(self):
        self.content = pn.Column("# Welcome to MetsTA!", acronym_section,
                                 "# The Go-To place for your extreme-resolution metabolomics data analysis need.",
                                 pn.Row(homepage_page, pn.pane.Image('Picture_Test.png', height=400)))
    def view(self):
        return self.content


class InstructionPage:
    def __init__(self):
        self.content = pn.Column("# Instructions and Considerations to keep in mind about MetsTA",
                                 instruction_page)
    def view(self):
        return self.content


class DataReading:
    def __init__(self):
        self.content = pn.Column("# Section 1: Data Input",
                                 pn.pane.HTML(desc_str.data_reading_opening_string),
                                section1page)

    def view(self):
        return self.content


class DataMetadata:
    def __init__(self):
        self.content = pn.Column("# Section 1.1: Selecting Metadata columns and defining your target",
                                 "Metadata columns are split into 4 categories: Formula, Annotated (compound name), Neutral Mass and Other columns. Sample columns are every column which was not selected in any metadata column.",
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
    """Perform Annotations based on available databases. Must provide the database filename and the name of the columns with the **ID**, the **Name** and **Formula** of the metabolites. **Cannot perform annotation without a selected Neutral Mass column.**
    You can annotate with multiple databases. However, each database is annotated individually.
    Annotation works by assigning to a m/z peak / feature all metabolites of a database that are within the provided error margin.
    **Currently, no adduct search is done to perform annotations (MetaboScape Data's Bucket Label should already take adducts into account).**
    Annotation from two different databases might annotate different metabolites for the same m/z peak / feature.
    Thus, each database annotation will generate **3 columns** added to the metadata: one with the **IDs** of the metabolites annotated, another with their **formula** and another with their **name**.""",
                                page2)

    def view(self):
        return self.content


class AnnDeDuplication:
    def __init__(self):
        self.content = pn.Column("# Section 2.1: Data Multiple Annotation De-Duplication - Metabolic Feature Merging",
                                page2_1)

    def view(self):
        return self.content


class DataPreTreatment:
    def __init__(self):
        
        self.content = pn.Column("# Section 3: Data Pre-Treatment",
                                 "Choose the pre-treatment to apply to the data. **Large Datasets might not show up right away on the right (refreshing may correct this).**",
                                page3)

    def view(self):
        return self.content


class ClassColours:
    def __init__(self):

        self.content = pn.Column("# Section 3.1: Select Colours for each Class",
                                 "These colours will be used in the different figures made hereafter (Venn diagrams, PCA, HCA, PLS-DA, Chemical Composition Series and plots in the Compound Finder search tool).",
                                 "If you are repeating analysis after having modified the dataset or the pre-treatment, a bunch of notifications may appear. Do not worry.",
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
                                 "This includes an overview analysis as well as Venn Diagrams and Intersection Plots of common and exclusive compounds of the different classes in the dataset.",
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
                                 ("In Univariate Analysis, each metabolite (variable) in the experimental dataset is tested individually to observe "
                                  'if there is a significant difference between a test and a control class. Thus, this does not take into account any '
                                  'interaction between metabolites as multivariate analysis does (and is expected in metabolites within a biological system).\n'
                                  'However, it provides the metabolites which are differentially expressed between 2 classes.'),
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


class PathwayAssignmentPage:
    def __init__(self):

        self.content = pn.Column("# Pathway Assignment (HMDB Annotation)",
                                 pn.pane.HTML(desc_str.path_assign_opening_string),
                                 path_assign_page)

    def view(self):
        return self.content


class BinSimPage:
    def __init__(self):

        self.content = pn.Column("# Performing BinSim Analysis",
                                 pn.pane.HTML(desc_str.BinSim_opening_string),
                                 binsim_analysis_page)

    def view(self):
        return self.content


class CompoundFinderPage:
    def __init__(self):

        self.content = pn.Column("# Find a specific compound",
                                 """Observe barplots and boxplots of a specific selected compound. The plots using averages and standard deviations of the classes **ignore missing values** of samples in those classes by default.
                                 See the sample bar plot to observe if there are low number of missing values in the classes for the searched compound.
                                 The class bar plots are only usable as good approximations if so. If not, they should be ignored, since missingness is an important part of the observed compound behaviour in the dataset.
                                 Missing values can be accounted for as 0 if that option is preferred however this is also not ideal. There are **2 checkboxes** for this parameter (2 plots) so the last one you have changed is the one that remains true for both plots.
                                 Check figure save name to confirm the exact parameters of how the figure was made if you do not remember the last one you have changed.""",
                                 comp_finder_page)

    def view(self):
        return self.content


class ReportGenerationPage:
    def __init__(self):

        self.content = pn.Column("# Metabolomics Data Analysis Report Generation",
                                 "## Initial Layout (Under Construction)",
                                 rep_gen_page)

    def view(self):
        return self.content




# Homepage

acronym_string = desc_str.acronym_string
homepage_string = desc_str.homepage_string

acronym_section = pn.pane.HTML(acronym_string, styles={'font-size': 'large'})
homepage_page = pn.pane.HTML(homepage_string)



# Instruction page
instruction_page_string = desc_str.instructions_page
instruction_page = pn.pane.HTML(instruction_page_string)


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


    def concat_annots(_, MS_df, annot_df):
        "Joins m/z peak data with annotation data in a single DataFrame."
        return pd.concat((MS_df, annot_df), axis=1)


    def reset(self):
        "Resets all relevant parameters."
        for param in self.param:
            if param not in ["name"]:
                setattr(self, param, self.param[param].default)


    def __init__(self, **params):

        super().__init__(**params)

        self.controls = pn.Param(self, parameters=['treated_df'], name='Pre-Treatment Selection')

# Initializing the Store
DataFrame_Store = DataFrame_Storage()


# Function to ensure statistical analysis buttons are disabled when data pre-treatment is changed

def _disabling_stat_analysis_buttons():
    "Disabling statistical analysis."
    # Enable all statistical analysis related buttons
    page5_button.disabled = True
    page6_button.disabled = True
    page7_button.disabled = True
    page8_button.disabled = True
    page9_button.disabled = True
    page10_button.disabled = True
    page11_button.disabled = True
    page12_button.disabled = True
    page13_button.disabled = True




# TODO: Data Visualization page does not reset figure parameters to default (Volcano plot does not reset colours but that is okay)

# Page 1 - Reading File
# TODO: Make it be able to read positive and negative ionization mode obtained data perhaps so it is not just suited to MetaboScape

class FileReading(param.Parameterized):
    """Class to store as attributes file read."""

    temp_target = param.Dict()
    read_df = param.DataFrame(pd.DataFrame())
    neutral_mass_column_inserted = param.Boolean(default=False)

    def reset(self):
        "Resets all relevant parameters."
        for param in self.param:
            if param not in ["name"]:
                setattr(self, param, self.param[param].default)

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
    "Controls the state of the button to confirm and read the file."
    if filename != '':
        confirm_button_filename.disabled = False
    else:
        confirm_button_filename.disabled = True

def _confirm_button_filename(event):
    "Reads the file given."

    # Read the file, updating widgets and parameters
    file.read_df, file.temp_target, file.neutral_mass_column_inserted = iaf.read_file(filename.filename, target_included_in_file.value)

    # Parameters to store for Report Generation
    RepGen.filename = filename.filename
    RepGen.target_included_in_file = target_included_in_file.value
    RepGen.neutral_mass_column = file.neutral_mass_column_inserted

    # Enabling button for next step
    section1page[2] = pn.widgets.DataFrame(file.read_df, disabled=True, sortable=False, reorderable=False)
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
    file.read_df, file.temp_target, file.neutral_mass_column_inserted = iaf.read_file('5yeasts_notnorm.csv', False)

    # Parameters to store for Report Generation
    RepGen.filename = 'Example Dataset (5yeasts_not_norm.csv)'
    RepGen.target_included_in_file = False
    RepGen.neutral_mass_column = file.neutral_mass_column_inserted

    # Enabling button for next step
    section1page[2] = pn.widgets.DataFrame(file.read_df, disabled=True, sortable=False, reorderable=False)
    confirm_button_step1.disabled = False
# Function happens when you press the button
load_example_df_button.on_click(_load_example_df_button)


# Confirm file, show next page, disable reading files, update columns of the dataset read
def _confirm_step1(event):
    "Performs actions to pass from step 1 page to step 1_1 page."
    # Enabling/Disabling appropriate Widgets
    page1_1_button.disabled = False
    confirm_button_filename.disabled = True
    load_example_df_button.disabled = True
    target_included_in_file.disabled = True
    filename.disabled = True
    confirm_button_target.disabled=True
    confirm_button_next_step_1_1.disabled=True

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
    "Get sample columns based on metadata columns selected, makes target editable wile providing an educated guess about it."
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
    confirm_button_target.disabled = False
    
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

# Disable further analysis when this is changed
@pn.depends(checkbox_formula.param.value, checkbox_annotation.param.value, radiobox_neutral_mass.param.value,
            checkbox_others.param.value, watch=True)
def _disable_remaining_analysis_from_metadata(a,b,c,d):
    "Disable analysis further on the pipeline"
    # Disable widgets on the current page
    confirm_button_target.disabled=True
    confirm_button_next_step_1_1.disabled=True

    # Disable sidebar buttons
    page1_2_button.disabled = True
    page2_button.disabled = True
    page2_1_button.disabled = True
    page3_button.disabled = True
    page4_button.disabled = True

    # Disable statistical analysis
    _disabling_stat_analysis_buttons()

# Make the target widget
target_placeholder = "A,A,A,A,A,B,B,B,B,B"
target_widget = pn.widgets.TextAreaInput(name='Define your Target (Class Labels)', placeholder=target_placeholder,
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
    "Controls the state of the button to confirm and read the target."
    if target_widget != '':
        confirm_button_target.disabled = False
    else:
        confirm_button_target.disabled = True


# Make sure that target makes sense in comparison to the number of sample columns
def _update_confirm_target(event):
    "Confirms the target selected."
    target = target_widget.value.split(',')
    # Filling the target storage with the correct target and default colours
    target_list.target = target_widget.value.split(',')
    target_list.classes = list(pd.unique(np.array(target_list.target)))
    target_list(target_widget.value.split(','), colours)
    sample_cols = target_list.sample_cols
    if len(sample_cols) != len(target):
        pn.state.notifications.error(
            f'Number of class labels ({len(target)}) is different than the number of sample columns ({len(sample_cols)}).')
        confirm_button_next_step_1_1.disabled = True
    else:
        confirm_button_next_step_1_1.disabled = False

    # Disable sidebar buttons
    page1_2_button.disabled = True
    page2_button.disabled = True
    page2_1_button.disabled = True
    page3_button.disabled = True
    page4_button.disabled = True

    # Disable statistical analysis
    _disabling_stat_analysis_buttons()


# Call the function
confirm_button_target.on_click(_update_confirm_target)

# Going to the next step function
def _confirm_button_next_step_1_1(event):
    "Performs actions to pass from step 1_1 page to step 1_2 page."
    page1_2_button.disabled = False
    confirm_button_next_step_2.disabled = True
    #page2_button.disabled = False
    # Assuring the initial layout of data filtering page
    while len(page1_2) > 3:
        page1_2.pop(-1)
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
    "Controls widget limits related to data filtering keyword based on the method chosen and the target."
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
filtered_df = pn.widgets.DataFrame(pd.DataFrame(), name='Filtered DataFrame', disabled=True, sortable=False, reorderable=False)
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

    # Disable posterior sidebar buttons
    page2_button.disabled = True
    page2_1_button.disabled = True
    page3_button.disabled = True
    page4_button.disabled = True
    # Disable statistical analysis
    _disabling_stat_analysis_buttons()

    confirm_button_next_step_2.disabled = False

    # Setup the page if not setup yet
    while len(page1_2) > 3:
        page1_2.pop(-1)
    page1_2.extend(['#### Characteristics of the Dataset',characteristics_df,'#### Filtered Dataset',
                    filtered_df,
                    confirm_button_next_step_2])

# Call the function
confirm_button_initial_filtering.on_click(_confirm_button_initial_filtering)

# Go to next step function and calling it
def _confirm_button_next_step_1_2(event):
    "Performs actions to pass from step 1_2 page to Data Annotation page."
    page2_button.disabled = False
    confirm_button_next_step_2_1.disabled = True

    while len(page2) > 1:
        page2.pop(-1)

    main_area.clear()
    show_page(pages["Data Annotation"])
confirm_button_next_step_2.on_click(_confirm_button_next_step_1_2)

# Initial page layout
page1_2 = pn.Column(pn.Row(filt_method, filt_method_tooltip), pn.Row(filt_kw, filt_kw_tooltip),
                    confirm_button_initial_filtering)




# Page 2 - Annotation of Metabolites
# TODO: Make it selectable if you want to search for possible adducts (if added possibility to read positive and negative m/z indexes)

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
        db_file_input = pn.widgets.TextInput(placeholder='placeholder.csv')
        
        static_db_abv = pn.widgets.StaticText(name='Abbreviation', value='')
        db_abv_input = pn.widgets.TextInput(placeholder='PLC')

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
            "Controls if the button to read database is pressable."
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
            "Reads Database inputted and stores parameters."
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
                #confirm_button_databases_read.disabled = False
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

def DB_dict_reset(DB_dict):
    "Resets the DB_dict section."
    for i in DB_dict.keys():
        DB_dict[i] = ['1']
        DB_dict[i] = DatabaseSection()
        
dbs_arrangement_all = pn.Row(DB_dict['1'].content, DB_dict['2'].content, DB_dict['3'].content,
                   DB_dict['4'].content, DB_dict['5'].content)
dbs_arrangement = pn.Row()

# Make the designated number of database sections appear
def _confirm_button_n_databases(event):
    "Updates layout making the database sections selected appear."
    # Updating the value
    n_databases.value = n_databases_show.value
    # Setting the databases read attribute to False
    for i in DB_dict:
        DB_dict[i].read.value = False

    # Disable posterior sidebar buttons
    page2_1_button.disabled = True
    page3_button.disabled = True
    page4_button.disabled = True
    # Disable statistical analysis
    _disabling_stat_analysis_buttons()

    # Setting up the layout to show
    titles = pn.Row()
    dbs_arrangement.clear()
    for i in range(n_databases.value):
        if len(dbs_arrangement_all[i]) == 7:
            dbs_arrangement_all[i].pop(-1)
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

    page2.append(confirm_button_databases_read)
    confirm_button_databases_read.disabled = True
    if n_databases.value == 0: # Case where no database is going to be used for annotation
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
    "Controls if the button to confirm the databases selected is pressable."
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
        confirm_button_databases_read.disabled=False

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
    "Controls the widget that appears as a parameter to complement the annotation margin method chosen."
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
    "Confirms the databases read and updates layout."
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

    performing_annotation_arrangement.clear()
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
        # Unnecessary I believe
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
    # Update layout of the page
    while len(performing_annotation_arrangement)>0:
        performing_annotation_arrangement.pop(-1)
    page2.append(performing_annotation_arrangement)
    confirm_button_annotation_perform.disabled=True

    # Perform metabolite annotation
    metabolite_annotation()

    # Store parameters used for annotation for Report Generation
    RepGen.annotation_margin_method = annotation_margin_method_radio.value
    RepGen.annotation_margin_ppm_deviation = annotation_ppm_deviation.value
    RepGen.annotation_margin_Da_deviation = annotation_Da_deviation.value

    # Update the information for annotation de-duplication
    data_ann_deduplicator.update_columns_with_annotations(annotated_df.value, checkbox_annotation, checkbox_formula)
    DataFrame_Store.original_df = DataFrame_Store.concat_annots(filtered_df.value, annotated_df.value)
    iaf.creating_has_match_column(DataFrame_Store, n_databases, checkbox_annotation)
    data_ann_deduplicator.annotated_df = DataFrame_Store.original_df
    data_ann_deduplicator.create_multiple_annotation_report()

    # Disable posterior sidebar buttons
    page2_1_button.disabled = True
    page3_button.disabled = True
    page4_button.disabled = True
    # Disable statistical analysis
    _disabling_stat_analysis_buttons()

    # Layout after annotation
    confirm_button_next_step_2_1.disabled = False
    page2.append(confirm_button_next_step_2_1)

confirm_button_annotation_perform.on_click(_press_confirm_annotation_perform)

# Button to next step
confirm_button_next_step_2_1 = pn.widgets.Button(icon=iaf.img_confirm_button, name='Next Step - Annotation De-Duplication',
                                                     button_type='success', disabled=False)

# Go to next step function and calling it
def _confirm_button_next_step_2_1(event):
    "Performs actions to pass from step 2_1 page to step 3 page."
    page2_1_button.disabled = False
    confirm_button_next_step_3.disabled = True

    # Reset the next page layout
    while len(page2_1)>3:
        page2_1.pop(-1)

    # Updating next page layout
    page2_1[1] = pn.pane.DataFrame(data_ann_deduplicator.mult_ann_report)

    # Update to show the Data Pre-Treatment page
    main_area.clear()
    show_page(pages["Annotation De-Duplication"])
confirm_button_next_step_2_1.on_click(_confirm_button_next_step_2_1)

#### This marks the separation to use mainly param instead of only panel




# Data Multiple Annotaton De-Duplication and Peak Merging Page

# Param Class to store parameters and data regarding Annotation De-Duplication
class AnnDeDuplication_Storage(param.Parameterized):
    "Class to store all information on Annotation De-Duplication."

    # Original DataFrame to work with
    annotated_df = param.DataFrame()

    # Columns with annotations
    mcid = param.List([])
    mcid_alt = param.Dict({}) # Alternate names including mentioning Previous Annotations

    # DataFrames to group information
    mult_ann_report = param.DataFrame(pd.DataFrame())
    mergings_performed = param.DataFrame()
    merge_description = param.DataFrame()
    merge_situations = param.DataFrame()
    full_merge_problems = param.DataFrame()
    merge_problems = param.DataFrame()
    merge_report = param.String()

    # Params used
    current_params = param.Dict()

    # Parameters for de-duplication
    consider_formula_cols = param.Boolean(default=True)
    text_problem_condition = param.String('How to Treat Scenarion 1 of Situation 4')
    problem_condition = param.String(default='Scenario 1 of Situation 4 like cases are not merged and are not shown.')


    def update_columns_with_annotations(self, ann_df, checkbox_annotation, checkbox_formula):
        "Creates the list with the annotation (compund and formula) columns in the data and update attributes."
        # Previous Annotations columns
        mcid_alt = {i:i+' - Prev. Ann.' for i in checkbox_annotation.value}
        mcid  = [i for i in checkbox_annotation.value]

        # Annotations columns made in the software
        for col in ann_df.columns:
            if col.startswith('Matched ') and col.endswith(' IDs'):
                mcid.append(col)
                mcid_alt[col] = col

        # Previous Formula Annotation columns
        mcid  = mcid + checkbox_formula.value
        for i in checkbox_formula.value:
            mcid_alt[i] = i+' - Prev. Ann.'

        # Updating attributes
        self.mcid = mcid
        self.mcid_alt = mcid_alt


    def create_multiple_annotation_report(self):
        "Builds DataFrame with report of possible multiple of the same annotations for each annotation used."
        # Set the DataFrame
        self.mult_ann_report = pd.DataFrame(index=['Nº of same annotations on multiple peaks',
                                                     'Total number of annotations for these cases',
                                                     'Maximum number of peaks with the same annotation'])

        try:
            for col in self.mcid:
                # See multiple annotations
                n_duplicates = self.annotated_df[col].value_counts()
                n_duplicates = n_duplicates[n_duplicates>1]
                if len(n_duplicates) == 0:
                    max_peaks = 0
                else:
                    max_peaks = n_duplicates.iloc[0]

                # Build the report
                self.mult_ann_report[self.mcid_alt[col]] = [n_duplicates.sum(), len(n_duplicates), max_peaks]

        except:
            pn.state.notifications.info('Multiple Annotation Report could not be compiled.')


    def perform_deduplication(self, ann_df, sample_cols, neutral_mass_col):
        "Performs Multiple Annotation peak merging, updates layout and returns results."

        # Performing deduplication
        # TODO: See deduplication function better
        annotated_data, mergings_performed, merging_situations, merge_description, mp = iaf.duplicate_disambiguator(
            self, ann_df, sample_cols, neutral_mass_col, mz_col=False)

        # Report
        self.merge_report = f'''Nº of Mergings: **{len(merge_description)}**
        Nº of Metabolic Features merged: **{pd.DataFrame(merge_description).T['Nº merged peaks'].sum()}**
        Nº of Metabolic Features dropped: **{len(ann_df) - len(annotated_data.index)}**
        Nº of Metabolic Features before merging: **{len(ann_df)}**
        Nº of Metabolic Features after merging: **{len(annotated_data.index)}**'''

        # Transforming and storing results
        self.annotated_df = annotated_data
        self.merge_description = pd.DataFrame(merge_description).T
        self.merge_description.index.name = 'New Indexes'
        self.mergings_performed = pd.DataFrame(pd.Series(mergings_performed), columns=['Nº of Mergings by Database'])
        self.merge_situations = pd.DataFrame(pd.Series(merging_situations), columns=['Nº of Mergings by Situation'])

        # Problems results
        if len(mp) > 0:
            self.merge_problems = iaf._merge_problems_creation(mp)
            self.full_merge_problems = self.merge_problems.copy()
        else:
            self.merge_problems = pd.DataFrame(columns=['DBs with same annotation', 'Annotation',
                                            'Nº of Peaks', 'Indexes', 'Contradiction With'])
            self.full_merge_problems = self.merge_problems.copy()
        self.merge_situations.loc['Possible Problem Cases'] = len(self.merge_problems)

        # Store parameters used
        self.current_params = {'consider_formula_cols': self.consider_formula_cols,
                              'problem_condition': self.problem_condition}


    def reset(self):
        "Reset parameters."
        for param in self.param:
            if param not in ["name"]:
                setattr(self, param, self.param[param].default)


    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets = {
            'consider_formula_cols': pn.widgets.Checkbox(name='Consider the formula columns for metabolic feature merging',
                                                        value=True, disabled=True),
            'text_problem_condition': pn.widgets.StaticText(name='', value='How to Treat Scenarion 1 of Situation 4',
                                                            styles={'font-weight': 'bold'}),
            'problem_condition': pn.widgets.RadioBoxGroup(name='',
                                value='Scenario 1 of Situation 4 like cases are not merged and are not shown.',
                                options=['Scenario 1 of Situation 4 like cases are not merged and are not shown.',
                                         'Scenario 1 of Situation 4 like cases are shown as problems to individually decide after merging.'],
                                inline=False, disabled=False),
        }

        # Control panel
        self.controls = pn.Param(self,
                                 parameters=['consider_formula_cols', 'text_problem_condition', 'problem_condition'],
                                 widgets=widgets, name='Parameters for Metabolic Feature Merging (De-Duplication)')


# Initialize metabolite annotation de-duplication storage class
data_ann_deduplicator = AnnDeDuplication_Storage()

# Widgets (mainly buttons) needed for the page layout

# To skip the section
skip_deduplication_button = pn.widgets.Button(name='Skip and Do Not Perform Annotation De-Duplication',
                                                     button_type='success', disabled=False)
# When pressing the button, skips this section and goes to next page - Data Pre-Treatment
def _skip_deduplication_button(event):
    "Goes to Data Pre-Treatment page without performing de-duplication."
    page3_button.disabled = False
    confirm_button_next_step_4.disabled = True
    save_data_dataframes_button.disabled = True
    # Join Filtered and Annotated dfs to make the original DataFrame for pre-treatment
    DataFrame_Store.original_df = DataFrame_Store.concat_annots(filtered_df.value, annotated_df.value)
    # Creates in the DataFrame a column ('Has Match?') that indicates if a feature was annotated either previously to being
    # inputted in this software or using the data annotation of this software
    iaf.creating_has_match_column(DataFrame_Store, n_databases, checkbox_annotation)

    # Disable posterior sidebar buttons
    page4_button.disabled = True
    # Disable statistical analysis
    _disabling_stat_analysis_buttons()

    # Informs De-Duplication was not performed for the report
    RepGen.deduplication_performed = 'Annotation De-Duplication was not performed'

    # Update to show the Data Pre-Treatment page
    main_area.clear()
    show_page(pages["Data Pre-Treatment"])
skip_deduplication_button.on_click(_skip_deduplication_button)

# To perform deduplication
perform_deduplication_button = pn.widgets.Button(name='Perform Annotation De-Duplication (Metabolic Feature Merging)',
                                                     button_type='success', disabled=False)
# When pressing the button, performs metabolic feature merging
def _perform_deduplication_button(event):
    "Performs multiple annotation metabolic feature merging and updates layout."

    # Disable posterior sidebar buttons
    page3_button.disabled = True
    page4_button.disabled = True
    see_merge_problems_button.disabled = False
    skip_merge_problem_cases_button.disabled = False
    # Disable statistical analysis
    _disabling_stat_analysis_buttons()

    # Updating layout
    while len(page2_1)>3:
        page2_1.pop(-1)

    # Loading Widget while de-duplication is happenning
    page2_1.append(pn.indicators.LoadingSpinner(value=True, size=90,
                                                        name='Perfoming Multiple Annotation Metabolic Feature Merging...'))

    # Perform Metabolic Feature Merging
    DataFrame_Store.original_df = DataFrame_Store.concat_annots(filtered_df.value, annotated_df.value)
    iaf.creating_has_match_column(DataFrame_Store, n_databases, checkbox_annotation)
    data_ann_deduplicator.perform_deduplication(DataFrame_Store.original_df, target_list.sample_cols,
                                            radiobox_neutral_mass.value)

    # Selecting merge problems to show
    if data_ann_deduplicator.problem_condition == 'Scenario 1 of Situation 4 like cases are not merged and are not shown.':
        data_ann_deduplicator.merge_problems = data_ann_deduplicator.full_merge_problems[
            data_ann_deduplicator.full_merge_problems['Nº of Peaks'] > 2]
    else:
        data_ann_deduplicator.merge_problems = data_ann_deduplicator.full_merge_problems
    data_ann_deduplicator.merge_situations.loc['Possible Problem Cases'] = len(data_ann_deduplicator.merge_problems)

    # Remove loading widget
    page2_1.pop(-1)

    # Middle section update
    middle_section_dedup[2] = data_ann_deduplicator.merge_report
    middle_section_dedup[3] = pn.Row(pn.pane.DataFrame(data_ann_deduplicator.merge_situations),
                                 pn.pane.DataFrame(data_ann_deduplicator.mergings_performed))
    if len(data_ann_deduplicator.merge_description) > 20:
        middle_section_dedup[5] = pn.pane.DataFrame(data_ann_deduplicator.merge_description, height=600)
    else:
        middle_section_dedup[5] = pn.pane.DataFrame(data_ann_deduplicator.merge_description)
    if len(data_ann_deduplicator.merge_problems) > 20:
        middle_section_dedup[7] = pn.pane.DataFrame(data_ann_deduplicator.merge_problems, height=600)
    else:
        middle_section_dedup[7] = pn.pane.DataFrame(data_ann_deduplicator.merge_problems)

    # Appending it
    page2_1.append(middle_section_dedup)
perform_deduplication_button.on_click(_perform_deduplication_button)

# Skip seeing problem cases
skip_merge_problem_cases_button = pn.widgets.Button(
    name='Skip Deciding on Individual Merging Problems (Will Remain Non-Merged)',
    button_type='success', disabled=False)
# When pressing the button, skips this section and goes to next page - Data Pre-Treatment
def _skip_merge_problem_cases_button(event):
    "Goes to Data Pre-Treatment page and saves metabolic feature merging performed."
    page3_button.disabled = False
    confirm_button_next_step_4.disabled = True
    save_data_dataframes_button.disabled = True
    # Pass merged DataFrame to DataFrame_Store
    DataFrame_Store.original_df = data_ann_deduplicator.annotated_df

    # Disable posterior sidebar buttons
    page4_button.disabled = True
    # Disable statistical analysis
    _disabling_stat_analysis_buttons()

    if RepGen.deduplication_performed != 'Annotation De-Duplication was performed and merge problems were individually observed/decided - Temp':
        # Informs De-Duplication was performed and merge problems were not looked at for the report
        RepGen.deduplication_performed = 'Annotation De-Duplication was performed but merge problems were not individually observed/decided'
    else:
        # Takes out the - Temp part of the string
        RepGen.deduplication_performed = 'Annotation De-Duplication was performed and merge problems were individually observed/decided'

    # Update to show the Data Pre-Treatment page
    main_area.clear()
    show_page(pages["Data Pre-Treatment"])
skip_merge_problem_cases_button.on_click(_skip_merge_problem_cases_button)

# See merging issues
see_merge_problems_button = pn.widgets.Button(name='Go Through Merge Problems Individually and Decide',
                                                     button_type='primary', disabled=False)
# When pressing the button, skips this section and goes to next page - Data Pre-Treatment
def _see_merge_problems_button(event):
    "Observe problem cases and updates layout."
    confirm_button_next_step_3.disabled = False

    # Disable posterior sidebar buttons
    page3_button.disabled = True
    page4_button.disabled = True
    see_merge_problems_button.disabled = True
    skip_merge_problem_cases_button.disabled = True
    # Disable statistical analysis
    _disabling_stat_analysis_buttons()

    # Updating layout
    while len(page2_1)>4:
        page2_1.pop(-1)

    # Merge Problems section update
    merge_problems_section_page.clear()

    # If there are Merge Problems
    if len(data_ann_deduplicator.merge_problems) != 0:
        merge_problems_section_page.append('### Merge Problem Cases')

        # Base for creating each merge problem part
        merge_problems_widgets.clear()
        merge_initial = pn.Column('#### Select Features to Merge',
                                        'Features with incompatible annotations cannot be merged (a notification will indicate this). **Note: MetaData columns labelled "Other" will be one of the values in the merged features randomly.**')

        # Creating the section for each merge problem
        for i in range(len(data_ann_deduplicator.merge_problems.index)):
            idx = data_ann_deduplicator.merge_problems.index[i]

            merge_problems_widgets[i] = pn.widgets.CheckBoxGroup(name='Peaks to Merge', value=[], options=idx,)
            ind_merge_problems = pn.GridSpec(mode='override')
            ind_merge_problems[0,0] = merge_initial
            ind_merge_problems[1,0] = merge_problems_widgets[i]
            ind_merge_problems[:, 1:4] = pn.pane.DataFrame(data_ann_deduplicator.annotated_df.loc[idx], height=400)

            merge_problems_section_page.append(pn.Column(f'### {data_ann_deduplicator.merge_problems.iloc[i, 1]}',
                                                    ind_merge_problems))

        confirm_mergeproblems_to_merge_button.disabled = False
        merge_problems_section_page.append(confirm_mergeproblems_to_merge_button)

        # Appending it
        page2_1.append(merge_problems_section_page)

    # If there are not merge problems
    else:
        merge_problems_section_page.append('### No Merge Problem Cases Detected')
        merge_problems_section_page.append(confirm_button_next_step_3)
        page2_1.append(merge_problems_section_page)
see_merge_problems_button.on_click(_see_merge_problems_button)

# Dictionary to store CheckBoxGroups Widgets for Problem Cases
merge_problems_widgets = {}

confirm_mergeproblems_to_merge_button = pn.widgets.Button(name='Perform Metabolic Feature Merging in Problem Cases of Selected Features (if 0 or only 1 feat. was chosen in a case, merging will not be attempted)',
                                description = 'If 0 or only 1 feature was selected in a case, merging will not be attempted. If merging cannot be performed with the selected features, a notification will appear.',
                                                     button_type='primary', disabled=False)
def _confirm_mergeproblems_to_merge_button(event):
    "Performs merging of selected peaks on merge problem cases, gives a report on the merging and updates the layout."
    # Obtaining the annotated columns
    ann_cols = list(checkbox_annotation.value)
    for key in data_ann_deduplicator.mcid:
        if key.startswith('Matched') and key.endswith(' IDs'):
            ann_cols.extend([key, 'Matched '+ key[8:-4] +' names', 'Matched '+ key[8:-4] +' formulas',
                                                key[8:-4] +' match count'])
    ann_cols.extend(checkbox_formula.value)
    merge_prob_description = pd.DataFrame(columns=data_ann_deduplicator.merge_description.columns)
    not_merged_desc = ['##### Metabolic Features not merged (from problem cases)']
    initial_len = len(data_ann_deduplicator.annotated_df)

    # Informs De-Duplication was performed and merge problems were looked at for the report
    RepGen.deduplication_performed = 'Annotation De-Duplication was performed and merge problems were individually observed/decided - Temp'

    # For each of the possible problems
    for i in merge_problems_widgets:
        idxs_to_merge = merge_problems_widgets[i].value
        # If more than one peak was selected try to merge them
        if len(idxs_to_merge) > 1:
            try:
                # If we can merge them
                df, m_d = iaf.individually_merging(data_ann_deduplicator, idxs_to_merge,
                             target_list.sample_cols,
                             ann_cols, radiobox_neutral_mass.value, mz_col=False)

                # Adjust everything
                data_ann_deduplicator.annotated_df = df
                m_d[list(m_d.keys())[0]]['DB'] = data_ann_deduplicator.merge_problems.iloc[i, 0]
                m_d[list(m_d.keys())[0]]['Repeating annotation'] = ' | '.join(data_ann_deduplicator.merge_problems.iloc[i, 1].values())
                merge_prob_description.loc[list(m_d.keys())[0]] = m_d[list(m_d.keys())[0]]

            except:
                m_d = f'{data_ann_deduplicator.merge_problems.iloc[i, 1]} with could not be merged: {idxs_to_merge} were incompatible.'
                not_merged_desc.append(m_d)

        else:
            m_d = f'{data_ann_deduplicator.merge_problems.iloc[i, 1]} was not merged.'
            not_merged_desc.append(m_d)

    merge_prob_description.index.name = 'New Indexes'

    # Updating the merging report
    if len(merge_prob_description) != 0:
        new_desc = data_ann_deduplicator.merge_report.split('\n')
        new_desc[0] = new_desc[0] + f' + **{len(merge_prob_description)}** (individually decided)'
        new_desc[1] = new_desc[1] + f' + **{merge_prob_description["Nº merged peaks"].sum()}** (individually decided)'
        new_desc[2] = new_desc[2] + f' + **{initial_len - len(data_ann_deduplicator.annotated_df)}** (individually decided)'
        new_desc[4] = new_desc[4].split('**')[0] + f'**{len(data_ann_deduplicator.annotated_df)}**'
        data_ann_deduplicator.merge_report = '\n'.join(new_desc)

    # Updating the merge situations and merge descriptions
    if len(merge_prob_description) != 0:
        m_p = pd.DataFrame(merge_prob_description['Situation'].value_counts())
        m_p.columns = ['Nº of Mergings by Situation']
        data_ann_deduplicator.merge_situations = pd.concat((data_ann_deduplicator.merge_situations, m_p))

        data_ann_deduplicator.merge_description = pd.concat((data_ann_deduplicator.merge_description, merge_prob_description))

    # Blocking the merge problem section so it cannot run without re-performing de-duplication
    for i in merge_problems_widgets:
        merge_problems_widgets[i].disabled = True
    confirm_mergeproblems_to_merge_button.disabled = True
    see_merge_problems_button.disabled = True
    skip_merge_problem_cases_button.disabled = True

    # Updating the layout of the page with the new merges performed
    while len(merge_problems_section_page) > (2+len(data_ann_deduplicator.merge_problems)):
        merge_problems_section_page.pop(-1)

    merge_problems_section_page.append('### Merge Problem Metabolic Feature Merging Description and Final Check')
    if len(merge_prob_description) != 0:
        merge_problems_section_page.append(pn.pane.DataFrame(merge_prob_description))
        merge_problems_section_page.append(data_ann_deduplicator.merge_report)
    if len(not_merged_desc) != 1:
        merge_problems_section_page.append('\n'.join(not_merged_desc))
    merge_problems_section_page.append(confirm_button_next_step_3)
confirm_mergeproblems_to_merge_button.on_click(_confirm_mergeproblems_to_merge_button)

# Setting up the general sections
# Middle section with details of metabolic peak merging performed
middle_section_dedup = pn.Column('### Details of Peak Merging Performed',
    '''Note: Metabolic features were merged first by the annotated columns indicated in the metadata, then by the annotation columns of the databases selected here and lastly by the formula columns (if selected).
    It is usual that the number of mergings decreases the further along the process we are since it is expected that many of the merged features to be merged in multiple databases.''',
                           data_ann_deduplicator.merge_report,
                          pn.Row(data_ann_deduplicator.merge_situations, data_ann_deduplicator.mergings_performed),
                        '#### Description of Every Metabolic Feature (Peak) Merged',
                          data_ann_deduplicator.merge_description,
                          '### Possible Problem Cases to see Individually',
                          data_ann_deduplicator.merge_problems,
                            pn.Row(skip_merge_problem_cases_button, see_merge_problems_button))

# End Section to see merge problems
merge_problems_section_page = pn.Column()

# Button to next step
confirm_button_next_step_3 = pn.widgets.Button(icon=iaf.img_confirm_button, name='Next Step - Data Pre-Treatment',
                                                     button_type='success', disabled=False)
confirm_button_next_step_3.on_click(_skip_merge_problem_cases_button)

page2_1 = pn.Column(pn.pane.HTML(desc_str.annotation_deduplication_opening_string),
                   data_ann_deduplicator.mult_ann_report,
                   pn.Row(skip_deduplication_button,
                          pn.Column(data_ann_deduplicator.controls, perform_deduplication_button)))




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
    norm_method = param.Selector(default="Total Intensity Sum")
    norm_kw = param.List(default=[])

    # Transformation
    tf_method = param.Selector(default="Generalized Logarithmic Transformation (glog)")
    tf_kw = param.Number(default=None)

    # Scaling
    scaling_method = param.Selector(default="Pareto Scaling")
    scaling_kw = param.String(default="Average")

    # Confirm Pre-Treatment Selection
    confirm_button = param.Boolean(default=False)


    # Function to confirm Pre-Treatment Selection and Updating DataFrames
    def _confirm_button_press(self, event):
        "Perform pre-treatment and update page layout."
        treat, proc, uni, meta, bin = iaf.performing_pretreatment(self, DataFrame_Store.original_df,
                                                                   target_list.target, target_list.sample_cols)
        DataFrame_Store.treated_df, DataFrame_Store.processed_df, DataFrame_Store.univariate_df = treat, proc, uni
        DataFrame_Store.metadata_df, DataFrame_Store.binsim_df = meta, bin

        # Locking in pre-treatment parameters chosen
        UnivarA_Store.locking_pretreatment_params(self)

        # Disable posterior sidebar buttons
        page4_button.disabled = True
        # Disable statistical analysis
        _disabling_stat_analysis_buttons()

        if len(DataFrame_Store.treated_df.columns) > 5000:
            page3[:2,2:5] = pn.Tabs(('Treated Data', pn.widgets.DataFrame(DataFrame_Store.treated_df, disabled=True,
                                                                        sortable=False, reorderable=False)),
                ('Metadata', pn.widgets.DataFrame(DataFrame_Store.metadata_df.T, disabled=True, sortable=False,
                                                reorderable=False)),
                ('BinSim Treated Data', pn.widgets.DataFrame(DataFrame_Store.binsim_df, disabled=True, sortable=False,
                                                reorderable=False)), height=600, dynamic=True)
        else:
            page3[:2,2:5] = pn.Tabs(('Treated Data', DataFrame_Store.treated_df),
                ('Metadata', DataFrame_Store.metadata_df.T),
                ('BinSim Treated Data', DataFrame_Store.binsim_df), height=600, dynamic=True)


        # See if there are nan values after pre-treatment
        n_nan_values = DataFrame_Store.treated_df.isnull().sum().sum()/(DataFrame_Store.treated_df.shape[0]*DataFrame_Store.treated_df.shape[1])
        # Remove extra parts of the page
        if len(page3[2, :].objects[(0, None, 1, None)]) > 4:
            page3[2, :].objects[(0, None, 1, None)].pop(0)

        # If there are no missing values as intended
        if n_nan_values == 0:
            confirm_button_next_step_4.disabled = False
            save_data_dataframes_button.disabled = False
            page13_button.disabled = False

        # If there are missing values, show an alert and inactivate further analysis
        else:
            alert_message = (f'The Dataset has **{n_nan_values*100:.2f} %** of NaN values after pre-treatment which conditions further statistical analysis.'
                '\nFor example, this may occur if you use **Zero** Missing Value Imputation with Generalized Logarithmic Transformation or Normalization by a Reference Feature that does not appear in every sample'
                ' (the latter is also not recommended even with other missing value imputations) but it is may not be limited to these combinations.'
                '\nTry to **change some parameters** of the Data Pre-Treament to solve this.')
            page3[2, :].objects[(0, None, 1, None)].insert(0, pn.pane.Alert(alert_message, alert_type='warning'))
            confirm_button_next_step_4.disabled = True
            save_data_dataframes_button.disabled = True
            page13_button.disabled = True


    def reset(self):
        "Resets all relevant parameters."
        for param in self.param:
            if param not in ["name"]:
                setattr(self, param, self.param[param].default)


    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets = {
            'mvi_method': pn.widgets.Select(name="Missing Value Imputation Method", value = 'Minimum of sample',
                            options=['Minimum of Sample', 'Minimum of Feature', 'Minimum of Data', 'Zero'],
                description='What Minimum to use as reference and what fraction of it for constant value imputation.'),
            'mvi_kw': pn.widgets.FloatSlider(name="Missing Value Imputation Fraction of Minimum Chosen",
                                             value=0.2, step=0.01),

            'norm_method': pn.widgets.Select(name="Normalization Method", value = 'Total Intensity Sum',
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
            'scaling_kw': pn.widgets.Select(name="Scaling Factor (for Level Scaling only)", value='Average', disabled=True,
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
    "Updates missing value imputation keyword widget limits based on the method chosen."
    if mvi_method == None:
        PreTreatment_Method.controls.widgets['mvi_kw'].value = 0
    PreTreatment_Method.controls.widgets['mvi_kw'].end = limits_mvi[
        PreTreatment_Method.controls.widgets['mvi_method'].value]

# Update options for normalization keyword based on normalization
options_norm = {"Reference Feature": None, "Total Intensity Sum": '', "PQN": ["mean", "median"],
                "Quantile": ["mean", "median"], 'None': ''}
@pn.depends(PreTreatment_Method.controls.widgets['norm_method'].param.value, watch=True)
def _update_options_norm(norm_method):
    "Updates normalization keyword widget limits based on the method chosen."
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
    "Updates scaling keyword widget limits based on the method chosen."
    if scaling_method == 'Level Scaling':
        PreTreatment_Method.controls.widgets['scaling_kw'].disabled = False
        PreTreatment_Method.controls.widgets['scaling_kw'].value = 'Average'

    else:
        PreTreatment_Method.controls.widgets['scaling_kw'].disabled = True
        PreTreatment_Method.controls.widgets['scaling_kw'].value = 'Average'

# Click button to confirm Pre-Treatment
PreTreatment_Method.controls.widgets['confirm_button'].on_click(PreTreatment_Method._confirm_button_press)


# Button to next step
confirm_button_next_step_4 = pn.widgets.Button(icon=iaf.img_confirm_button, name='Next Step - Class Colours',
                                                     button_type='success', disabled=True)

# Go to next step function and calling it
def _confirm_button_next_step_4(event):
    "Performs actions to pass from step 3 page to step 4 page."
    page4_button.disabled = False

    # Setting up the layout of page 4
    page4.clear()
    n_classes = len(target_list.color_classes)
    for row in range(0, n_classes, 5):
        new_row = pn.Row()
        n_max = min([5, n_classes-row])
        for col in range(n_max):
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
    "Saves Data main DataFrames."
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

    def reset(self):
        "Reset parameters."
        for param in self.param:
            if param not in ["name", "updater"]:
                setattr(self, param, self.param[param].default)

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
    "Performs actions to pass from step 4 page to transitional page, while preparing all statistical analysis (and making HCA and PCA analysis)."

    # Pass the colours picked to the target storage
    n_classes = len(target_list.color_classes)
    for row in range(0, n_classes, 5):
        n_max = min([5, n_classes-row])
        for col in range(n_max):
            key = list(target_list.color_classes.keys())[row+col]
            target_list.color_classes[key] = page4[row//5][col].value

    # Resetting all statistical analysis performed before (if repeating analysis)
    if reset_time.value == 0:
        if confirm_button_next_step_transitionalpage.clicks > 1:
            # Common and Exclusive Compound page . Has to be after Data diversity visualization reset
            com_exc_compounds.reset()
            while len(end_page_comexc) > 0:
                end_page_comexc.pop(-1)

            # Unsupervised Analysis page
            n_components_compute.value = 10
            PCA_params.reset()
            HCA_params.reset()

            # Supervised Analysis page
            plsda_feat_imp_show_annots_only.value = False
            plsda_feat_imp_show_annots_only.disabled = True
            rec_comp_indicator_widget.value = None
            PLSDA_store.soft_reset()
            pls_optim_section[1] = PLSDA_store.optim_figure[0]
            save_plsda_feat_imp_button.disabled = True
            while len(pls_results_section) > 5:
                pls_results_section.pop(-1)
            pls_results_section[0][1][1] = PLSDA_store.n_results
            pls_results_section[3] = pn.pane.DataFrame(PLSDA_store.feat_impor)

            rf_feat_imp_show_annots_only.value = False
            rf_feat_imp_show_annots_only.disabled = True
            RF_store.reset()
            rf_optim_section[1] = RF_store.optim_figure[0]
            save_rf_feat_imp_button.disabled = True
            while len(rf_results_section) > 5:
                rf_results_section.pop(-1)
            rf_results_section[0][1][1] = RF_store.n_results
            rf_results_section[3] = pn.pane.DataFrame(RF_store.feat_impor)

            # Univariate Analysis page
            UnivarA_Store.reset()
            while len(univar_analysis_page) > 2:
                univar_analysis_page.pop(-1)

            # Data Diversity Visualization page
            iaf._group_compounds_per_class(com_exc_compounds, target_list, DataFrame_Store) # Add temporary compounds per class dfs
            dataviz_store.reset()
            vk_plots.clear()
            kmd_plots.clear()
            while len(vk_page) > 1:
                vk_page.pop(-1)
            while len(kmd_page) > 1:
                kmd_page.pop(-1)
            while len(ccs_page) > 1:
                ccs_page.pop(-1)

            # Pathway Assignment page
            PathAssign_store.reset()
            while len(path_matching_section) > 3:
                path_matching_section.pop(-1)
            while len(hmdb_id_searching_section) > 3:
                hmdb_id_searching_section.pop(-1)

            # BinSim Analysis Page
            # PCA
            n_components_compute_binsim.value = 10
            PCA_params_binsim.reset()
            PCA_params_binsim.binsim_flag = True
            for _, w in PCA_params_binsim.controls.widgets.items():
                w.disabled = True
            middle_page_PCA_binsim[0,1:3] = 'To plot a PCA'
            end_page_PCA_binsim[0] = 'To plot explained variance figure'
            end_page_PCA_binsim[1] = 'To plot matrices of PCA projections'
            # HCA
            HCA_params_binsim.reset()
            HCA_params_binsim.binsim_flag = True
            # PLS-DA
            plsda_feat_imp_show_annots_only_binsim.value = False
            plsda_feat_imp_show_annots_only_binsim.disabled = True
            rec_comp_indicator_binsim_widget.value = None
            PLSDA_store_binsim.soft_reset()
            pls_optim_section_binsim[1] = PLSDA_store_binsim.optim_figure[0]
            save_plsda_feat_imp_binsim_button.disabled = True
            while len(pls_results_section_binsim) > 5:
                pls_results_section_binsim.pop(-1)
            pls_results_section_binsim[0][1][1] = PLSDA_store_binsim.n_results
            pls_results_section_binsim[3] = pn.pane.DataFrame(PLSDA_store_binsim.feat_impor)
            PLSDA_store_binsim.binsim_flag = True
            PLSDA_store_binsim.controls.widgets['scale'].disabled = True
            PLSDA_store_binsim.controls_optim.widgets['scale'].disabled = True
            # RF
            rf_feat_imp_show_annots_only_binsim.value = False
            rf_feat_imp_show_annots_only_binsim.disabled = True
            RF_store_binsim.reset()
            rf_optim_section_binsim[1] = RF_store_binsim.optim_figure[0]
            save_rf_feat_imp_binsim_button.disabled = True
            while len(rf_results_section_binsim) > 5:
                rf_results_section_binsim.pop(-1)
            rf_results_section_binsim[0][1][1] = RF_store_binsim.n_results
            rf_results_section_binsim[3] = pn.pane.DataFrame(RF_store_binsim.feat_impor)
            RF_store_binsim.binsim_flag = True

            # Compound Finder search tool page
            comp_finder.reset()
            comp_finder_page[1] = 'DataFrame of the Searched Compound'
            comp_finder_page[3] = 'Sample Bar Plot of the Searched Compound'
            comp_finder_page[5] = 'Class Bar Plot of the Searched Compound'
            comp_finder_page[7] = 'Class Boxplot of the Searched Compound'
            com_exc_compounds.reset()


    reset_time.value = 0

    # Enable all statistical analysis related buttons
    page5_button.disabled = False
    page6_button.disabled = False
    page7_button.disabled = False
    page8_button.disabled = False
    page9_button.disabled = False
    page10_button.disabled = False
    page11_button.disabled = False
    page12_button.disabled = False
    page13_button.disabled = False

    # Initial Calculations for PCA and storing initial plots
    principaldf, var, loadings = metsta.compute_df_with_PCs_VE_loadings(DataFrame_Store.treated_df,
                                       n_components=n_components_compute.value,
                                       whiten=True, labels=target_list.target, return_var_ratios_and_loadings=True)
    PCA_params.pca_scores = principaldf
    PCA_params.explained_variance = var
    PCA_params.pca_loadings = pd.DataFrame(loadings)
    PCA_params.PCA_plot[0] = iaf._plot_PCA(PCA_params, target_list)
    PCA_filename_string = 'PCA_plot'
    if PCA_params.ellipse_draw:
        if PCA_params.confidence != 0:
            PCA_filename_string = PCA_filename_string + f'_ellipse({PCA_params.confidence*100}%confidence)'
        else:
            PCA_filename_string = PCA_filename_string + f'_ellipse({PCA_params.confidence_std}std)'
    middle_page_PCA[0,1:3] = pn.pane.Plotly(PCA_params.PCA_plot[0], config = {'toImageButtonOptions': {'filename': PCA_filename_string, 'scale':4,}})
    PCA_params.exp_var_fig_plot[0] = iaf._plot_PCA_explained_variance(PCA_params)
    end_page_PCA[0] = pn.pane.Plotly(PCA_params.exp_var_fig_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_exp_var_plot', 'scale':4,}})
    PCA_params.scatter_PCA_plot[0] = iaf._scatter_PCA_plot(PCA_params, target_list)
    end_page_PCA[1] = pn.pane.Plotly(PCA_params.scatter_PCA_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_scatter_plot', 'scale':4,}})

    # Initial calculations for HCA and storing initial plots
    HCA_params.controls.widgets['dist_metric'].options = ['euclidean', 'cityblock', 'minkowski', 'seuclidean',
                    'sqeuclidean', 'braycurtis', 'canberra', 'chebyshev', 'correlation', 'cosine']
    HCA_params.dists = dist.pdist(DataFrame_Store.treated_df, metric=HCA_params.dist_metric)
    HCA_params.Z = hier.linkage(HCA_params.dists, method=HCA_params.link_metric)
    HCA_params.HCA_plot[0] = iaf._plot_HCA(HCA_params, target_list)
    page_HCA[0:6,1:4] = HCA_params.HCA_plot[0]

    # Updating Widgets for Univariate Analysis
    UnivarA_Store._update_widgets()

    # Updating Widgets for Supervised Analysis
    PLSDA_store.n_folds_limits_and_class_update(target_list)
    RF_store.n_folds_limits_and_class_update(target_list)

    # Updating Widgets for Pathway Assignment
    PathAssign_store._update_widgets(DataFrame_Store.metadata_df, pathway_db)

    # BinSim page updating widgets
    HCA_params_binsim.controls.widgets['dist_metric'].options = ['dice', 'hamming', 'jaccard', 'kulczynski1', 'rogerstanimoto',
                                                             'russellrao', 'sokalmichener','sokalsneath', 'yule']
    HCA_params_binsim.dist_metric = 'jaccard' # This change will also trigger the computation of the base Dendrogram
    PLSDA_store_binsim.n_folds_limits_and_class_update(target_list)
    RF_store_binsim.n_folds_limits_and_class_update(target_list)

    # Obtaining possibilities to search for based on type of identifier for the compound finder tool
    comp_finder._calculate_possible_options(DataFrame_Store.metadata_df, checkbox_annotation,
                                        checkbox_formula, radiobox_neutral_mass)

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
PathAssign_A = pn.widgets.Button(name='Pathway Assignment Analysis', button_type='default')
BinSim_A = pn.widgets.Button(name='BinSim Specific Analysis', button_type='default')
CompFinder_A = pn.widgets.Button(name='Compound Finder', button_type='default')
ToBeAdded_A = pn.widgets.Button(name='More to be added', button_type='danger', disabled=True)

# Functions for pressing each button
def _confirm_button_ComExc_A(event):
    "Button to go to the Common and Exclusive Compound page."
    main_area.clear()
    show_page(pages["Common and Exclusive Compounds"])
ComExc_A.on_click(_confirm_button_ComExc_A)

def _confirm_button_Unsup_A(event):
    "Button to go to the Unsupervised Analysis page."
    main_area.clear()
    show_page(pages["Unsupervised Analysis"])
Unsup_A.on_click(_confirm_button_Unsup_A)

def _confirm_button_Sup_A(event):
    "Button to go to the Supervised Analysis page."
    main_area.clear()
    show_page(pages["Supervised Analysis"])
Sup_A.on_click(_confirm_button_Sup_A)

def _confirm_button_Univariate_A(event):
    "Button to go to the Univariate Analysis page."
    main_area.clear()
    show_page(pages["Univariate Analysis"])
Univariate_A.on_click(_confirm_button_Univariate_A)

def _confirm_button_DataViz_A(event):
    "Button to go to the Data Diversity Visualization Analysis page."
    main_area.clear()
    show_page(pages["Data Visualization"])
DataViz_A.on_click(_confirm_button_DataViz_A)

def _confirm_button_PathAssign_A(event):
    "Button to go to the HMDB Pathway Matching/Assignment page."
    main_area.clear()
    show_page(pages["Pathway Assignment"])
PathAssign_A.on_click(_confirm_button_PathAssign_A)

def _confirm_button_BinSim_A(event):
    "Button to go to the BinSim Analysis page."
    main_area.clear()
    show_page(pages["BinSim Analysis"])
BinSim_A.on_click(_confirm_button_BinSim_A)

def _confirm_button_CompFinder_A(event):
    "Button to go to the Compound Finder Search Tool page."
    main_area.clear()
    show_page(pages["Compound Finder"])
CompFinder_A.on_click(_confirm_button_CompFinder_A)

# Transitional page Layout
transitional_page = pn.Column(pn.Row(ComExc_A, Unsup_A, Sup_A, Univariate_A, DataViz_A),
                             '## Other Options:',
                             pn.Row(PathAssign_A, BinSim_A, CompFinder_A, ToBeAdded_A))




# Page for Common and Exclusive Compounds

# Param Class to store parameters and data regarding Common and Exclusive Compounds
class ComExc_Storage(param.Parameterized):
    "Class to store all information on common and exclusive compounds and to plot Venn diagrams and Intersection Plots."

    # Dictionaries to group information
    groups = param.Dict()
    group_dfs = param.Dict()
    group_dfs_ids = param.Dict()
    groups_description = param.String(default='')

    # DataFrame with the common metabolites
    common_all = param.DataFrame()
    common_all_id = param.DataFrame()

    # Dictionary with DataFrames for exclusive compounds
    exclusives = param.Dict()
    exclusives_id = param.Dict()

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

    # Intersection Plot parameters
    inter_class_subset = param.List(default=list())
    inter_include_counts_percentages = param.String(default='Show Nº of metabolites')
    dpi_inter = param.Number(default=200)

    # Storing figure
    Venn_plot = param.List(default=['Pane for Venn Diagram'])
    IntersectionPlot = param.List(default=['Pane for Intersection Plot1', 'Pane for Intersection Plot2'])


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
        subsetdf_comexc_section_page[0:4,1:3] = pn.pane.DataFrame(self.specific_cl_df, height=500)
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


    # Update the Intersection Plot based on specifications chosen
    @param.depends('inter_class_subset', 'inter_include_counts_percentages', 'dpi_inter', watch=True)
    def _update_Intersectionplot(self):
        "Update the Intersection Plots based on parameters given."
        if 1 < len(self.inter_class_subset):
            self.IntersectionPlot = []

            # Select the relevant data - all metabolites and annotated metabolites
            groups_dict = {}
            groups_dict_ids = {}
            for df in self.inter_class_subset:
                groups_dict[df] = self.group_dfs[df].index
                groups_dict_ids[df] = self.group_dfs_ids[df].index

            ups = from_contents(groups_dict)
            self.IntersectionPlot.append(iaf._plot_intersection_plots(self, groups_dict, ups))
            interplot_page[2] = pn.pane.Matplotlib(self.IntersectionPlot[0], height=300)
            plt.close()

            ups_ids = from_contents(groups_dict_ids)

            if sum(DataFrame_Store.metadata_df['Has Match?']) != 0:
                self.IntersectionPlot.append(iaf._plot_intersection_plots(self, groups_dict_ids, ups_ids))

                interplot_page[5] = pn.pane.Matplotlib(self.IntersectionPlot[1], height=300)
            else:
                interplot_page[5] = 'No annotations have been found in the dataset to compute this Intersection Plot.'
            plt.close()
        else:
            pn.state.notifications.info(f'Intersection Plot can only be made with 2 or more classes. You currently have {len(self.inter_class_subset)} classes.')


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
            'inter_class_subset': pn.widgets.CheckBoxGroup(name='Classes', value=target_list.classes,
                                options=target_list.classes, inline=False, disabled=False),
            'inter_include_counts_percentages': pn.widgets.RadioBoxGroup(name='Show:',
                                value='Show Nº of metabolites', options=['Show Nº and % of metabolites',
                                'Show Nº of metabolites', 'Show neither'], inline=False, disabled=False),
            'dpi_inter': pn.widgets.IntInput(name="DPI (Resolution)", value=200, step=10, start=100, disabled=False,
                                    description='Set the resolution of Intersection Plots'),
        }
        self.class_subset = []
        #self.venn_class_subset = target_list.classes
        #self.inter_class_subset = target_list.classes

        # Control panel for the overview section, for the Venn diagram section and for the Intersection plot section
        return (pn.Param(self, parameters=['class_subset', 'df_type', 'annot'], widgets=widgets, name='Subset of Data to See'),
                pn.Param(self, parameters=['venn_class_subset', 'venn_alpha', 'type_of_venn', 'dpi_venn'], widgets=widgets,
                         name='Parameters to draw Venn Diagram'),
                pn.Param(self, parameters=['inter_class_subset', 'inter_include_counts_percentages','dpi_inter'],
                         widgets=widgets, name='Parameters to draw Intersection Plots', default_layout=pn.Row))


    def reset(self):
        "Reset parameters."
        for param in self.param:
            if param not in ["name"]:
                setattr(self, param, self.param[param].default)


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
            'inter_class_subset': pn.widgets.CheckBoxGroup(name='Classes', value=target_list.classes,
                                options=target_list.classes, inline=False, disabled=False),
            'inter_include_counts_percentages': pn.widgets.RadioBoxGroup(name='Show:',
                                value='Show Nº of metabolites', options=['Show Nº and % of metabolites',
                                'Show Nº of metabolites', 'Show neither'], inline=False, disabled=False),
            'dpi_inter': pn.widgets.IntInput(name="DPI (Resolution)", value=200, step=10, start=100, disabled=False,
                                    description='Set the resolution of Intersection Plots'),
        }

        # Control panel for the overview section
        self.controls = pn.Param(self,
                                 parameters=['class_subset', 'df_type', 'annot'],
                                 widgets=widgets, name='Subset of Data to See')
        # Control panel for the Venn diagram section
        self.venn_controls = pn.Param(self,
                                 parameters=['venn_class_subset', 'venn_alpha', 'type_of_venn', 'dpi_venn'],
                                 widgets=widgets, name='Parameters to draw Venn Diagram')
        # Control panel for the Intersection Plots section
        self.interplot_controls = pn.Param(self, parameters=['inter_class_subset', 'inter_include_counts_percentages',
                                                             'dpi_inter'],
                                 widgets=widgets, name='Parameters to draw Intersection Plots', default_layout=pn.Row)


# Initialize common and exclusive compound storage
com_exc_compounds = ComExc_Storage()

# Widgets for the page
compute_ComExc_button = pn.widgets.Button(name='Compute', button_type='success', height=50)
checkbox_com_exc = pn.widgets.CheckBoxGroup(name='Include:', value=['Venn Diagram', 'Intersection Plot'],
                                            options=['Venn Diagram', 'Intersection Plot'], inline=True)

# When pressing the button, performs common and exclusive compound calculations, and sets up the different pages in the tabs section
def _compute_ComExc_button(event):
    "Actions after pressing the button"

    iaf._group_compounds_per_class(com_exc_compounds, target_list, DataFrame_Store) # Add compounds per class dfs
    iaf._compute_com_exc_compounds(com_exc_compounds) # Add common to all and exclusive to each class compounds dfs

    new_controls = com_exc_compounds._update_widgets() # Update Widgets
    com_exc_compounds.controls = new_controls[0]
    com_exc_compounds.venn_controls = new_controls[1]
    com_exc_compounds.interplot_controls = new_controls[2]

    subsetdf_comexc_section_page[0:2,0] = com_exc_compounds.controls
    venn_page[0,0][0] = com_exc_compounds.venn_controls
    interplot_page[0] = com_exc_compounds.interplot_controls

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
        com_exc_compounds.venn_class_subset = target_list.classes
        end_page_comexc.append(('Venn Diagram', venn_page))
    if 'Intersection Plot' in checkbox_com_exc.value:
        com_exc_compounds.inter_class_subset = target_list.classes
        end_page_comexc.append(('Intersection Plot', interplot_page))

    com_exc_compounds._update_specific_cl_df() # Update the DataFrame shown in the overview tab
    #com_exc_compounds.venn_class_subset = com_exc_compounds.venn_controls[1].value # Updating starting subset value for Venn diagram
    #com_exc_compounds.inter_class_subset = com_exc_compounds.interplot_controls[0][1].value # Updating starting subset value for Intersection Plots


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
    "Downloads common (or exclusive) compounds of selected classes with selected characteristics DataFrame."
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
subsetdf_comexc_section_page[0:4,1:3] = pn.pane.DataFrame(com_exc_compounds.specific_cl_df, height=500)


# Widget to save Venn Diagram (needed since it is a matplotlib plot instead of a plotly plot)
save_Venn_diag_button = pn.widgets.Button(name='Save as a png (in current folder)', button_type='success',
                                         icon=iaf.download_icon)
# When pressing the button, downloads the figure
def _save_Venn_diag_button(event):
    "Save Venn Diagram plot."
    filename_string = f'Venn_diagram_{com_exc_compounds.type_of_venn}_classes'
    for cl in com_exc_compounds.venn_class_subset:
        filename_string = filename_string + '_'+cl
    com_exc_compounds.Venn_plot[0].savefig(filename_string, dpi=com_exc_compounds.dpi_venn)
    pn.state.notifications.success(f'Venn Diagram successfully saved.')
save_Venn_diag_button.on_click(_save_Venn_diag_button)

# Create specific section of the page for the Venn Diagram tab
venn_page = pn.GridSpec(mode='override')
venn_page[0,0] = pn.Column(com_exc_compounds.venn_controls, save_Venn_diag_button)
venn_page[0,1:3] = com_exc_compounds.Venn_plot[0]


# Widgets to save Intersection Plots (needed since they are matplotlib plots instead of a plotly plots)
save_IntersectionPlot_allmets_button = pn.widgets.Button(name='Save as a png (in current folder)', button_type='success',
                                         icon=iaf.download_icon)
save_IntersectionPlot_annotatedmets_button = pn.widgets.Button(name='Save as a png (in current folder)', button_type='success',
                                         icon=iaf.download_icon)
# When pressing the button, downloads the figures
def _save_IntersectionPlot_allmets_button(event):
    "Save Intersection Plot with all metabolites."
    filename_string = 'IntersectionPlot_all_metabolites_classes'
    for cl in com_exc_compounds.inter_class_subset:
        filename_string = filename_string + '_'+cl
    com_exc_compounds.IntersectionPlot[0].savefig(filename_string + '.png', dpi=com_exc_compounds.dpi_inter)
    pn.state.notifications.success(f'Intersection Plot (all metabolites) successfully saved.')
save_IntersectionPlot_allmets_button.on_click(_save_IntersectionPlot_allmets_button)

def _save_IntersectionPlot_annotatedmets_button(event):
    "Save Intersection Plot with annotated metabolites."
    filename_string = 'IntersectionPlot_annotated_metabolites_classes'
    for cl in com_exc_compounds.inter_class_subset:
        filename_string = filename_string + '_'+cl
    com_exc_compounds.IntersectionPlot[1].savefig(filename_string + '.png', dpi=com_exc_compounds.dpi_inter)
    pn.state.notifications.success(f'Intersection Plot (annotated metabolites) successfully saved.')
save_IntersectionPlot_annotatedmets_button.on_click(_save_IntersectionPlot_annotatedmets_button)

# Create specific section of the page for the Intersection Plot tab
interplot_page = pn.Column(com_exc_compounds.interplot_controls, # Control parameters for Intersection Plot
                      '## Intersection Plot with all metabolites of the dataset',
                       com_exc_compounds.IntersectionPlot[0],
                       save_IntersectionPlot_allmets_button,
                       '## Intersection Plot with only annotated metabolites in the dataset',
                       com_exc_compounds.IntersectionPlot[1],
                       save_IntersectionPlot_annotatedmets_button)


end_page_comexc = pn.Accordion(toggle=True)

comexc_page = pn.Column(initial_page_comexc, end_page_comexc)




# Page for Unsupervised Analysis

# Tab for PCA analysis
# TODO: PCA is straight away computed with 10 components. If your data has less than 10 features, this will lead to an error and everything will fail.
# Should this possibility be taken into account?

# Param Class to store parameters and data regarding PCA
class PCA_Storage(param.Parameterized):
    "Class use to store PCA and PCA related plots parameters, dataframes and figures."

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

    # BinSim Flag
    binsim_flag = param.Boolean(default=False)
    current_pages_associated = param.List(default=[])

    # Update the PCA plot
    @param.depends('n_dimensions', 'PCx', 'PCy', 'PCz', 'ellipse_draw', 'confidence', 'confidence_std', 'dot_size', watch=True)
    def _update_PCA_plot(self):
        "Plots PCA with the set parameters."
        self.PCA_plot[0] = iaf._plot_PCA(self, target_list)
        filename_string = 'PCA_plot'
        if self.ellipse_draw:
            if self.confidence != 0:
                filename_string = filename_string + f'_ellipse({self.confidence*100}%confidence)'
            else:
                filename_string = filename_string + f'_ellipse({self.confidence_std}std)'
        if self.binsim_flag:
            filename_string = filename_string + f'_BinSim'
        self.current_pages_associated[0][0,1:3] = pn.pane.Plotly(self.PCA_plot[0], config = {'toImageButtonOptions': {'filename': filename_string, 'scale':4}})


    # Function to see if Z-axis can be edited (3D) or not (2D)
    @param.depends('n_dimensions', watch=True)
    def _update_PCz_disabled(self):
        "Controls PCz widget in PCA plots."
        if self.n_dimensions == '2 Components':
            self.controls.widgets['PCz'].disabled = True
        elif self.n_dimensions == '3 Components':
            self.controls.widgets['PCz'].disabled = False


    # Function updating the possible components to choose based on number of components of the PCA and updating the explained variance plot
    @param.depends('n_components', watch=True)
    def _update_PC_options(self):
        "Controlling options and plots based on the number of components PCA was computed with."
        # Controlling widget options for components to show in the projection plot
        if int(self.PCx[3:]) > self.n_components:
            self.PCx = 'PC 1'
        if int(self.PCy[3:]) > self.n_components:
            self.PCy = 'PC 1'
        if int(self.PCz[3:]) > self.n_components:
            self.PCz = 'PC 1'
        self.controls.widgets['PCx'].options = ['PC '+str(i+1) for i in range(self.n_components)]
        self.controls.widgets['PCy'].options = ['PC '+str(i+1) for i in range(self.n_components)]
        self.controls.widgets['PCz'].options = ['PC '+str(i+1) for i in range(self.n_components)]

        # Updating the explained variance figure with the number of components PCA was computed
        self.exp_var_fig_plot[0] = iaf._plot_PCA_explained_variance(self)
        filename_string = 'PCA_exp_var_plot'
        if self.binsim_flag:
            filename_string = filename_string + f'_BinSim'
        self.current_pages_associated[1][0] = pn.pane.Plotly(self.exp_var_fig_plot[0],
                                            config = {'toImageButtonOptions': {'filename': filename_string, 'scale':4,}})

        # Control the many PCA scatter plots to include less than 6 components if PCA was computed with less than 6 components
        if type(self.scatter_PCA_plot[0]) != str:
            if len(self.scatter_PCA_plot[0].data[0]['dimensions']) < 6:
                self.scatter_PCA_plot[0] = iaf._scatter_PCA_plot(self, target_list)
                filename_string = 'PCA_scatter_plot'
                if self.binsim_flag:
                    filename_string = filename_string + f'_BinSim'
                self.current_pages_associated[1][1] = pn.pane.Plotly(self.scatter_PCA_plot[0],
                                            config = {'toImageButtonOptions': {'filename': filename_string, 'scale':4,}})
            else:
                if self.n_components < 6:
                    self.scatter_PCA_plot[0] = iaf._scatter_PCA_plot(self, target_list)
                    filename_string = 'PCA_scatter_plot'
                    if self.binsim_flag:
                        filename_string = filename_string + f'_BinSim'
                    self.current_pages_associated[1][1] = pn.pane.Plotly(self.scatter_PCA_plot[0], config = {'toImageButtonOptions': {'filename': filename_string, 'scale':4,}})


    # Function enabling/disabling the confidence level parameters for ellipses
    @param.depends('ellipse_draw', watch=True)
    def _update_ellipse_options(self):
        "Controls ellipse widgets based on if they are drawn or not."
        if self.ellipse_draw:
            self.controls.widgets['confidence'].disabled = False
            if self.confidence == 0:
                self.controls.widgets['confidence_std'].disabled = False
        else:
            self.controls.widgets['confidence'].disabled = True
            self.controls.widgets['confidence_std'].disabled = True


    # Function enabling/disabling the confidence level based on std parameter for ellipses
    @param.depends('confidence', watch=True)
    def _update_ellipse_std_options(self):
        "Controls ellipse based on std. confidence widget based on the usual confidence widget."
        if self.confidence == 0:
            self.controls.widgets['confidence_std'].disabled = False
        else:
            self.controls.widgets['confidence_std'].disabled = True


    def reset(self):
        "Reset parameters."
        for param in self.param:
            if param not in ["name", "current_pages_associated"]:
                setattr(self, param, self.param[param].default)


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
    "Computes PCA."
    if PCA_params.binsim_flag:
        df = DataFrame_Store.binsim_df
    else:
        df = DataFrame_Store.treated_df
    principaldf, var, loadings = metsta.compute_df_with_PCs_VE_loadings(df,
                                       n_components=n_components_compute.value,
                                       whiten=True, labels=target_list.target, return_var_ratios_and_loadings=True)
    PCA_params.pca_scores = principaldf
    PCA_params.explained_variance = var
    PCA_params.pca_loadings = pd.DataFrame(loadings)
    PCA_params.n_components = n_components_compute.value
    PCA_params.controls.widgets['n_components'].value = n_components_compute.value

compute_PCA_button.on_click(_compute_PCA_button)



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

# Page sections that will change based on HCA_params
PCA_params.current_pages_associated.append(middle_page_PCA)
PCA_params.current_pages_associated.append(end_page_PCA)


page_PCA = pn.Column(initial_page_PCA, middle_page_PCA, end_page_PCA)



# Tab for HCA analysis
# Param Class to store parameters and data regarding HCA
# TODO: Each time you modify the HCA, a new matplotlib plot is created, thus this could lead to memory issues with many manipulations
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

    # BinSim Flag
    binsim_flag = param.Boolean(default=False)
    current_pages_associated = param.List(default=[])


    # Update the HCA plot
    @param.depends('dist_metric', 'link_metric', watch=True)
    def _update_HCA_compute_plot(self):
        "Computes and updates HCA based on distance and linkage metrics."
        # Recompute Distance and Linkage matrices
        if self.binsim_flag:
            self.dists = dist.pdist(DataFrame_Store.binsim_df, metric=self.dist_metric)
        else:
            self.dists = dist.pdist(DataFrame_Store.treated_df, metric=self.dist_metric)
        self.Z = hier.linkage(self.dists, method=self.link_metric)
        # Plot the new HCA
        self.HCA_plot[0] = iaf._plot_HCA(self, target_list)
        self.current_pages_associated[0][0:6,1:4] = pn.pane.Matplotlib(self.HCA_plot[0], dpi=self.dpi)


    @param.depends('fig_text', 'fig_x', 'fig_y', 'col_threshold', 'dpi', watch=True)
    def _update_HCA_plot(self):
        "Updates HCA plot based on figure parameters."
        if type(self.Z) == np.ndarray:
            self.HCA_plot[0] = iaf._plot_HCA(self, target_list)
            self.current_pages_associated[0][0:6,1:4] = pn.pane.Matplotlib(self.HCA_plot[0], dpi=self.dpi)


    def reset(self):
        "Reset parameters."
        for param in self.param:
            if param not in ["name", "current_pages_associated"]:
                setattr(self, param, self.param[param].default)


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
            'col_threshold': pn.widgets.FloatInput(name="Color Threshold",
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


# Widget to save HCA plot (needed since it is a matplotlib plot instead of a plotly plot)
save_HCA_plot_button = pn.widgets.Button(name='Save as a png (in current folder)', button_type='success',
                                         icon=iaf.download_icon)
# When pressing the button, downloads the figure
def _save_HCA_plot_button(event):
    "Saves HCA plot."
    filename = f'HCA_plot_{HCA_params.dist_metric}Dist_{HCA_params.link_metric}Linkage.png'
    HCA_params.HCA_plot[0].savefig(filename, dpi=HCA_params.dpi)
    pn.state.notifications.success(f'Dendrogram successfully saved.')
save_HCA_plot_button.on_click(_save_HCA_plot_button)

# Organization for the HCA page
page_HCA = pn.GridSpec(mode='override')
page_HCA[0:5,0] = HCA_params.controls
page_HCA[0:6,1:4] = HCA_params.HCA_plot[0]
page_HCA[5,0] = save_HCA_plot_button

# Page sections that will change based on HCA_params
HCA_params.current_pages_associated.append(page_HCA)

# Page with PCA and HCA analysis
unsup_analysis_page = pn.Tabs(('PCA', page_PCA), ('HCA', page_HCA))




# Page for Supervised Analysis

# Section for the PLS-DA analysis
plsda_opening_string = desc_str.plsda_opening_string


# Param Class to store parameters and data regarding PLS-DA
class PLSDA_Storage(param.Parameterized):
    "Class to store PLS-DA models, parameters and results."

    # PLS-DA flag and dataset to treat
    binsim_flag = param.Boolean(default=False)
    current_pages_associated = param.List()

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
    current_other_plsda_params = param.Dict({})

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
    positive_class = param.String(default='')
    roc_n_iter = param.Number(default=10)
    confirm_button_roc = param.Boolean(default=False)

    # Storing figures
    optim_figure = param.List(default=['To Plot the Optimization PLS Figure'])
    PLS_plot = param.List(default=['To Plot the PLS Projection Figure'])
    ROC_figure = param.List(default=['To Plot the ROC Curve(s)'])
    perm_figure = param.List(default=['To Plot the Permutation Test Figure'])

    # Update number of folds limits
    def n_folds_limits_and_class_update(self, target_list):
        "Updates the limits of the CV fold number based on the number of samples per class and classes to choose for ROC curve."
        # Updating the end value limits
        min_samples_in_class = pd.Series(target_list.target).value_counts().min()
        self.controls_optim.widgets['n_fold'].end = min_samples_in_class
        self.controls.widgets['n_fold'].end = min_samples_in_class

        # Updating the proper value if it is above the highest possible number
        if min_samples_in_class < self.controls_optim.widgets['n_fold'].value:
            self.controls_optim.widgets['n_fold'].value = min_samples_in_class
            self.controls.widgets['n_fold'].value = min_samples_in_class
            self.n_fold = min_samples_in_class

        # Updating the classes available to choose as the positive for ROC curve
        if len(target_list.classes) == 2:
            self.controls_roc.widgets['positive_class'].options = target_list.classes
            self.controls_roc.widgets['positive_class'].value = target_list.classes[0]
            self.controls_roc.widgets['positive_class'].disabled = False
            self.positive_class = target_list.classes[0]
        else:
            self.controls_roc.widgets['positive_class'].options = ['Not Available']
            self.controls_roc.widgets['positive_class'].value = 'Not Available'
            self.controls_roc.widgets['positive_class'].disabled = True
            self.positive_class = 'Not Available'


    # Function to confirm the optimization
    def _confirm_optim_button(self, event):
        "Performs optimization, plots the optimization figure and updated corresponding layout."

        # Loading Widget while the PLS-DA optimization is being performed
        self.current_pages_associated[0][1] = pn.indicators.LoadingSpinner(value=True, size=90,
                                                            name='Performing Optimization of PLS-DA Components...')

        # Round to integer the values in the range
        self.n_min_max_components = int(np.round(self.n_min_max_components[0])), int(np.round(self.n_min_max_components[1]))

        # See what dataset to use
        if self.binsim_flag:
            plsda_df = DataFrame_Store.binsim_df
        else:
            plsda_df = DataFrame_Store.treated_df

        # Performs optimization
        PLS_optim = metsta.optim_PLSDA_n_components(plsda_df, target_list.target, # Data and target
                                encode2as1vector=True,
                                max_comp=self.n_min_max_components[1], # Max. number of components to search
                                min_comp=self.n_min_max_components[0], # Min. number of components to search
                                kf=None, n_fold=self.n_fold, # Nº of folds
                                scale=self.scale)

        # Saving parameters for optimization
        self.current_other_plsda_params.update({'n_fold_optim': self.n_fold, 'scale_optim': self.scale,
                                                'n_min_max_components': self.n_min_max_components})

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
        self.optim_figure = [fig,]
        self.rec_components = rec_comp
        # Update the layout
        filename_string = f'PLS_optim_plot_{self.n_min_max_components[0]}to{self.n_min_max_components[1]}components_'
        filename_string = filename_string + f'{self.n_fold}-foldstratCV_scale{self.scale}'
        if self.binsim_flag:
            filename_string = filename_string + '_BinSim'
        self.current_pages_associated[0][1] = pn.pane.Plotly(self.optim_figure[0],
                                              config = {'toImageButtonOptions': {'filename': filename_string, 'scale':4,}})
        self.current_pages_associated[0][0][1].value = self.rec_components


    # Function to fit the PLS-DA model and obtain results
    def _confirm_plsda_button(self, event):
        "Fits a PLS-DA model, gives model performance estimations and feature importance lists, updating the layout."

        # Loading Widget while the PLS-DA model is being fitted and assessed
        self.current_pages_associated[1][0][1][1] = pn.indicators.LoadingSpinner(value=True, size=90,
                                                            name='Fitting and Assessing PLS-DA model...')

        # Sees what type of feature importance metric will be used
        if self.imp_feature_metric == 'VIP':
            feat_type = 'VIP'
        elif self.imp_feature_metric == 'Coefficients':
            feat_type = 'Coef'
        else:
            feat_type = 'Weights'

        # See what dataset to use
        if self.binsim_flag:
            plsda_df = DataFrame_Store.binsim_df
        else:
            plsda_df = DataFrame_Store.treated_df

        # Fits a PLS-DA model under a stratified cross validation scheme
        PLSDA_results = metsta.PLSDA_model_CV(plsda_df, target_list.target, # Data and target
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
        self.current_pages_associated[1][0][1][1] = pn.pane.DataFrame(self.n_results)
        self.current_pages_associated[1][3] = pn.pane.DataFrame(self.feat_impor, height=600)
        self.current_pages_associated[3].value = False
        self.current_pages_associated[3].disabled = False
        self.current_pages_associated[4].disabled = False

        # Fit a PLS-DA model with all samples and store model and x_scores
        self.models[0], self.x_scores = ma.fit_PLSDA_model(
            plsda_df, target_list.target, # Data and target
            n_comp=self.n_components, return_scores=True,
            scale=self.scale,
            encode2as1vector=True, lv_prefix='LV ', label_name='Label')
        self.complete_soft_reset()

        # Name of the file
        filename_string = f'PLS_plot_({self.n_components}comp)'
        if self.n_dimensions == '2 Components':
            if self.ellipse_draw:
                if self.confidence != 0:
                    filename_string = filename_string + f'_ellipse({self.confidence*100}%confidence)'
                else:
                    filename_string = filename_string + f'_ellipse({self.confidence_std}std)'
        if self.binsim_flag:
            filename_string = filename_string + '_BinSim'

        self.PLS_plot[0] = iaf._plot_PLS(self, target_list)
        self.current_pages_associated[2][0,1:3] = pn.pane.Plotly(self.PLS_plot[0], config={'toImageButtonOptions': {'filename': filename_string, 'scale':4}})

        if len(self.current_pages_associated[1]) == 5:
            self.current_pages_associated[1].append('### PLS Projection Section')
            self.current_pages_associated[1].append(self.current_pages_associated[2])
            self.current_pages_associated[1].append('### PLS-DA Permutation Test Section')
            self.current_pages_associated[1].append(pn.Column(pn.pane.HTML(permutation_test_description),
                                                pn.pane.LaTeX(pvalue_equation_string, styles={'font-size': '14pt'})))
            self.current_pages_associated[1].append(pn.Row(self.controls_permutation,
                                        self.perm_figure[0]))
            self.current_pages_associated[1].append('### PLS-DA Receiver Operating Characteristic (ROC) Curve Section')
            self.current_pages_associated[1].append(pn.pane.HTML(ROC_curve_description))
            self.current_pages_associated[1].append(pn.Row(self.controls_roc, self.ROC_figure[0]))


    # Update the PLS Projection plot
    @param.depends('n_dimensions', 'LVx', 'LVy', 'LVz', 'ellipse_draw', 'confidence', 'confidence_std', 'dot_size', watch=True)
    def _update_PLS_plot(self):
        "Updates PLS plot based on figure parameters."
        # Name of the file
        filename_string = f'PLS_plot_({self.n_components}comp)'
        if self.n_dimensions == '2 Components':
            if self.ellipse_draw:
                if self.confidence != 0:
                    filename_string = filename_string + f'_ellipse({self.confidence*100}%confidence)'
                else:
                    filename_string = filename_string + f'_ellipse({self.confidence_std}std)'
        if self.binsim_flag:
            filename_string = filename_string + '_BinSim'

        self.PLS_plot[0] = iaf._plot_PLS(self, target_list)
        self.current_pages_associated[2][0,1:3] = pn.pane.Plotly(self.PLS_plot[0], config={'toImageButtonOptions': {'filename': filename_string, 'scale':4}})


    # Set of Functions controlling parameters for PLS Projection
    # Function to see if Z-axis can be edited (3D) or not (2D)
    @param.depends('n_dimensions', watch=True)
    def _update_LVz_disabled(self):
        "Controls LVz widget in PLS plots."
        if self.n_dimensions == '2 Components':
            self.controls_projection.widgets['LVz'].disabled = True
        elif self.n_dimensions == '3 Components':
            self.controls_projection.widgets['LVz'].disabled = False


    # Function updating the possible Latent variables to choose based on number of LVs of the PLS-DA
    @param.depends('n_components', watch=True)
    def _update_LV_options(self):
        "Updates LV options to choose."
        if int(self.LVx[3:]) > self.n_components:
            self.LVx = 'LV 1'
        if int(self.LVy[3:]) > self.n_components:
            self.LVy = 'LV 1'
        if int(self.LVz[3:]) > self.n_components:
            self.LVz = 'LV 1'
        self.controls_projection.widgets['LVx'].options = ['LV '+str(i+1) for i in range(self.n_components)]
        self.controls_projection.widgets['LVy'].options = ['LV '+str(i+1) for i in range(self.n_components)]
        self.controls_projection.widgets['LVz'].options = ['LV '+str(i+1) for i in range(self.n_components)]


    # Function enabling/disabling the confidence level parameters for ellipses
    @param.depends('ellipse_draw', watch=True)
    def _update_plsda_ellipse_options(self):
        "Controls ellipse widgets based on if they are drawn or not."
        if self.ellipse_draw:
            self.controls_projection.widgets['confidence'].disabled = False
            if self.confidence == 0:
                self.controls_projection.widgets['confidence_std'].disabled = False
        else:
            self.controls_projection.widgets['confidence'].disabled = True
            self.controls_projection.widgets['confidence_std'].disabled = True


    # Function enabling/disabling the confidence level based on std parameter for ellipses
    @param.depends('confidence', watch=True)
    def _update_plsda_ellipse_std_options(self):
        "Controls ellipse based on std. confidence widget based on the usual confidence widget."
        if self.confidence == 0:
            self.controls_projection.widgets['confidence_std'].disabled = False
        else:
            self.controls_projection.widgets['confidence_std'].disabled = True


    # Function to confirm the permutation test to perform
    def _confirm_button_permutation(self, event):
        "Performs a permutation test for the PLS-DA model and updates the layout."

        # Loading Widget while the Permutation Test is being performed
        self.current_pages_associated[1][9][1] = pn.indicators.LoadingSpinner(value=True, size=90,
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

        # See what dataset to use
        if self.binsim_flag:
            plsda_df = DataFrame_Store.binsim_df
        else:
            plsda_df = DataFrame_Store.treated_df

        # Performs the Permutation Test on the PLS-DA model under a stratified cross validation scheme
        perm_results_PLSDA = metsta.permutation_PLSDA(
            plsda_df, target_list.target, # Data and target
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
        perm_figure = iaf._plot_permutation_test(perm_results_PLSDA, plsda_df, self.n_fold,
                                                     self.perm_metric, 'PLS-DA Permutation Test')
        self.perm_figure =[perm_figure,]

        # Update the layout
        self.current_pages_associated[1][9][1] = pn.pane.Matplotlib(self.perm_figure[0], height=600)
        self.controls_permutation.widgets['save_figure_button_permutation'].disabled = False


    # Function to confirm the ROC Curves to plot
    def _confirm_button_roc(self, event):
        "Plots a ROC Curve to assess PLS-DA model performance and updates the layout."

        # Loading Widget while the ROC Curves are being computed
        self.current_pages_associated[1][12][1] = pn.indicators.LoadingSpinner(value=True, size=90,
                                                            name='Computing Receiver Operating Characteristic Curves...')

        # See what dataset to use
        if self.binsim_flag:
            plsda_df = DataFrame_Store.binsim_df
        else:
            plsda_df = DataFrame_Store.treated_df

        # Computes the ROC Curve (whether you have 2 or more classes) and returns the plots and corresponding filenames
        roc_fig, filename = iaf._plot_PLSDA_ROC_curve(self, plsda_df, target_list)
        if self.binsim_flag:
            filename = filename + '_BinSim'
        self.ROC_figure = [roc_fig,]

        # Saving filename
        self.current_other_plsda_params.update({'ROC_filename':filename})

        # Update the layouts
        self.current_pages_associated[1][12][1] = pn.pane.Plotly(self.ROC_figure[0],
                                                config={'toImageButtonOptions': {'filename': filename, 'scale':4}})


    def soft_reset(self):
        "Reset parameters."
        for param in self.param:
            if param not in ["name", 'n_dimensions', 'LVx', 'LVy', 'LVz', 'ellipse_draw', 'confidence', 'confidence_std', 'dot_size',
                             'current_pages_associated']:
                setattr(self, param, self.param[param].default)
        self.current_other_plsda_params.clear()
        self.current_other_plsda_params = {}


    def complete_soft_reset(self):
        "Resets figure parameters."
        for param in ['n_dimensions', 'LVx', 'LVy', 'LVz', 'ellipse_draw', 'confidence', 'confidence_std', 'dot_size']:
            setattr(self, param, self.param[param].default)


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
                start=3, end=40, value=5, step=1, styles={'font-weight': 'bold'}),
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
                start=3, end=40, value=5, step=1, disabled=True,
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
                start=3, end=40, value=5, step=1, disabled=True,
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

        widgets_PLS_ROC = {'n_components': pn.widgets.IntInput(name='Number of Components for PLS-DA model',
                start=1, end=40, value=5, step=1, disabled=True,
                description='This parameter was chosen when fitting the PLS-DA model and assessing its performance'),
            'n_fold': pn.widgets.IntInput(name="Number of folds for stratified cross-validation",
                value=5, start=2, end=20, disabled=True,
                description='This parameter was chosen when fitting the PLS-DA model and assessing its performance'),
            'scale': pn.widgets.Checkbox(name='Check if you did not perform scaling in the pre-treatment section. Performs auto scale.',
                value=False, disabled=True),
            'positive_class': pn.widgets.Select(name='Choose the positive class:', value='', options=['']),
            'roc_n_iter': pn.widgets.IntSlider(name="Nº of Iterations to perform", start=1, value=10, end=50,
                step=1),
            'confirm_button_roc': pn.widgets.Button(name="Compute ROC Curve", button_type='primary'),
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

        self.controls_roc = pn.Param(self, parameters=['n_components', 'n_fold', 'scale', 'positive_class',
                                                        'roc_n_iter', 'confirm_button_roc'],
                                 widgets=widgets_PLS_ROC, name='Parameters for PLS-DA ROC Curve')

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
    "Save PLS-DA permutation figure."
    filename_string = f'PLS-DA_permutation_test_{PLSDA_store.current_plsda_params_permutation["n_permutations"]}perm_'
    filename_string = filename_string + f'{PLSDA_store.current_plsda_params_permutation["n_components"]}comp_'
    filename_string = filename_string + f'{PLSDA_store.current_plsda_params_permutation["n_folds"]}-foldstratCV_scale'
    filename_string = filename_string + f'{PLSDA_store.current_plsda_params_permutation["scale"]}_metric'
    filename_string = filename_string + f'{PLSDA_store.current_plsda_params_permutation["perm_metric"]}'
    if PLSDA_store.binsim_flag:
        filename_string = filename_string + '_BinSim'
    PLSDA_store.perm_figure[0].savefig(filename_string+'.png', dpi=PLSDA_store.dpi)
    pn.state.notifications.success(f'Figure {filename_string} successfully saved.')

# Click button to save the aforementioned figure
PLSDA_store.controls_permutation.widgets['save_figure_button_permutation'].on_click(_save_figure_button_permutation_PLSDA)

# Click button to compute PLS-DA model ROC curves and obtain the corresponding plots
PLSDA_store.controls_roc.widgets['confirm_button_roc'].on_click(PLSDA_store._confirm_button_roc)

# Widget to add recommended number of components to page
rec_comp_indicator_widget = pn.indicators.Number(name='Recommended Components (based on max. Q2)', font_size='14pt', title_size='14pt',
                                    value=PLSDA_store.rec_components)

# Organizing the optimization section of the page
pls_optim_section = pn.Row(pn.Column(PLSDA_store.controls_optim, rec_comp_indicator_widget,
                                    'A lower number of components with a similar Q2 may be preferable than the number shown'),
                           PLSDA_store.optim_figure[0])


# Results Section of the PLS-DA page
# Specific Widget for PLS results section of the page, shows DataFrame with only annotated metabolites or all metabolites
plsda_feat_imp_show_annots_only = pn.widgets.Checkbox(name='Only show annotated metabolites in feature importance table',
                                                       value=False, disabled=True)

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
        plsda_params = PLSDA_store.current_plsda_params
        # Building the datafile name
        filename_string = f'PLS-DA_FeatImp_{plsda_params["feat_imp"]}_model_params_components{plsda_params["n_components"]}'
        filename_string = filename_string + f'_{plsda_params["n_folds"]}-foldstratCV_iterations{plsda_params["n_iterations"]}'
        filename_string = filename_string + f'_scale{plsda_params["scale"]}'
        if PLSDA_store.binsim_flag:
            filename_string = filename_string + '_BinSim'
        filename_string = filename_string + f'.xlsx'

        # Saving the file
        PLSDA_store.feat_impor.to_excel(filename_string)
        pn.state.notifications.success(f'{filename_string} successfully saved.')
    except:
        pn.state.notifications.error(f'File could not be saved.')

save_plsda_feat_imp_button.on_click(_save_plsda_feat_imp_button)



# Permutation Test associated widgets
# This description and equation will also be used for the Random Forest section
permutation_test_description = desc_str.permutation_test_description
pvalue_equation_string = r'p-value = \(\frac{1 + \text{times permutated model has better performance than non-permutated model}}{\text{number of permutations}}\)'


# ROC Curves Description
# This description will also be used for the Random Forest section
ROC_curve_description = desc_str.ROC_curve_description


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

# Page sections / Widgets that will change based on PLSDA_store
PLSDA_store.current_pages_associated.append(pls_optim_section)
PLSDA_store.current_pages_associated.append(pls_results_section)
PLSDA_store.current_pages_associated.append(pls_proj_page)
PLSDA_store.current_pages_associated.append(plsda_feat_imp_show_annots_only)
PLSDA_store.current_pages_associated.append(save_plsda_feat_imp_button)

# Layout of the PLS-DA page
page_PLSDA = pn.Column(pn.pane.HTML(plsda_opening_string), pls_optim_section, pls_results_section)



# Section for the Random Forest analysis
rf_opening_string = desc_str.rf_opening_string


# Param Class to store parameters and data regarding Random Forests
class RF_Storage(param.Parameterized):
    "Class to store Random Forest models, parameters and results."

    # PLS-DA flag and dataset to treat
    binsim_flag = param.Boolean(default=False)
    current_pages_associated = param.List()

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
    current_other_rf_params = param.Dict({})

    # Permutation Test
    n_perm = param.Number(default=500)
    perm_metric = param.String(default='Accuracy')
    dpi = param.Number(default=200)
    confirm_button_permutation = param.Boolean(default=False)
    save_figure_button_permutation = param.Boolean(default=False)
    current_rf_params_permutation = param.Dict()

    # ROC Curve
    positive_class = param.String(default='')
    roc_n_iter = param.Number(default=10)
    confirm_button_roc = param.Boolean(default=False)

    # Storing figures
    optim_figure = param.List(default=['To Plot the Optimization RF Figure'])
    ROC_figure = param.List(default=['To Plot the ROC Curve(s)'])
    perm_figure = param.List(default=['To Plot the Permutation Test Figure'])

    # Update number of folds limits
    def n_folds_limits_and_class_update(self, target_list):
        "Updates the limits of the CV fold number based on the number of samples per class and classes to choose for ROC curve."
        # Updating the end value limits
        min_samples_in_class = pd.Series(target_list.target).value_counts().min()
        self.controls_optim.widgets['n_fold'].end = min_samples_in_class
        self.controls.widgets['n_fold'].end = min_samples_in_class

        # Updating the proper value if it is above the highest possible number
        if min_samples_in_class < self.controls_optim.widgets['n_fold'].value:
            self.controls_optim.widgets['n_fold'].value = min_samples_in_class
            self.controls.widgets['n_fold'].value = min_samples_in_class
            self.n_fold = min_samples_in_class

        # Updating the classes available to choose as the positive for ROC curve
        if len(target_list.classes) == 2:
            self.controls_roc.widgets['positive_class'].options = target_list.classes
            self.controls_roc.widgets['positive_class'].value = target_list.classes[0]
            self.controls_roc.widgets['positive_class'].disabled = False
            self.positive_class = target_list.classes[0]
        else:
            self.controls_roc.widgets['positive_class'].options = ['Not Available']
            self.controls_roc.widgets['positive_class'].value = 'Not Available'
            self.controls_roc.widgets['positive_class'].disabled = True
            self.positive_class = 'Not Available'


    # Function to confirm the optimization
    def _confirm_optim_button(self, event):
        "Performs optimization, plots the optimization figure and updated corresponding layout."

        # Loading Widget while the Random Forst optimization is being performed
        self.current_pages_associated[0][1] = pn.indicators.LoadingSpinner(value=True, size=90,
                                                            name='Performing Optimization of RF nº of trees...')

        # Round to integer the values in the range
        self.n_min_max_trees = int(np.round(self.n_min_max_trees[0])), int(np.round(self.n_min_max_trees[1]))

        # See what dataset to use
        if self.binsim_flag:
            rf_df = DataFrame_Store.binsim_df
        else:
            rf_df = DataFrame_Store.treated_df

        # Perform Optimization
        rf_optim = iaf._optimization_n_trees_rf(self, rf_df, target_list.target)
        self.optim_scores = list(rf_optim['mean_test_score'])
        self.optim_ntrees = list(rf_optim['param_n_estimators'])

        # Saving parameters for optimization
        self.current_other_rf_params.update({'n_fold_optim': self.n_fold, 'n_min_max_trees': self.n_min_max_trees,
                                            'n_interval': self.n_interval})

        # RF optimization Figure
        rf_optim_results = pd.DataFrame([self.optim_ntrees, self.optim_scores],
                                       index=['Number of Trees', 'Model Accuracy (estimated by CV)']).T
        self.optim_figure = [px.line(rf_optim_results, x='Number of Trees', y='Model Accuracy (estimated by CV)',
                title='Random Forest Optimization Plot', range_y=(0, 1.05),
                range_x=(self.n_min_max_trees[0] - 10, self.n_min_max_trees[1] + 10)
                 ),]

        # Update the layout
        filename_string = f'RF_optim_plot_{self.n_fold}-foldStratCV_{self.n_min_max_trees[0]}to{self.n_min_max_trees[1]}trees'
        filename_string = filename_string + f'({self.n_interval}interval)'
        if self.binsim_flag:
            filename_string = filename_string + '_BinSim'
        self.current_pages_associated[0][1] = pn.pane.Plotly(self.optim_figure[0],
                                              config = {'toImageButtonOptions': {'filename': filename_string, 'scale':4,}})


    # Function to fit the RF model and obtain results
    def _confirm_rf_button(self, event):
        "Fits a Random Forest model, gives model performance estimations and feature importance lists, updating the layout."

        # Loading Widget while the Random Forest model is being fitted and assessed
        self.current_pages_associated[1][0][1][1] = pn.indicators.LoadingSpinner(value=True, size=90,
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

        # See what dataset to use
        if self.binsim_flag:
            rf_df = DataFrame_Store.binsim_df
        else:
            rf_df = DataFrame_Store.treated_df

        # Fits a RF model under a stratified cross validation scheme
        RF_results = metsta.RF_model(rf_df, target_list.target, # Data and labels
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
        self.current_pages_associated[1][0][1][1] = pn.pane.DataFrame(self.n_results)
        self.current_pages_associated[1][3] = pn.pane.DataFrame(self.feat_impor, height=600)
        self.current_pages_associated[2].value = False
        self.current_pages_associated[2].disabled = False
        self.current_pages_associated[3].disabled = False

        if len(self.current_pages_associated[1]) == 5:
            self.current_pages_associated[1].append('### Random Forest Permutation Test Section')
            self.current_pages_associated[1].append(pn.Column(pn.pane.HTML(permutation_test_description),
                                                 pn.pane.LaTeX(pvalue_equation_string, styles={'font-size': '14pt'})))
            self.current_pages_associated[1].append(pn.Row(self.controls_permutation,
                                          self.perm_figure[0]))
            self.current_pages_associated[1].append('### Random Forest Receiver Operating Characteristic (ROC) Curve Section')
            self.current_pages_associated[1].append(pn.pane.HTML(ROC_curve_description))
            self.current_pages_associated[1].append(pn.Row(self.controls_roc, self.ROC_figure[0]))


    # Function to confirm the permutation test to perform
    def _confirm_button_permutation(self, event):
        "Performs a permutation test for the Random Forest model and updates the layout."

        # Loading Widget while the Permutation Test is being performed
        self.current_pages_associated[1][7][1] = pn.indicators.LoadingSpinner(value=True, size=90,
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

        # See what dataset to use
        if self.binsim_flag:
            rf_df = DataFrame_Store.binsim_df
        else:
            rf_df = DataFrame_Store.treated_df

        # Performs the Permutation Test on the Random Forest model under a stratified cross validation scheme
        perm_results_RF = metsta.permutation_RF(
            rf_df, target_list.target, # Data and target
            iter_num=self.n_perm, # Nº of permutations to do in your test - around 500 should be enough
            n_trees=self.n_trees, # Number of trees in the model
            cv=None, n_fold=self.n_fold, # Choose the number of folds
            random_state=None, # Random seed given to make the permutations rng class labels
            metric=perm_metric) # Choose a metric to use to evaluate if the model is significant

        # Store Parameters used for the model
        self.current_rf_params_permutation = {'n_trees': self.n_trees, 'n_folds': self.n_fold,
                                     'n_permutations': self.n_perm, 'perm_metric': self.perm_metric}

        # Plot the permutation test
        perm_figure = iaf._plot_permutation_test(perm_results_RF, rf_df, self.n_fold,
                                                     self.perm_metric, 'Random Forest Permutation Test')
        self.perm_figure = [perm_figure,]

        # Update the layout
        self.current_pages_associated[1][7][1] = pn.pane.Matplotlib(self.perm_figure[0], height=600)
        self.controls_permutation.widgets['save_figure_button_permutation'].disabled = False


    # Function to confirm the ROC Curves to plot
    def _confirm_button_roc(self, event):
        "Plots a ROC Curve to assess Random Forest model performance and updates the layout."

        # Loading Widget while the ROC Curves are being computed
        self.current_pages_associated[1][10][1] = pn.indicators.LoadingSpinner(value=True, size=90,
                                                            name='Computing Receiver Operating Characteristic Curves...')

        # See what dataset to use
        if self.binsim_flag:
            rf_df = DataFrame_Store.binsim_df
        else:
            rf_df = DataFrame_Store.treated_df

        # Computes the ROC Curve (whether you have 2 or more classes) and returns the plots and corresponding filenames
        roc_fig, filename = iaf._plot_RF_ROC_curve(self, rf_df, target_list)
        if self.binsim_flag:
            filename = filename + '_BinSim'
        self.ROC_figure = [roc_fig,]

        # Saving filename
        self.current_other_rf_params.update({'ROC_filename':filename})

        # Update the layouts
        self.current_pages_associated[1][10][1] = pn.pane.Plotly(self.ROC_figure[0],
                                                config={'toImageButtonOptions': {'filename': filename, 'scale':4}})


    def reset(self):
        "Reset parameters."
        for param in self.param:
            if param not in ["name", "current_pages_associated"]:
                setattr(self, param, self.param[param].default)
        self.current_other_rf_params.clear()
        self.current_other_rf_params = {}


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

        widgets_RF_ROC = {'n_trees': pn.widgets.IntInput(name='Number of Trees for Random Forest model',
                start=1, end=600, value=200, step=1, disabled=True,
                description='This parameter was chosen when fitting the Random Forest model and assessing its performance'),
            'n_fold': pn.widgets.IntInput(name="Number of folds for stratified cross-validation",
                value=5, start=2, end=20, disabled=True,
                description='This parameter was chosen when fitting the Random Forest model and assessing its performance'),
            'positive_class': pn.widgets.Select(name='Choose the positive class:', value='', options=['']),
            'roc_n_iter': pn.widgets.IntSlider(name="Nº of Iterations to perform", start=1, value=10, end=50,
                step=1),
            'confirm_button_roc': pn.widgets.Button(name="Compute ROC Curve", button_type='primary'),
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

        self.controls_roc = pn.Param(self, parameters=['n_trees', 'n_fold', 'positive_class',
                                                        'roc_n_iter', 'confirm_button_roc'],
                                 widgets=widgets_RF_ROC, name='Parameters for Random Forest ROC Curve')

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
    "Saves Random Forest permutation figure."
    filename_string = f'RF_permutation_test_{RF_store.current_rf_params_permutation["n_permutations"]}perm_'
    filename_string = filename_string + f'{RF_store.current_rf_params_permutation["n_trees"]}trees_'
    filename_string = filename_string + f'{RF_store.current_rf_params_permutation["n_folds"]}-foldstratCV_metric'
    filename_string = filename_string + f'{RF_store.current_rf_params_permutation["perm_metric"]}'
    if RF_store.binsim_flag:
        filename_string = filename_string + '_BinSim'
    RF_store.perm_figure[0].savefig(filename_string+'.png', dpi=RF_store.dpi)
    pn.state.notifications.success(f'Figure {filename_string} successfully saved.')

# Click button to save the aforementioned figure
RF_store.controls_permutation.widgets['save_figure_button_permutation'].on_click(_save_figure_button_permutation_RF)

# Click button to compute Random Forest model ROC curves and obtain the corresponding plots
RF_store.controls_roc.widgets['confirm_button_roc'].on_click(RF_store._confirm_button_roc)


# Optimization section of the page
# Organizing the optimization section of the page
rf_optim_section = pn.Row(RF_store.controls_optim, RF_store.optim_figure[0])


# Results Section of the Random Forest page
# Specific Widget for Random Forest results section, shows DataFrame with only annotated metabolites or all metabolites
rf_feat_imp_show_annots_only = pn.widgets.Checkbox(name='Only show annotated metabolites in feature importance table',
                                                       value=False, disabled=True)

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
        rf_params = RF_store.current_rf_params
        # Building the datafile name
        filename_string = f'RF_FeatImp_Gini_model_params_{rf_params["n_trees"]}trees_{rf_params["n_folds"]}-foldstratCV_'
        filename_string = filename_string + f'iterations{rf_params["n_iterations"]}'
        if RF_store.binsim_flag:
            filename_string = filename_string + '_BinSim'
        filename_string = filename_string + '.xlsx'

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

# Page sections / Widgets that will change based on RF_store
RF_store.current_pages_associated.append(rf_optim_section)
RF_store.current_pages_associated.append(rf_results_section)
RF_store.current_pages_associated.append(rf_feat_imp_show_annots_only)
RF_store.current_pages_associated.append(save_rf_feat_imp_button)


# Layout of the RF page
page_RF = pn.Column(pn.pane.HTML(rf_opening_string), rf_optim_section, rf_results_section)


# Layout of the full supervised analysis page
sup_analysis_page = pn.Tabs(('PLS-DA', page_PLSDA), ('RF ', page_RF))




# Page for Univariate Analysis
univ_opening_string = desc_str.univ_opening_string

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
    scaling_kw = param.String(default="Average")

    # Unsupervised Analysis Main Parameters
    control_class = param.Selector(default='')
    test_class = param.Selector(default='')
    univariate_test_str = param.String('What Univariate Test to perform')
    univariate_test = param.Selector(default='T-Test (Parametric)')
    expected_equal_variance = param.Boolean(default=True)
    p_value_threshold = param.Number(default=0.05)
    fold_change_threshold = param.Number(default=2)
    univariate_df = param.DataFrame()
    univariate_df_set = param.Dict()

    # Store Univariate Test Parameters
    current_univ_params = param.Dict({})

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

        # Saving parameters
        self.current_univ_params = {'Control Class': self.control_class, 'Test Class': self.test_class,
                                    'Test': self.univariate_test, 'Expected Equal Var.': self.expected_equal_variance,
                                    'p-value': self.p_value_threshold, 'Fold Change Threshold': self.fold_change_threshold}

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
        results_df = self.univariate_results_non_filt.copy()
        results_df['-log10(Adj. p-value)'] = -np.log10(results_df['FDR adjusted p-value'])

        expression = []
        for i in results_df.index:
            if i not in self.univariate_results.index:
                expression.append('Non-Significant')
            elif results_df.loc[i, 'log2FC'] < 0:
                expression.append('Downregulated')
            else:
                expression.append('Upregulated')
        results_df['Expression'] = expression

        Volcano_fig = iaf._plot_Volcano_plot(results_df, self,)
        self.Volcano_fig = [Volcano_fig,]
        # Update the layout
        univ_parameters = self.current_univ_params
        filename_string = f'VolcanoPlot - {univ_parameters["Test Class"]}_vs_{univ_parameters["Control Class"]}'
        filename_string = filename_string + f'_{univ_parameters["Test"]}_pvalue{univ_parameters["p-value"]}_FC'
        filename_string = filename_string + f'{univ_parameters["Fold Change Threshold"]}'
        layout_volcano[0,1:3] = pn.pane.Plotly(self.Volcano_fig[0], config={'toImageButtonOptions': {
                   'filename': filename_string, 'scale':4}})


    @param.depends('test_class_subset', 'show_annots_only', watch=True)
    def _update_intersections(self):
        "Update the intersections DataFrame and description."
        iaf._univariate_intersections(self, DataFrame_Store)
        layout_final_section[0][1] = self.inter_description
        layout_final_section[1] = pn.pane.DataFrame(self.specific_cl_df, height=600)


    def reset(self):
        "Reset parameters."
        for param in self.param:
            if param not in ["name", 'color_non_sig', 'color_down_sig', 'color_up_sig', 'filt_method', 'filt_kw', 'mvi_method',
                             'mvi_kw', 'norm_method', 'norm_kw', 'tf_method', 'tf_kw', 'scaling_method', 'scaling_kw']:
                setattr(self, param, self.param[param].default)


    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets = {
            'control_class': pn.widgets.Select(name='Control Class', options=target_list.classes),
            'test_class': pn.widgets.Select(name='Test Class', options=target_list.classes),
            'univariate_test_str': pn.widgets.StaticText(value='What Univariate Test to perform',
                                                         styles={'font-size': 'medium', 'font-weight': 'bold'}),
            'univariate_test': pn.widgets.RadioBoxGroup(name='What Univariate Test to perform',
                    value = 'T-Test (Parametric)',
                    options=['T-Test (Parametric)', 'Mann-Whitney Test (Non-Parametric)']),
            'expected_equal_variance': pn.widgets.Checkbox(name='Consider Variance between classes as equal (only affects T-Test)', value=True),
            'p_value_threshold': pn.widgets.EditableFloatSlider(name='P-value threshold', start=0, end=1, value=0.05,
                    step=0.001),
            'fold_change_threshold': pn.widgets.FloatInput(name='Fold Change Threshold', value=2, step=0.1, start=1,
                    description='''The threshold is set so as only selecting features as significant if either the control class or the test class average is "chosen" threhold times higher than the opposing class.
                    1 will essentially skip this threshold and use only the p-value threshold. Setting 1 for p-value threshold has the same effect for the p-value step.
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
    univ_parameters = UnivarA_Store.current_univ_params
    test_performed = univ_parameters["Test"].split(' ')[0]
    filename_string = f'Univar_res_{univ_parameters["Test Class"]}_vs_{univ_parameters["Control Class"]}_{test_performed}'
    filename_string = filename_string + f'_pvalue{univ_parameters["p-value"]}_FC{univ_parameters["Fold Change Threshold"]}'
    if test_performed == 'T-Test':
        if univ_parameters["Expected Equal Var."]:
            filename_string = filename_string + f'_equalvariance.xlsx'
        else:
            filename_string = filename_string + f'_notequalvariance.xlsx'
    else:
        filename_string = filename_string + f'.csv'

    # Saving the file
    UnivarA_Store.univariate_results.to_excel(filename_string)
    pn.state.notifications.success(f'{filename_string} successfully saved.')

save_univariate_results_button.on_click(_save_univariate_results_button)

# Widget button to save dataframe shown when seeing intersections of multiple univariate analysis of
# different chosen test classes against a chosen control class
save_multiple_univariate_button = pn.widgets.Button(name='Save shown Dataframe with int. values as .csv file (in current folder)',
                                                button_type='warning', icon=iaf.download_icon)

# When pressing the button, downloads the dataframe (filename quite big)
def _save_multiple_univariate_button(event):
    "Saves multiple univariate (features significant in the control class versus the selected test classes) DataFrame created."
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
    univ_parameters = UnivarA_Store.current_univ_params
    test_performed = univ_parameters["Test"].split(' ')[0]
    filename_string = annot_string + '_significant_in_testing_each_chosen_test_class' + class_subset_string + 'against'
    filename_string = filename_string + f'_control_{univ_parameters["Control Class"]}_{test_performed}'
    filename_string = filename_string + f'_pvalue{univ_parameters["p-value"]}_FC{univ_parameters["Fold Change Threshold"]}'
    if test_performed == 'T-Test':
        if univ_parameters["Expected Equal Var."]:
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
    univ_parameters = UnivarA_Store.current_univ_params
    filename_string = f'VolcanoPlot - {univ_parameters["Test Class"]}_vs_{univ_parameters["Control Class"]}'
    filename_string = filename_string + f'_{univ_parameters["Test"]}_pvalue{univ_parameters["p-value"]}_FC'
    filename_string = filename_string + f'{univ_parameters["Fold Change Threshold"]}'
    layout_volcano[0,1:3] = pn.pane.Plotly(UnivarA_Store.Volcano_fig[0], config={'toImageButtonOptions': {
                'filename': filename_string, 'scale':4}})

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
# TODO: (PROBLEM) Legend does not appear in Van Krevelen Plot - make it appear
# TODO: Make Kendrick Mass Defect Plots have different sizes based on avg. intensity like in VK plots?
# TODO: KMD Plots of different classes have different dot sizes - why?? What can even cause this? - No idea

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
    ccs_df = param.DataFrame()

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

        if len(self.vk_formula_to_consider) == 0:
            pn.state.notifications.error('At least 1 Formula Annotation column must be provided for VK plot.')
            vk_plots.clear()
            self.VanKrevelen_plot = ['Pane for Van Krevelen Plot']
            self.vanKrevelen_filenames = []
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
            kmd_plots.clear()
            self.KendrickMD_plot = ['Pane for Chemical Composition Series']
            self.KendrickMD_filenames = []
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
    def _compute_CCS_plot(self,
                          series_order=('CHO', 'CHOS', 'CHON', 'CHNS', 'CHONS', 'CHOP', 'CHONP','CHONSP', 'other')):
        "Computes and plots the chemical compostion series plot, updating the corresponding layout."

        if len(self.ccs_formula_to_consider) == 0:
            pn.state.notifications.error('At least 1 Formula Annotation column must be provided for CCS plot.')
            # Clear the page
            while len(ccs_page) > 2:
                ccs_page.pop(-1)
            self.CCS_plot = ['Pane for Chemical Composition Series']
            raise ValueError('At least 1 Formula Annotation column must be provided for CCS plot.')

        # Initialize figure
        fig = go.Figure()

        # Initialize the DataFrame
        series_df = pd.DataFrame()

        # Description of the number of formulas assigned considered for each class and from how many peaks they came from.
        desc_string = [f'**Description with number of formulas assigned considered**', '']

        # For each class
        for g in com_exc_compounds.group_dfs:
            # Calculating H/C and O/C ratios and series classes for data diversity plots.
            forms = com_exc_compounds.group_dfs[g].dropna(subset=self.ccs_formula_to_consider, how='all')
            elems = iaf.create_element_counts(forms, formula_subset=self.ccs_formula_to_consider)

            # Calculating counts of each series for each class
            counts = elems['Series'].value_counts().reindex(series_order)
            series_df[g] = counts # Store the counts in the dataframe

            # Adding the class' bars to the plot
            if self.ccs_bar_plot_type == 'Horizontal': # Horizontal bar plot
                fig.add_trace(go.Bar(name=g, x=counts, y=series_order, orientation='h', marker_color=target_list.color_classes[g]))
            else: # Vertical bar plot
                fig.add_trace(go.Bar(name=g, x=series_order, y=counts, orientation='v', marker_color=target_list.color_classes[g]))

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
        self.ccs_df = series_df.T
        filename = 'CCS_Plot_formulacolumns'
        # Create appropriate filename
        for cl in self.ccs_formula_to_consider:
            filename = filename + f'_{cl}'

        # Updating chemical composition series layout
        if len(ccs_page) == 2:
            ccs_page.append(self.ccs_desc)
            ccs_page.append(pn.pane.Plotly(self.CCS_plot[0], height=600,
                            config = {'toImageButtonOptions': {'filename': filename, 'scale':4,}}))
            ccs_page.append(pn.pane.DataFrame(series_df.T, height=400))
            ccs_page.append(save_ccs_table_button)
        else:
            ccs_page[2] = self.ccs_desc
            ccs_page[3] = pn.pane.Plotly(self.CCS_plot[0], height=600,
                            config = {'toImageButtonOptions': {'filename': filename, 'scale':4,}})
            ccs_page[4] = pn.pane.DataFrame(series_df.T, height=400)


    def update_widgets(self):
        "Update the needed widget values."
        # Calculate specific class DataFrames in case it has not been calculated before
        try:
            for g in com_exc_compounds.group_dfs:
                g
        except:
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
        self._compute_CCS_plot()


    def reset(self):
        "Reset parameters."
        self.vk_highlight_by = 'None'
        for param in self.param:
            if param not in ["name", 'vk_formula_to_consider', 'kmd_formula_to_consider', 'ccs_formula_to_consider', 'vk_highlight_by']:
                setattr(self, param, self.param[param].default)
        self.vk_highlight_by = 'Rank'
        for param in ['vk_formula_to_consider', 'kmd_formula_to_consider', 'ccs_formula_to_consider']:
            setattr(self, param, self.param[param].default)


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
            'vk_draw_class_rectangle': pn.widgets.Checkbox(name='Draw Peptide, Lignins, Tannins, Nucleotides and Phytochemical area rectangles',
                                                            value=False, disabled=True),
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
    "Computes and updates layout of page."
    if len(vk_page) == 1:
        vk_page.append(dataviz_store.controls_vk)
        vk_page.append(vk_plots)
    if len(kmd_page) == 1:
        kmd_page.append(dataviz_store.controls_kmd)
        kmd_page.append(kmd_plots)
    if len(ccs_page) == 1:
        ccs_page.append(dataviz_store.controls_ccs)
    dataviz_store.update_widgets()
compute_vk_kmd_ccs_button.on_click(_compute_vk_kmd_ccs_button)


# Van Krevelen Plot Section
vk_opening_string = desc_str.vk_opening_string

vk_plots = pn.Column()

vk_page = pn.Column(pn.pane.HTML(vk_opening_string))


# Kendrick Mass Defect Plot Section
kmd_opening_string = desc_str.kmd_opening_string

kmd_plots = pn.Column()

kmd_page = pn.Column(pn.pane.HTML(kmd_opening_string))


# Chemical Composition Series Section
ccs_opening_string = desc_str.ccs_opening_string

# Widget to save dataframe of chemical compositions series
save_ccs_table_button = pn.widgets.Button(name='Save Chemical Composition Table obtained as .xlsx (in current folder)',
                                                button_type='warning', icon=iaf.download_icon, disabled=False)

# When pressing the button, downloads the dataframe (builds the appropriate filename)
def _save_ccs_table_button(event):
    "Save Chemical Composition Series Table as an Excel."
    try:
        # Building the datafile name
        filename = 'CCS_Table_formulacolumns'
        # Create appropriate filename
        for cl in dataviz_store.ccs_formula_to_consider:
            filename = filename + f'_{cl}'

        # Saving the file
        dataviz_store.ccs_df.to_excel(filename + '.xlsx')
        pn.state.notifications.success(f'{filename}.xlsx successfully saved.')
    except:
        pn.state.notifications.error(f'File could not be saved.')

save_ccs_table_button.on_click(_save_ccs_table_button)

ccs_page = pn.Column(pn.pane.HTML(ccs_opening_string))


data_viz_page = pn.Column(compute_vk_kmd_ccs_button, pn.Tabs(('Van Krevelen Plot', vk_page), ('Kendrick Mass Defect Plot', kmd_page),
                        ('Chemical Composition Series', ccs_page)))




# Page for Pathway Assignment (based on HMDB)
# Param Class to store parameters and data regarding the Path Assignment Page
class PathAssignment_Storage(param.Parameterized):
    "Class to contain parameters and relevant DataFrames to the pathway matching to HMDB IDs."

    # HMDB ID columns
    hmdb_id_cols = param.List(default=[])
    current_hmdb_id_cols = param.List(default=[])

    # DataFrame to store pathway information
    pathway_db = param.DataFrame()

    # DataFrame to store pathway assignments
    pathway_assignments = param.DataFrame(pd.DataFrame(columns=pathway_db.columns))
    assign_desc = param.String()

    # In case specific HMDB are selected to see their attributions
    chosen_hmdb_ids = param.String()
    chosen_hmdb_ids_assigns = param.Dict({})


    def _pathway_assignments(self, metadata_df, pathway_db):
        "Performs pathway matching between HMDB IDs found and the pathway database."

        # Obtain the list of every HMDB ID annotated in the columns given
        HMDB_list = [] # To store
        error_cols = []

        # Go through every hmdb_id_cols provided
        for id_cols in self.hmdb_id_cols:
            # Get the column and dropna null values
            filt_df = metadata_df[id_cols].dropna()

            # Every idx
            for idx in filt_df.index:
                value = filt_df.loc[idx]
                # If a single string HMDB idx
                if type(value) == str:
                    try:
                        if value not in HMDB_list: # Add if not already present
                            HMDB_list.append(value)
                    except: # If not a string
                        if id_cols not in error_cols:
                            error_cols.append(id_cols)

                # If a list of HMDB idxs
                else:
                    try:
                        for hmdb_idx in range(len(filt_df.loc[idx])):
                            curr_id = filt_df.loc[idx][hmdb_idx]
                            if curr_id not in HMDB_list: # Add if not already present
                                HMDB_list.append(curr_id)
                    except: # If not a list-like or has non string values in the list
                        if id_cols not in error_cols:
                            error_cols.append(id_cols)

        for col in error_cols:
            pn.state.notifications.info(f'{col} provided contained values that were not strings, list-like iterable of strings or nulls as values. Or it had non-string elements in the list-like iterable.')

        # Build the pathways assignment DataFrame
        pathways_assignment = pd.DataFrame(index=HMDB_list, columns=pathway_db.columns)
        # Match the HMDB IDs to their information in the pathway db (if it exists)
        for idx in pathways_assignment.index:
            if idx in pathway_db.index:
                pathways_assignment.loc[idx] = pathway_db.loc[idx]
        pathways_assignment.index.name = 'HMDB IDs'

        # Assign to attribute and store search made
        self.pathway_assignments = pathways_assignment
        self.current_hmdb_id_cols = self.hmdb_id_cols


    def _path_assign_describer(self):
        "Generates Description of assignment made."

        hmdb_ids = [] # Get actual HMDB IDs found
        complete_len = len(self.pathway_assignments.index)
        for idx in self.pathway_assignments.index:
            # Has to fulfill these conditions
            if len(idx) == 11:
                if idx.startswith('HMDB'):
                    hmdb_ids.append(idx)

        path_filt = self.pathway_assignments.dropna()
        self.assign_desc = f'From the **{complete_len}** IDs provided, **{len(hmdb_ids)}** were HMDB-like IDs. From those, **{len(path_filt)}** were matched to the pathway database.'


    def _hmdb_id_search(self, pathway_db):
        "Searches provdied list of HMDB IDs in the pathway db and in the dataset and returns results based on what was found."

        # Obtain the IDs inputted in the text box
        chosen_ids = self.chosen_hmdb_ids.split('\n')

        # Clear the dictionary to start with
        self.chosen_hmdb_ids_assigns.clear()

        for hmdb_id in chosen_ids:
            # See if the hmdb_id is in the database and see if it is in your data
            hmdb_id_in_db = False
            hmdb_id_in_dataset = False
            if hmdb_id in pathway_db.index:
                hmdb_id_in_db = True
            if hmdb_id in PathAssign_store.pathway_assignments.index:
                hmdb_id_in_dataset = True

            # If no pathway matching is found in our pathway database for the hmdb id we are looking for
            if not hmdb_id_in_db:
                # But it is found annotated in the dataset
                if hmdb_id_in_dataset:
                    desc = f'**{hmdb_id}** was annotated in the current dataset but no pathway matching was found in our database.'
                    self.chosen_hmdb_ids_assigns[hmdb_id] = desc

                # And also it is not found annotated in the dataset
                else:
                    desc = f'**{hmdb_id}** neither was annotated in the current dataset or was associated with pathways in our database.'
                    self.chosen_hmdb_ids_assigns[hmdb_id] = desc

            # If there is pathway matching found in our pathway database
            if hmdb_id_in_db:
                # And it is also found annotated in the dataset
                if hmdb_id_in_dataset:
                    df = pathway_db.loc[[hmdb_id]].explode(column=['Pathway Name', 'Pathway ID'])
                    PathAssign_store.chosen_hmdb_ids_assigns[hmdb_id] = df

                # But it is not found annotated in the dataset
                else:
                    df = pathway_db.loc[[hmdb_id]].explode(column=['Pathway Name', 'Pathway ID'])
                    desc = f'**{hmdb_id}** ({df["HMDB Name"].iloc[0]}) was associated with pathways in our database (see below) but was **NOT found in the current dataset**.'
                    self.chosen_hmdb_ids_assigns[hmdb_id] = [desc, df]


    def _update_widgets(self, metadata_df, pathway_db):
        "Updates Widgets based on dataset information."

        self.controls.widgets['hmdb_id_cols'].options = list(metadata_df.columns)

        self.pathway_db = pathway_db
        self.controls_hmdb.widgets['chosen_hmdb_ids'].options = list(pathway_db.index)


    def reset(self):
        "Resets all relevant parameters."
        for param in self.param:
            if param not in ["name"]:
                setattr(self, param, self.param[param].default)


    def __init__(self, **params):

        super().__init__(**params)
        widgets_initial = {
            'hmdb_id_cols': pn.widgets.CrossSelector(name='Choose data columns with HMDB IDs from metadata columns:',
                value=[], options=[])
        }

        widgets_final = {
            'chosen_hmdb_ids':pn.widgets.TextAreaInput(
                name='Enter specific HMDB IDs to search for (follow placeholder formatting):',
                rows=10, placeholder='''HMDB0000001\nHMDB0001045\nHMDB0000125\nHMDB0123678''', max_length=20000)
        }

        self.controls = pn.Param(self, parameters=['hmdb_id_cols'], widgets=widgets_initial,
                                 name='Choose data columns with HMDB IDs (in string or list format) from metadata columns:')

        self.controls_hmdb = pn.Param(self, parameters=['chosen_hmdb_ids'], widgets=widgets_final,
                                 name='Specific search of HMDBs in Database')


# Load pathway database into a DataFrame
with open('RAMP_ID_pathways_improved.pickle', 'rb') as handle:
    pathway_db = pickle.load(handle)
pathway_db = pathway_db[['HMDB Name', 'Pathway Name', 'Pathway ID']]

# Initialize store for PathwayAssignment parameters
PathAssign_store = PathAssignment_Storage(pathway_db=pathway_db)

# Pathway matching section widgets
# Widget to confirm the HMDB ID columns
confirm_hmdb_cols_button = pn.widgets.Button(name='Confirm HMDB ID columns selected and perform Pathway Assignment',
                                             button_type='success', icon=iaf.img_confirm_button)
def _confirm_hmdb_cols_button(event):
    "Perform the pathway assignment matching to the HMDB IDs."
    confirm_hmdb_ids_to_search_button.disabled = False
    # Perform the assignment
    PathAssign_store._pathway_assignments(DataFrame_Store.metadata_df, pathway_db)
    PathAssign_store._path_assign_describer() # Get the description of the assignment

    # Update page layout
    while len(path_matching_section) > 3:
        path_matching_section.pop(-1)
    path_matching_section.append(PathAssign_store.assign_desc)
    path_matching_section.append(pathassign_show_matched_only)
    if len(PathAssign_store.pathway_assignments) > 5:
        path_matching_section.append(pn.pane.DataFrame(PathAssign_store.pathway_assignments, height=800))
    else:
        path_matching_section.append(pn.pane.DataFrame(PathAssign_store.pathway_assignments))
    path_matching_section.append(save_pathway_assignments_button)
# Calling the button function
confirm_hmdb_cols_button.on_click(_confirm_hmdb_cols_button)


# Specific Widget for middle section of the page, shows DataFrame with all metabolites annotated or only those with matched pathways
pathassign_show_matched_only = pn.widgets.Checkbox(name='Only show HMDB ID metabolites with matched pathways', value=False)

# Change the DataFrame shown based on checkbox
@pn.depends(pathassign_show_matched_only.param.value, watch=True)
def _pathway_assignments_df(pathassign_show_matched_only):
    "Update the layout based on if we are showing all annotated metabolites or only ones with matching pathways."
    # Select DataFrame
    if pathassign_show_matched_only:
        df_to_show = PathAssign_store.pathway_assignments.dropna()
    else:
        df_to_show = PathAssign_store.pathway_assignments

    # Update the layout
    if len(path_matching_section) > 5:
        if len(df_to_show) > 5:
            path_matching_section[5] = pn.pane.DataFrame(df_to_show, height=800)
        else:
            path_matching_section[5] = pn.pane.DataFrame(df_to_show, height=800)

# Widget to save dataframe of univariate analysis performed in .csv format
save_pathway_assignments_button = pn.widgets.Button(name='Save pathway matching DataFrame as .xlsx (in current folder)',
                                                button_type='warning', icon=iaf.download_icon)

# When pressing the button, downloads the dataframe (builds the appropriate filename)
def _save_pathway_assignments_button(event):
    "Save pathway assignments made."
    # Building the datafile name
    filename_string = f'HMDB_IDs_PathwaysMatching_to_IDs_in'
    for col in PathAssign_store.current_hmdb_id_cols:
        filename_string = filename_string + f'_{col}'
    filename_string = filename_string + '.xlsx'

    # Saving the file
    PathAssign_store.pathway_assignments.to_excel(filename_string)
    pn.state.notifications.success(f'{filename_string} successfully saved.')
# Calling the function
save_pathway_assignments_button.on_click(_save_pathway_assignments_button)


# Pathway matching section page layout
path_matching_section = pn.Column('## Pathway Assignment and Matching Section',
                                  PathAssign_store.controls, confirm_hmdb_cols_button)


# HMDB ID searching section page widgets
# Widget to confirm the HMDB ID columns
confirm_hmdb_ids_to_search_button = pn.widgets.Button(name='Confirm HMDB IDs selected to search for associated pathways',
                                             button_type='success', icon=iaf.img_confirm_button, disabled=True)
def _confirm_hmdb_ids_to_search_button(event):
    "Search for the HMDB IDs provided on the pathway database and onthe dataset analysed and update page layout accordingly."

    # Perform the search
    PathAssign_store._hmdb_id_search(pathway_db)

    # Reset the layout for each search
    while len(hmdb_id_searching_section) > 3:
        hmdb_id_searching_section.pop(-1)

    # Update the page layout according to the HMDB searched and the results found
    for key, value in PathAssign_store.chosen_hmdb_ids_assigns.items():

        # If it is a string, this means the HMDB identifier was not found in the database
        if type(value) == str:
            compound_section = pn.Column(f'### {key}', value)
            hmdb_id_searching_section.append(compound_section)

        # If it is a list, this means the HMDB identifier was found in the database but was not in the dataset
        elif type(value) == list:
            if len(value[1]) > 20:
                compound_section = pn.Column(f'### {key} - {value[1]["HMDB Name"].iloc[0]}',
                                            pn.pane.Alert(value[0], alert_type='warning'),
                                            pn.pane.DataFrame(value[1], height=600))
            else:
                compound_section = pn.Column(f'### {key} - {value[1]["HMDB Name"].iloc[0]}',
                                            pn.pane.Alert(value[0], alert_type='warning'),
                                            pn.pane.DataFrame(value[1]))
            hmdb_id_searching_section.append(compound_section)

        # If it is a DataFrame, this means the HMDB identifier was found in the database and in the dataset
        else:
            if len(value) > 20:
                compound_section = pn.Column(f'### {key} - {value["HMDB Name"].iloc[0]}',
                                            pn.pane.DataFrame(value, height=600))
            else:
                compound_section = pn.Column(f'### {key} - {value["HMDB Name"].iloc[0]}',
                                            pn.pane.DataFrame(value))
            hmdb_id_searching_section.append(compound_section)

# Calling the button function
confirm_hmdb_ids_to_search_button.on_click(_confirm_hmdb_ids_to_search_button)

# HMDB ID searching page layout
hmdb_id_searching_section = pn.Column('## HMDB ID Pathway Searching Section',
                                      PathAssign_store.controls_hmdb,
                                      confirm_hmdb_ids_to_search_button)


# Overall Page layout
path_assign_page = pn.Column(path_matching_section, hmdb_id_searching_section)




# Page for BinSim Analysis

# Page will have 4 tabs, one for each of the main statistical analysis: PCA, HCA, PLS-DA and Random Forest

# Starting with the PCA section

# Running initial param to store PCA details - BinSim version
PCA_params_binsim = PCA_Storage()
PCA_params_binsim.binsim_flag = True # Raising BinSim flag
# Making PCA Projection Controlling widgets disabled
for n, w in PCA_params_binsim.controls.widgets.items():
    w.disabled = True

# Extra widgets for the page
compute_PCA_binsim_button = pn.widgets.Button(name='Compute', button_type='success')
n_components_compute_binsim = pn.widgets.IntInput(name='Number of Components to Compute:',
                                           value=10, step=1, start=2, end=20,
                                          description='Select 2-20 components.')


# Modified version of this function so PCA will be computed and images will be plotted when it is a new analysis
# When pressing the button, runs the PCA with the designated number of components and the binsim treated data
def _compute_PCA_binsim_button(event):
    "Computes PCA, plots the figures in the respective pages if not present."

    # Select DataFrame
    if PCA_params_binsim.binsim_flag:
        df = DataFrame_Store.binsim_df
    else:
        df = DataFrame_Store.treated_df

    # Calculate PCA and store results
    principaldf, var, loadings = metsta.compute_df_with_PCs_VE_loadings(df,
                                       n_components=n_components_compute_binsim.value,
                                       whiten=True, labels=target_list.target, return_var_ratios_and_loadings=True)
    PCA_params_binsim.pca_scores = principaldf
    PCA_params_binsim.explained_variance = var
    PCA_params_binsim.pca_loadings = pd.DataFrame(loadings)
    PCA_params_binsim.n_components = n_components_compute_binsim.value # This may catalyse the 'change in n_components' react function
    PCA_params_binsim.controls.widgets['n_components'].value = n_components_compute_binsim.value

    # Plot the figures - PCA Projection
    PCA_params_binsim.PCA_plot[0] = iaf._plot_PCA(PCA_params_binsim, target_list)
    PCA_filename_string = 'PCA_plot'
    if PCA_params_binsim.ellipse_draw:
        if PCA_params_binsim.confidence != 0:
            PCA_filename_string = PCA_filename_string + f'_ellipse({PCA_params_binsim.confidence*100}%confidence)'
        else:
            PCA_filename_string = PCA_filename_string + f'_ellipse({PCA_params_binsim.confidence_std}std)'
    PCA_filename_string = PCA_filename_string + '_BinSim'
    middle_page_PCA_binsim[0,1:3] = pn.pane.Plotly(PCA_params_binsim.PCA_plot[0],
                                    config = {'toImageButtonOptions': {'filename': PCA_filename_string, 'scale':4,}})

    # Expected Variance Plot
    if type(PCA_params_binsim.scatter_PCA_plot[0]) == str:
        PCA_params_binsim.exp_var_fig_plot[0] = iaf._plot_PCA_explained_variance(PCA_params_binsim)
    end_page_PCA_binsim[0] = pn.pane.Plotly(PCA_params_binsim.exp_var_fig_plot[0],
                                    config = {'toImageButtonOptions': {'filename': 'PCA_exp_var_plot_BinSim', 'scale':4,}})

    # PCA Scatter Plot
    if type(PCA_params_binsim.scatter_PCA_plot[0]) == str:
        PCA_params_binsim.scatter_PCA_plot[0] = iaf._scatter_PCA_plot(PCA_params_binsim, target_list)
    end_page_PCA_binsim[1] = pn.pane.Plotly(PCA_params_binsim.scatter_PCA_plot[0],
                                    config = {'toImageButtonOptions': {'filename': 'PCA_scatter_plot_BinSim', 'scale':4,}})

    # Enabling Widgets to control the PCA Projection
    for n, w in PCA_params_binsim.controls.widgets.items():
        if n not in ['PCz', 'confidence_std']:
            w.disabled = False

    PCA_params_binsim._update_ellipse_options()
    PCA_params_binsim._update_PCz_disabled()

compute_PCA_binsim_button.on_click(_compute_PCA_binsim_button)


# Initializing layout of the PCA section
# First section of the page
initial_page_PCA_binsim = pn.GridSpec(mode='override')
initial_page_PCA_binsim[0,:2] = n_components_compute_binsim
initial_page_PCA_binsim[0,2] = compute_PCA_binsim_button

# Middle section of the page
middle_page_PCA_binsim = pn.GridSpec(mode='override')
middle_page_PCA_binsim[0,0] = PCA_params_binsim.controls
middle_page_PCA_binsim[0,1:3] = 'To plot a PCA'

# Final section of the page
end_page_PCA_binsim = pn.Column('To plot explained variance figure',
                                'To plot matrices of PCA projections')

# Page sections that will change based on HCA_params
PCA_params_binsim.current_pages_associated.append(middle_page_PCA_binsim)
PCA_params_binsim.current_pages_associated.append(end_page_PCA_binsim)


# Complete PCA page
page_binsim_PCA = pn.Column(pn.pane.HTML("<strong>Principal Component Analysis (PCA)</strong>"),
                            initial_page_PCA_binsim,
                            middle_page_PCA_binsim,
                            end_page_PCA_binsim)



# HCA section of the BinSim analysis
HCA_params_binsim = HCA_Storage()
HCA_params_binsim.binsim_flag = True # Setting the BinSim flag

# Widget to save HCA plot binsim treated (needed since it is a matplotlib plot instead of a plotly plot)
save_HCA_plot_binsim_button = pn.widgets.Button(name='Save as a png (in current folder)', button_type='success',
                                         icon=iaf.download_icon)

# When pressing the button, downloads the figure
def _save_HCA_plot_binsim_button(event):
    "Saves HCA plot (BinSim version)."
    filename = f'HCA_plot_{HCA_params_binsim.dist_metric}Dist_{HCA_params_binsim.link_metric}Linkage_BinSim.png'
    HCA_params_binsim.HCA_plot[0].savefig(filename, dpi=HCA_params_binsim.dpi)
    pn.state.notifications.success(f'Dendrogram successfully saved.')

save_HCA_plot_binsim_button.on_click(_save_HCA_plot_binsim_button)

# Organization for the HCA page section in BinSim Analysis
HCA_binsim = pn.GridSpec(mode='override')
HCA_binsim[0:5,0] = HCA_params_binsim.controls
HCA_binsim[0:6,1:4] = HCA_params_binsim.HCA_plot[0]
HCA_binsim[5,0] = save_HCA_plot_binsim_button

# Page sections that will change based on HCA_params_binsim
HCA_params_binsim.current_pages_associated.append(HCA_binsim)

# Complete HCA page
page_binsim_HCA = pn.Column("<strong>Hierarchical Clustering Analysis (HCA)</strong>",
                           pn.pane.HTML("""Since BinSim treated data is binary (either 0 or 1s), only distance metrics
                                        suited to this type of boolean data will be available."""),
                           HCA_binsim)



# PLS-DA section of the BinSim analysis

# Running initial param to store PLSDA details - BinSim version
PLSDA_store_binsim = PLSDA_Storage()
PLSDA_store_binsim.binsim_flag = True # Setting the BinSim flag
# Disabling the scale attribute since it should not be used with BinSim
PLSDA_store_binsim.controls.widgets['scale'].disabled = True
PLSDA_store_binsim.controls_optim.widgets['scale'].disabled = True

# Click button to confirm PLS Optimization
PLSDA_store_binsim.controls_optim.widgets['confirm_optim_button'].on_click(PLSDA_store_binsim._confirm_optim_button)

# Click button to fit the PLS-DA model and obtain model performance metrics
PLSDA_store_binsim.controls.widgets['confirm_plsda_button'].on_click(PLSDA_store_binsim._confirm_plsda_button)

# Click button to perform PLS-DA Permutation Test
PLSDA_store_binsim.controls_permutation.widgets['confirm_button_permutation'].on_click(
    PLSDA_store_binsim._confirm_button_permutation)

# Function to save the Permutation Test figure as png
def _save_figure_button_permutation_PLSDA_binsim(event):
    "Save PLS-DA permutation figure - BinSim version."
    filename_string = f'PLS-DA_permutation_test_{PLSDA_store_binsim.current_plsda_params_permutation["n_permutations"]}perm_'
    filename_string = filename_string + f'{PLSDA_store_binsim.current_plsda_params_permutation["n_components"]}comp_'
    filename_string = filename_string + f'{PLSDA_store_binsim.current_plsda_params_permutation["n_folds"]}-foldstratCV_scale'
    filename_string = filename_string + f'{PLSDA_store_binsim.current_plsda_params_permutation["scale"]}_metric'
    filename_string = filename_string + f'{PLSDA_store_binsim.current_plsda_params_permutation["perm_metric"]}'
    if PLSDA_store_binsim.binsim_flag:
        filename_string = filename_string + '_BinSim'
    PLSDA_store_binsim.perm_figure[0].savefig(filename_string+'.png', dpi=PLSDA_store_binsim.dpi)
    pn.state.notifications.success(f'Figure {filename_string} successfully saved.')

# Click button to save the aforementioned figure
PLSDA_store_binsim.controls_permutation.widgets['save_figure_button_permutation'].on_click(
    _save_figure_button_permutation_PLSDA_binsim)

# Click button to compute PLS-DA model ROC curves and obtain the corresponding plots
PLSDA_store_binsim.controls_roc.widgets['confirm_button_roc'].on_click(PLSDA_store_binsim._confirm_button_roc)

# Widget to add recommended number of components to page
rec_comp_indicator_binsim_widget = pn.indicators.Number(name='Recommended Components (based on max. Q2)',
                                    font_size='14pt', title_size='14pt', value=PLSDA_store_binsim.rec_components)


# Initial layout of the optimization section of PLS-DA BinSim analysis
pls_optim_section_binsim = pn.Row(pn.Column(PLSDA_store_binsim.controls_optim, rec_comp_indicator_binsim_widget,
                                    'A lower number of components with a similar Q2 may be preferable than the number shown'),
                           PLSDA_store_binsim.optim_figure[0])


# Results Section of the PLS-DA section of BinSim analysis
# Specific Widget for PLS results section of the page, shows DataFrame with only annotated metabolites or all metabolites
plsda_feat_imp_show_annots_only_binsim = pn.widgets.Checkbox(
    name='Only show annotated metabolites in feature importance table', value=False, disabled=True)

# Change the DataFrame shown based on checkbox
@pn.depends(plsda_feat_imp_show_annots_only_binsim.param.value, watch=True)
def _layout_plsda_feat_import_dataframe_binsim(plsda_feat_imp_show_annots_only_binsim):
    "Update the layout based on if we are showing all metabolites or only annotated ones."
    # Select DataFrame
    if plsda_feat_imp_show_annots_only_binsim:
        df_to_show = PLSDA_store_binsim.feat_impor[PLSDA_store_binsim.feat_impor['Has Match?']]
    else:
        df_to_show = PLSDA_store_binsim.feat_impor

    # Update the layout
    pls_results_section_binsim[3] = pn.pane.DataFrame(df_to_show, height=600)


# Widget to save dataframe with features ordered by importance
save_plsda_feat_imp_binsim_button = pn.widgets.Button(
    name='Save PLS-DA Feature Importance table (BinSim version) obtained as .xlsx (in current folder)',
    button_type='warning', icon=iaf.download_icon, disabled=True)

# When pressing the button, downloads the dataframe (builds the appropriate filename)
def _save_plsda_feat_imp_binsim_button(event):
    "Save PLS-DA Feature Importance results as an Excel."
    try:
        plsda_params = PLSDA_store_binsim.current_plsda_params
        # Building the datafile name
        filename_string = f'PLS-DA_FeatImp_{plsda_params["feat_imp"]}_model_params_components{plsda_params["n_components"]}'
        filename_string = filename_string + f'_{plsda_params["n_folds"]}-foldstratCV_iterations{plsda_params["n_iterations"]}'
        filename_string = filename_string + f'_scale{plsda_params["scale"]}'
        if PLSDA_store_binsim.binsim_flag:
            filename_string = filename_string + '_BinSim'
        filename_string = filename_string + f'.xlsx'

        # Saving the file
        PLSDA_store_binsim.feat_impor.to_excel(filename_string)
        pn.state.notifications.success(f'{filename_string} successfully saved.')
    except:
        pn.state.notifications.error(f'File could not be saved.')

save_plsda_feat_imp_binsim_button.on_click(_save_plsda_feat_imp_binsim_button)


# PLS projection section of the page - BinSim Version
pls_proj_page_binsim = pn.GridSpec(mode='override')
pls_proj_page_binsim[0,0] = PLSDA_store_binsim.controls_projection
pls_proj_page_binsim[0,1:3] = 'To plot a PLS'

# Layout of the full results section (Partial, more is added when fitting PLS-DA model)
pls_results_section_binsim = pn.Column(pn.Row(PLSDA_store_binsim.controls,
                                       pn.Column('### Model Performance Metrics', PLSDA_store_binsim.n_results)),
                                '### Feature Importance Table',
                                plsda_feat_imp_show_annots_only_binsim,
                               pn.pane.DataFrame(PLSDA_store_binsim.feat_impor),
                               save_plsda_feat_imp_binsim_button,)

# Page sections / Widgets that will change based on PLSDA_store_binsim
PLSDA_store_binsim.current_pages_associated.append(pls_optim_section_binsim)
PLSDA_store_binsim.current_pages_associated.append(pls_results_section_binsim)
PLSDA_store_binsim.current_pages_associated.append(pls_proj_page_binsim)
PLSDA_store_binsim.current_pages_associated.append(plsda_feat_imp_show_annots_only_binsim)
PLSDA_store_binsim.current_pages_associated.append(save_plsda_feat_imp_binsim_button)


# Complete layout of the PLS-DA section of BinSim analysis
page_binsim_PLSDA = pn.Column(pn.pane.HTML(plsda_opening_string),
                              pls_optim_section_binsim,
                              pls_results_section_binsim)



# Random Forest section of the BinSim Analysis

# Running initial param to store RF details - BinSim version
RF_store_binsim = RF_Storage()
RF_store_binsim.binsim_flag = True # Setting the BinSim Flag

# Click button to confirm Random Forest Optimization
RF_store_binsim.controls_optim.widgets['confirm_optim_button'].on_click(RF_store_binsim._confirm_optim_button)

# Click button to fit the Random Forest model and obtain model performance metrics
RF_store_binsim.controls.widgets['confirm_rf_button'].on_click(RF_store_binsim._confirm_rf_button)

# Click button to perform Random Forest Permutation Test
RF_store_binsim.controls_permutation.widgets['confirm_button_permutation'].on_click(
    RF_store_binsim._confirm_button_permutation)

# Function to save the Permutation Test figure as png
def _save_figure_button_permutation_RF_binsim(event):
    "Saves Random Forest permutation figure - BinSim version."
    filename_string = f'RF_permutation_test_{RF_store_binsim.current_rf_params_permutation["n_permutations"]}perm_'
    filename_string = filename_string + f'{RF_store_binsim.current_rf_params_permutation["n_trees"]}trees_'
    filename_string = filename_string + f'{RF_store_binsim.current_rf_params_permutation["n_folds"]}-foldstratCV_metric'
    filename_string = filename_string + f'{RF_store_binsim.current_rf_params_permutation["perm_metric"]}'
    if RF_store_binsim.binsim_flag:
        filename_string = filename_string + '_BinSim'
    RF_store_binsim.perm_figure[0].savefig(filename_string+'.png', dpi=RF_store_binsim.dpi)
    pn.state.notifications.success(f'Figure {filename_string} successfully saved.')


# Click button to save the aforementioned figure
RF_store_binsim.controls_permutation.widgets['save_figure_button_permutation'].on_click(
    _save_figure_button_permutation_RF_binsim)

# Click button to compute Random Forest model ROC curves and obtain the corresponding plots
RF_store_binsim.controls_roc.widgets['confirm_button_roc'].on_click(RF_store_binsim._confirm_button_roc)


# Initial layout of the optimization section of Random Forest BinSim analysis
rf_optim_section_binsim = pn.Row(RF_store_binsim.controls_optim, RF_store_binsim.optim_figure[0])


# Results Section of the Random Forest section of BinSim analysis
# Specific Widget for Random Forest results section of the page, shows DataFrame with only annotated metabolites or all metabolites
rf_feat_imp_show_annots_only_binsim = pn.widgets.Checkbox(
    name='Only show annotated metabolites in feature importance table', value=False, disabled=True)

# Change the DataFrame shown based on checkbox
@pn.depends(rf_feat_imp_show_annots_only_binsim.param.value, watch=True)
def _layout_rf_feat_import_dataframe_binsim(rf_feat_imp_show_annots_only_binsim):
    "Update the layout based on if we are showing all metabolites or only annotated ones."
    # Select DataFrame
    if rf_feat_imp_show_annots_only_binsim:
        df_to_show = RF_store_binsim.feat_impor[RF_store_binsim.feat_impor['Has Match?']]
    else:
        df_to_show = RF_store_binsim.feat_impor

    # Update the layout
    rf_results_section_binsim[3] = pn.pane.DataFrame(df_to_show, height=600)


# Widget to save dataframe with features ordered by importance
save_rf_feat_imp_binsim_button = pn.widgets.Button(
    name='Save Random Forest Feature Importance table (BinSim version) obtained as .xlsx (in current folder)',
    button_type='warning', icon=iaf.download_icon, disabled=True)

# When pressing the button, downloads the dataframe (builds the appropriate filename)
def _save_rf_feat_imp_binsim_button(event):
    "Save Random Forest Feature Importance results as an Excel."
    try:
        rf_params = RF_store_binsim.current_rf_params
        # Building the datafile name
        filename_string = f'RF_FeatImp_Gini_model_params_{rf_params["n_trees"]}trees_{rf_params["n_folds"]}-foldstratCV_'
        filename_string = filename_string + f'iterations{rf_params["n_iterations"]}'
        if RF_store_binsim.binsim_flag:
            filename_string = filename_string + '_BinSim'
        filename_string = filename_string + '.xlsx'

        # Saving the file
        RF_store_binsim.feat_impor.to_excel(filename_string)
        pn.state.notifications.success(f'{filename_string} successfully saved.')
    except:
        pn.state.notifications.error(f'File could not be saved.')

save_rf_feat_imp_binsim_button.on_click(_save_rf_feat_imp_binsim_button)


# Layout of the full results section (Partial, more is added when fitting RF model)
rf_results_section_binsim = pn.Column(pn.Row(RF_store_binsim.controls,
                                       pn.Column('### Model Performance Metrics', RF_store_binsim.n_results)),
                                '### Feature Importance Table',
                                rf_feat_imp_show_annots_only_binsim,
                               pn.pane.DataFrame(RF_store_binsim.feat_impor),
                               save_rf_feat_imp_binsim_button,)

# Page sections / Widgets that will change based on RF_store_binsim
RF_store_binsim.current_pages_associated.append(rf_optim_section_binsim)
RF_store_binsim.current_pages_associated.append(rf_results_section_binsim)
RF_store_binsim.current_pages_associated.append(rf_feat_imp_show_annots_only_binsim)
RF_store_binsim.current_pages_associated.append(save_rf_feat_imp_binsim_button)


# Complete layout of the PLS-DA section of BinSim analysis
page_binsim_RF = pn.Column(pn.pane.HTML(rf_opening_string),
                           rf_optim_section_binsim,
                           rf_results_section_binsim)



# Complete BinSim Analysis Page Layout
binsim_analysis_page = pn.Column(pn.Tabs(('PCA', page_binsim_PCA), ('HCA', page_binsim_HCA),
                                        ('PLS-DA', page_binsim_PLSDA), ('RF ', page_binsim_RF)))




# Page for Compound Finder
class CompoundFinder(param.Parameterized):
    "Class to contain parameters and figures related to a compound finder tool."

    # Parameters to detect which compound to see
    id_type = param.String(default='Metabolite Bucket Label')
    id_comp = param.String()
    confirm_id_button = param.Boolean(default=False)

    # Store Possibilities
    bucket_to_idxs = param.Dict()
    name_to_idxs = param.Dict()
    formula_to_idxs = param.Dict()
    neutral_mass_to_idxs = param.Dict()

    # DataFrame to store found ID
    id_df = param.DataFrame()
    current_params = param.Dict()

    # Figure parameters
    sample_bar_plot_type = param.String(default='Horizontal')
    class_bar_plot_type = param.String(default='Horizontal')
    class_boxplot_points = param.String(default='Only outliers')
    ignore_missing_values = param.Boolean(default=True)

    # Storing Figures
    sample_bar_plot = param.List(default=['To plot a bar plot with sample intensities'])
    class_bar_plot = param.List(default=['To plot a bar plot with class avg. intensities'])
    class_boxplot = param.List(default=['To plot a boxplot with class avg. intensities'])


    def _calculate_possible_options(self, metadata_df, checkbox_annotation, checkbox_formula, radiobox_neutral_mass):
        "Obtain a dictionary detailing all the possible peak/annotation options and corresponding idxs where they appear."

        # Get columns list for annotation possibilities
        name_cols = checkbox_annotation.value + [
            i for i in metadata_df.columns if i.startswith('Matched') and i.endswith('names')]
        formula_cols = checkbox_formula.value + [
            i for i in metadata_df.columns if i.startswith('Matched') and i.endswith('formulas')]

        # Obtain the dictionaries (keys are type of identifiers and values are the idxs where they appear)
        self.name_to_idxs = iaf.build_annotation_to_idx_dict(metadata_df, name_cols)
        self.formula_to_idxs = iaf.build_annotation_to_idx_dict(metadata_df, formula_cols)
        self.bucket_to_idxs = {str(idx): [idx] for idx in metadata_df.index}

        # If there is a neutral mass column
        if radiobox_neutral_mass.value != 'None':
            self.neutral_mass_to_idxs = {
                str(metadata_df.loc[idx, radiobox_neutral_mass.value]): [idx] for idx in metadata_df.index}
            self.controls.widgets['id_type'].options = ['Metabolite Bucket Label', 'Metabolite Name', 'Metabolite Formula',
                                                       'Metabolite Neutral Mass']

        # If there isn't remove it as an option
        else:
            self.neutral_mass_to_idxs = {}
            self.controls.widgets['id_type'].options = ['Metabolite Bucket Label', 'Metabolite Name', 'Metabolite Formula']

        # Default option for id_type - setting up widgets
        self.controls.widgets['id_comp'].options = list(metadata_df.index)
        self.controls.widgets['id_comp'].placeholder = metadata_df.index[0]


    def find_id_df(self, processed_df, groups):
        "Find the DataFrame subsection correponding to the identifier given."

        # Finding the DataFrame
        if self.id_type == 'Metabolite Bucket Label':
            finder = processed_df.loc[self.bucket_to_idxs[self.id_comp]].copy()
        elif self.id_type == 'Metabolite Formula':
            finder = processed_df.loc[self.formula_to_idxs[self.id_comp]].copy()
        elif self.id_type == 'Metabolite Name':
            finder = processed_df.loc[self.name_to_idxs[self.id_comp]].copy()
        elif self.id_type == 'Metabolite Neutral Mass':
            finder = processed_df.loc[self.neutral_mass_to_idxs[self.id_comp]].copy()
        else:
            pn.state.notifications('Type of identifier provided not recognized.')

        # Adding Average and Standard Deviation columns for each class
        for g in groups:
            finder[g+' Average'] = finder[finder.columns.intersection(groups[g])].mean(axis=1)
            finder[g+' std'] = finder[finder.columns.intersection(groups[g])].std(axis=1)

        return finder # Return DataFrame


    def _confirm_id_button(self, event):
        "Actions to find the identifier provided and plot the plots."

        # Calculate specific class DataFrames in case it has not been calculated before
        if type(com_exc_compounds.group_dfs) != dict:
            iaf._group_compounds_per_class(com_exc_compounds, target_list, DataFrame_Store) # Add compounds per class dfs

        # Saving parameters
        self.current_params = {'id_type': self.id_type, 'id_comp': self.id_comp}

        # Obtaining the DataFrame from the processed data
        self.id_df = self.find_id_df(DataFrame_Store.processed_df, com_exc_compounds.groups)

        # Plotting the sample bar plot
        self.sample_bar_plot[0] = iaf.plot_sample_bar_plot(self, target_list, com_exc_compounds)

        # Plotting the class bar plot
        self.class_bar_plot[0] = iaf.plot_class_bar_plot(self, target_list, com_exc_compounds)

        # Plotting the class boxplot
        self.class_boxplot[0] = iaf.plot_class_boxplot(self, target_list, com_exc_compounds)

        # Update the layout
        comp_finder_page[1] = pn.pane.DataFrame(comp_finder.id_df)
        comp_finder_page[3] = pn.pane.Plotly(self.sample_bar_plot[0], config={'toImageButtonOptions': {
                   'filename': f'NormInt_SampleBarPlot_{self.id_type}_{self.id_comp}', 'scale':4}})
        if self.ignore_missing_values:
            mv = 'Ignored'
        else:
            mv = 'As0'
        comp_finder_page[5] = pn.pane.Plotly(self.class_bar_plot[0], config={'toImageButtonOptions': {
                   'filename': f'NormInt_ClassBarPlot_{self.id_type}_{self.id_comp}_MissValues{mv}', 'scale':4}})
        comp_finder_page[7] = pn.pane.Plotly(self.class_boxplot[0], config={'toImageButtonOptions': {
                   'filename': f'NormInt_ClassBoxplot_{self.id_type}_{self.id_comp}_MissValues{mv}', 'scale':4}})


    def reset(self):
        "Resets all relevant parameters and Widgets."
        self.name_to_idxs = None
        self.formula_to_idxs = None
        self.neutral_mass_to_idxs = None
        self.id_df = None
        self.current_params = None
        self.sample_bar_plot = ['To plot a bar plot with sample intensities']
        self.class_bar_plot = ['To plot a bar plot with class avg. intensities']
        self.class_boxplot = ['To plot a boxplot with class avg. intensities']
        self.id_type = 'Metabolite Bucket Label'
        self.controls.widgets['id_type'].value = 'Metabolite Bucket Label'
        self.controls.widgets['id_type'].options = ['Metabolite Bucket Label', 'Metabolite Name', 'Metabolite Formula', 'Metabolite Neutral Mass']
        self.id_comp = self.param['id_comp'].default
        self.controls.widgets['id_comp'].options = ['']
        self.controls.widgets['id_comp'].placeholder = ''
        self.bucket_to_idxs = None


    # Update the Widget options
    @param.depends('id_type', watch=True)
    def _update_widgets(self):
        "Update Widgets (specifically id_comp) based on the id_type chosen."

        self.controls.widgets['id_comp'].value = '' # Reset the value
        self.id_comp = ''

        # Update the options based on id_type chosen
        if self.id_type == 'Metabolite Bucket Label':
            self.controls.widgets['id_comp'].options = list(self.bucket_to_idxs.keys())
            self.controls.widgets['id_comp'].placeholder = list(self.bucket_to_idxs.keys())[0]

        elif self.id_type == 'Metabolite Name':
            self.controls.widgets['id_comp'].options = list(self.name_to_idxs.keys())
            self.controls.widgets['id_comp'].placeholder = list(self.name_to_idxs.keys())[0]

        elif self.id_type == 'Metabolite Formula':
            self.controls.widgets['id_comp'].options = list(self.formula_to_idxs.keys())
            self.controls.widgets['id_comp'].placeholder = list(self.formula_to_idxs.keys())[0]

        else:
            self.controls.widgets['id_comp'].options = list(self.neutral_mass_to_idxs.keys())
            self.controls.widgets['id_comp'].placeholder = list(self.neutral_mass_to_idxs.keys())[0]


    @param.depends('sample_bar_plot_type', watch=True)
    def _update_sample_bar_plot(self):
        "Update the normalized sample bar plot (horizontal or vertical)."
        # Plotting the sample bar plot
        self.sample_bar_plot[0] = iaf.plot_sample_bar_plot(self, target_list, com_exc_compounds)
        # Update the layout
        comp_finder_page[3] = pn.pane.Plotly(self.sample_bar_plot[0], config={'toImageButtonOptions': {
                   'filename': f'NormInt_SampleBarPlot_{self.id_type}_{self.id_comp}', 'scale':4}})


    @param.depends('class_bar_plot_type', 'ignore_missing_values', watch=True)
    def _update_class_bar_plot(self):
        "Update the normalized class bar plot."
        # Plotting the class bar plot
        self.class_bar_plot[0] = iaf.plot_class_bar_plot(self, target_list, com_exc_compounds)
        # Update the layout
        if self.ignore_missing_values:
            mv = 'Ignored'
        else:
            mv = 'As0'
        comp_finder_page[5] = pn.pane.Plotly(self.class_bar_plot[0], config={'toImageButtonOptions': {
                   'filename': f'NormInt_ClassBarPlot_{self.id_type}_{self.id_comp}_MissValues{mv}', 'scale':4}})


    @param.depends('class_boxplot_points', 'ignore_missing_values', watch=True)
    def _update_class_boxplot(self):
        "Update the normalized class boxplot."
        # Plotting the class boxplot
        self.class_boxplot[0] = iaf.plot_class_boxplot(self, target_list, com_exc_compounds)
        # Update the layout
        if self.ignore_missing_values:
            mv = 'Ignored'
        else:
            mv = 'As0'
        comp_finder_page[7] = pn.pane.Plotly(self.class_boxplot[0], config={'toImageButtonOptions': {
                   'filename': f'NormInt_ClassBoxplot_{self.id_type}_{self.id_comp}_MissValues{mv}', 'scale':4}})


    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets = {
            'id_type': pn.widgets.RadioButtonGroup(name='Choose which type of identifier you want to use for your compound',
                    value='Metabolite Bucket Label',
                    options=['Metabolite Bucket Label', 'Metabolite Name', 'Metabolite Formula', 'Metabolite Neutral Mass'],
                    button_type='success', description='Choose which type of identifier you want to use for your compound.'),
            'id_comp': pn.widgets.AutocompleteInput(name="Type the compound identifier you want to see",
                    value = '', options=[''], search_strategy='includes', case_sensitive=False,
                    placeholder='Glutathione'),
            'confirm_id_button': pn.widgets.Button(name="Find Compound", button_type='primary'),
        }

        # Widget Figure params
        widgets_fig_param = {
            'sample_bar_plot_type': pn.widgets.Select(name='Type of bar plot:',
                                value='Horizontal', options=['Horizontal', 'Vertical']),
            'class_bar_plot_type': pn.widgets.Select(name='Type of bar plot:',
                                value='Horizontal', options=['Horizontal', 'Vertical']),
            'class_boxplot_points': pn.widgets.Select(name='What points appear in the Boxplot:',
                                value='Only Outliers', options=['Only outliers', 'All points']),
            'ignore_missing_values': pn.widgets.Checkbox(value=True,
                name='Class based plots ignore missing values (if not selected missing values are treated as 0).')
        }

        self.controls = pn.Param(self, parameters=['id_type', 'id_comp', 'confirm_id_button'],
                                 widgets=widgets, name='Compound Identifier to Find')

        self.controls_fig = pn.Param(self, parameters=['sample_bar_plot_type', 'class_bar_plot_type',
                                                       'class_boxplot_points', 'ignore_missing_values'],
                                 widgets=widgets_fig_param, name='Figure Parameters')

# Initialize Store
comp_finder = CompoundFinder()

# Calling the function when the button is pressed
comp_finder.controls.widgets['confirm_id_button'].on_click(comp_finder._confirm_id_button)


# Initializing the layout of the page
comp_finder_page = pn.Column(comp_finder.controls,
                            'DataFrame of the Searched Compound',
                             pn.Column('### Sample Bar Plot',
                                       comp_finder.controls_fig.widgets['sample_bar_plot_type']),
                            'Sample Bar Plot of the Searched Compound',
                             pn.Column('### Class Bar Plot',
                                    comp_finder.controls_fig.widgets['class_bar_plot_type'],
                                    comp_finder.controls_fig.widgets['ignore_missing_values']),
                            'Class Bar Plot of the Searched Compound',
                            pn.Column('### Class Boxlot',
                                      comp_finder.controls_fig.widgets['class_boxplot_points'],
                                    comp_finder.controls_fig.widgets['ignore_missing_values']),
                            'Class Boxplot of the Searched Compound',)




# Report Generation Page
class ReportGeneration(param.Parameterized):
    "Class use to store parameters useful for report generation."

    # Important attributes to build the report

    # Related to File Reading
    filename = param.String('')
    target_included_in_file = param.Boolean(default=False)
    neutral_mass_column = param.Boolean(default=False)

    # Related to Annotation
    annotation_margin_method = param.String('')
    annotation_margin_ppm_deviation = param.Number()
    annotation_margin_Da_deviation = param.Number()

    # Related to Annotation De-Duplication
    deduplication_performed = param.String('Annotation De-Duplication was not performed')

    # Checkboxes to select which analysis will be shown
    com_exc_analysis = param.List(default=[])
    unsup_analysis = param.List(default=['PCA', 'HCA',])
    sup_analysis = param.List(default=[])
    univ_analysis = param.Boolean(default=False)
    dataviz_analysis =param.List(default=[])
    pathassign_analysis = param.Boolean(default=False)
    BinSim_analysis = param.List(default=[])


    def reset(self):
        "Resets all relevant parameters."
        for param in self.param:
            if param not in ["name"]:
                setattr(self, param, self.param[param].default)


    def __init__(self, **params):

        super().__init__(**params)
        # Base Widgets
        widgets = {
            'com_exc_analysis': pn.widgets.CheckBoxGroup(name='Common and Exclusive Compounds Analysis',
                    value=[], options=['Overview', 'Venn Diagram', 'Intersection Plot']),
            'unsup_analysis': pn.widgets.CheckBoxGroup(name='Unsupervised Analysis',
                    value=['PCA', 'HCA'], options=['PCA', 'HCA']),
            'sup_analysis': pn.widgets.CheckBoxGroup(name='Supervised Analysis',
                    value=[], options=['PLS-DA', 'Random Forest']),
            'univ_analysis': pn.widgets.Checkbox(name='Univariate Analysis', value=False),
            'dataviz_analysis': pn.widgets.CheckBoxGroup(name='Data Diversity Visualization Analysis',
                    value=[], options=['Van Krevelen', 'Kendrick Mass Defect', 'Chem. Comp. Series']),
            'pathassign_analysis': pn.widgets.Checkbox(name='HMDB Pathway Assignment', value=False),
            'BinSim_analysis': pn.widgets.CheckBoxGroup(name='BinSim Analysis',
                    value=[], options=['PCA', 'HCA', 'PLS-DA', 'Random Forest']),
        }
        self.controls = pn.Param(self, parameters=['com_exc_analysis', 'unsup_analysis', 'sup_analysis', 'univ_analysis',
                                                   'dataviz_analysis', 'pathassign_analysis', 'BinSim_analysis'],
                                widgets=widgets, name='', default_layout=pn.Row)


# Initializing the store for Report Generation
RepGen = ReportGeneration()

# Setting up the widgets for the page layout
desc_repgen = pn.Row('#### Common and Exclusive Compound Analysis',
              '#### Unsupervised Analysis',
              '#### Supervised Analysis',
              '#### Univariate Analysis',
              '#### Data Diversity Visualization Analysis',
              '#### HMDB Pathways Assignment',
              '#### BinSim Analysis',)

# Widget to select folder where report figures and tables will be created in
folder_selection = pn.widgets.TextAreaInput(name='Folder Name where Report and associated Figures and Tables will be Downloaded to (Do Not put a name of a pre-existing folder)',
                                           placeholder='Report_folder', value='Report', rows=1)

# Widgets to select what type of analysis will be in the report
widgets_repgen = pn.Row()
widgets_repgen.extend([RepGen.controls.widgets[i] for i in RepGen.controls.widgets])

# Widget to perform the generation of the report
report_generation_button = pn.widgets.Button(
    name='Generate Report File and Folder (Report + Figures + Tables) of Analysis Performed',
    button_type='warning', icon=iaf.download_icon)
def _report_generation_button(event):
    "Calls the function that generates the report while pressed."

    while len(rep_gen_page) > 5:
        rep_gen_page.pop(-1)

    # Loading Widget while report is being generated
    rep_gen_page.append(pn.indicators.LoadingSpinner(value=True, size=100,
                                                        name='Report is currently being generated with selected statistical analysis...'))

    # Perform Report Generation
    try:
        ReportGenerator(folder_selection.value, RepGen, file, checkbox_annotation, checkbox_formula, radiobox_neutral_mass, checkbox_others,
                        target_list, UnivarA_Store, characteristics_df, DataFrame_Store, n_databases, DB_dict, verbose_annotated_compounds,
                        data_ann_deduplicator, com_exc_compounds, PCA_params, HCA_params, PLSDA_store, RF_store, dataviz_store, PathAssign_store,
                        PCA_params_binsim, HCA_params_binsim, PLSDA_store_binsim, RF_store_binsim, rep_gen_page)
    except:
        while len(rep_gen_page) > 5:
            rep_gen_page.pop(-1)
        pn.state.notifications.error('Report was not successfully generated.')
        raise ValueError('Report was not successfully generated.')

    # When finished
    pn.state.notifications.success(f'Report successfully generated in {folder_selection.value} folder.')
    rep_gen_page[5] = f'### Report successfully generated in {folder_selection.value} folder'
report_generation_button.on_click(_report_generation_button)


# Initializing the page
rep_gen_page = pn.Column(pn.pane.HTML(desc_str.report_opening_string),
                         folder_selection,
                         desc_repgen,
                         widgets_repgen,
                         report_generation_button)




# RESET 'page'

# Reset panel widgets and associated functions
reset_panel_float_text = pn.widgets.StaticText(name='', value='Do you want to reset all parameters of the software? (This is a faux Reset with a few tricks but it seems to do the job)')
reset_panel_float_yes_button = pn.widgets.Button(name='Yes', button_type='danger')
reset_panel_float_no_button = pn.widgets.Button(name='No', button_type='default')
reset_floatpanel = pn.layout.FloatPanel(reset_panel_float_text, reset_panel_float_yes_button, reset_panel_float_no_button,
                                      name='Reset', margin=20, position='center', config={"headerControls": {"close": "remove"}},
                                      contained=False)


def RESET(event):
    "Open the Reset floatpanel for confirmation (or not) that you want to perform a reset."
    reset_panel = pn.Column(reset_floatpanel, height=200)
    #main_area.clear()
    if len(main_area) == 1:
        main_area.append(reset_panel)
        reset_floatpanel.status = 'normalized'
    else:
        main_area.append(reset_panel)
        reset_floatpanel.status = 'normalized'


def No_Reset(event):
    "Closing the Reset floatpanel."
    reset_floatpanel.status = 'closed'
    main_area.pop(-1)


reset_time = pn.widgets.IntInput(value=0)
def Yes_Reset(event):
    "Resetting the program and all variables."
    # Resetting pre-treatment steps

    # Resetting Data Reading page speficically (parameters and widgets)
    file.reset()
    filename.value = None
    filename.disabled = False
    target_included_in_file.value = False
    target_included_in_file.disabled = False
    load_example_df_button.disabled = False
    confirm_button_filename.disabled = True
    confirm_button_step1.disabled = True
    section1page[2] = file.read_df

    # Resetting other Pre-Treatment related Param classes
    PreTreatment_Method.reset()

    # Disable all sidebar buttons except data reading
    page1_1_button.disabled = True
    page1_2_button.disabled = True
    page2_button.disabled = True
    page2_1_button.disabled = True
    page3_button.disabled = True
    page4_button.disabled = True
    _disabling_stat_analysis_buttons()

    # Resetting all statistical analysis performed

    # Common and Exclusive Compound page
    com_exc_compounds.reset()
    while len(end_page_comexc) > 0:
        end_page_comexc.pop(-1)

    # Unsupervised Analysis page
    n_components_compute.value = 10
    PCA_params.reset()
    HCA_params.reset()

    # Supervised Analysis page
    plsda_feat_imp_show_annots_only.value = False
    plsda_feat_imp_show_annots_only.disabled = True
    rec_comp_indicator_widget.value = None
    PLSDA_store.soft_reset()
    pls_optim_section[1] = PLSDA_store.optim_figure[0]
    save_plsda_feat_imp_button.disabled = True
    while len(pls_results_section) > 5:
        pls_results_section.pop(-1)
    pls_results_section[0][1][1] = PLSDA_store.n_results
    pls_results_section[3] = pn.pane.DataFrame(PLSDA_store.feat_impor)

    rf_feat_imp_show_annots_only.value = False
    rf_feat_imp_show_annots_only.disabled = True
    RF_store.reset()
    rf_optim_section[1] = RF_store.optim_figure[0]
    save_rf_feat_imp_button.disabled = True
    while len(rf_results_section) > 5:
        rf_results_section.pop(-1)
    rf_results_section[0][1][1] = RF_store.n_results
    rf_results_section[3] = pn.pane.DataFrame(RF_store.feat_impor)

    # Univariate Analysis page
    UnivarA_Store.reset()
    while len(univar_analysis_page) > 2:
        univar_analysis_page.pop(-1)

    # Data Diversity Visualization page
    iaf._group_compounds_per_class(com_exc_compounds, target_list, DataFrame_Store) # Add compounds per class dfs
    dataviz_store.reset()
    vk_plots.clear()
    kmd_plots.clear()
    while len(vk_page) > 1:
        vk_page.pop(-1)
    while len(kmd_page) > 1:
        kmd_page.pop(-1)
    while len(ccs_page) > 1:
        ccs_page.pop(-1)

    # Pathway Assignment page
    PathAssign_store.reset()
    while len(path_matching_section) > 3:
        path_matching_section.pop(-1)
    while len(hmdb_id_searching_section) > 3:
        hmdb_id_searching_section.pop(-1)

    # Binary Simplification (BinSim) analysis page
    # PCA
    n_components_compute_binsim.value = 10
    PCA_params_binsim.reset()
    PCA_params_binsim.binsim_flag = True
    for _, w in PCA_params_binsim.controls.widgets.items():
        w.disabled = True
    middle_page_PCA_binsim[0,1:3] = 'To plot a PCA'
    end_page_PCA_binsim[0] = 'To plot explained variance figure'
    end_page_PCA_binsim[1] = 'To plot matrices of PCA projections'
    # HCA
    HCA_params_binsim.reset()
    HCA_params_binsim.binsim_flag = True
    # PLS-DA
    plsda_feat_imp_show_annots_only_binsim.value = False
    plsda_feat_imp_show_annots_only_binsim.disabled = True
    rec_comp_indicator_binsim_widget.value = None
    PLSDA_store_binsim.soft_reset()
    pls_optim_section_binsim[1] = PLSDA_store_binsim.optim_figure[0]
    save_plsda_feat_imp_binsim_button.disabled = True
    while len(pls_results_section_binsim) > 5:
        pls_results_section_binsim.pop(-1)
    pls_results_section_binsim[0][1][1] = PLSDA_store_binsim.n_results
    pls_results_section_binsim[3] = pn.pane.DataFrame(PLSDA_store_binsim.feat_impor)
    PLSDA_store_binsim.binsim_flag = True
    PLSDA_store_binsim.controls.widgets['scale'].disabled = True
    PLSDA_store_binsim.controls_optim.widgets['scale'].disabled = True
    # RF
    rf_feat_imp_show_annots_only_binsim.value = False
    rf_feat_imp_show_annots_only_binsim.disabled = True
    RF_store_binsim.reset()
    rf_optim_section_binsim[1] = RF_store_binsim.optim_figure[0]
    save_rf_feat_imp_binsim_button.disabled = True
    while len(rf_results_section_binsim) > 5:
        rf_results_section_binsim.pop(-1)
    rf_results_section_binsim[0][1][1] = RF_store_binsim.n_results
    rf_results_section_binsim[3] = pn.pane.DataFrame(RF_store_binsim.feat_impor)
    RF_store_binsim.binsim_flag = True

    # Compound Finder search tool page
    comp_finder.reset()
    comp_finder_page[1] = 'DataFrame of the Searched Compound'
    comp_finder_page[3] = 'Sample Bar Plot of the Searched Compound'
    comp_finder_page[5] = 'Class Bar Plot of the Searched Compound'
    comp_finder_page[7] = 'Class Boxplot of the Searched Compound'
    com_exc_compounds.reset()

    # Report Generation page
    RepGen.reset()

    # Target reset
    target_list.reset()
    reset_time.value = 1

    # Resetting page layouts of the pre-treatment section
    # Data Metadata page
    checkbox_formula.value = ['Formula',]
    checkbox_annotation.value = ['Name',]
    radiobox_neutral_mass.value = 'Neutral Mass'
    checkbox_others.value = []
    checkbox_samples.value = []
    target_widget.disabled = True
    target_widget.value = ''

    # Data Filtering page
    filt_method.value = "Total Samples"
    filt_kw.value = 2
    while len(page1_2) > 3:
        page1_2.pop(-1)

    # Data Annotation page
    n_databases_show.value = 1
    n_databases.value = 1
    annotation_margin_method_radio.value = 'PPM Deviation'
    annotation_ppm_deviation.value = 1
    annotation_Da_deviation.value = 0.001
    annotated_df.value = pd.DataFrame(index=filtered_df.value.index)
    #DB_dict_reset(DB_dict)
    dbs_arrangement.clear()
    dbs_arrangement_all.clear()
    dbs_arrangement_all.extend([DB_dict['1'].content, DB_dict['2'].content, DB_dict['3'].content,
                   DB_dict['4'].content, DB_dict['5'].content])
    while len(page2) > 1:
        page2.pop(-1)

    # Annotation De-Duplication page
    data_ann_deduplicator.reset()
    merge_problems_widgets.clear()
    while len(page2_1)>3:
        page2_1.pop(-1)
    page2_1[1] = pn.pane.DataFrame(data_ann_deduplicator.mult_ann_report)
    middle_section_dedup[2] = data_ann_deduplicator.merge_report
    middle_section_dedup[3] = pn.Row(data_ann_deduplicator.merge_situations, data_ann_deduplicator.mergings_performed)
    middle_section_dedup[5] = data_ann_deduplicator.merge_description
    middle_section_dedup[7] = data_ann_deduplicator.merge_problems
    merge_problems_section_page.clear()

    # Close the floatpanel
    reset_floatpanel.status = 'closed'
    main_area.pop(-1)

    # Show the homepage
    # Update the Main page
    main_area.clear()
    show_page(pages["Index"])


# Functions to apply on each main button of the Reset floatpanel
reset_panel_float_no_button.on_click(No_Reset)
reset_panel_float_yes_button.on_click(Yes_Reset)




# Overall layout of the program and initialization
# The pages we have
pages = {
    "Index": OpeningPage(),
    "Instructions": InstructionPage(),
    "Data Reading": DataReading(),
    "Data Metadata": DataMetadata(),
    "Data Filtering": DataFiltering(),
    "Data Annotation": DataAnnotation(),
    "Annotation De-Duplication": AnnDeDuplication(),
    "Data Pre-Treatment": DataPreTreatment(),
    "Class Colours": ClassColours(),
    "Transitional Page": TransitionalPage(),
    "Common and Exclusive Compounds": CommonExclusivePage(),
    "Unsupervised Analysis": UnsupervisedAnalysisPage(),
    "Supervised Analysis": SupervisedAnalysisPage(),
    "Univariate Analysis": UnivariateAnalysisPage(),
    "Data Visualization": DataVisualizationPage(),
    "Pathway Assignment": PathwayAssignmentPage(),
    "BinSim Analysis": BinSimPage(),
    "Compound Finder": CompoundFinderPage(),
    "Report Generation": ReportGenerationPage(),
}

# Function to show the selected page - needs update
def show_page(page_instance):
    "Shows page chosen in the main area."
    main_area.clear()
    def clear_again():
        main_area.clear()
    clear_again()
    #main_area.get_param_values()
    main_area.append(page_instance.view())


# Create and organize the sidebar
# The sidebar will include buttons to every main page and will be used as a tool for navigation

# Define buttons to navigate between pages
index_button = pn.widgets.Button(name="Home", button_type="success", icon=iaf.home_icon, disabled=False)
instruction_button = pn.widgets.Button(name="Instructions and Considerations", button_type="warning", icon=iaf.instruction_icon, disabled=False)
page1_button = pn.widgets.Button(name="Data Reading", button_type="primary")
page1_1_button = pn.widgets.Button(name="Data Metadata", button_type="primary", disabled=True)
page1_2_button = pn.widgets.Button(name="Data Filtering", button_type="primary", disabled=True)
page2_button = pn.widgets.Button(name="Data Annotation", button_type="primary", disabled=True)
page2_1_button = pn.widgets.Button(name="Annotation De-Duplication", button_type="primary", disabled=True)
page3_button = pn.widgets.Button(name="Data Pre-Treatment", button_type="primary", disabled=True)
page4_button = pn.widgets.Button(name="Class Colours", button_type="primary", disabled=True)
page5_button = pn.widgets.Button(name="Common/Exclusive Comp.", button_type="default", disabled=True)
page6_button = pn.widgets.Button(name="Unsupervised Analysis", button_type="default", disabled=True)
page7_button = pn.widgets.Button(name="Supervised Analysis", button_type="default", disabled=True)
page8_button = pn.widgets.Button(name="Univariate Analysis", button_type="default", disabled=True)
page9_button = pn.widgets.Button(name="Data Visualization", button_type="default", disabled=True)
page10_button = pn.widgets.Button(name="Pathway Assignment", button_type="default", disabled=True)
page11_button = pn.widgets.Button(name="BinSim Analysis", button_type="default", disabled=True)
page12_button = pn.widgets.Button(name="Compound Finder", button_type="default", disabled=True)
page13_button = pn.widgets.Button(name="Report Generation", button_type="default", disabled=True)
RESET_button = pn.widgets.Button(name="RESET", button_type="danger", disabled=False)

# Set up button click callbacks
index_button.on_click(lambda event: show_page(pages["Index"]))
instruction_button.on_click(lambda event: show_page(pages["Instructions"]))
page1_button.on_click(lambda event: show_page(pages["Data Reading"]))
page1_1_button.on_click(lambda event: show_page(pages["Data Metadata"]))
page1_2_button.on_click(lambda event: show_page(pages["Data Filtering"]))
page2_button.on_click(lambda event: show_page(pages["Data Annotation"]))
page2_1_button.on_click(lambda event: show_page(pages["Annotation De-Duplication"]))
page3_button.on_click(lambda event: show_page(pages["Data Pre-Treatment"]))
page4_button.on_click(lambda event: show_page(pages["Class Colours"]))
page5_button.on_click(lambda event: show_page(pages["Common and Exclusive Compounds"]))
page6_button.on_click(lambda event: show_page(pages["Unsupervised Analysis"]))
page7_button.on_click(lambda event: show_page(pages["Supervised Analysis"]))
page8_button.on_click(lambda event: show_page(pages["Univariate Analysis"]))
page9_button.on_click(lambda event: show_page(pages["Data Visualization"]))
page10_button.on_click(lambda event: show_page(pages["Pathway Assignment"]))
page11_button.on_click(lambda event: show_page(pages["BinSim Analysis"]))
page12_button.on_click(lambda event: show_page(pages["Compound Finder"]))
page13_button.on_click(lambda event: show_page(pages["Report Generation"]))
RESET_button.on_click(RESET)


sidebar = pn.Column(index_button, instruction_button, '## Data Pre-Processing and Pre-Treatment', page1_button, page1_1_button, page1_2_button, page2_button, page2_1_button,
                    page3_button, page4_button, '## Statistical Analysis', page5_button, page6_button, page7_button, page8_button, page9_button, page10_button, page11_button,
                    page12_button, '## Report Generation', page13_button, '## To Reset', RESET_button)


app = pn.template.BootstrapTemplate(title='Testing MetsTA', sidebar=[sidebar], main=[main_area])

app.show(websocket_max_message_size=100*1024*1014, http_server_kwargs={'max_buffer_size': 100*1024*1014})