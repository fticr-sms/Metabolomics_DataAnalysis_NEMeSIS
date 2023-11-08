
## Needed imports
import pandas as pd
import numpy as np
import panel as pn
import param
import random
import seaborn as sns
import holoviews as hv
import plotly.express as px
import plotly.graph_objects as go
from aux_metsta_functions import plot_confidence_ellipse
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hier
import matplotlib.pyplot as plt
import matplotlib as mpl
import interface_aux_functions as iaf

# metanalysis_standard.py file
import metanalysis_standard as metsta

# The initial pages, especially the read file one does not have the nomenclature that I started using later on
# for the different widgets as well as organization
pn.extension('plotly', 'floatpanel', notifications=True)
hv.extension('plotly')
pn.config.sizing_mode="stretch_width"

# Define pages as classes
# Initial Pages class building
class OpeningPage:
    def __init__(self):
        self.content = pn.Column("# Welcome to MetSta!",
                                 "# The Go-To place for your extreme-resolution metabolomics data analysis need.",
                                pn.pane.Image('Picture_Test.png'))
    def view(self):
        return self.content

class DataReading:
    def __init__(self):
        self.content = pn.Column("# Section 1: Data Input", "Inputting your Excel or csv MetaboScape file.",
                                section1page)

    def view(self):
        return self.content

class DataMetadata:
    def __init__(self):
        self.content = pn.Column("# Section 1.1: Selecting MetaData",
                                page1_1)

    def view(self):
        return self.content
    
class DataFiltering:
    def __init__(self):
        self.content = pn.Column("# Section 1.2: Selecting Data Filtering Method",
                                page1_2)

    def view(self):
        return self.content


# Define pages as classes
class DataAnnotation:
    def __init__(self):
        self.content = pn.Column("# Section 2: Data Annotation", 
    """Perform Annotations based on available databases. 
    You can annotate with multiple databases. However, each database is annotated individually.
    Annotation works by assigning to a m/z peak / feature all metabolites of a database that are within the provided error margin.
    Annotation from two different databases might annotate different metabolites for the same m/z peak / feature.
    Thus, each database annotation will generate 3 columns added to the metadata: one with the IDs of the metabolites annotated, another with their formula and another with their name.""",
                                page2)

    def view(self):
        return self.content

# Define pages as classes
class DataPreTreatment:
    def __init__(self):
        
        self.content = pn.Column("# Section 3: Data Pre-Treatment", "Choose the pre-treatment to apply to the data.",
                                page3)

    def view(self):
        return self.content

# Define pages as classes
class ClassColours:
    def __init__(self):

        self.content = pn.Column("# Section 3.1: Select Colours for each Class", "To maintain consistency in the plots.",
                                 "### Choose the colors for each class", page4)

    def view(self):
        return self.content

# Define pages as classes
class TransitionalPage:
    def __init__(self):

        self.content = pn.Column("# Select which Analysis you would like to do:", transitional_page)

    def view(self):
        return self.content

class CommonExclusivePage:
    def __init__(self):

        self.content = pn.Column("# Seeing Common and Exclusive Compounds Between Biological Classes",
                                 comexc_page)

    def view(self):
        return self.content

class UnsupervisedAnalysisPage:
    def __init__(self):

        self.content = pn.Column("# Performing Unsupervised Analysis", "Including PCA and Hierarchical Clustering",
                                 "##### Known problem: Changing HCA plot characteristics many many times can lead to memory issues.",
                                 unsup_analysis_page)

    def view(self):
        return self.content

class SupervisedAnalysisPage:
    def __init__(self):

        self.content = pn.Column("# Performing Supervised Analysis", "Including Random Forest and PLS-DA",
                                 sup_analysis_page)

    def view(self):
        return self.content

class UnivariateAnalysisPage:
    def __init__(self):

        self.content = pn.Column("# Performing Univariate Analysis", "Including Univariate and Fold-Change Analysis",
                                 univar_analysis_page)

    def view(self):
        return self.content

class DataVisualizationPage:
    def __init__(self):

        self.content = pn.Column("# Plots to Visualize Your Data", "Van Krevelen and Chemical Composition Series",
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

# Create the main area and display the first page
main_area = OpeningPage().content



# Param Class to store all DataFrames
# TODO: Put previous dataframes - filtered_df and annotated_df here as well

# Contains before treatment data, treated_data, processed_data, univariate_data, meta_data, bin_data
class DataFrame_Storage(param.Parameterized):

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
        return pd.concat((MS_df, annot_df), axis=1)

# Initializing the Store
DataFrame_Store = DataFrame_Storage()


# Page 1 - Reading File
# TODO: Make Reset button to read other datasets.

# Page 1 - Reading File

# Widgets and reacting functions of page 1
filename = pn.widgets.FileInput(name='Choose file', accept='.csv,.xlsx,.xls')

confirm_button_filename = pn.widgets.Button(name='Read File', button_type='primary', disabled=True)
tooltip_file = pn.widgets.TooltipIcon(value="Provided file must come from MetaboScape.")
dataframe_to_show = pn.widgets.DataFrame(pd.DataFrame(), name='Data')

def read_file(event):
    "Function to read the file given."
    if filename.value == '':
        dataframe_to_show.value = pd.DataFrame()
    elif filename.filename.endswith('.csv'):
        file = pd.read_csv(filename.filename)
        file.insert(1, 'Neutral Mass', file['Bucket label'].str.replace('Da', '').astype('float'))
        #important for database match
        file = file.set_index('Bucket label')

        # Replaces zeros with numpy nans. Essential for data processing
        file = file.replace({0:np.nan})

        # Samples names frequently have 00000. Use this code to make them 'cleaner'.
        # If not needed just skip this part (use #)
        def renamer(colname):
            # Util to optionally remove all those 00000 from sample names
            return ''.join(colname.split('00000'))
        file.columns = [renamer(i) for i in file.columns]
        dataframe_to_show.value = file

    elif filename.filename.endswith('.xlsx') or filename.filename.endswith('.xls'):
        file = pd.read_excel(filename.filename)
        file = file.set_index('Bucket label')
        dataframe_to_show.value = file

    else:
        pn.state.notifications.error('Provided file is not an Excel or a csv file.')


# Update button so it can be pressed after you put something in the filename
@pn.depends(filename.param.filename, watch=True)
def _update_confirm_button_filename(filename):
    if filename != '':
        confirm_button_filename.disabled = False
    else:
        confirm_button_filename.disabled = True

# Function happens when you press the button        
confirm_button_filename.on_click(read_file)

img_confirm_button = '''<svg xmlns="http://www.w3.org/2000/svg" 
    class="icon icon-tabler icon-tabler-check" width="24" height="24" viewBox="0 0 24 24" stroke-width="2"
    stroke="currentColor" fill="none" stroke-linecap="round" stroke-linejoin="round">
    <path stroke="none" d="M0 0h24v24H0z" fill="none"></path>
    <path d="M5 12l5 5l10 -10"></path>
    </svg>'''

# Enabling button for next step
def _update_confirm_step1(event):
    confirm_button_step1.disabled = False

# Confirm file, show next page, disable reading files, update columns of the dataset read
def _confirm_step1(event):
    page1_1_button.disabled = False
    confirm_button_filename.disabled = True
    filename.disabled = True
    DataFrame_Store.read_df = dataframe_to_show.value
    checkbox_formula.options = list(DataFrame_Store.read_df.columns)
    checkbox_annotation.options = list(DataFrame_Store.read_df.columns)
    radiobox_neutral_mass.options = ['None'] + list(DataFrame_Store.read_df.columns)
    checkbox_others.options = list(DataFrame_Store.read_df.columns)
    checkbox_samples.options = list(DataFrame_Store.read_df.columns)
    main_area.clear()
    show_page(pages["Data Metadata"])
    # reset button

# Make button and call the appropriate functions when the buttons are pressed
confirm_button_step1 = pn.widgets.Button(icon=img_confirm_button, name='Confirm - Next Step', button_type='success',
                                         disabled=True)
confirm_button_filename.on_click(_update_confirm_step1)
confirm_button_step1.on_click(_confirm_step1)

# Setting up the page layout
section1page = pn.Column(filename, pn.Row(confirm_button_filename, tooltip_file),
                        dataframe_to_show, confirm_button_step1)




# Page 1-1 - Selecting sample columns and Target

# Making checkbox widgets for each category
checkbox_formula = pn.widgets.CheckBoxGroup(
    name='Formula', value=['Formula'], options=list(dataframe_to_show.value.columns),
    inline=False, disabled=False)

checkbox_annotation = pn.widgets.CheckBoxGroup(
    name='Annotation', value=['Name'], options=list(dataframe_to_show.value.columns),
    inline=False, disabled=False)

# Select only one instead of multiple
radiobox_neutral_mass = pn.widgets.RadioBoxGroup(
    name='Neutral Mass', value='Neutral Mass', options=['None'] + list(dataframe_to_show.value.columns),
    inline=False, disabled=False)

checkbox_others = pn.widgets.CheckBoxGroup(
    name='Others', options=list(dataframe_to_show.value.columns),
    inline=False, disabled=False)

checkbox_samples = pn.widgets.CheckBoxGroup(
    name='Samples', options=list(dataframe_to_show.value.columns),
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
confirm_button_column_selection = pn.widgets.Button(icon=img_confirm_button, name='Confirm Columns', button_type='success')

# Confirm column selection, update sample columns, make target editable while providing a possible target
def _update_confirm_column_selection(event):
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
    # Attempt to make 'an automatic target'
    if checkbox_samples.value != '':
        tg = [i.split('_')[0] for i in checkbox_samples.value]
        target_widget.placeholder = ','.join(tg)
        target_widget.value = ','.join(tg)
    else:
        target_widget.placeholder = target_placeholder
    
confirm_button_column_selection.on_click(_update_confirm_column_selection)

# Make the target widget
target_placeholder = "A,A,A,A,A,B,B,B,B,B"
target_widget = pn.widgets.TextAreaInput(name='Target - Class Labels', placeholder=target_placeholder, 
                                         max_length=5000, height=100,
                                        disabled=True)
target_tooltip = pn.widgets.TooltipIcon(value="Provide the class labels of your samples. No spaces between labels.")

# Make button be pressable when you have a target
@pn.depends(target_widget.param.value, watch=True)
def _update_read_target_button(target_widget):
    if target_widget != '':
        confirm_button_target.disabled = False
    else:
        confirm_button_target.disabled = True

confirm_button_target = pn.widgets.Button(icon=img_confirm_button, name='Confirm Target', button_type='success',
                                         disabled=True)
confirm_button_next_step_1_1 = pn.widgets.Button(icon=img_confirm_button,
                                                 name='Next Step - Data Filtering and Characteristics',
                                                 button_type='success', disabled=True)

# Make sure that target makes sense in comparison to the number of sample columns
def _update_confirm_target(event):
    target = target_widget.value.split(',')
    sample_cols = checkbox_samples.value
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
filt_method = pn.widgets.Select(name="Feature Filter Method", value="total_samples", 
                                       options=['total_samples', 'class_samples', None])
filt_kw = pn.widgets.IntSlider(name="Feature Filter Keyword", value=2, start=0, end=len(target_widget.value.split(',')))
limits_filt = {"total_samples": len(target_widget.value.split(',')), 
               "class_samples": pd.Series(target_widget.value.split(',')).value_counts().min(), None: 1}
# Change the IntSlider limits based on the number of class labels
@pn.depends(target = target_widget.param.value,
            method = filt_method.param.value, watch=True)
def _update_filt_kw_limits(target, method):
    if method == 'total_samples':
        filt_kw.end = len(target.split(','))
    elif method == 'class_samples':
        filt_kw.end = pd.Series(target.split(',')).value_counts().min()
    else:
        filt_kw.end = 1

filt_method_tooltip = pn.widgets.TooltipIcon(value=
    """'total_samples' requires a feature to appear in at least x samples in the whole dataset to be retained.
    'class_samples' requires a feature to appear in at least x samples of at least 1 class to be retained.""")
filt_kw_tooltip = pn.widgets.TooltipIcon(value=
    """How many samples a feature has to appear to be retained based on the method chosen before.""")

# Preparing DataFrames
filtered_df = pn.widgets.DataFrame(pd.DataFrame(), name='Filtered DataFrame')
characteristics_df = pn.widgets.DataFrame(pd.DataFrame(), name='Characteristics DataFrame')

# Button to perform filtering
confirm_button_initial_filtering = pn.widgets.Button(icon=img_confirm_button, name='Perform Filtering',
                                                     button_type='success', disabled=False)

# Button to next step
confirm_button_next_step_2 = pn.widgets.Button(icon=img_confirm_button, name='Next Step - Annotation',
                                                     button_type='success', disabled=False)

# Perform filtering
def initial_filtering(df, sample_cols, target=None, filt_method='total_samples', filt_kw=2):
    meta_cols = [i for i in df.columns if i not in sample_cols]
    temp = metsta.basic_feat_filtering(df[sample_cols].T, target=target,
                                filt_method=filt_method, # Method
                                filt_kw=filt_kw) # Make a sample have to appear
    df = pd.concat((df[meta_cols].reindex(temp.columns), temp.T), axis=1)
    data_characteristics = metsta.characterize_data(df[sample_cols].T, target=target)
    
    filtered_df.value = df
    characteristics_df.value = pd.DataFrame(pd.Series(data_characteristics)).iloc[1:]

# Call filtering function and extend the page with data characteristics and results of filtering
def _confirm_button_initial_filtering(event):
    sample_cols = checkbox_samples.value
    target = target_widget.value.split(',')
    initial_filtering(DataFrame_Store.read_df,
                      sample_cols, target=target, filt_method=filt_method.value, filt_kw=filt_kw.value)
    annotated_df.value = pd.DataFrame(index=filtered_df.value.index)
    if len(page1_2) == 3:
        page1_2.extend(['#### Characteristics of the Dataset',characteristics_df,'#### Filtered Dataset',filtered_df,
                       confirm_button_next_step_2])

# Call the function
confirm_button_initial_filtering.on_click(_confirm_button_initial_filtering)

# Go to next step function and calling it
def _confirm_button_next_step_1_2(event):
    page2_button.disabled = False
    main_area.clear()
    show_page(pages["Data Annotation"])
confirm_button_next_step_2.on_click(_confirm_button_next_step_1_2)

# Initial page layout
page1_2 = pn.Column(pn.Row(filt_method, filt_method_tooltip), pn.Row(filt_kw, filt_kw_tooltip),
                    confirm_button_initial_filtering)




# Page 2 - Annotation of Metabolites
# TODO: Call notification error when database reading does not go well - highlight which filled space 
# caused it to not go well (should be on the easy side)

# Read the database function
def read_database(filename, abv, ID_col, name_col, formula_col):
    try:
        if filename.endswith('.csv'):
            db = pd.read_csv(filename).set_index(ID_col)
            if abv == 'HMDB':
                db[name_col] = db[name_col].str.replace("b'", "")
                db[name_col] = db[name_col].str.replace("'", "")
        elif filename.endswith('.xlsx'):
            db = pd.read_excel(filename).set_index(ID_col)
            if abv == 'HMDB':
                db[name_col] = db[name_col].str.replace("b'", "")
                db[name_col] = db[name_col].str.replace("'", "")
        else:
            raise ValueError('File Format not accepted. Only csv and xlsx files are accepted.')
        ##
        db['Mass'] = db[formula_col].dropna().apply(metsta.calculate_monoisotopic_mass)
        return db
    except:
        pn.state.notifications.error('Database could not be read. Check if all parameters given are correct.',
                                     duration=5000)

# Widgets for selecting number of databases
n_databases_show = pn.widgets.IntInput(name='Nº of Databases to annotate', value=1, step=1, start=0, end=5)
n_databases = pn.widgets.IntInput(name='Nº of Databases to annotate', value=1, step=1, start=0, end=5)
tooltip_n_databases = pn.widgets.TooltipIcon(value="Select how many (0-5) databases you want to use for annotation.")
# Button to perform filtering
confirm_button_n_databases = pn.widgets.Button(icon=img_confirm_button, name='Select Databases',
                                                     button_type='success', disabled=False)

# Class to organize how a single Database section will be shown and that can be repeated
class DatabaseSection():
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
        
        db = pn.widgets.DataFrame(pd.DataFrame(), name='Database')
        
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
                db.value = read_database(filename, abv, ID_col, name_col, formula_col)
                confirm_button_db.param.clicks = 0
                self.db = db
                if len(self.content) == 6:
                    self.content.append(f'Database {self.abv} has {len(self.db.value)} metabolites.')
                else:
                    self.content[6] = f'Database {self.abv} has {len(self.db.value)} metabolites.'
                self.read.value = True
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

    if n_databases.value == 0:
        page2.append(confirm_button_databases_read)
        confirm_button_databases_read.disabled = False

confirm_button_n_databases.on_click(_confirm_button_n_databases)

# Initial page layout
page2 = pn.Column(pn.Row(n_databases_show, tooltip_n_databases, confirm_button_n_databases))
confirm_button_databases_read = pn.widgets.Button(icon=img_confirm_button, name='Confirm Databases',
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
        
annotation_param_selection = pn.Column(annotation_margin_method_radio, 
                                       pn.Row(annotation_ppm_deviation, tooltip_annotation),
                                       confirm_button_annotation_perform)

# Make the annotation part of the layout appear and disable button to confirm databases
def _confirm_button_databases_read(event):
    confirm_button_databases_read.disabled = True
    confirm_button_annotation_perform.disabled = False
    annotated_df.value = pd.DataFrame(index=filtered_df.value.index)
    
    page2.append(annotation_param_selection)
confirm_button_databases_read.on_click(_confirm_button_databases_read)

# Widgets for annotation part
performing_annotation_arrangement = pn.Column() # Start with empty page for the annotation widgets
tqdm_database = {i+1:pn.widgets.Tqdm() for i in range(n_databases.end)}
verbose_annotated_compounds = {i+1:pn.widgets.StaticText(name='', value=f'') for i in range(n_databases.end)}
annotated_df = pn.widgets.DataFrame(pd.DataFrame(index=filtered_df.value.index))

# Function to perform metabolite annotation (also contributes to updating the page)
def metabolite_annotation():
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
                    radiobox_neutral_mass.value][a])/filtered_df.value[radiobox_neutral_mass.value][a])*10**6
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
    while len(performing_annotation_arrangement)>0:
        performing_annotation_arrangement.pop(-1)
    page2.append(performing_annotation_arrangement)
    confirm_button_annotation_perform.disabled=True
    metabolite_annotation()
    page2.append(confirm_button_next_step_3)

confirm_button_annotation_perform.on_click(_press_confirm_annotation_perform)

# Button to next step
confirm_button_next_step_3 = pn.widgets.Button(icon=img_confirm_button, name='Next Step - Data Pre-Treatment',
                                                     button_type='success', disabled=False)

# Go to next step function and calling it
def _confirm_button_next_step_3(event):
    page3_button.disabled = False
    DataFrame_Store.original_df = DataFrame_Store.concat_annots(filtered_df.value, annotated_df.value)
    creating_has_match_column()
    main_area.clear()
    show_page(pages["Data Pre-Treatment"])
confirm_button_next_step_3.on_click(_confirm_button_next_step_3)

def creating_has_match_column():
    "Creating the has match column to be used later (common/exclusive compounds)."
    if n_databases.value == 0:
        for i in DataFrame_Store.original_df.index:
            df = DataFrame_Store.original_df.loc[[i]]
            hasmatch = df[checkbox_annotation.value].notnull().values.any()
            DataFrame_Store.original_df.at[i, 'Has Match?'] = hasmatch

    else:
        cols_to_see = checkbox_annotation.value + list(DataFrame_Store.original_df.columns[-n_databases.value*4:])

        DataFrame_Store.original_df['Has Match?'] = np.nan
        for i in DataFrame_Store.original_df.index:
            df = DataFrame_Store.original_df.loc[[i]]
            hasmatch = df[cols_to_see].notnull().values.any()
            DataFrame_Store.original_df.at[i, 'Has Match?'] = hasmatch




#### This marks the separation to use mainly param instead of only panel

# Page 3 - Data Pre-Treatment
# TODO: There is currently an error when running the pre-treatment due to the layout, it does not seem to affect functionality though

# Param to encompass the choice of Pre-Treatment methods
# TODO: Introduce way to read reference sample for PQN normalization
class PreTreatment(param.Parameterized):
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
    # TODO: Make Metadata appear in a more clean way
    def _confirm_button_press(self, event):
        performing_pretreatment(self, DataFrame_Store, target_widget, checkbox_samples)
        page3[:5,2:5] = pn.Tabs(('Treated Data', DataFrame_Store.treated_df),
            ('Metadata', DataFrame_Store.metadata_df.T),
            ('BinSim Treated Data', DataFrame_Store.binsim_df), height=600, dynamic=True)
        confirm_button_next_step_4.disabled = False
        #page3[5, :] = confirm_button_next_step_4

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
            'norm_kw': pn.widgets.MultiChoice(name="Normalization Keyword", max_items=1,
                                              option_limit=8, search_option_limit=8, disabled=True,
                        description="""If you do not find your reference feature (Norm. by Ref. Feat.) or your sample (if you are choosing a sample for PQN), start writing it in the box below until it appears."""),

            'tf_method': pn.widgets.Select(name="Transformation Method",
                                           value = 'Generalized Logarithmic Transformation (glog)',
                            options=['Generalized Logarithmic Transformation (glog)', None]),
            'tf_kw': pn.widgets.FloatInput(name="Transformation Keyword", value=None,
                    description="""Lambda value in glog transformation. Leave empty or 0 for usual log transformation.
                    Ignored if Transformation Method is None.
            glog equation is: log\N{SUBSCRIPT TWO}((y + \u221A(y\N{SUPERSCRIPT TWO} + lambda\N{SUPERSCRIPT TWO}) ) /2)
                    """),

            'scaling_method': pn.widgets.Select(name="Scaling Method", value = 'Pareto Scaling',
                            options=['Pareto Scaling', 'Auto Scaling', 'Mean Centering', 'Range Scaling',
                                     'Vast Scaling', 'Level Scaling', None]),
            'scaling_kw': pn.widgets.Select(name="Scaling Keyword for Level Scaling only", value='', disabled=True,
                                           options=['Average', 'Median']),

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
        PreTreatment_Method.controls.widgets['norm_kw'].options = options_norm[norm_method] + list(checkbox_samples.value)
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

def performing_pretreatment(PreTreatment_Method, DataFrame_Store, target_widget, checkbox_samples):
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

    # Data to pass, target and sample columns
    data = DataFrame_Store.original_df
    target = target_widget.value.split(',')
    sample_cols = checkbox_samples.value

    a,b,c,d,e = metsta.filtering_pretreatment(data, target, sample_cols,
                      filt_method=None, filt_kw=2, # Filtering based on number of times features appear
                      extra_filt=None, # Filtering based on annotation of features ('Formula' or 'Name')
                      mvi=mvi, mvi_kw=mvi_kw, # Missing value imputation
                      norm=norm, norm_kw=norm_kw, # Normalization
                      tf=tf, tf_kw=tf_kw, # Transformation
                      scaling=scaling, scaling_kw=scaling_kw) # Scaling

    DataFrame_Store.treated_df, DataFrame_Store.processed_df, DataFrame_Store.univariate_df = a, b, c
    DataFrame_Store.metadata_df, DataFrame_Store.binsim_df = d, e

# Button to next step
confirm_button_next_step_4 = pn.widgets.Button(icon=img_confirm_button, name='Next Step - Class Colours',
                                                     button_type='success', disabled=True)

# Go to next step function and calling it
def _confirm_button_next_step_4(event):
    page4_button.disabled = False

    # Filling the target storage with the correct target and default colours
    target_list.target = target_widget.value.split(',')
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

    main_area.clear()
    show_page(pages["Class Colours"])
confirm_button_next_step_4.on_click(_confirm_button_next_step_4)

# Organize Page Layout
page3 = pn.GridSpec(mode='override')
page3[:5,0:2] = PreTreatment_Method.controls
page3[:5,2:5] = pn.Tabs(('Treated Data', DataFrame_Store.treated_df),
        ('Metadata', DataFrame_Store.metadata_df),
       ('BinSim Treated Data', DataFrame_Store.binsim_df), height=600, dynamic=True)
page3[5, :] = confirm_button_next_step_4




# Page 4 - Choosing Class Colours

# Default colours - 10 possible colours
colours = sns.color_palette('tab10', 10)

# Function to fill the target and colours of the above class
def TargetStorage_filling(target_list, colours):

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

# From Stack Overflow
def RGB(col): return '#%02x%02x%02x' % (int(col[0]), int(col[1]), int(col[2]))

# Class to store target and to store target colours
class TargetStorage(param.Parameterized):

    # Function to do when calling the function
    updater = param.Callable(TargetStorage_filling)

    # Setting up parameters
    target = param.List(default=target_widget.value.split(','))
    color_classes = param.Dict(default={})

    def __init__(self, **params):
        super().__init__(**params)

    def __call__(self, target, colours):
        self.color_classes = self.updater(target, colours)

        return self

target_list = TargetStorage()

# Setting up empty page
page4 = pn.Column()

# Button to next step and to confirm colours
confirm_button_next_step_transitionalpage = pn.widgets.Button(icon=img_confirm_button, name='Next Step - Analysis',
                                                     button_type='success', disabled=False)
# Confirm colours, go to next step function and calling it, enabling all buttons for analysis and performing initial computations for each analysis
def _confirm_button_next_step_5(event):
    n_classes = len(target_list.color_classes)
    for row in range(0, n_classes, 5):
        for col in range(5):
            key = list(target_list.color_classes.keys())[row+col]
            target_list.color_classes[key] = page4[row//5][col].value
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
    PCA_params.PCA_plot[0] = _plot_PCA()
    middle_page_PCA[0,1:3] = pn.pane.Plotly(PCA_params.PCA_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_plot',}})
    PCA_params.exp_var_fig_plot[0] = _plot_PCA_explained_variance()
    end_page_PCA[0] = pn.pane.Plotly(PCA_params.exp_var_fig_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_exp_var_plot',}})
    PCA_params.scatter_PCA_plot[0] = _scatter_PCA_plot()
    end_page_PCA[1] = pn.pane.Plotly(PCA_params.scatter_PCA_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_scatter_plot',}})

    # Initial calculations for HCA and storing initial plots
    HCA_params.dists = dist.pdist(DataFrame_Store.treated_df, metric=HCA_params.dist_metric)
    HCA_params.Z = hier.linkage(HCA_params.dists, method=HCA_params.link_metric)
    HCA_params.HCA_plot[0] = _plot_HCA()
    page_HCA[0:6,1:4] = HCA_params.HCA_plot[0]

    main_area.clear()
    show_page(pages["Transitional Page"])
confirm_button_next_step_transitionalpage.on_click(_confirm_button_next_step_5)




ComExc_A = pn.widgets.Button(name='Common/Exclusive Comp.', button_type='primary')
Unsup_A = pn.widgets.Button(name='Unsupervised Analysis', button_type='default')
Sup_A = pn.widgets.Button(name='Supervised Analysis', button_type='success')
Univariate_A = pn.widgets.Button(name='Univariate Analysis', button_type='warning')
DataViz_A = pn.widgets.Button(name='Data Visualization', button_type='danger')

pn.Row(ComExc_A, Unsup_A, Sup_A, Univariate_A, DataViz_A)

BinSim_A = pn.widgets.Button(name='BinSim Specific Analysis', button_type='success')
CompFinder_A = pn.widgets.Button(name='Compound Finder', button_type='warning')
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

transitional_page = pn.Column(pn.Row(ComExc_A, Unsup_A, Sup_A, Univariate_A, DataViz_A),
                             '### Other Options:',
                             pn.Row(BinSim_A, CompFinder_A, ToBeAdded_A))




# Page for Common and Exclusive Compounds
comexc_page = pn.Column()

# Page for Unsupervised Analysis

# Tab for PCA analysis
# TODO: PCA is straight away computed with 10 components. If your data has less than 10 features, this will lead to an error and everything will fail.
# Should this possibility be taken into account?

# Param Class to store parameters and data regarding PCA
class PCA_Storage(param.Parameterized):
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

    # Computed PCA df
    pca_scores = param.DataFrame()
    explained_variance = param.Array()
    pca_loadings = param.DataFrame()

    # Storing figures
    PCA_plot = param.List(default=['a'])
    exp_var_fig_plot = param.List(default=['a'])
    scatter_PCA_plot = param.List(default=['a'])

    # Update the PCA plot
    @param.depends('n_dimensions', 'PCx', 'PCy', 'PCz', 'ellipse_draw', 'confidence', 'confidence_std', watch=True)
    def _update_PCA_plot(self):
        self.PCA_plot[0] = _plot_PCA()
        middle_page_PCA[0,1:3] = pn.pane.Plotly(PCA_params.PCA_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_plot',}})

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
        }
        self.controls = pn.Param(self, parameters=['n_dimensions', 'PCx', 'PCy', 'PCz',
                                                   'ellipse_draw', 'confidence', 'confidence_std'],
                                 widgets=widgets, name='Parameters for PCA Projection plot')

# Running initial param to store PCA details
PCA_params = PCA_Storage()

def _plot_PCA():
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

    return final_PCA_plot

def _plot_PCA_explained_variance():
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

def _scatter_PCA_plot():
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
    PCA_params.controls.widgets['PCx'].options = ['PC '+str(i+1) for i in range(components)]
    PCA_params.controls.widgets['PCy'].options = ['PC '+str(i+1) for i in range(components)]
    PCA_params.controls.widgets['PCz'].options = ['PC '+str(i+1) for i in range(components)]

    PCA_params.exp_var_fig_plot[0] = _plot_PCA_explained_variance()
    end_page_PCA[0] = pn.pane.Plotly(PCA_params.exp_var_fig_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_exp_var_plot',}})

    if len(PCA_params.scatter_PCA_plot[0].data[0]['dimensions']) < 4:
        PCA_params.exp_var_fig_plot[0] = _plot_PCA_explained_variance()
        end_page_PCA[1] = pn.pane.Plotly(PCA_params.exp_var_fig_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_scatter_plot',}})
    else:
        if components < 4:
            PCA_params.scatter_PCA_plot[0] = _scatter_PCA_plot()
            end_page_PCA[1] = pn.pane.Plotly(PCA_params.scatter_PCA_plot[0], config = {'toImageButtonOptions': {'filename': 'PCA_scatter_plot',}})

# Function enabling/disabling the confidence level parameters for ellipses
@pn.depends(PCA_params.controls.widgets['ellipse_draw'].param.value,
            watch=True)
def _update_ellipse_options(ellipse):
    if ellipse:
        PCA_params.controls.widgets['confidence'].disabled = False
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
download_icon = '''<svg xmlns="http://www.w3.org/2000/svg" class="icon icon-tabler icon-tabler-download" width="24" height="24" viewBox="0 0 24 24" stroke-width="2" stroke="currentColor" fill="none" stroke-linecap="round" stroke-linejoin="round">
   <path stroke="none" d="M0 0h24v24H0z" fill="none"></path>
   <path d="M4 17v2a2 2 0 0 0 2 2h12a2 2 0 0 0 2 -2v-2"></path>
   <path d="M7 11l5 5l5 -5"></path>
   <path d="M12 4l0 12"></path>
</svg>'''
save_HCA_plot_button = pn.widgets.Button(name='Save as a png (in current folder)', button_type='success',
                                         icon=download_icon)
# When pressing the button, downloads the figure
def _save_HCA_plot_button(event):
    HCA_params.HCA_plot[0].savefig('HCA_plot.png', dpi=HCA_params.dpi)
save_HCA_plot_button.on_click(_save_HCA_plot_button)

# Organization for the HCA page
page_HCA = pn.GridSpec(mode='override')
page_HCA[0:5,0] = HCA_params.controls
page_HCA[0:6,1:4] = HCA_params.HCA_plot[0]
page_HCA[5,0] = save_HCA_plot_button

# Page with PCA and HCA analysis
unsup_analysis_page = pn.Tabs(('PCA', page_PCA), ('HCA', page_HCA))




# Page for Supervised Analysis
sup_analysis_page = pn.Column()

# Page for Univariate Analysis
univar_analysis_page = pn.Column()

# Page for Data Visualization
data_viz_page = pn.Column()

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

# Define buttons to navigate between pages
index_button = pn.widgets.Button(name="Index", button_type="primary", disabled=True)
page1_button = pn.widgets.Button(name="Data Reading", button_type="primary")
page1_1_button = pn.widgets.Button(name="Data Metadata", button_type="default", disabled=True)
page1_2_button = pn.widgets.Button(name="Data Filtering", button_type="default", disabled=True)
page2_button = pn.widgets.Button(name="Data Annotation", button_type="primary", disabled=True)
page3_button = pn.widgets.Button(name="Data Pre-Treatment", button_type="primary", disabled=True)
page4_button = pn.widgets.Button(name="Class Colours", button_type="primary", disabled=True)
page5_button = pn.widgets.Button(name="Common/Exclusive Comp.", button_type="primary", disabled=True)
page6_button = pn.widgets.Button(name="Unsupervised Analysis", button_type="default", disabled=True)
page7_button = pn.widgets.Button(name="Supervised Analysis", button_type="success", disabled=True)
page8_button = pn.widgets.Button(name="Univariate Analysis", button_type="warning", disabled=True)
page9_button = pn.widgets.Button(name="Data Visualization", button_type="danger", disabled=True)
page10_button = pn.widgets.Button(name="BinSim Analysis", button_type="success", disabled=True)
page11_button = pn.widgets.Button(name="Compound Finder", button_type="warning", disabled=True)
RESET_button = pn.widgets.Button(name="RESET", button_type="danger", disabled=False)

# Set up button click callbacks
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


# Create the sidebar
sidebar = pn.Column(index_button, page1_button, page1_1_button, page1_2_button, page2_button, page3_button, page4_button,
                   page5_button, page6_button, page7_button, page8_button, page9_button, page10_button,
                   page11_button, 'To Reset (TODO)', RESET_button)

app = pn.template.FastListTemplate(title='Testing MetSta', sidebar=[sidebar], main=[main_area])

app.show()