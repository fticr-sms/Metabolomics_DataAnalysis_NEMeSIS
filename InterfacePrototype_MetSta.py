
## Needed imports
import pandas as pd
import numpy as np
import panel as pn
import param
import random
import seaborn as sns

# metanalysis_standard.py file
import metanalysis_standard as metsta

# The initial pages, especially the read file one does not have the nomenclature that I started using later on
# for the different widgets as well as organization
pn.extension(notifications=True)
pn.config.sizing_mode="stretch_width"

# Define pages as classes
# Initial Pages class building
class OpeningPage:
    def __init__(self):
        self.content = pn.Column("# Welcome to MetSta!",
                                 "# The Go-To place for your extreme-resolution metabolomics data analysis need.",
    pn.pane.Image('Picture_Test.png'))

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




# Page 1 - Reading File
# TODO: Make Reset button to read other datasets.

# Page 1 - Reading File
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


# Widgets and reacting functions of page 1
filename = pn.widgets.FileInput(name='Choose file', accept='.csv,.xlsx,.xls')

confirm_button_filename = pn.widgets.Button(name='Read File', button_type='primary', disabled=True)
tooltip_file = pn.widgets.TooltipIcon(value="Provided file must come from MetaboScape.")
file_extension = pn.widgets.RadioBoxGroup(name='File type', options=['CSV','Excel'], inline=False)
dataframe_to_show = pn.widgets.DataFrame(pd.DataFrame(), name='Data')

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

# Enanling button for next step
def _update_confirm_step1(event):
    confirm_button_step1.disabled = False

# Confirm file, show next page, disable reading files, update columns of the dataset read
def _confirm_step1(event):
    page1_1_button.disabled = False
    confirm_button_filename.disabled = True
    filename.disabled = True
    file_extension.disabled = True
    checkbox_formula.options = list(dataframe_to_show.value.columns)
    checkbox_annotation.options = list(dataframe_to_show.value.columns)
    radiobox_neutral_mass.options = ['None'] + list(dataframe_to_show.value.columns)
    checkbox_others.options = list(dataframe_to_show.value.columns)
    checkbox_samples.options = list(dataframe_to_show.value.columns)
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
    cols = list(dataframe_to_show.value.columns)
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
    initial_filtering(dataframe_to_show.value, 
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

# Widgets fot annotation part
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
                annotated_df.value.at[annotated_df.value.index[a], match_count_col] = 0

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
    main_area.clear()
    show_page(pages["Data Pre-Treatment"])
confirm_button_next_step_3.on_click(_confirm_button_next_step_3)




#### This marks the separation to use mainly param instead of only panel
#### I believe most areas before do not need param, except for the TODO marked in the next part

# Param Class to store all DataFrame
# TODO: Put previous dataframes - dataframe_to_show, filtered_df and annotated_df here as well and put this class at the begginning of program

# Contains before treatment data, treated_data, processed_data, univariate_data, meta_data, bin_data
class DataFrame_Storage(param.Parameterized):

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

# Page 3 - Data Pre-Treatment
# TODO: There is currently an error when running the pre-treatment due to the layout, it does not seem to affect functionality though
# TODO: Normalization by a Reference Feature has to be fixed in metabolinks

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
        page3[:4,2:5] = pn.Tabs(('Treated Data', DataFrame_Store.treated_df),
            ('Metadata', DataFrame_Store.metadata_df.T),
            ('BinSim Treated Data', DataFrame_Store.binsim_df), height=600, dynamic=True)
        confirm_button_next_step_4.disabled = False
        page3[5, :] = confirm_button_next_step_4

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
page3 = pn.GridSpec()
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
# Confirm colours, go to next step function and calling it
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
unsup_analysis_page = pn.Column()

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

# Create the sidebar
sidebar = pn.Column(index_button, page1_button, page1_1_button, page1_2_button, page2_button, page3_button, page4_button,
                   page5_button, page6_button, page7_button, page8_button, page9_button, page10_button,
                   page11_button)
# Create the main area and display the first page
main_area = OpeningPage().content

app = pn.template.FastListTemplate(title='Testing MetSta', sidebar=[sidebar], main=[main_area])

app.show()