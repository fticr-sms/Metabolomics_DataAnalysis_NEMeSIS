
## Needed imports
import pandas as pd
import numpy as np
import panel as pn

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
        self.params  = {'filt_method': filt_method.value, 'filt_kw': filt_kw.value,
                        'mvi_method': mvi_method.value, 'mvi_kw': mvi_kw.value,
                        'norm_method': norm_method.value}

    def view(self):
        return self.content
    
    def retrieve_params(self):
        self.params  = {'filt_method': filt_method.value, 'filt_kw': filt_kw.value,
                        'mvi_method': mvi_method.value, 'mvi_kw': mvi_kw.value,
                        'norm_method': norm_method.value}
        return self.params
    

# Page 1 - Reading File
# TODO: Substitute file reading by fileinput widget.
# TODO: Make Reset button to read other datasets.

def read_file(event):
    "Function to read the file given."
    if filename.value == '':
        dataframe_to_show.value = pd.DataFrame()
        return '', ''
    elif file_extension.value == 'CSV':
        if filename.value.endswith('.csv'):
            file = pd.read_csv(filename.value)
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
        else:
            pn.state.notifications.error('Provided file is not a csv file.')
    else:
        if filename.value.endswith('.xlsx') or filename.value.endswith('.xls'):
            file = pd.read_excel(filename.value)
            file = file.set_index('Bucket label')
        else:
            pn.state.notifications.error('Provided file is not an Excel file.')
    
    dataframe_to_show.value = file


# Widgets and reacting functions of page 1
filename = pn.widgets.TextInput(name='Filename', placeholder='Dataset.csv', max_length=100)

confirm_button_filename = pn.widgets.Button(name='Read File', button_type='primary', disabled=True)
tooltip_file = pn.widgets.TooltipIcon(value="Provided file must come from MetaboScape.")
file_extension = pn.widgets.RadioBoxGroup(name='File type', options=['CSV','Excel'], inline=False)
#file_input = pn.widgets.FileInput(name='Choose file', accept='.csv,.xlsx,.xls') # Substitute read_filename and file_extension in the future.
dataframe_to_show = pn.widgets.DataFrame(pd.DataFrame(), name='Data')

# Update button so it can be pressed after you put something in the filename
@pn.depends(filename.param.value, watch=True)
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
section1page = pn.GridSpec()
section1page[0:4,   0  ] = pn.Column(filename, pn.Row(file_extension, tooltip_file), confirm_button_filename)
section1page[0:4,   1:5] = dataframe_to_show
section1page[5,   :] = confirm_button_step1




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

# Change the IntSlider limits based on method chosen and dataset
#@pn.depends(filt_method.param.value, watch=True)
#def _update_limit_filt(filt_method):
#    filt_kw.end = limits_filt[filt_method]

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
# TODO: Perform annotation (funny, taking into account the name of the page)
# TODO: Call notification error when database reading does not go well - highlight which filled space 
# caused it to not go well (should be on the easy side)

# Read the database function
def read_database(filename, abv, ID_col, name_col, formula_col):
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

# Widgets for selecting number of databases
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
                confirm_button_databases_read.disabled = False

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
    titles = pn.Row()
    dbs_arrangement = pn.Row()
    for i in range(n_databases.value):
        dbs_arrangement.append(dbs_arrangement_all[i])
        titles.append(f'#### Database {i+1}')
    
    # Keep the page layout organized
    if len(page2) == 1:
        page2.append(titles)
        page2.append(pn.Row(dbs_arrangement))
    else:
        page2[1] = pn.Row(titles)
        page2[2] = pn.Row(dbs_arrangement)

confirm_button_n_databases.on_click(_confirm_button_n_databases)

# Initial page layout
page2 = pn.Column(pn.Row(n_databases, tooltip_n_databases, confirm_button_n_databases))
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
    
    if len(page2) == 4:
        page2.append(annotation_param_selection)
confirm_button_databases_read.on_click(_confirm_button_databases_read)




# Page 3 - Data Pre-Treatment (very incomplete)
# TODO: Most everything

# Missing Value Imputations
mvi_method = pn.widgets.Select(name="Missing Value Imputation Method", value="min_sample", 
                                       options=['min_sample', 'min_feat', 'min_data', 'zero'])
mvi_kw = pn.widgets.FloatSlider(name="Missing Value Imputation Keyword", value=0.2, start=0, end=1, step=0.01)
limits_mvi = {"min_sample": 1, "min_feat": 1, "min_data": 1, "zero": 0}
# Update the keyword slider based on methodology chosen
@pn.depends(mvi_method.param.value, watch=True)
def _update_limit_mvi(mvi_method):
    mvi_kw.end = limits_mvi[mvi_method]
    mvi_kw.value = 0.2

# Normalization
norm_method = pn.widgets.Select(name="Normalization Method", value="ref_feat", 
                                       options=['ref_feat', 'total_sum', 'PQN', 'Quantile', None])
#norm_kw = pn.widgets.FloatSlider(name="Missing Value Imputation Keyword", value=0.2, start=0, end=1, step=0.01)

# Layout of the INCOMPLETE page
page3 = pn.Column(mvi_method, mvi_kw, norm_method)



# Overall layout of the program and initialization
# The pages we have
pages = {
    "Index": OpeningPage(),
    "Data Reading": DataReading(),
    "Data Metadata": DataMetadata(),
    "Data Filtering": DataFiltering(),
    "Data Annotation": DataAnnotation(),
    "Data Pre-Treatment": DataPreTreatment()
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

# Set up button click callbacks
page1_button.on_click(lambda event: show_page(pages["Data Reading"]))
page1_1_button.on_click(lambda event: show_page(pages["Data Metadata"]))
page1_2_button.on_click(lambda event: show_page(pages["Data Filtering"]))
page2_button.on_click(lambda event: show_page(pages["Data Annotation"]))
page3_button.on_click(lambda event: show_page(pages["Data Pre-Treatment"]))

# Create the sidebar
sidebar = pn.Column(index_button, page1_button, page1_1_button, page1_2_button, page2_button, page3_button)

# Create the main area and display the first page
main_area = OpeningPage().content

app = pn.template.FastListTemplate(title='Testing MetSta', sidebar=[sidebar], main=[main_area])

app.show()