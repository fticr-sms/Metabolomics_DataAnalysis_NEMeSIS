
## Needed imports
import pandas as pd
import numpy as np
import panel as pn
import param
import holoviews as hv
from io import BytesIO
import plotly.express as px
from metabolinks import read_data_from_xcel, align
from metabolinks.peak_alignment import alignment_summary

# Activating extensions
pn.extension('plotly', 'floatpanel', 'katex', notifications=True)
hv.extension('plotly')
pn.config.sizing_mode="stretch_width"

# Important Text
alignment_description_string = """
<p>This Data Alignment software is geared towards <strong>Direct Infusion mass spectrometry data</strong>. That is,
metabolic features will be aligned based on their mass values alone, and it does not accept other parameters such as
retention time or collision cross-section.</p>
<p>To perform data alignment, you must choose an <strong>Excel</strong> file (.xlsx or .xls) file. This file formatting
should be as follow: it should have <strong>one sample per Excel sheet</strong>; the <strong>name</strong> of the Excel
sheet should correspond to the <strong>sample name</strong>; each sheet should have in its <strong>first column the mass
values</strong> (<em>m/z</em>, neutral mass or equivalent) and in its <strong>second column the&nbsp;</strong>corresponding
<strong>intensity values</strong>; and finally, the <strong>first row</strong> should have the <strong>name of the two
columns</strong>, for example, &apos;<em>m/z</em>&apos; and &apos;I&apos; (this name should be consistent between samples).
The example file &apos;example_samples_to_align.xlsx&apos; is available in the interface folder as guidance.</p>
<p>The data alignment perform is made using the <strong><em>align</em></strong> function from the Metabolinks Python package.
It uses 2 parameters: the <strong>PPM Deviation Tolerance</strong> that defines the maximum tolerance (in parts per million)
to group metabolic features from different samples together and the <strong>Minimum Sample Number Appearance</strong> that
defines the minimum number of samples a metabolic feature must appear in to be kept in the final aligned Dataset.</p>
<p>After alignment, a short description of the alignment process will be shown as well as the aligned dataset (if this
aligned dataset is large, only the first 10,000 metabolic features will be shown). The dataset and an abridged description
of the metabolic feature alignments can be saved. The aligned dataset can be directly used in our MetsTA data analysis
software.</p>"""


class AlignmentPage:
    def __init__(self):
        self.content = pn.Column("# Data Alignment Section", alignment_description_string, alignment_page)
    def view(self):
        return self.content


# Contains relevant parameters regarding Data Alignment
class DataAlignmentObject(param.Parameterized):
    "Class to contain all relevant parameters regarding Data Alignment."

    # Filename
    file = param.String()

    # Samples in the list
    samples = param.List()

    # Aligned DataFrame
    aligned_df = param.DataFrame()

    # Parameters for Alignment
    ppmtol = param.Number(default=1.0) # PPM Deviation Tolerance to make the aligned buckets
    min_samples = param.Number(default=2) # Min. number of samples an aligned feature must appear to not be discarded
    confirm_button_alignment = param.Boolean(default=False)

    # Alignment Description
    desc = param.DataFrame()
    aligned_desc = param.String(default='')


    def update_widgets(self, filename):
        "Updated min_samples widget based on number of samples in read file."

        self.file = filename
        self.controls.widgets['min_samples'].end = len(self.samples)
        if len(self.samples) < self.min_samples:
            self.controls.widgets['min_samples'].value = len(self.samples)
            self.min_samples = len(self.samples)


    def _confirm_button_alignment(self, event):
        "Performs Data Alignment and updates layout."

        # Updating layout
        while len(alignment_page)>2:
            alignment_page.pop(-1)

        # Loading Widget while de-duplication is happenning
        alignment_page.append(pn.indicators.LoadingSpinner(value=True, size=120,
                                                        name='Performing Data Alignment...'))

        try:
            # Perform Alignment
            aligned, desc = align(self.samples,
                          min_samples=1, 
                          ppmtol=self.ppmtol,
                          return_alignment_desc=True,
                          verbose=False)

            n_groups = len(aligned)

            # Get only featuress that appeared in at least min_samples samples
            desc = desc[desc['# features'] >= self.min_samples]
            aligned = aligned.iloc[desc.index]

            # Storing files
            self.aligned_df = aligned
            self.desc = desc

            # Building the description
            n_total_features = 0
            for samp in self.samples:
                n_total_features += len(samp)
            desc_list = [f'A total of **{n_total_features} features** were extracted in **{len(self.samples)} samples**. Data Alignment was performed with PPM Deviation tolerance of {self.ppmtol}.',
                         f'**{n_groups}** total groupings were found. From these, **{n_groups-len(self.aligned_df)}** groups were discarded (appeared in less than **{self.min_samples} samples**).',
                         f'The aligned dataset had the remaining **{len(self.aligned_df)} metabolic features**.',
                         '']
            summary_string = alignment_summary(self.desc, self.ppmtol).split('\n')
            desc_list.extend(summary_string)
            self.aligned_desc = '<br />'.join(desc_list)

            # Updating Layout
            alignment_page.pop(-1)
            alignment_page.append(self.aligned_desc)
            if len(self.aligned_df) < 8000:
                alignment_page.append(pn.pane.DataFrame(self.aligned_df, height=600))
            else:
                alignment_page.append(pn.pane.DataFrame(self.aligned_df.iloc[:8000], height=600))
            alignment_page.append(saving_df_widgets)

            pn.state.notifications.success(f'Data sucessfully aligned.')

        except:
            while len(alignment_page)>2:
                alignment_page.pop(-1)
            pn.state.notifications.error(f'Data could not be aligned.')


    def reset(self):
        "Resets all relevant parameters."
        for param in self.param:
            if param not in ["name"]:
                setattr(self, param, self.param[param].default)


    def __init__(self, **params):

        super().__init__(**params)

        widgets = {
            'ppmtol': pn.widgets.FloatInput(name='PPM Deviation Tolerance (for metabolic feature grouping)',
                start=0.1, end=1000, value=1, step=0.1, styles={'font-weight': 'bold'}),
            'min_samples': pn.widgets.IntSlider(name='Minimum Sample Number Appearance (for a metabolic feature)',
                start=1, value=2, step=1, end=3, styles={'font-weight': 'bold'}),
            'confirm_button_alignment': pn.widgets.Button(name="Perform Data Alignment", button_type='warning')
        }

        self.controls = pn.Param(self, parameters=['ppmtol', 'min_samples', 'confirm_button_alignment'],
                                 widgets=widgets, name='Alignment Parameters')

# Initialize Object
alignment_storage = DataAlignmentObject()

# Data Alignment Function happens when you press the corresponding button
alignment_storage.controls.widgets['confirm_button_alignment'].on_click(alignment_storage._confirm_button_alignment)


# Needed Widgets
filename = pn.widgets.FileInput(name='Choose file with samples', accept='.xlsx,.xls')

confirm_button_filename = pn.widgets.Button(name='Read Samples', button_type='primary', disabled=True)

# Update button so it can be pressed after you put something in the filename
@pn.depends(filename.param.filename, watch=True)
def _update_confirm_button_filename(filename):
    "Controls the state of the button to confirm and read the file."
    if filename != '':
        confirm_button_filename.disabled = False
    else:
        confirm_button_filename.disabled = True


def _confirm_button_filename(event):
    "Reads the samples in the provided file."

    # Updating layout
    while len(alignment_page)>1:
        alignment_page.pop(-1)
    while len(data_reading_section)>2:
        data_reading_section.pop(-1)

    try:
        # Reading the file
        data = read_data_from_xcel(BytesIO(filename.value), header=[0])

        # Editing the file format a bit
        full_data = []
        for i in range(len(data)):
            single_sample = list(data.values())[i][0]
            single_sample.columns = [list(data.keys())[i],]
            full_data.append(single_sample)

        # Storing the samples
        alignment_storage.samples = full_data
        data_reading_section.append(f'**{len(full_data)} samples** were read from {filename.filename}.')

        # Plotting a bar plot
        samp_bar_plot = px.bar(x=[len(i) for i in alignment_storage.samples],
             y=[i.columns[0] for i in alignment_storage.samples], orientation='h',
            title='Metabolic Features per Sample')
        samp_bar_plot.update_layout(
            xaxis_title="Number of Metabolic Features",
            yaxis_title="Samples",
            height=25*len(alignment_storage.samples))
        data_reading_section.append(pn.pane.Plotly(samp_bar_plot, config={'toImageButtonOptions': {
                   'filename': f'FeaturePerSamplePlot_{filename.filename}', 'scale':4}}))

        # Update the possible minimum number of samples and the layout
        alignment_storage.update_widgets(filename.filename)
        alignment_page.append(data_alignment_section)

        pn.state.notifications.success(f'Samples successfully read from {filename.filename}.')

    except:
        pn.state.notifications.error(f'Could not successfully read samples from {filename.filename}. Confirm formatting.')


# Function happens when you press the button        
confirm_button_filename.on_click(_confirm_button_filename)


# Widget to select name of the aligned dataset
dataset_download_name = pn.widgets.TextAreaInput(name='Name to save the dataset as:',
                                           placeholder='aligned_df', value='aligned_df', rows=1)

# Button to save files
save_alignment_files_button = pn.widgets.Button(name='Save Aligned Dataset and Alignment Description',
                                                button_type='success', disabled=False)

def _save_alignment_files_button(event):
    "Saves aligned DataFrame and description."

    try:
        # Saving the file
        alignment_storage.aligned_df.to_csv(f'{dataset_download_name.value}.csv')
        alignment_storage.desc.to_excel(f'{dataset_download_name.value}_desc.xlsx')

        pn.state.notifications.success(f'Aligned dataset successfully saved in {dataset_download_name.value}.csv.')

    except:
        pn.state.notifications.error(f'Dataset could not be saved in {dataset_download_name.value}.csv.')

# Function happens when you press the button        
save_alignment_files_button.on_click(_save_alignment_files_button)


# Section of the interface regarding reading the file
data_reading_section = pn.Column(filename, confirm_button_filename)

# Section of the interface regarding 
data_alignment_section = pn.Column(alignment_storage.controls)

# Section to save the aligned dataset
saving_df_widgets = pn.Column(dataset_download_name, save_alignment_files_button)

alignment_page = pn.Column(data_reading_section)

# Create the main area and display the first page
main_area = AlignmentPage().content

app = pn.template.BootstrapTemplate(title='Data Alignment Interface', sidebar=[], main=[main_area])

app.show(websocket_max_message_size=100*1024*1014, http_server_kwargs={'max_buffer_size': 100*1024*1014})
