
## Introduction to the illustrative jupyter notebook containing the workflow in NEMeSIS

These jupyter notebooks contain a lot of descriptions to help guide the user already. Thus, to avoid repetition we will focus on their general aspects instead of detailing each step like it is done for the graphical interface.

The larger jupyter notebook version of NEMeSIS is `NEMeSIS_Data_Analysis_V2.8.ipynb`. This notebook includes every step of the main steps of data analysis from the alignment of _m/z_ peak lists representing samples to biological interpretation. Furthermore, an extra notebook regarding the conversion of raw spectral data in open mzML format to _m/z_ peak lists is available in `Raw_Spectra_Processing.ipynb` (see here for details [mzML Spectra Conversion](jupyter_docs.md#mzml-spectra-conversion)). This was placed as a separate noteboook to lower memory necessity per individual notebook making the experience smoother. The interactive graphs with a very high number of points in `Raw_Spectra_Processing.ipynb` could affect performance in downstream steps.

The main notebook is split into many different steps from step 0 to 11 which are described in the table of contents. Using this table is the fastest way to navigate the (quite large) notebook.

Step 1 is the equivalent to Stage 1 of the graphical interface: **Data Reading**.

Steps 1.1, 1.2, 1.3 and 2 are the equivalent to Stage 2 of the graphical interface: **Data Pre-Processing and Pre-Treatment**

Steps 3 to 11 are the equivalent of Stage 3 of the graphical interface: **Data Analysis and Biological Interpretation**

Finally, the program is capable of outputting trated data ([Data Exporting](jupyter_docs.md#data-exporting)) that can be used as input for other statistical analysis not present in the main NEMeSIS program in independent side modules ([Independent Side Modules](jupyter_docs.md#independent-side-modules)). These independent side modules are standalone jupyter notebooks that perform statistical analysis not present in the main notebook and can be adapted to specific statistical analysis desired by the researched. We present two independent side modules [here](jupyter_docs.md#independent-side-modules) that are, simultaneously, example side modules to be adapted by researchers to create others and that contain very useful analysis based on graph representations of samples built to consider possible biochemical transformations between metabolites - sMDiNs and FDiGNNs - to be performed based on the type of analyses desired - more information on these side modules [here](jupyter_docs.md#independent-side-modules).

![GUI Structure](img/GUIStructure.png)

## Definition of Parameters

Looking at the first cell after `Step 1` shows the idea on how the jupyter notebook works.

Each cell may have at its beginning a set of parameter that can be edited as shown below. These include the parameter name and the current parameter that can be edited. Usually, the parameter and what should be put there is explained in a comment after the parameter (or in the lines before if the explanation is large). In some cases, if the parameter only has a set number of options allowed, they are also enumerated in a comment.

```
    filename = '5yeasts_notnorm.csv' # Name of your file
    ## Indicate if the file read includes the target (sample classes) in its first row
    target_in_file = False
    idx_masses = 'Neutral' # 'Neutral' (neutral masses), 'Positive' (Obtained from ESI+), 'Negative' (Obtained from ESI-)
    # or 'None' (Neutral mass column cannot be inferred)
```

The rest of the cell then contributes to the current step of the analysis being performed. Most functions used in the data analysis are stored in the `metanalysis_standard.py` file while functions and code related to rendering figures or showing results are in the notebook itself to be more easily edited.

## Data Reading

!!! info

    The first part of this section shares high similarities with the [Stage 1: Data Input](GUI_docs.md#stage-1-data-input) section of the [Get started with GUI](GUI_docs.md) page.

NEMeSIS accepts previosuly aligned 2D tables with samples on one axis and metabolic features on the other, samples represented as lists of _m/z_ peaks or spectral raw data sample in mzML format. The first steps of the notebook is reading your input data, which requires strict formatting from the files. It can accept either **data tables** with samples having already been previously aligned or samples represented as **lists of _m/z_ peaks**. The type of data to be used is chosen by choosing the `aligned_samples` parameter (first cell of Step 1) to `True` if data to be used is in form of **data tables** or `False` if it is **lists of _m/z_ peaks**. For the former, you can then skip to the `Data Matrix (Aligned Samples) Reading Section` to define the filename and reading specific parameters while, for the latter, you go to the `Data Alignment Section` to define the filename and alignment specific parameters.

For using **spectral raw data sample in mzML format**, the mzML files must be converted in a separate notebook capable of transforming these spectra into **data tables** and/or **lists of _m/z_ peaks** to use as input for the main notebook. See [mzML Spectra Conversion](jupyter_docs.md#mzml-spectra-conversion) for more information.

The allowed formatting for **data tables** is exemplified in the figure below with the addition that files can either be `.csv` or `.xlsx` files. Metabolic Features should be represented in the rows, while samples and metadata should be represented in columns. 

The first column in your data must correspond to the identifier of your metabolic features, more specifically, to a mass value if possible - black box - and they should be unique (non-repeating). The name of the column will be overwritten to be `Bucket label`. If the columns are mass values capable of being intepreted as floats (numbers), then a `Mass` column in your data will be added. This column is later necessary to perform `Data Annotation` (but nothing else). You have 3 possibilities to choose from: `Neutral` where the masses in your index are taken as the neutral masses creating a `Neutral Mass` column, `m/z (Positive)` or `m/z (Negative)` where the masses are assumed to be _m/z_ values obtained in positive or negative (respectively) ionization mode creating a `Probable m/z` column.

!!! danger

    Avoid having a column named `Neutral Mass` or `m/z` in your data, since it will be overwritten.

The remaining columns are free to be any metadata you have in your data (green box) and the samples with the corresponding intensity values (purple box). If you are using LC-MS or GC-MS or any other kind of hyphenated method and have another column other than mass to characterize your features such as retention time, you can still analyse them and identify this columns as metadata but they will not be used for data annotation or any other step in the analysis.

Optionally, the first row (after the column names) of your data may include the class labels of each of your samples (orange box) that will be used to automatically generate the target (class labels) later on.

![Data Format Example](img/ExampleData.png)

The allowed formatting for **lists of _m/z_ peaks** to be aligned is exemplified in the figure below. The software will expect a **single Excel (`.xlsx` or `.xls`) file**. It should have **one sample per Excel sheet**; the **name of the Excel sheet should correspond to the sample name** (orange box); each sheet should have in its **first column the mass values** (m/z, neutral mass or equivalent) and in its **second column the corresponding intensity values**; and finally, the first row should have the name of the two columns, for example, 'm/z' and 'I' (this name should be consistent between samples if possible). The example file `example_samples_to_align.xlsx` is available in the `Files_To_Align` folder as guidance.

![Formatting of Mass Peak Lists](img/SampleToBeAligned.png)

!!! info

    The nature of the jupyter notebooks allows the researcher to tweak the data reading functions so they accept data in different formats as well.

## Figure and Table Saving

The jupyter notebook allows saving every table and figure made in the analysis. For most of the more relevant cases, this is already prepared to be performed, although the default option is turned off.

For Figures, most have at the end of the code cell that produces something like the small code example shown below for the PCA Figure. By removing the '**#**' (comment), the figure will be saved when the cell is ran again (name can be adjusted).

``` { .yaml .no-copy }
#f.savefig('Name_PCAplot.png', dpi=600) # Save the figure
```

For the most significant Tables, they have usually a whole code cell to format and save the table like in the example shown below to save the list of important features to build Random Forest models. Here, at the beginning of the cell, there is usually a parameter than can be changed between `False` and `True` to save or not save the table. In the example shown, the parameter name is `SAVE_IMP_FEAT`. By changing it to `True`, the table would be saved. Another example is the `GENERATE_Excel_file` parameter in the `Common and Exclusive Compound Analysis` step.

``` { .yaml .no-copy }
# Saving Important feature dataset in an excel
SAVE_IMP_FEAT = False

# Saving the most important features by their fraction 'frac_feat_impor'.
# If None, saving the most important features based on a threshold 'VIP_Score_threshold'.
# If also None, save the full dataset of all features
frac_feat_impor = 0.02 # Fraction of features to save, If None the variable in the next line is used.
score_threshold = None # Only used if variable above is None, threshold of score to consider a feature important.

if SAVE_IMP_FEAT:
    if frac_feat_impor:
        max_idx = int(frac_feat_impor*len(imp_feats_rf))
        filt_imp_feats_rf = imp_feats_rf.iloc[:max_idx]
        filt_imp_feats_rf.to_excel(f'RF_ImpFeat_{frac_feat_impor*100}%.xlsx')
    elif score_threshold:
        filt_imp_feats_rf = imp_feats_rf[imp_feats_rf['Gini Importance'] > score_threshold]
        filt_imp_feats_rf.to_excel(f'RF_ImpFeat_GiniImpgreater{score_threshold}.xlsx')
    else:
        imp_feats_rf.to_excel(f'RF_FeatByImportance.xlsx')
```

## Data Exporting

For exporting the data tables obtained after pre-processing and pre-treatment, a cell is available for saving different formats of the file in the NEMeSIS folder. This has the purpose to save data for the researcher to use somewhere else but also to be used as input for the [Independent Side Modules](jupyter_docs.md#independent-side-modules). The multiple files are made because they can be useful in different contexts. Here, we show these files with their default name but this name can be changed in the same cell code mentioned. These are:

- `Export_TreatedData.xlsx` - Excel with 4 sheets, one containing the normalized data (without missing value imputation, transformation or scaling) with metadata, another the treated intensity data without metadata, another with the dataset treated with Binary Simplification (or Spectral Digitalization) method and the last with the treated values after missing value imputation and normalization (but before transformation and scaling).
- `Export_Target.txt` - File containing the class labels (target) of the dataset's samples in the same order as they appear in the dataset.
- `Export_TreatedData.pickle` - Treated intensity data without metadata in pickle format to guarantee the 'Bucket Labels' suffer no changes (roundings).
- `Export_ProcData.pickle` - Normalized data (without missing value imputation, transformation or scaling) with metadata in pickle format to guarantee the 'Bucket Labels' suffer no changes (roundings).

!!! info

    The pickle versions of the dataset exist to eliminate the possiblity of rounding errors made to the index of the data tables, keeping them consistent between analyses.

## mzML Spectra Conversion

The mzML Spectra Conversion can be performed in an independent side jupytre notebook available in the NEMeSIS installation called `Raw_Spectra_Processing.ipynb`. This file can converts the open format mzML spectra into **lists of _m/z_ peaks** that can be aligned into 2D **data tables**. For this purpose, the [pyopenms](https://pyopenms.readthedocs.io/en/latest/) Python package is used.

The names of the file should be inputted in the 2nd code cell of the notebook in a list as exemplified in the notebook itself. The most critical parameter in `signal_to_noise` ratio can also be adjusted in the following cell. Furthermore, other parameters that can be adjusted are in the cell where the conversion is performed and with a comment as seen in code section below. To use these parameters, the # before each line can be removed and the parameter adjusted.

``` { .yaml .no-copy }
    param = cnt.getParameters()
    param.setValue("signal_to_noise", signal_to_noise)
    #param.setValue("spacing_difference_gap", 7.0)
    #param.setValue("spacing_difference", 5.0)
    #param.setValue("SignalToNoise:win_len", 200.0) # def 200
    #param.setValue("SignalToNoise:bin_count", 200) #def 200
    #param.setValue("SignalToNoise:min_required_elements", 100)
    #param.setValue("SignalToNoise:max_intensity", 1000000)
    #param.setValue("SignalToNoise:auto_mode", -1)
```

Explanations of each parameter are shown with extensive documentation is present in [https://openms.de/current_doxygen/html/classOpenMS_1_1PeakPickerHiRes.html](https://openms.de/current_doxygen/html/classOpenMS_1_1PeakPickerHiRes.html). Credit to pyopenms devs.

The **lists of _m/z_ peaks** can be directly exported to an excel in the format appropriate for the Data Alignment section of the main jupyter notebook. Furthermore, the next section can also performed the data alignment of the **lists of _m/z_ peaks** equivalent to the main jupyter notebooks, that can also be exported to a `.csv` file that can be the input of the main jupyter notebook or the NEMeSIS GUI as well.

##### Spectra Visualization

!!! info

    This section shares high similarities with the first part of the [Spectra Visualization](GUI_docs.md#spectra-visualization) section of the [Get started with GUI](GUI_docs.md) page.

If files were succesfully converted and data alignment was performed, the last section allows to visualize and assess the quality of the spectral data processing by observing the raw, centroided (after spectral processing) and aligned (after data alignment) spectra in interactive graphs of a few samples simultaneously. They are interactive so close ups can be selected for each graph. This allows to observe the genral shape of peaks that were kept through the process and which were discarded and to see if the quality of the aligned spectra is as desired or expected.

The image below shows an example. You can choose which samples to see, if the spectra values should be normalized (1), and which mass ranges to plot initially.
{ .annotate }

1. Normalization in the raw, centroided and aligned spectra is made by dividing their intensities the sum of all intensities in the **raw data** specifically.

!!! warning

    Spectra contain a very large number of points which combined to their interactive nature, make it heavy resource demanding for the software and the computer. Hence, although it depends on the computer available, we do not recommend selecting more than 3-4 samples simultaneously. Moreover, it may take a few seconds until the spectra are fully operational.

![Mass Spectra Example](img/MSpectra.png)

## Independent Side Modules

The modular nature of the statistical analysis together with the Python implementation of NEMeSIS allows the addition of extra statistical analysis tailored to the specific dataset analysed by using independent side modules. These are illustrated how they can be applied by resorting to jupyter notebooks. To this end, NEMeSIS allows exporting data after the Data Pre-Processing and Pre-Treatment steps to be used by the independent modules. In order to facilitate the creation of these side modules, we present two side modules named `sMDiN_analysis_module.ipynb` and `FDiGNN_analysis_module.ipynb` that introduce how the data saved by the main workflows of NEMeSIS can be inputted and used.

Moreover, these notebooks apply methodologies developed in the FT-ICR-MS-Lisboa laboratory group that we believe could improve and complement most analysis performed using the NEMeSIS software. The `sMDiN_analysis_module.ipynb` notebook applies the sample Mass-Difference Network (sMDiN) methodology as described in its paper [here](https://doi.org/10.3389/fmolb.2022.917911) and the `FDiGNN_analysis_module.ipynb` notebook applied the Formula-Difference Graph Neural Network (FDiGNN) methodology as described here (**put paper DOI when accepted**) - currently submitted paper under review. Both methodologies focus on the representation of metabolomics as graph (non-tabular) representations instead of the 2D Data Tables usually used. These graph representations are built by taking into account possible chemical transformations between the detected metabolites, thus inserting into the data structural information that represents these possible transformations. The analysis performed will then include this implicit information stored in the graphs and use it to inform supervised models and metabolite importance, which can complement the conventional workflow. These methodolgoies are not included in the main NEMeSIS workflow since they are not a part of the conventional workflow most researchers start with to analyse data quality and because their use on the graphical interface would be more clunky as a consequence of the graph based approaches which could lead to a loss of nuance. Despite this, their complementarity and usefulness makes them indispensable to add as a side analysis to be performed.

sMDiNs besides using a graph representation also use feature occurrence data instead of intensity data. The built a graph representing each sample that should be characteristic of each biological class and perform network analyis on the set of graphs constructed. The network analysis made is tied to the type of information that is desired to be extracted from the graphs (observe the performance of supervised models to see if the extracted data is reliable, that is, models with low performance cannot provide reliable information). For example, Weighted Mass-Difference based Building blocks Impact can be used to observe the difference in prevalence of different biochemical transformations used to build the graphs between classes - see paper [here](https://doi.org/10.3389/fmolb.2022.917911) for details.

FDiGNN builds a graph representation for each sample, adding attributes to the nodes (each representing a metabolite) such as intensity or feature occurrence. It then trains a Graph Neural Network model (with many parameters that can be adjusted and optimized) that can accept those graphs directly as input. The architecture of the network itself is also available to be altered and modified to the researcher's desires. By using PINNI (Probability-Impact based Network Node Importance), a measure of importance to each metabolite/node can be assigned. This prioritizes nodes near other nodes (in the network) that are important for discrimination. Thus, the methodology prioritizes the highlighting of network sections (associated metabolites by possible chemical transformations) over individual metabolite importance. Pathway analysis can also be performed to overlap this important sections to known pathways in order to improve biological interpretation. This unique emphasis on possible 'interactions' between metabolites is a great advantage of FDiGNNs. However, fitting neural networks can be computational heavy especially with larger networks which may lead to the model taking a few hours to be fit and information extracted.

As for structure, the first few code cells in each case show how the exported data can be read back into the jupyter notebook obtaining all the main data tables (treated data, annotated data, BinSim treated data) as well as the target with the class labels of each sample. They also show the overall details that can be adjusted at the beginning, in these cases, the number of folds to perform stratified cross-validation and the number of iterations to repeat the cross-validation when fitting and estimating performance of supervised models. These parameters can change based on the type of analysis performed.

After this section, each notebook contemplates the necessary procedures to apply their respective methodologies with detailed informations on the different steps present in the notebooks themselves. The `FDiGNN_analysis_module.ipynb` notebook has, at the end, a mini Dash application that will open on the browser to browse and analyse results from the FDiGNN analysis and metabolite importances extracted (as well as pathway analysis) facilitating interpretation.

