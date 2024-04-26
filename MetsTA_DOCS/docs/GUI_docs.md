
## Starting the Graphical Interface of MetsTA

I am still unsure of what will be the final version to open the graphical interface. 

**CHANGE THIS WHEN THE STARTING METHOD IS DECIDED.**

## Introduction to the Graphical Interface of MetsTA

When opening the software, you are greeted with a page such as the one presented below that should open in your browser of choice. 

Here, you can see on the on the left-hand side of the page (in a black box), you get an index-like column with buttons for every main page of the software. Most will be locked off in the beginning since they require other steps to be completed first before accessing. This is the main way to access different pages in the software. The first button (`Home`) and initial homepage of the software leads you to an introduction of what this software is and what analysis steps it includes. The second button (`Considerations and Instructions`) takes you to a few explanations and clarification on the workings of the software, which are also present here at the end of the page in [Tips and Precautions](GUI_docs.md#tips-and-precautions). From then on, we go onto the data analysis pages.

The green box marks the main section where the page currently selected appears. The red box in the upper right corner has a small ring in it. This ring will indicate to you if the program is currently computing anything. If it is, there will be a small black section going around the ring. While the program is computing, you cannot update anything else in it.

![GUI Opening Page](img/OpeningPage_GUI.png)

The data analysis pipeline is divided into 3 stages:

1. Data Reading (a single page so it is in the Data Pre-Processing and Pre-Treatment section)
2. Data Pre-Processing and Pre-Treatment
3. Data Analysis and Biological Intepretation (Statistical Analysis)

![GUI Structure](img/GUIStructure.png)

## Stage 1: Data Reading

The first section of the data analysis pipeline is data reading, that is, the data you will input to the program and thus requires a strict formatting from the files. The software accepts already data tables with the raw data of the sample having already been previously aligned. From then on, the software is capable of performing the remaining data analysis.

This formatting is exemplified in the figure below with the addition that files can either be `.csv` or `.xlsx` files. Metabolic Features should be represented in the rows, while samples and metadata should be represented in columns. 

The first column in your data must correspond to the identifier of your metabolic features, more specifically, to a mass value if possible - black box - and they should be unique (non-repeating). The name of the column will be overwritten to be `Bucket label`. If the columns are mass values capable of being intepreted as floats (numbers), then a `Mass` column in your data will be added. This column is later necessary to perform `Data Annotation` (but nothing else). You have 3 possibilities to choose from: `Neutral` where the masses in your index are taken as the neutral masses creating a `Neutral Mass` column, `m/z (Positive)` or `m/z (Negative)` where the masses are assumed to be _m/z_ values obtained in positive or negative (respectively) ionization mode creating a `Probable m/z` column.

!!! danger

    Avoid having a column named `Neutral Mass` in your data, since it will be overwritten.

The remaining columns are free to be any metadata you have in your data (green box) and the samples with the corresponding intensity values (purple box). If you are using LC-MS or GC-MS or any other kind of hyphenated method and have another column other than mass to characterize your features such as retention time, you can still analyse them and identify this columns as metadata but they will not be used for data annotation or any other step in the analysis.

Optionally, the first row (after the column names) of your data may include the class labels of each of your samples (orange box) that will be used to automatically generate the target (class labels) later on.

After data is read, the next step of the analysis will be unlocked (loading an example dataset to test the software is possible).

![Data Format Example](img/ExampleData.png)

!!! info

    If you want, you can load the parameters used in a previous analysis at this stage as well. Further details about this in [Report Generation and Parameter Saving](GUI_docs.md#report-generation-and-parameter-saving).

## Stage 2: Data Pre-Processing and Pre-Treatment

The second stage of the software is the Data Pre-Processing and Pre-Treatment that includes many sub-sections: [Metadata Selection](GUI_docs.md#metadata-selection), [Data Filtering](GUI_docs.md#data-filtering), [Data Annotation](GUI_docs.md#data-annotation), [De-Duplication of Annotated Features](GUI_docs.md#de-duplication-of-ann-features) and [Data Pre-Treatment](GUI_docs.md#data-pre-treatment).

All these sections are mandatory and are progressed and unlocked one by one in a linear fashion. You can redo previous parts of the pre-treatment. If you **change any part of the pre-treatment**, posterior pre-treatment steps and statistical analysis will immediately **be locked** and will **erase previous statistical analysis** if any has been performed. This is to avoid situations where the user wanted to change the pre-treatment but forgot to re-apply posterior pre-treatment steps. Thus, it is necessary to keep coherency.

#### Metadata Selection

Here, you will be asked to select which columns of your data are metadata (non-sample columns). If this metadata includes formula assignment or metabolite compound annotation columns, you can select them in their corresponding spots so they can be used downstream as that. Furthermore, you can select a column that contains the mass values of your metabolic features, whether it is the `Neutral Mass` or `Probable m/z` columns created in the previous section or another prepared previously that will be used for `Data Annotation` (if `None` is chosen, `Data Annotation` cannot be performed). Other metadata will not be used in the analysis.

After selecting it, either a target will try to be inferred from the name of your samples or it will be taken from your data (if it was present). In either case, you can edit your target to your needs, before confirming it and moving to the next section.

![Metadata Selection Example](img/MetaDataSelection.png)

#### Data Filtering

Data Filtering is used to _clean_ the dataset of features that appear in very few samples. You can select to perform Data Filtering either by keeping features if they appear at least in **n** samples from the total number of samples in your dataset (`Total Samples` method) or in **n** samples of at least one specific class in your data (`Class Samples` method). This **n** is referred to as the _feature filter keyword_ and can go from 1 to the number of samples in your data in the first case or to the number of samples of your least populated class in the second. After Data Filtering, a small table detailing some characteristics of your data as well as your dataset will be shown.

#### Data Annotation

Next up, we have Data Annotation. This step is only possible if you have a Mass column as selected in the earlier section (either created in the software based on mass values provided or previously obtained) that represents the masses of the metabolites detected as numbers. 

You can skip this step by annotation with _0_ databases. Or you can select 1 to 5 databases to perform independent annotations. A database must have a compound ID column, a compound name column and a compound formula column (that is used to calculate the compound theoretical mass used for annotation). To use a database, it should be on the directory where MetsTA is located and you must provide a series of informations such as the file name and the column identifier of the ID, name and formula of the compounds. You will also be asked for an abbreviated name of the database to use (for example, HMDB for Human Metabolome Database).

After loading your desired databases, you can choose a maximum threshold of deviation between the theoretical masses of the compounds in the database and the neutral masses in the dataset. This can be a flat threshold (`Absolute Dalton Deviation`) or a parts per million based threshold (`PPM Deviation`) (1).
{ .annotate }

1. A usual value for extreme-resolution data could be a ppm deviation between 0.5 and 1 ppm.

Finally, you can decide which adducts to search by inputting the adduct name and consequent mass shift they cause on a neutral mass. This should be inputted in the following format "Adduct_Name : Mass_Shift_to_neutral_mass". As an example: "[M+H]+ : 1.007276451988935". The software automatically provides an example of common adducts to search based on if `Positive (m/z)` ([M+H]+, [M+Na]+, [M+K]+), `Negative (m/z)` ([M-H]-, [M+Cl]-) or `Neutral` ([M]) were chosen in Step 1.

With these parameters selected, you can perform the annotation. The annotation is made by comparing the dataset's masses with the databases. Every compound in the database that fall within the error margin threshold chosen to a given metabolite feature is annotated to that feature. Thus, very often features have multiple annotations associated (for example, every isomer in a database will be annotated if one is). This annotation is based on the fact that the mass value of a feature is not enough to differentiate and choose between different isomers. The results of the annotation of our data based on an individual database is the creation of 3 metadata columns: for each metabolic feature, one has a list of compound IDs annotated, another of compound names and another of compound formulas (1). This is made independently for each database selected, that is, every single database selected will generate 3 new columns with their own annotations (2).
{ .annotate }

1. The lists are in the same order in the 3 columns. Thus, the 3rd compound ID in its list corresponds to the 3rd compound name in its list.
2. This is not ideal but avoids the issue of database merging leading to the same metabolite being represented multiple times due to a lack of a standardized metabolite identifier.

#### De-Duplication of Ann. Features

Blah blah, Lorem Ipsum even maybe a bit of dolor.

#### Data Pre-Treatment

Data Pre-Treatment is a series of successive operations that include: missing value imputation, normalizations, transformations and scaling. The only mandatory operation is missing value imputation since the existence of missing values does not allow the application of many statistical methods downstream. These procedures have the aim of highlighting relevant biological variation while reducing the effect of undesired variation ([van den Berg et al., 2006](https://doi.org/10.1186/1471-2164-7-142)) and make the values in dataset more amenable to posterior downstream analysis.

Each category has multiple different approaches available which we will present next and, since data pre-treatment highly impacts posterior data analysis, careful deliberation on the combination of pre-treatments to use is advised.

!!! danger

    Not all possible combinations of pre-treatments are viable. There are some combinations which included incompatible methods between them or with your data. As two examples, **Zero** missing value imputation cannot be used with a generalized logarithmic transformation and normalization by a reference feature cannot be used if your reference feature does not appear in **all the dataset samples**. This incompatibility will generate missing values in the data which the software will warn you about.

***Missing Value Imputation***

Missing Value Imputation aims to replace missing values in the dataset (metabolic features that do not appear in certain samples despite appearing in others) with "probable" intensity values to make possible many different statistical analysis. There are 3 types of missing values: Missed At Random (MAR), Missed Not At Random (MNAR) or Missed Completely At Random (MCAR). The methods available in this software assume the missing values are MNAR, that is, that they were missed not due to instrumental error but due to the metabolic feature being absent or in concentrations / intensities below the detection limit in that sample. Thus, these methods replace this missing values with low values below most other intensity values in the data.

| Missing Value Imputation Methods Available | Description                          |
| ------------------------------------------ | ------------------------------------ |
| Minimum of Sample                          | Replace missing values with a fraction (0 to 1) of the minimum intensity value in the corresponding sample. |
| Minimum of Feature                         | Replace missing values with a fraction (0 to 1) of the minimum intensity value in the corresponding metabolic feature. |
| Minimum of Data                            | Replace missing values with a fraction (0 to 1) of the minimum intensity value in the dataset. |
| Zero Imputation                            | Replace missing values with 0. |

***Normalization***

| Normalization Methods Available            | Description                          |
| ------------------------------------------ | ------------------------------------ |
| By Total Sum of Intensities                | Dividing each intensity by the total sum of intensities of the corresponding sample (all values between 0 and 1). |
| By a Reference Feature                     | Dividing each intensity by the intensity of the reference feature in the corresponding sample. |
| Probabilistic Quotient Normalization (PQN) | Normalization method proposed by Dieterle et al. in 2006 - see [paper here](https://doi.org/10.1021/ac051632c) for details. |
| Quantile Normalization                     | As explained by Bolstad et al. in 2003 - see [paper here](https://doi.org/10.1093/bioinformatics/19.2.185) for details. |
| None                                       | No Normalization Performed. |

!!! tip

    Besides the incompatibilities mentioned earlier, we recommend not to use Quantile Normalization when your data has a high percentage of missing values, which is more likely with direct infusion extreme-resolution data that this software aims to analyse.

***Transformation***

| Transformation Methods Available              | Description                          |
| --------------------------------------------- | ------------------------------------ |
| Generalized Logarithmic Transformation (glog) | Apply the equation: log<sub>2</sub>(y + (y<sup>2</sup> + lambda<sup>2</sup>) /2). When lambda = 0, it equivalent to a logarithmic transformation |
| None                                          | No Transformation Performed. |

***Scaling***

| Scaling Methods Available | Description                          |
| ------------------------- | ------------------------------------ |
| Mean Centering            | $\widetilde{x}_{ij} = x_{ij} - \overline{x}_{i}$ |
| Pareto Scaling            | $\widetilde{x}_{ij} = \dfrac{x_{ij} - \overline{x}_{i}}{\sqrt{s_{i}}}$ |
| Auto Scaling              | $\widetilde{x}_{ij} = \dfrac{x_{ij} - \overline{x}_{i}}{s_{i}}$ |
| Range Scaling             | $\widetilde{x}_{ij} = \dfrac{x_{ij} - \overline{x}_{i}}{x_{i_{max}} - x_{i_{min}}}$ |
| Vast Scaling              | $\widetilde{x}_{ij} = \dfrac{x_{ij} - \overline{x}_{i}}{s_{i}}$ . $\dfrac{\overline{x}_{i}}{s_{i}}$  (associated [paper](https://doi.org/10.1016/S0003-2670(03)00094-1)) |
| Level Scaling             | $\widetilde{x}_{ij} = \dfrac{x_{ij} - \overline{x}_{i}}{\overline{x}_{i}}$ (the scaling factor can be the median of ${x}_{i}$ as an alternative)|
| None                      | No Scaling Performed. |

!!! info

    Scaling equations obtained from van den Berg et al., 2006, paper available [here](https://doi.org/10.1186/1471-2164-7-142).


## Stage 3: Data Analysis and Interpretation

The third and final stage of the software is the Data Analysis and Interpretation (Statistical Analysis) that also includes many sub-sections: [Common and Exclusive Compound Analysis](GUI_docs.md#com-and-exc-compound-analysis), [Unsupervised Analysis](GUI_docs.md#unsupervised-analysis), [Supervised Analysis](GUI_docs.md#supervised-analysis), [Univariate Analysis](GUI_docs.md#univariate-analysis), [Data Diversity Visualization](GUI_docs.md#data-diversity-visualization), [HMDB IDs to Pathways Assignment](GUI_docs.md#hmdb-ids-to-pathways-assignment), [BinSim Treated Data Analysis](GUI_docs.md#binsim-treated-data-analysis) and a [Compound Search Tool](GUI_docs.md#compound-search-tool).

#### Com. and Exc. Compound Analysis

Blah blah, Lorem Ipsum even maybe a bit of dolor.

#### Unsupervised Analysis

Blah blah, Lorem Ipsum even maybe a bit of dolor.

#### Supervised Analysis

Blah blah, Lorem Ipsum even maybe a bit of dolor.

#### Univariate Analysis

Blah blah, Lorem Ipsum even maybe a bit of dolor.

#### Data Diversity Visualization

Blah blah, Lorem Ipsum even maybe a bit of dolor.

#### HMDB IDs to Pathways Assignment

Blah blah, Lorem Ipsum even maybe a bit of dolor.

#### BinSim Treated Data Analysis

Blah blah, Lorem Ipsum even maybe a bit of dolor.

#### Compound Search Tool

With this search tool, you can localize a metabolic feature of interest in your dataset by its index value (`Bucket label`), its name (by searching all annotation columns either selected as annotation or made during analysis), formula (by searching all formula assignment columns) or neutral mass (by searching the `Neutral Mass` column if it exists). A series of bar plots and box plots as shown below will then be computed to show the variation of the metabolic feature between samples and classes.

!!! info

    The values shown are obtained from the original data inputted with only normalization (the one chosen in Data Pre-Treatment being applied). No Missing Value Imputation or posterior parts of the Data Pre-Treatment are applied.

![Compound Search Tool Example](img/SearchToolExample.png)

## Report Generation and Parameter Saving

The graphical interface offers two different approaches to save your analysis made besides the tables and figures you can download during the analysis. These are the report generation function and the parameter saving and loading function.

The report generation creates a folder with a user chosen name. This folder contains a main read-only word .docx file and a myriad of images, interactive figures and tables that support the report. The report includes all of the data reading and data pre-processing and pre-treatment stages and includes the statistical analyses chosen by the user in the corresponding page. When generating the report, a pop-up will appear confirming if there were issues or not in creating the report. The report attempts to describe the analysis made including all the relevant parameters used for the current analysis (1). Furthermore, these parameters are also used for the name of the different figures and tables generated.
{ .annotate }

1. If the parameters selected were changed after the analysis is made and the analysis was not redone, then the parameters in the report should still be the ones originally used in the analysis.

!!! note

    If you find a mistake in the generated report, such as a parameter used in the analysis not being the one that is in the report or in the associated figures, please inform us thorugh [https://github.com/fticr-sms/Metabolomics_DataAnalysis_Pipeline/issues](https://github.com/fticr-sms/Metabolomics_DataAnalysis_Pipeline/issues). We would greatly appreciate so we can iron out any loose ends and problems with the software.

Parameter saving is a function that specifically aims to save the **currently used parameters** in the analysis for future use on another analysis. It is available at two different locations in the interface: at the end of the data pre-processing and pre-treatment stages and in the report generation page. The former only saves parameters regarding the data pre-processing and pre-treatment steps while the latter saves all currently used parameters in the software. Although most more significant parameters are saved, there are some that cannot be saved since they are very specific to the dataset currently under analysis. A few examples are the metadata columns selected or the control/test class used in Univariate Analysis.

The idea behind parameter saving is that you can load in these parameters and the `Data Reading` page of the software to use in other posterior analysis. Furthermore, when loading in there are some parameters that might not be applicable to the current case such as the feature used in Normalization by a Reference Feature or the minimum number of samples used in data filtering. For these cases, a check will be made to see if they can be used and if not, the default value will remain.

!!! warning

    Some parameters are saved but should still always be adapted based on the dataset analyse. A case and point is that of the number of components used for building PLS-DA models. This number highly affects the performance of the models and is characteristic of the dataset used. Thus it should always be adapted to the current dataset. Other parameters such as pre-treatment related can be used on a standard preferred pipeline independent of the dataset analysed.

!!! info

    Within the saved parameters, there are some that are saved based on what is currently present in the program and not what was originally used for the analysis: these are the methods used for the common and exclusive compound analysis, the minimum and maximum number of components used for optimization of the PLS-DA, the model performance evaluation metrice used, the dpi of the permutation figure and number of iterations used in the ROC analysis of the PLS-DA section (for both intensity treated data and BinSim treated data) and the model performance evaluation metrice used and the dpi of the permutation figure of the Random Forest section (for both intensity treated data and BinSim treated data).

## Closing and Resetting the Software

To close the program, we recommend going to the **terminal / command line** tab that opened when you opened the software and closing it, which will terminate the program. Closing the software tab in your browser without closing the terminal will keep the program running in the background.

If you want to redo your analysis with another dataset, you may use the `RESET` button available at the end of the index on the left-hand side. This will reset all your current analysis and parameters chosen. A similar effect but with a softer reset happens when you change a parameter in any of the pre-processing or pre-treatment stages, however you will not be able to read another dataset until you perform a full reset of the software.

!!! bug

    After resetting the software once using the `RESET` button, when clicking it again, the window pane to confirm the reset might not appear. To fix this issue, refresh the image and it should appear once again.

## Tips and Precautions

Blah blah, Lorem Ipsum even maybe a bit of dolor.

Image and table files names as well as in report generations

Refresh the page if something seems out of place