# Metabolomics_DataAnalysis_Pipeline
Repository for the metabolomics data analysis pipeline

## Repository by:

- FT-ICR-MS-Lisboa Laboratory Group

# Required Python Packages and Database File Examples

**Anaconda and Anaconda related packages**

- Most packages needed are automatically installed when you install Anaconda (see tutorial here https://docs.anaconda.com/anaconda/install/index.html).

**Metabolinks and UpSetPlot**

- Open 'Command Line' or 'Linha de comandos' on your pc.
- Run the line 'pip install metabolinks'.
- Run the line 'pip install UpSetPlot'.
- Restart jupyter.

**pyvenn**

- Go to https://github.com/tctianchi/pyvenn/blob/master/venn.py.
- Over 'Raw' (see where 'Watch' is on the upper left and move look down), click 'Save link as...'.
- Save it on the same folder you have the rest of the files.

**Databases**

Right now, the HMDB, LOTUS and DrugBank databases that are used in the Jupyter notebook as examples are present in this repository.

However, in the future, these databases will have corresponding python files that you can run to extract, parse and automatically treat them so you can always have the most current version of them.

# Patch Notes

#### V1.0 - Release Version (03/10/2022)

- First release version of the metabolomics data analysis notebooks.

#### V1.1 - Feature Importance Correction (14/10/2022)

- Correcting bug in Random Forest and PLS-DA feature importance workflow that would not give the correct names from the dataframe to the most important features (it would give the first names in the dataset instead of the important ones).
- Updating instructions for Metabolinks installation.

#### V1.2 - Quality of Life Improvements and Minor Corrections (03/11/2022)

- Added Table of Contents for quick navigation of the different steps in the notebook (also allows further analysis to be added to notebooks without making it clunkier and more difficult to navigate).
- Added Patch Notes to end of notebook for user to follow changes on notebook.
- Correcting bug in Random Forest and PLS-DA feature importance workflow mentioned in V1.1 for the PythonAnnotation Notebook.
- Added feature to step 5 by adding a way to see common and exclusive names/formulas as well as overall compounds/masses. Use id_selection to choose either 'Name' ('Matched names') or 'Formula' ('Matched formulas').
- Other minor changes and corrections.

#### V1.3 - Filtering_pretreatment function and PLS-DA analysis metrics bug fixes (06/12/2022)

- Fixed multiple bugs in `filtering_pretreatment` function where args for filt_method, filt_kw, extra_filt, extra_filt_kw, scaling and scaling_kw would only use the default values even if the user changed it in the function.
- Fixed issue where calculation of f1-score, precision and recall metrics of PLS-DA would fail if the number of samples per class was not a multiple of the number of folds in cross validation.

#### V1.4 - Improvements to common and exclusive features analysis (03/01/2023)

- Added automatic Venn Diagram plotting for common and exclusive features analysis (venn.py must now be in the same folder)
- Added automatic UpSet plot construction (an alternative to Venn Diagram) which may be preferable for larger numbers of classes (may require installing a new package - pip install upsetplot)

#### V2.0 - Reorganization and merging of notebooks, addition of multiple quality of life features, univariate analysis, etc. (08/02/2023)

- Reorganization of notebook for a more streamlined and easy-to-understand experience.
- Merged the two notebooks based on if you had performed metabolite annotation on MetaboScape into one.
- Put all essential functions in a separate file - `metanalysis_standard.py` - to de-clutter notebook. Recommended to check this file too see better the options available for each.
- Updated Step 0 instructions to install metabolinks and other packages.
- Added way to annotate your data with multiple databases in the notebook (a set of 4 columns for each database) and remade annotation section.
- Added way to check for duplicate (or more) annotations and remove those when possible by specific procedures - step 1.3.
- Changed default imputation method to 'min_sample' from 'min_feat'.
- Fixed bug in 'extra_filt' that could be made during `basic_feat_filtering`.
- Venn Diagram now follows each classes' colours as set up during step 1.1.
- Increased the number of options available to plot UpSetPlot including showing the number of annotated compounds.
- Slight change in ROC curve calculation and plotting and bugfix.
- Added more meta_data to important features tables from Random Forest and PLS-DA and allowed output of important feature tables.
- Added Univariate Analysis and Fold Change analysis - Step 6.
- Made points colored in the Van Krevelen based on either the log of their average intensity or the rank of their average intensity (more useful).
- Made Chemical Composition Series appear for all classes in the same plot (more convenient to compare).
- Added Kendrick Mass Defect plots to step 7.
- Added a BinSim specific step for Unsupervised and Supervised analysis.
- Made the compound finder more universal accepting formulas, _m/z_ or Names. Added boxplot and average bar graph to the usual bar graph in this section.
- Minor typo fixes and added explanations.

#### V2.01 - Inserted references (09/02/2023)

- Put clear references to the sources used for the different functions and steps.
- Minor changes.

#### V2.02 - Git-hub Repository (17/02/2023)

- Changes to put the project in a git-hub repository (patch notes now in Readme.md file of repository).
- Added Readme.md and requirements.txt file.
- Updated package installation instructions.
- Changed HMDB database to be a xlsx file (smaller file) and minor bugfix in database reading.

#### V2.1 - Permutation Tests, Volcano Plots and Loadings in PCA (06/06/2023)

- Updated and added PCA function to metanalysis.py (now retrieves loadings).
- Added Permutation Testing for Random Forest and PLS-DA supervised analysis that can accept one of multiple metrics.
- Added ROC Curve for PLS-DA and updated ROC section for Random Forest in notebook.
- Added Volcano Plot prototype in visualization of unsupervised analysis.
- Added PCA and Loading plot (biplot) based on the counts of each chemical composition series in each sample.
- Some minor changes and additions in comments.

#### V2.2 - HMDB Compound Pathway Assignment (17/11/2023)

- Added possibility of reading target directly from the excel or csv file with your data (changed initial filtering place to account for this).
- Added section that allows to see pathways related to HMDB compounds annotated.
- Fixed 'Normalized Intensity' position (y axis) in compound finder plots.
- Fixed deprecation issue in Kendrick Mass Defect plots.

#### Current Patch modifications (29/11/2023)

- Bugfix in chemical composition series plots where using multiple database formula annotations would only include formulas assigned in all chosen databases instead of in at least one of them.
- Small change in Kendrick Mass Defect plots allowing points to be coloured by chemical composition series based on database formula annotation performed in the notebook instead of only being able if formulas came from MetaboScape 'SmartFormula'.

#### Current Patch modifications (04/12/2023)

- Small bugfix in compound finder search tool when looking for m/z values, finder DataFrame would not update.

#### In Future Patches:

- Re-update all section 7 of notebook for better functionality.
- Add sMDiN analysis to notebooks.
- Add XGBoost analysis.
- Add customization of fonts in figures.
- Grant possiblity of extracting more dataframes.
- (Long term) provide graphical user interface as notebook alternative.