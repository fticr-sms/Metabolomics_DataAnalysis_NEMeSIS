# Welcome to NEMeSIS

NEMeSIS (**N**ovel **E**xtreme-resolution **Me**tabolomics **S**treamlined data analysis **I**ntegrated **S**oftware) is a software geared towards FT-ICR-MS extreme-resolution mass spectrometry data. For this objective, it offers a complete pipeline from **raw data (in mzML format) or peak lists to data interpretation** in two different formats (both built with the Python programming language) based on the users preference.

- Jupyter Notebook for a more hands on approach for intermediate (or higher) experienced Python users to take advantage of the implemented metabolomics data analysis workflow and customize it to their liking. The file to start to experience NEMeSIS in this way is `Met_Data_Analysis_Standard.ipynb` - see details on [Get started with jupyter](jupyter_docs.md) page.
- Graphical Interface for a more streamlined experience for non-experienced Python users where no coding knowledge is necessary and results are shown in an interactive format. The file to start to experience MetsTA in this way is `InterfacePrototype_Standard.py`. Detailed instructions and ways to start up the program are in the [Get started with GUI](GUI_docs.md) page.

Here, on this page, we will detail on how to install this software and dependencies to run them. Then, we will have a separate page for each main way to use NEMeSIS with a tutorial on how to use them.

Software made by: **FT-ICR and Structural Mass Spectrometry Laboratory, BioISI, Faculdade de Ciências da Universidade de Lisboa** ([http://ft-icr.rd.ciencias.ulisboa.pt/](http://ft-icr.rd.ciencias.ulisboa.pt/)).

If you use this software in your untargeted metabolomcis data analysis, we would be grateful if you cite our paper introducing it:

- Paper Link: <strong>(LINK)</strong><br>
- Citation: <strong>(CITATION)</strong>

## Installing NEMeSIS

NEMeSIS is currently available as an open source software that can be downloaded or cloned from [https://github.com/fticr-sms/Metabolomics_DataAnalysis_Pipeline](https://github.com/fticr-sms/Metabolomics_DataAnalysis_Pipeline). This includes a `venn.py` file originally provided in [https://github.com/tctianchi/pyvenn/blob/master/venn.py](https://github.com/tctianchi/pyvenn/blob/master/venn.py) and used according to its license. It is a Python-based software and does requires the Python programming language and has many dependencies on different Python packages.

### Installing Python and other packages with Anaconda

We recommend the installation of Anaconda to use the software, with a tutorial available at [https://docs.anaconda.com/anaconda/install/index.html](https://docs.anaconda.com/anaconda/install/index.html).

From here, we present two ways to install the remaining packages needed, the first one package-by-package and the second which will install every package rapidly (requirements.txt based).

##### Package-by-package

On the command line (or the AnacondaPrompt) in your pc, run these code lines (one at a time):

``` { .yaml .no-copy }
pip install metabolinks
pip install UpSetPlot
conda install -c conda-forge py-xgboost
conda install -c pyviz panel
conda install -c anaconda param
conda install -c conda-forge holoviews
conda install -c plotly plotly
pip install python-docx==1.1.0
pip install kaleido==0.1.0post1
pip install pyopenms
```

!!! tip

    If you cannot use `conda` on your command line, use AnacondaPrompt.

##### Installing with requirements.txt

On the command line (or the AnacondaPrompt) in your pc, go inside the directory where you have installed NEMeSIS or open the command line directly in this directory. Then run the following line:

```
pip install -r requirements.txt
```

This will install using `pip` all the packages used in the software (including some already installed with Anaconda).

### Databases for Annotation

The Git-Hub repository includes currently 2 databases (as `.csv` or `.xlsx` files) last updated in February 2023 that you may use for annotation in the software. These are from the Human Metabolome Database ([HMDB](10.1093/nar/gkab1062), [https://hmdb.ca/](https://hmdb.ca/)) for human related metabolites and the [LOTUS Database](https://doi.org/10.7554/eLife.70780) ([https://lotus.naturalproducts.net/](https://lotus.naturalproducts.net/)) for natural products.

