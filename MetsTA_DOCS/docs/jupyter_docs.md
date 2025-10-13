
## Introduction to the jupyter notebook version of NEMeSIS

The jupyter notebook version of NEMeSIS contains a lot of descriptions in the notebook itself to help guide the user already. Thus, to avoid repetition we will focus on general aspects of it instead of detailing each step like it is done for the graphical interface.

It is split into many different steps from step 0 to 11 which are described in the table of contents. Using this table is the fastest way to navigate the (quite large) notebook.

Step 1 is the equivalent to Stage 1 of the graphical interface: **Data Reading**.

Steps 1.1, 1.2, 1.3 and 2 are the equivalent to Stage 2 of the graphical interface: **Data Pre-Processing and Pre-Treatment**

Steps 3 to 11 are the equivalent of Stage 3 of the graphical interface: **Data Analysis and Biological Interpretation**

Finally, independent side modules are standalone jupyter notebooks

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