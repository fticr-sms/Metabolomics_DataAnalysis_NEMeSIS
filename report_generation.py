
# Necessary imports
import docx
from docx.shared import Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
import os
from stat import S_IREAD, S_IRGRP, S_IROTH
import panel as pn
import pandas as pd
import numpy as np

def ReportGenerator(folder, RepGen, file, checkbox_annotation, checkbox_formula, radiobox_neutral_mass, checkbox_others,
                     target_list, UnivarA_Store, characteristics_df, DataFrame_Store, n_databases, DB_dict, verbose_annotated_compounds,
                     data_ann_deduplicator, com_exc_compounds, PCA_params, HCA_params, PLSDA_store, RF_store):
    "Makes a read-only Word file with the metabolomics data analysis performed of selected statistical analysis."

    # Create Folder to put the report in
    if not os.path.exists(folder):
        os.mkdir(folder)

    # Creating the Word Document
    document = docx.Document()

    # Adding the title
    document.add_heading('Generated Data Analysis Report', 0)

    # Adding the disclaimer
    disclaimer = document.add_paragraph('This report is an attempt to describe the metabolomics data analysis performed using ')
    disclaimer.add_run('the MetsTA software. We cannot assure that the descriptions shown here are 100% correct, therefore ')
    disclaimer.add_run('revision of the report to confirm the information shown is advised. Any found inconsistencies are ')
    disclaimer.add_run('encouraged to be reported so the reason can be found and fixed.')

    # Adding the first major section
    document.add_heading('Data Pre-Processing and Pre-Treatment', level=1)

    # Chapter of section
    document.add_heading('Data Reading and Partition', level=2)

    # 1st paragraph - Dataset Read
    p1 = document.add_paragraph('Dataset analysed was present in the file: ')
    p1.add_run(RepGen.filename).bold = True
    p1.add_run('. This dataset had a total of ')

    dataset_cols = len(file.read_df.columns)
    # If a neutral mass column was added to the dataset
    if RepGen.neutral_mass_column:
        nm_string = f' Furthermore, a "Neutral Mass" column was added to the data by interpreting the bucket label given as floats (numbers).'
        dataset_cols = dataset_cols - 1
    else:
        nm_string = ''

    # If the target was in the first row of the dataset read or not
    if RepGen.target_included_in_file:
        dataset_size = len(file.read_df) + 1
        p1.add_run(f'{dataset_size} rows by {dataset_cols} columns including metadata and the first row (has dataset target).')
        p1.add_run(f' That is, the dataset has {dataset_size - 1} detected features.')
        
    else:
        dataset_size = len(file.read_df)
        p1.add_run(f'{dataset_size} rows by {dataset_cols} columns including metadata.')

    # Neutral mass section
    p1.add_run(nm_string)

    # 2nd paragraph - Metdata Columns
    # Putting metaadata columns in string format
    ann_cols = ', '.join(checkbox_annotation.value)
    form_cols = ', '.join(checkbox_formula.value)
    other_cols = ', '.join(checkbox_others.value)

    # Annotation and Formula columns section
    p2 = document.add_paragraph(f'From the {dataset_cols} columns, {len(checkbox_annotation.value)} (')
    p2.add_run(f'{ann_cols}').italic = True
    p2.add_run(f') was/were selected as columns with data annotations; {len(checkbox_formula.value)} (')
    p2.add_run(f'{form_cols}').italic = True
    p2.add_run(f') was/were selected as columns with formula annotations; ')

    # Neutral Mass and Other Metadata columns section
    if radiobox_neutral_mass.value != 'None':
        p2.add_run(f'{radiobox_neutral_mass.value}').italic = True
        p2.add_run(f' was selected as the column representing neutral mass values; and {len(checkbox_others.value)} (')
    else:
        p2.add_run(f'no column was selected to represent neutral mass values; and {len(checkbox_others.value)} (')
    p2.add_run(f'{other_cols}').italic = True
    p2.add_run(f') was/were selected as metadata columns with other information (which were not used).')

    # 3rd paragraph - Target and Sample Columns
    class_str = ', '.join(target_list.classes)

    # Text on the target
    p3 = document.add_paragraph(f'From the {dataset_cols} columns, the remaining {len(target_list.sample_cols)} were samples.')
    p3.add_run(f'These samples belonged to {len(target_list.classes)} classes: ')
    p3.add_run(f'{class_str}').bold = True
    p3.add_run(f'. It had an average of {(len(target_list.sample_cols)/len(target_list.classes)):.2f} samples per class. ')
    p3.add_run(f'The target was the following:')

    # Presenting target as classes
    target_line = document.add_paragraph('')
    target_line.add_run('Target (Class Labels)').bold = True

    for t in range(0, len(target_list.sample_cols), 5):
        table_target = document.add_table(rows=2, cols=6, style='Table Grid')
        sample_cells = table_target.rows[0].cells
        target_cells = table_target.rows[1].cells

        sample_cells[0].text = 'Samples'
        target_cells[0].text = 'Target'
        min_r = min(6, len(target_list.sample_cols)-t+1)
        for col in range(1,min_r):
            sample_cells[col].text = str(target_list.sample_cols[col-1+t])
            target_cells[col].text = str(target_list.target[col-1+t])

        # Give space between tables
        document.add_paragraph()


    # Chapter of section
    document.add_heading('Data Filtering and Pre-Treatment', level=2)

    disclaimer2 = document.add_paragraph('')
    disclaimer2.add_run('Annotation and Annotation De-Duplication (if performed) happenned after data ').bold = True
    disclaimer2.add_run('filtering but before data pre-treatment, that is, metabolic feature merging happenned before data pre-treatment').bold = True
    disclaimer2.add_run(' despite being later in the report.').bold = True

    p4 = document.add_paragraph('Data Filtering and Pre-Treatment parameters used were as following:')

    # Data Filtering
    if UnivarA_Store.filt_method == 'None':
        document.add_paragraph('Filtering Method: None', style='List Bullet'
        )
    else:
        pre_treat_bullet = document.add_paragraph(
            f'Filtering Method: Features that appeared in at least ', style='List Bullet'
        )
        pre_treat_bullet.add_run(f'{UnivarA_Store.filt_kw} samples').bold = True
        if UnivarA_Store.filt_method == 'Total Samples':
            pre_treat_bullet.add_run(' of Total (All) Samples in the dataset')
        elif UnivarA_Store.filt_method == 'Class Samples':
            pre_treat_bullet.add_run(' of at least one Class in the dataset')

    # Missing Value Imputation
    if UnivarA_Store.mvi_method == 'Zero':
        document.add_paragraph(f'Missing Value Imputation Method: Zero Imputation', style='List Bullet'
            )
    else:
        document.add_paragraph(f'Missing Value Imputation Method: {UnivarA_Store.mvi_kw*100}% of {UnivarA_Store.mvi_method}',
                            style='List Bullet')

    # Normalization Method
    norm_paragraph = document.add_paragraph(f'Normalization Method:', style='List Bullet'
            )
    if UnivarA_Store.norm_method == 'Reference Feature':
        norm_paragraph.add_run(f' Normalization by a Reference Feature - {UnivarA_Store.norm_kw}')
    elif UnivarA_Store.norm_method == 'Total Intensity Sum':
        norm_paragraph.add_run(f' Normalization by Total Sum of Intensities')
    elif UnivarA_Store.norm_method == 'PQN':
        if UnivarA_Store.norm_kw in ['mean', 'median']:
            n_kw = 'the ' + str(UnivarA_Store.norm_kw[0]) + ' of all samples as reference'
        else:
            n_kw = str(UnivarA_Store.norm_kw[0]) + ' sample as reference'
        norm_paragraph.add_run(f' PQN (Probabilistic Quotient Normalization) using {n_kw}')
    elif UnivarA_Store.norm_method == 'Quantile':
        norm_paragraph.add_run(f' Quantile Normalization using the {str(UnivarA_Store.norm_kw[0])} of all samples as reference')
    elif UnivarA_Store.norm_method == 'None':
        norm_paragraph.add_run(f' None')
    else:
        pn.state.notifications.error('Normalization method used not recognized.')
        raise ValueError('Normalization method used not recognized for Report Generation.')

    # Transformation Method
    if UnivarA_Store.tf_method == 'Generalized Logarithmic Transformation (glog)':
        if UnivarA_Store.tf_kw:
            document.add_paragraph(
                f'Transformation Method: Generalized Logarithmic Transformation (glog) with lambda equal to {UnivarA_Store.tf_kw}',
                                style='List Bullet'
                )
        else:
            document.add_paragraph(
                f'Transformation Method: Generalized Logarithmic Transformation (glog) with lambda equal to 0',
                                style='List Bullet'
                )
    elif not UnivarA_Store.tf_method:
        document.add_paragraph(f'Transformation Method: None', style='List Bullet'
                )
    else:
        pn.state.notifications.error('Transformation method used not recognized.')
        raise ValueError('Transformation method used not recognized for Report Generation.')

    # Scaling Method
    if UnivarA_Store.scaling_method == 'Level Scaling':
        document.add_paragraph(
                f"Scaling Method: Level Scaling (using feature's {UnivarA_Store.scaling_kw} as scaling factors)",
                                style='List Bullet'
                )
    else:
        document.add_paragraph(f"Scaling Method: {UnivarA_Store.scaling_method}", style='List Bullet'
                )

    # Dataset extra characteristics
    range_min, range_max = characteristics_df.value.loc['feature value ranges'][0].split(' - ')
    p5 = document.add_paragraph(f"After pre-treatment, the dataset had {len(DataFrame_Store.treated_df.columns)} features. ")
    p5.add_run(f'Absolute intensity values (without pre-treatment) ranged between {range_min[1:]} to {range_max[:-1]}, ')
    p5.add_run(f'averaging at {characteristics_df.value.loc["feature value average (std)"][0]}.')
        
    document.add_paragraph(f"Treated DataFrame was downloaded as 'Report_treated_df.csv'.")


    # Chapter of section
    document.add_heading('Data Annotation and De-Duplication', level=2)
    # TODO: Add PPM DEVIATION / ABSOLUTE DALTON DEVIATION used to the report

    ann_pg = document.add_paragraph(f'Data Annotation in this software was made using {n_databases.value} database(s).')

    if n_databases.value != 0:
        # Database Annotation based on PPM Deviation
        if RepGen.annotation_margin_method == 'PPM Deviation':
            ann_pg.add_run(f' Database matching was performed using a maximum of {RepGen.annotation_margin_ppm_deviation} ppm')
            ann_pg.add_run(f' (parts per million) deviation of the theoretical compound mass (based on formula) to the ')
            ann_pg.add_run(f'experimental mass in the Neutral Mass Column ({radiobox_neutral_mass.value}) chosen.')

        # Database Annotation based on flat Da deviation
        else:
            ann_pg.add_run(f' Database matching was performed using a maximum of {RepGen.annotation_margin_ppm_deviation} Da')
            ann_pg.add_run(f' deviation of the theoretical compound mass (based on formula) to the ')
            ann_pg.add_run(f'experimental mass in the Neutral Mass Column ({radiobox_neutral_mass.value}) chosen.')

    # Information on Databases used
    for db in range(n_databases.value):
        key = str(db + 1)
        db_info = DB_dict[key]
        document.add_paragraph(f"{db_info.abv}: Database read from '{db_info.file}' using {db_info.IDcol} column as IDs. It had {len(db_info.db.value)} metabolites.",
                            style='List Bullet')

    # Information on the annotations made with each Database used
    if n_databases.value != 0:
        document.add_paragraph(f'Annotation with these databases led to the annotation of :')
        for db in range(n_databases.value):
            db_info = DB_dict[str(db + 1)]
            document.add_paragraph(f'{db_info.abv}: {verbose_annotated_compounds[db+1].value.split( )[1]} compounds',
                                style='List Bullet')

    # Annotation De-Duplication section
    p6 = document.add_paragraph(f'After the annotations performed, the annotations (including previous annotations made) were ')
    p6.add_run('checked. The following table was obtained from observing the multiplicity of the annotations.')

    # Table with the report of the multiple annotations
    table_ann = document.add_table(rows=4, cols=len(data_ann_deduplicator.mult_ann_report.columns)+1, style='Light Grid')

    def fill_word_table(table, df):
        "Function to convert a pandas DataFrame to a table in the Word Report."
        
        # Get the Index column and fill it
        hdr_cells = table.columns[0].cells
        idx_pos = 1
        for idx in df.index:
            hdr_cells[idx_pos].text = idx
            idx_pos += 1

        # Go through every other column and fill it
        n_col = 1
        for col in df.columns:
            # Get the col
            new_col = table.columns[n_col].cells
            n_col += 1
            new_col[0].text = col # Name of the column
            a = 1
            # Fill remaining values
            for value in df[col]:
                new_col[a].text = str(value)
                a += 1
        table.autofit = True
        return table

    table_ann = fill_word_table(table_ann, data_ann_deduplicator.mult_ann_report)

    # Give space after table
    document.add_paragraph()

    # Annotation De-duplication section
    if RepGen.deduplication_performed == 'Annotation De-Duplication was not performed':
        document.add_paragraph(f'Multiple Annotation De-Duplication and Metabolic Feature merging was not performed.')

    else:
        # If annotation De-Duplication was made

        # Type of metabolic merging performed
        if RepGen.deduplication_performed == 'Annotation De-Duplication was performed but merge problems were not individually observed/decided':
            p7 = document.add_paragraph(f'Multiple Annotation De-Duplication and Metabolic Feature merging was performed but')
            p7.add_run(' possible merge problems were not individually observed. ')
        else:
            p7 = document.add_paragraph(f'Multiple Annotation De-Duplication and Metabolic Feature merging was performed and')
            p7.add_run(' merge problems were individually observed and decided whether they were merged. ')

        # If Formula columns were considered
        if data_ann_deduplicator.current_params['consider_formula_cols']:
            p7.add_run('Previously Formula annotated columns were also considered in metabolic feature merging. ')
        else:
            p7.add_run('Previously Formula annotated columns were not considered in metabolic feature merging. ')

        # If Scenario 1 of Case 4 of mergings was considered as potential problem
        if data_ann_deduplicator.current_params['problem_condition'] == 'Scenario 1 of Situation 4 like cases are not merged and are not shown.':
            p7.add_run('Furthermore, Scenario 1 of Situation 4 like cases (see description in MetsTA software for further ')
            p7.add_run('details) were not merged and are not shown (not considered possible problems).')
        else:
            p7.add_run('Furthermore, Scenario 1 of Situation 4 like cases (see description in MetsTA software for further ')
            p7.add_run('details) are shown as possible problems to individually decide.')

        # Describing Mergings
        target_line = document.add_paragraph('')
        target_line.add_run('Describing Metabolic Feature Mergings').bold = True

        # Table with descriptions of types of merging made
        table_merge_situations = document.add_table(rows=len(data_ann_deduplicator.merge_situations.index)+1,
                                                    cols=len(data_ann_deduplicator.merge_situations.columns)+1,
                                                    style='Light Grid')
        table_merge_situations = fill_word_table(table_merge_situations, data_ann_deduplicator.merge_situations)

        # Give space after table
        document.add_paragraph()

        # Merge Report - Treating the string first
        m_report = data_ann_deduplicator.merge_report
        m_report = m_report.replace('**', '')
        m_report = m_report.replace('  ', '')
        document.add_paragraph(m_report)

        # Saving merge description as Excel
        data_ann_deduplicator.merge_description.to_excel(folder+'/Report_Metabolic_Feature_Peak_Merging_Description.xlsx')
        document.add_paragraph("Detailed information of metabolic feature peak merging made saved in 'Report_Metabolic_Feature_Peak_Merging_Description.xlsx'.")


    # Final Description
    p8 = document.add_paragraph("")
    p8.add_run("Finished Data Pre-Treatment").bold = True
    p8.add_run(f". Dataset had {len(target_list.sample_cols)} samples and {len(DataFrame_Store.treated_df.columns)} metabolic features,")
    p8.add_run(f" from which {sum(DataFrame_Store.metadata_df['Has Match?'])} are annotated.")
    p8.add_run(f" Dataset after annotation and metabolic feature merging (if performed) were saved. 'Report_Complete_non_Treated_df.csv'")
    p8.add_run(f" has original intensity values and metadata. 'Report_Complete_Treated_df.csv' has treated intensity values and ")
    p8.add_run(f"metadata. 'Report_Treated_df.csv' only has the treated intensity values.")

    # Saving the respective DataFrames
    DataFrame_Store.original_df.to_csv(folder+'/Report_Complete_non_Treated_df.csv')
    DataFrame_Store.treated_df.to_csv(folder+'/Report_Treated_df.csv')
    pd.concat((DataFrame_Store.metadata_df, DataFrame_Store.treated_df.T), axis=1).to_csv(folder+'/Report_Complete_Treated_df.csv')

    # End of section
    document.add_page_break()




    # Adding the second major section
    document.add_heading('Statistical Analysis', level=1)

    # Getting the Statistical Methods that will be in 
    stat_methods = []
    for i in RepGen.com_exc_analysis:
        if i == 'Overview':
            stat_methods.append('Common and Exclusive Compound Overview')
        else:
            stat_methods.append(i)
    for i in RepGen.unsup_analysis:
        stat_methods.append(i)
    for i in RepGen.sup_analysis:
        stat_methods.append(i)
    if RepGen.univ_analysis:
        stat_methods.append('Univariate Analysis')
    for i in RepGen.dataviz_analysis:
        stat_methods.append(i)
    if RepGen.pathassign_analysis:
        stat_methods.append('HMDB Pathway Assignment')
    for i in RepGen.BinSim_analysis:
        stat_methods.append('BinSim - ' + i)

    document.add_paragraph(f'This section includes the statistical analysis chosen to be in this report. These are: {", ".join(stat_methods)}.')

    # Common and Exclusive Compound Section
    if len(np.intersect1d(['Common and Exclusive Compound Overview', 'Venn Diagram', 'Intersection Plot'], stat_methods)) != 0:
        document.add_heading('Common and Exclusive Compound Analysis', level=2)

        # Overview Section
        if 'Common and Exclusive Compound Overview' in stat_methods:
            # Heading
            document.add_heading('Common and Exclusive Compound Overview', level=3)

            # In case analysis was not performed
            if len(com_exc_compounds.com_exc_desc) == 0:
                document.add_paragraph('Common and Exclusive Compound Analysis was not performed. Thus, this section will be skipped.')

            else:
                # Descriptions of peaks per class and exclusive peaks per class
                document.add_paragraph(com_exc_compounds.groups_description.replace('<br />', '\n').replace('**', ''))
                document.add_paragraph()
                document.add_paragraph(com_exc_compounds.com_exc_desc.replace('<br />', '\n').replace('**', ''))

        # Venn Diagram Section
        if 'Venn Diagram' in stat_methods:
            # Heading
            document.add_heading('Venn Diagram', level=3)

            # If there is a Venn Diagram
            if type(com_exc_compounds.Venn_plot[0]) != str:
                # Description
                venn_pg = document.add_paragraph('Venn Diagram downloaded as currently present in the program. Venn Diagram ')
                venn_pg.add_run(f"shown includes the following classes: ")
                venn_pg.add_run(f"{', '.join(com_exc_compounds.venn_class_subset)}.").bold = True

                # Creating filename for the Venn Figure and saving it
                filename_string = f'/Report_Venn_diagram_{com_exc_compounds.type_of_venn}_classes'
                for cl in com_exc_compounds.venn_class_subset:
                    filename_string = filename_string + '_'+cl
                com_exc_compounds.Venn_plot[0].savefig(folder+filename_string, dpi=com_exc_compounds.dpi_venn)

                # Adding figure
                document.add_picture(folder+filename_string+'.png', width=Cm(15))

            # If there is no Venn Diagram
            else:
                venn_pg = document.add_paragraph('Venn Diagram was not created during the Data Analysis. ')
                venn_pg.add_run(f"Thus, this section will be skipped.")

        # Intersection Plot Section
        if 'Intersection Plot' in stat_methods:
            # Heading
            document.add_heading('Intersection Plots', level=3)

            # If there are Intersection Plots
            if type(com_exc_compounds.IntersectionPlot[0]) != str:
                # Description
                inter_pg = document.add_paragraph('Intersection Plots downloaded as currently present in the program. ')
                inter_pg.add_run(f"Intersection Plots shown include the following classes: ")
                inter_pg.add_run(f"{', '.join(com_exc_compounds.inter_class_subset)}.").bold = True

                # Creating filename for the Intersection Plots and saving them
                filename_string_all_mets = f'/Report_IntersectionPlot_all_metabolites_classes'
                filename_string_ann_mets = f'/Report_IntersectionPlot_annotated_metabolites_classes'
                for cl in com_exc_compounds.inter_class_subset:
                    filename_string_all_mets = filename_string_all_mets + '_'+cl
                    filename_string_ann_mets = filename_string_ann_mets + '_'+cl
                com_exc_compounds.IntersectionPlot[0].savefig(folder + filename_string_all_mets + '.png',
                                                            dpi=com_exc_compounds.dpi_inter)
                if type(com_exc_compounds.IntersectionPlot[1]) != str:
                    com_exc_compounds.IntersectionPlot[1].savefig(folder + filename_string_ann_mets + '.png',
                                                                dpi=com_exc_compounds.dpi_inter)

                # Adding description and figure - all metabolites
                inter_pg2 = document.add_paragraph()
                inter_pg2.add_run('Intersection Plot with all metabolites of the dataset').bold = True
                inter_pg2.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
                document.add_picture(folder+filename_string_all_mets + '.png', width=Cm(14))

                # Adding description and figure - annotated metabolites
                if type(com_exc_compounds.IntersectionPlot[1]) != str:
                    inter_pg3 = document.add_paragraph()
                    inter_pg3.add_run('Intersection Plot with only annotated metabolites in the dataset').bold = True
                    inter_pg3.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
                    document.add_picture(folder+filename_string_ann_mets + '.png', width=Cm(14))
                else:
                    inter_pg3 = document.add_paragraph('Intersection Plot with only annotated metabolites in the dataset')
                    inter_pg3.add_run(' was not created since annotations were not found in the dataset.')

            # If there is no Intersection Plots
            else:
                inter_pg = document.add_paragraph('Intersection Plots were not created during the Data Analysis. ')
                inter_pg.add_run(f"Thus, this section will be skipped.")

        # End of section
        document.add_page_break()


    # Unsupervised Analysis Section
    if len(np.intersect1d(['PCA', 'HCA'], stat_methods)) != 0:
        document.add_heading('Unsupervised Analysis', level=2)

        # PCA Section
        if 'PCA' in stat_methods:
            # Heading
            document.add_heading('Principal Component Analysis (PCA)', level=3)

            # If there is a PCA Projection
            if type(PCA_params.PCA_plot[0]) != str:
                # Description
                pca_pg = document.add_paragraph(f'Principal Component Analysis was made with ')
                pca_pg.add_run(f'{len(PCA_params.controls.widgets["PCx"].options)} components.')

                pca_pg2 = document.add_paragraph()
                if PCA_params.n_dimensions == '2 Components':
                    pca_pg2.add_run(f'2-D PCA Projection Plot of Components: {PCA_params.PCx} and {PCA_params.PCy}.').bold = True
                else:
                    pca_pg2.add_run(f'3-D PCA Projection Plot of Components: {PCA_params.PCx}, {PCA_params.PCy} and {PCA_params.PCz}.').bold = True
                pca_pg2.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER

                # Creating filename for the PCA Projection and saving it
                filename_string = '/Report_PCA_plot'
                if PCA_params.ellipse_draw:
                    if PCA_params.confidence != 0:
                        filename_string = filename_string + f'_ellipse({PCA_params.confidence*100}%confidence)'
                    else:
                        filename_string = filename_string + f'_ellipse({PCA_params.confidence_std}std)'
                PCA_params.PCA_plot[0].write_image(folder+filename_string+'.png', scale=4)
                PCA_params.PCA_plot[0].write_html(folder+filename_string+".html")

                # Adding figure
                document.add_picture(folder+filename_string+'.png', width=Cm(12.5))

                # If there is an explained variance plot
                if type(PCA_params.PCA_plot[0]) != str:
                    # Description
                    pca_pg = document.add_paragraph()
                    pca_pg.add_run('Cumulative Explained Variance (by Principal Component) Plot').bold = True
                    pca_pg.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER

                    # Creating filename for the PCA explained variance plot and saving it
                    filename_string = '/Report_PCA_exp_var_plot'
                    PCA_params.exp_var_fig_plot[0].write_image(folder+filename_string+'.png', scale=4)
                    PCA_params.exp_var_fig_plot[0].write_html(folder+filename_string+".html")

                    # Adding figure
                    document.add_picture(folder+filename_string+'.png', width=Cm(12.5))

                # If there is a PCA Projection (Scatter) Plot of the Principal Components
                if type(PCA_params.scatter_PCA_plot[0]) != str:
                    # Description
                    pca_pg = document.add_paragraph()
                    pca_pg.add_run('2-D PCA Projection (Scatter) Plot of the Principal Components').bold = True
                    pca_pg.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER

                    # Creating filename for the PCA Projection (Scatter) Plot of the Principal Components and saving it
                    filename_string = '/Report_PCA_scatter_plot'
                    PCA_params.scatter_PCA_plot[0].write_image(folder+filename_string+'.png', scale=4)
                    PCA_params.scatter_PCA_plot[0].write_html(folder+filename_string+".html")

                    # Adding figure
                    document.add_picture(folder+filename_string+'.png', width=Cm(15))

            # If there is no Venn Diagram
            else:
                pca_pg = document.add_paragraph('Principal Component Analysis was not performed during the Data Analysis. ')
                pca_pg.add_run(f"Thus, this section will be skipped.")


        # HCA Section
        if 'HCA' in stat_methods:
            # Heading
            document.add_heading('Hierarchical Clustering Analysis (PCA)', level=3)

            # If there is a dendrogram
            if type(HCA_params.HCA_plot[0]) != str:
                # Description
                hca_pg = document.add_paragraph(f'Hierarchical Clustering Analysis was made with ')
                hca_pg.add_run(f'{HCA_params.dist_metric} distance metric and {HCA_params.link_metric} linkage metric.')

                # Creating filename for the PCA Projection and saving it
                filename = f'/Report_HCA_plot_{HCA_params.dist_metric}Dist_{HCA_params.link_metric}Linkage'
                HCA_params.HCA_plot[0].savefig(folder + filename_string + '.png', dpi=HCA_params.dpi)

                # Adding figure
                document.add_picture(folder+filename_string+'.png', width=Cm(12))

        # End of section
        document.add_page_break()



    # Supervised Analysis Section
    if len(np.intersect1d(['PLS-DA', 'Random Forest'], stat_methods)) != 0:
        document.add_heading('Supervised Analysis', level=2)

        # PLS-DA Section
        if 'PLS-DA' in stat_methods:
            # Heading
            document.add_heading('Partial Least Squares - Discriminant Analysis (PLS-DA)', level=3)

            # If there is a PLS optimization
            if type(PLSDA_store.optim_figure[0]) != str:
                # Minor Heading
                document.add_heading('Number of Components Optimization', level=4)

                # Description
                n_min_components, n_max_components = PLSDA_store.current_other_plsda_params["n_min_max_components"]
                plsda_pg = document.add_paragraph(f'First, an optimization of the number of components to use to fit a ')
                plsda_pg.add_run(f'PLS-DA model was made by fitting PLS models from {n_min_components}')
                plsda_pg.add_run(f' to {n_max_components} components and evaluating their Q2 (mainly) and R2 scores estimated')
                plsda_pg.add_run(f' by {PLSDA_store.current_other_plsda_params["n_fold_optim"]}-fold stratified')
                plsda_pg.add_run(f' cross-validation (scale = {PLSDA_store.current_other_plsda_params["scale_optim"]}). Results')
                plsda_pg.add_run(f' are shown in the figure below with a maximum Q2 value with {PLSDA_store.rec_components} ')
                plsda_pg.add_run(f'components.')

                # Plot title
                plsda_pg2 = document.add_paragraph()
                plsda_pg2.add_run('PLS-DA Component Optimization Plot').bold = True
                plsda_pg2.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER

                # Creating filename for the PLS optimization plot and saving it
                filename_string = f'/Report_PLS_optim_plot_{n_min_components}to{n_max_components}'
                filename_string = filename_string + f'components_{PLSDA_store.current_other_plsda_params["n_fold_optim"]}-'
                filename_string = filename_string + f'foldstratCV_scale{PLSDA_store.current_other_plsda_params["scale_optim"]}'
                PLSDA_store.optim_figure[0].write_image(folder+filename_string+'.png', scale=4)
                PLSDA_store.optim_figure[0].write_html(folder+filename_string+".html")

                # Adding figure
                document.add_picture(folder+filename_string+'.png', width=Cm(15))

            else:
                document.add_paragraph(f'Optimization of the number of components to fit the PLS-DA model was not performed.')


            # If the PLS-DA model was fitted
            if type(PLSDA_store.models[0]) != str:
                # Minor Heading
                document.add_heading('PLS-DA Model Fitting and Results', level=4)

                # Description
                plsda_params = PLSDA_store.current_plsda_params
                plsda_pg3 = document.add_paragraph(f'PLS-DA model was fitted using {plsda_params["n_components"]} components')
                plsda_pg3.add_run(f' estimating model performance and feature importance ({plsda_params["feat_imp"]}')
                plsda_pg3.add_run(f' score method) by {plsda_params["n_folds"]}-fold stratified cross-validation (scale = ')
                plsda_pg3.add_run(f'{plsda_params["scale"]}) repeated {plsda_params["n_iterations"]} times (randomized folds ')
                plsda_pg3.add_run(f'in cross-validation). Model Performance results are shown in the table below (by the ')

                # Building the filename and saving PLS-DA Feature Importance Table
                filename_string = f'/Report_PLS-DA_FeatImp_{plsda_params["feat_imp"]}_model_params_components'
                filename_string = filename_string + f'{plsda_params["n_components"]}_{plsda_params["n_folds"]}-foldstratCV_'
                filename_string = filename_string + f'iterations{plsda_params["n_iterations"]}_scale{plsda_params["scale"]}.xlsx'
                PLSDA_store.feat_impor.to_excel(folder+filename_string)

                # Finishing up the description
                plsda_pg3.add_run(f"metrics chosen) and Feature Importance Table was saved as '{filename_string[1:]}'. ")
                plsda_pg3.add_run(f'Furthermore, below we also show a PLS Projection of chosen Latent Variables.')

                # Table with results of PLS-DA
                table_plsda_results = document.add_table(rows=len(PLSDA_store.n_results.index)+1,
                                                            cols=len(PLSDA_store.n_results.columns)+1,
                                                            style='Light Grid')
                table_plsda_results = fill_word_table(table_plsda_results, PLSDA_store.n_results)

                # Give space after table
                document.add_paragraph()

                plsda_pg4 = document.add_paragraph()
                if PLSDA_store.n_dimensions == '2 Components':
                    plsda_pg4.add_run(f'2-D PLS Projection Plot of Latent Variables: {PLSDA_store.LVx} and {PLSDA_store.LVy}.').bold = True
                else:
                    plsda_pg4.add_run(f'3-D PCA Projection Plot of Latent Variables: {PLSDA_store.LVx}, {PLSDA_store.LVy}').bold = True
                    plsda_pg4.add_run(f' and {PLSDA_store.LVz}.').bold = True
                plsda_pg4.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER

                # Creating filename for the PLS Projection and saving it
                filename_string = '/Report_PLS_projection_plot'
                if PLSDA_store.n_dimensions == '2 Components':
                    if PLSDA_store.ellipse_draw:
                        if PLSDA_store.confidence != 0:
                            filename_string = filename_string + f'_ellipse({PLSDA_store.confidence*100}%confidence)'
                        else:
                            filename_string = filename_string + f'_ellipse({PLSDA_store.confidence_std}std)'
                PLSDA_store.PLS_plot[0].write_image(folder+filename_string+'.png', scale=4)
                PLSDA_store.PLS_plot[0].write_html(folder+filename_string+".html")

                # Adding figure
                document.add_picture(folder+filename_string+'.png', width=Cm(15))


                # If there is Permutation Test plot
                if type(PLSDA_store.perm_figure[0]) != str:
                    # Minor Heading
                    document.add_heading('PLS-DA Permutation Test', level=4)

                    # Description
                    plsda_pg5 = document.add_paragraph('Permutation Test was performed with the same number of components as before - ')
                    plsda_pg5.add_run(f"{PLSDA_store.current_plsda_params_permutation['n_components']} - and with ")
                    plsda_pg5.add_run(f"{PLSDA_store.current_plsda_params_permutation['n_permutations']} permutations. Model ")
                    plsda_pg5.add_run(f"performance was evaluated by {PLSDA_store.current_plsda_params_permutation['perm_metric']}")
                    plsda_pg5.add_run(f" estimated with {PLSDA_store.current_plsda_params_permutation['n_folds']}-fold stratified ")
                    plsda_pg5.add_run(f"cross-validation (scale = {PLSDA_store.current_plsda_params_permutation['scale']}).")

                    # Creating filename for the PLS-DA Permutation Test and saving it
                    filename_string = f'/Report_PLS-DA_permutation_test_{PLSDA_store.current_plsda_params_permutation["n_permutations"]}perm_'
                    filename_string = filename_string + f'{PLSDA_store.current_plsda_params_permutation["n_components"]}comp_'
                    filename_string = filename_string + f'{PLSDA_store.current_plsda_params_permutation["n_folds"]}-foldstratCV_scale'
                    filename_string = filename_string + f'{PLSDA_store.current_plsda_params_permutation["scale"]}_metric'
                    filename_string = filename_string + f'{PLSDA_store.current_plsda_params_permutation["perm_metric"]}'
                    PLSDA_store.perm_figure[0].savefig(folder+filename_string+'.png', dpi=PLSDA_store.dpi)

                    # Adding figure
                    document.add_picture(folder+filename_string+'.png', width=Cm(15))


                # If there is ROC Curve plot
                if type(PLSDA_store.ROC_figure[0]) != str:
                    # Minor Heading
                    document.add_heading('PLS-DA ROC Curve', level=4)

                    # Description
                    roc_params = PLSDA_store.current_other_plsda_params['ROC_filename'].split('_')
                    if len(target_list.classes) > 2:
                        plsda_pg6 = document.add_paragraph('Since there are more than 2 classes, ROC curves were computed with a ')
                        plsda_pg6.add_run(f"1vsAll scheme for each of the {len(target_list.classes)} classes.")
                    else:
                        plsda_pg6 = document.add_paragraph('Since there are only 2 classes, ROC curves were computed considering ')
                        plsda_pg6.add_run(f"{roc_params[2][:-8]}").bold = True
                        plsda_pg6.add_run(f" as the positive class.")
                    plsda_pg6.add_run(f' Other parameters were maintained: number of components - {roc_params[3][:-10]}; ')
                    plsda_pg6.add_run(f"{roc_params[4][:-7]} stratified cross-validation; scale - {roc_params[6][5:]}. ")
                    plsda_pg6.add_run(f"Finally, ROC Curve estimation was repeated {roc_params[5][:-10]} times (randomized ")
                    plsda_pg6.add_run(f"cross-validation folds).")

                    # Creating filename for the PLS-DA ROC Curve and saving it
                    filename_string = '/Report_' + PLSDA_store.current_other_plsda_params['ROC_filename']
                    PLSDA_store.ROC_figure[0].write_image(folder+filename_string+'.png', scale=4)
                    PLSDA_store.ROC_figure[0].write_html(folder+filename_string+".html")

                    # Adding figure
                    document.add_picture(folder+filename_string+'.png', width=Cm(15))

            # If there is no PLS-DA model fitted
            else:
                plsda_pg3 = document.add_paragraph('PLS-DA model was not fitted during the Data Analysis. ')
                plsda_pg3.add_run(f"Thus, this section will be skipped.")


        # Random Forest Section
        if 'Random Forest' in stat_methods:
            # Heading
            document.add_heading('Random Forest (RF)', level=3)

            # If there is a Random Forest Tree Number optimization
            if type(RF_store.optim_figure[0]) != str:
                # Minor Heading
                document.add_heading('Number of Trees Optimization', level=4)

                # Description
                n_min_trees, n_max_trees = RF_store.current_other_rf_params["n_min_max_trees"]
                rf_pg = document.add_paragraph(f'First, an optimization of the number of trees to use to fit a ')
                rf_pg.add_run(f'Random Forest model was made by fitting RF models from {n_min_trees} to {n_max_trees} trees in ')
                rf_pg.add_run(f'{RF_store.current_other_rf_params["n_interval"]} tree steps, evaluating their model accuracy')
                rf_pg.add_run(f' estimated by {RF_store.current_other_rf_params["n_fold_optim"]}-fold stratified')
                rf_pg.add_run(f' cross-validation. Results are shown in the figure below.')

                # Plot title
                rf_pg2 = document.add_paragraph()
                rf_pg2.add_run('Random Forest Number of Trees Optimization Plot').bold = True
                rf_pg2.paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER

                # Creating filename for the Random Forest optimization plot and saving it
                filename_string = f'/Report_RF_optim_plot_{RF_store.current_other_rf_params["n_fold_optim"]}-foldStratCV_'
                filename_string = filename_string + f'{n_min_trees}to{n_max_trees}trees_'
                filename_string = filename_string + f'({RF_store.current_other_rf_params["n_interval"]}interval)'
                RF_store.optim_figure[0].write_image(folder+filename_string+'.png', scale=4)
                RF_store.optim_figure[0].write_html(folder+filename_string+".html")

                # Adding figure
                document.add_picture(folder+filename_string+'.png', width=Cm(15))

            else:
                document.add_paragraph(f'Optimization of the number of trees to fit the Random Forest model was not performed.')


            # If the Random Forest model was fitted
            if type(RF_store.models[0]) != str:
                # Minor Heading
                document.add_heading('Random Forest Model Fitting and Results', level=4)

                # Description
                rf_params = RF_store.current_rf_params
                rf_pg3 = document.add_paragraph(f'Random Forest model was fitted with {rf_params["n_trees"]}')
                rf_pg3.add_run(f' trees estimating model performance and feature importance by {rf_params["n_folds"]}')
                rf_pg3.add_run(f'-fold stratified cross-validation repeated {rf_params["n_iterations"]} times (randomized ')
                rf_pg3.add_run(f'folds in cross-validation). Model Performance results are shown in the table below (by the ')

                # Building the filename and saving Random Forest Feature Importance Table
                filename_string = f'/Report_RF_FeatImp_Gini_model_params_{rf_params["n_trees"]}trees_{rf_params["n_folds"]}'
                filename_string = filename_string + f'-foldstratCV_iterations{rf_params["n_iterations"]}.xlsx'
                RF_store.feat_impor.to_excel(folder+filename_string)

                # Finishing up the description
                rf_pg3.add_run(f"metrics chosen) and Feature Importance Table was saved as '{filename_string[1:]}'.")

                # Table with results of Random Forest
                table_rf_results = document.add_table(rows=len(RF_store.n_results.index)+1,
                                                            cols=len(RF_store.n_results.columns)+1,
                                                            style='Light Grid')
                table_rf_results = fill_word_table(table_rf_results, RF_store.n_results)

                # Give space after table
                document.add_paragraph()


                # If there is Permutation Test plot
                if type(RF_store.perm_figure[0]) != str:
                    # Minor Heading
                    document.add_heading('Random Forest Permutation Test', level=4)

                    # Description
                    rf_params_perm = RF_store.current_rf_params_permutation
                    rf_pg4 = document.add_paragraph('Permutation Test was performed with the same number of trees as before - ')
                    rf_pg4.add_run(f"{rf_params_perm['n_trees']} - and with {rf_params_perm['n_permutations']} permutations. ")
                    rf_pg4.add_run(f"Model performance was evaluated by {rf_params_perm['perm_metric']} estimated with ")
                    rf_pg4.add_run(f"{rf_params_perm['n_folds']}-fold stratified cross-validation.")

                    # Creating filename for the Random Forest Permutation Test and saving it
                    filename_string = f'/Report_RF_permutation_test_{rf_params_perm["n_permutations"]}perm_'
                    filename_string = filename_string + f'{rf_params_perm["n_trees"]}comp_'
                    filename_string = filename_string + f'{rf_params_perm["n_folds"]}-foldstratCV_metric'
                    filename_string = filename_string + f'{rf_params_perm["perm_metric"]}'
                    RF_store.perm_figure[0].savefig(folder+filename_string+'.png', dpi=RF_store.dpi)

                    # Adding figure
                    document.add_picture(folder+filename_string+'.png', width=Cm(15))


                # If there is ROC Curve plot
                if type(RF_store.ROC_figure[0]) != str:
                    # Minor Heading
                    document.add_heading('Random Forest ROC Curve', level=4)

                    # Description
                    roc_params = RF_store.current_other_rf_params['ROC_filename'].split('_')
                    if len(target_list.classes) > 2:
                        rf_pg5 = document.add_paragraph('Since there are more than 2 classes, ROC curves were computed with a ')
                        rf_pg5.add_run(f"1vsAll scheme for each of the {len(target_list.classes)} classes.")
                    else:
                        rf_pg5 = document.add_paragraph('Since there are only 2 classes, ROC curves were computed considering ')
                        rf_pg5.add_run(f"{roc_params[2][:-8]}").bold = True
                        rf_pg5.add_run(f" as the positive class.")
                    rf_pg5.add_run(f' Other parameters were maintained: number of trees - {roc_params[3][:-5]}; ')
                    rf_pg5.add_run(f"{roc_params[4][:-7]} stratified cross-validation. Finally, ROC Curve estimation was ")
                    rf_pg5.add_run(f"repeated {roc_params[5][:-10]} times (randomized cross-validation folds).")

                    # Creating filename for the Random Forest ROC Curve and saving it
                    filename_string = '/Report_' + RF_store.current_other_rf_params['ROC_filename']
                    RF_store.ROC_figure[0].write_image(folder+filename_string+'.png', scale=4)
                    RF_store.ROC_figure[0].write_html(folder+filename_string+".html")

                    # Adding figure
                    document.add_picture(folder+filename_string+'.png', width=Cm(15))

            # If there is no Random Forest model fitted
            else:
                rf_pg3 = document.add_paragraph('Random Forest model was not fitted during the Data Analysis. ')
                rf_pg3.add_run(f"Thus, this section will be skipped.")

        # End of section
        document.add_page_break()




    document.save(folder+'/Report.docx')
    os.chmod(folder+'/Report.docx', S_IREAD|S_IRGRP|S_IROTH)