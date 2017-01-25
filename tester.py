from stepone import StepOne as so
x = so('C:/Users/ljaworski/Dropbox/Lab/PCR/Script')
x.get_eds_files()
x.get_sample_ids()
# for file in x.eds_files:
#     print(file)
#     x.extract_ct_values(file, contam_params=('H2O', 35, 10), f_contam=False)
    # x.extract_ct_values(file, contam_params=('H2O', 35, 10), linked_operation=(True, False))
    # x.add_melt_temps(file, linked_operation=(False, True))
x.extract_ct_values(x.eds_files[0])
# x.add_melt_temps(x.eds_files[0], linked_operation=(False, True))
# y = so()
# y.extract_ct_values(x.eds_files[1])
x.df, repeats = x.scrub_contamination(x.df)
x.df, clean_controls = x.scrub_controls(x.df)
repeats, dirty_controls = so.scrub_controls(repeats)
repeats.sort_values(by=['Detector', 'Sample'], inplace=True)
# repeats[['Detector', 'ID']].to_excel(x.dir+'/repeats.xlsx')

# import openpyxl
# import pandas as pd
#
# wb = openpyxl.load_workbook('C:/Users/Student2/Dropbox/Lab/PCR/Control Genes/PCR Analysis for Genorme.xlsx')
# df = pd.DataFrame()
# genes = ['RPL4', 'HPRT I', 'ACTB', 'YWHAZ', 'TBP', 'HMBS', '18s', 'GAPDH']
# for ws in wb:
#     for col in ws.iter_cols():
#         col2 = [val.value for val in col]
#         imp_col_val = col2[1:10], col2[11:20], col2[21:30]
#         if isinstance(imp_col_val[0][0], int):
#             temp_df = pd.DataFrame({'Detector': genes, 'Cq': imp_col_val[0][1:]})
#             temp_df['Sample'] = imp_col_val[0][0]
#             temp_df['Exp_Grp'] = ws.title
#             temp_df['Day'] = 1
#             df = df.append(temp_df, ignore_index=True)
#         if isinstance(imp_col_val[1][0], int):
#             temp_df = pd.DataFrame({'Detector': genes, 'Cq': imp_col_val[1][1:]})
#             temp_df['Sample'] = imp_col_val[1][0]
#             temp_df['Exp_Grp'] = ws.title
#             temp_df['Day'] = 5
#             df = df.append(temp_df, ignore_index=True)
#         if isinstance(imp_col_val[2][0], int):
#             temp_df = pd.DataFrame({'Detector': genes, 'Cq': imp_col_val[2][1:]})
#             temp_df['Sample'] = imp_col_val[2][0]
#             temp_df['Exp_Grp'] = ws.title
#             temp_df['Day'] = 10
#             df = df.append(temp_df, ignore_index=True)
# df.to_excel('C:/Users/Student2/Dropbox/Lab/PCR/Control Genes/GeneNorm_df.xlsx')
# print(df)
