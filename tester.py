from stepone import StepOne as so
x = so('C:/Users/Lukas/Dropbox/Lab/PCR/Script')
x.get_eds_files()
x.get_sample_ids()
# for file in x.eds_files:
#     print(file)
    # x.extract_ct_values(file, contam_params=('H2O', 35, 10))
    # x.extract_cq_values(file, contam_params=('H2O', 35, 10))
    # x.add_melt_temps(file)
x.extract_multicomponent(x.eds_files[0])
# x.add_melt_temps(x.eds_files[0], linked_operation=(False, True))
# y = so()
# y.extract_ct_values(x.eds_files[1])
# x.df, repeats = x.scrub_contamination(x.df)
# x.df, clean_controls = x.scrub_controls(x.df)
# repeats, dirty_controls = so.scrub_controls(repeats)
# repeats.sort_index(inplace=True)
# repeats[['Detector', 'ID']].to_excel(x.dir+'/repeats.xlsx')


