from stepone import StepOne as so
import pandas as pd
x = so('C:/Users/Student2/Dropbox/Lab/Experiments/PCR/a3-gcr(ox)-NP_AF/Script')
x.get_sample_ids()
x.dir = 'C:/Users/Student2/Dropbox/Lab/Experiments/PCR/a3-gcr(ox)-NP_AF/eds_files'
x.get_eds_files()
for file in x.eds_files:
    print(file)
    x.extract_cq_values(file, contam_params=('H2O', 35, 10))
    x.add_melt_temps(file)
x.dir = 'C:/Users/Student2/Dropbox/Lab/Experiments/PCR/a3-gcr(ox)-NP_AF/Script'
# x.extract_multicomponent(x.eds_files[0])
# x.df, repeats = x.scrub_contamination(x.df)
x.df, clean_controls = x.scrub_controls(x.df)
# repeats, dirty_controls = so.scrub_controls(repeats)
x.df.sort_index(inplace=True)
idx = pd.IndexSlice
a,b = x.df.index.levels
for id_num in a:
    x.df.loc[idx[id_num, :], 'SampleName'] = x.sample_ids[id_num]
file_want = x.df[['Cq', 'SampleName']]
x.df.to_excel(x.dir+'/debug.xlsx')
file_want.to_csv(x.dir+'/pcrfile.txt', sep='\t', na_rep='40.0')
# repeats[['Detector', 'ID']].to_excel(x.dir+'/repeats.xlsx')


