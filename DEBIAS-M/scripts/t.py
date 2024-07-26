from biom_df import load_table
import pandas as pd

data = load_table("./data/10317_190329_all.biom")

meta = pd.read_table("./metadata/10317_20240607-135622.txt", sep='\t')


meta_asd = meta[['sample_name', 'asd']]
data_asd = list(data.columns[1:])



meta_asd_list = meta_asd.loc[meta_asd["asd"] == 'Diagnosed by a medical professional (doctor, physician assistant)']
meta_asd_list = list(meta_asd_list["sample_name"])


not_in = []
m_list = list(meta_asd['sample_name'])
for value in data_asd:
    if value not in m_list:
        not_in.append(value)

print("not in metadata frame: ", not_in,
      meta_asd.loc[meta_asd["sample_name"]=="10317.000115111"])

asd = []
for value in data_asd:
    if value in meta_asd_list:
        asd.append(value)

print(asd)
print(len(meta_asd), len(meta_asd_list), len(not_in), len(data_asd))