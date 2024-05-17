import pandas as pd
import InputOutput_class as io

"""
#input
    raw_data_path = '/Users/lucas/16S_amplicon_analysis/phyloseq/final_table.csv'
#output
    Xlsx_output = '/Users/lucas/16S_amplicon_analysis/OUT_16S_Amplicon_Analysis.xlsx'
"""
def main(Xlsx_output:str, raw_data_path:str):
    DF_raw_data = io.ReadFiles.readToDataframe(raw_data_path, None)
    DF_taxonomic_data = DF_raw_data.iloc[:, 1:7] # Extracting the taxonomic data (Kingdom,Phylum,Class,Order,Family,Genus)
    DF_ASV_data = DF_raw_data.iloc[:, 7:] # Extracting the ASV data
    DF_ASV_data = DF_ASV_data.drop(columns='V23')
    for x in range(0, len(DF_ASV_data.columns)):
        DF_ASV_data.rename(columns={DF_ASV_data.columns[x]: 'Abundance of ' + DF_ASV_data.columns[x]}, inplace=True)


    Header = ['Combined abundance', 'Min', 'Max', 'Mean', 'Median', 'Std']
    DF_ASV_data_describe = pd.DataFrame(columns=Header)

    for x in range(0, len(DF_ASV_data)):
        DF_ASV_data_describe.loc[x, 'Combined abundance'] = DF_ASV_data.iloc[x,:].sum()
        DF_ASV_data_describe.loc[x, 'Min'] = DF_ASV_data.iloc[x,:].min()
        DF_ASV_data_describe.loc[x, 'Max'] = DF_ASV_data.iloc[x,:].max()
        DF_ASV_data_describe.loc[x, 'Mean'] = DF_ASV_data.iloc[x,:].mean()
        DF_ASV_data_describe.loc[x, 'Median'] = DF_ASV_data.iloc[x,:].median()
        DF_ASV_data_describe.loc[x, 'Std'] = DF_ASV_data.iloc[x,:].std()
        
    TABLE_raw_data = pd.DataFrame(columns = ['ID','Name'])
    TABLE_raw_data['ID'] = DF_raw_data['ASV']
    TABLE_raw_data['Name'] = DF_raw_data['ASV']
    TABLE_raw_data = pd.concat([TABLE_raw_data, DF_ASV_data_describe], axis=1)
    TABLE_raw_data = pd.concat([TABLE_raw_data, DF_ASV_data], axis=1)
    TABLE_raw_data = pd.concat([TABLE_raw_data, DF_taxonomic_data], axis=1)

    TABLE_raw_data_out = io.WriteFiles.writeToXlsx(Xlsx_output, 'Raw_Data', TABLE_raw_data, 0, 0)
    if TABLE_raw_data_out == True:
        return True
    else:
        return False