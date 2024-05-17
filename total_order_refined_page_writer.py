import pandas as pd
import numpy as np
import InputOutput_class as io

from  matplotlib.colors import LinearSegmentedColormap  # genera il codice colore per la colormap del log2
colormap=LinearSegmentedColormap.from_list('rg',["g", "w", "r"], N=64)

"""
# paths
    Xlsx_output = "/Users/lucas/16S_amplicon_analysis/OUT_16S_Amplicon_Analysis.xlsx"
    Legend_csv_path = '/Users/lucas/16S_amplicon_analysis/Legend.csv'
    TotalCSV_path = '/Users/lucas/16S_amplicon_analysis/filtered_rows_total_page.csv'
"""
def data_table(total_page:pd.DataFrame, legend:pd.DataFrame):   
    NGSID_ASVs = []
    ngs_id = legend['NGS_ID'].tolist()
    ngs_id_perc = [x + '%' for x in ngs_id] # genera la lista con i nomi delle colonne percentuali es S23%, S24%... per recuperare i dati da TotalCSV

    for pos in range(0, len(legend)):
        NGSID_ASVs.append(total_page.loc[:,legend.loc[pos, 'NGS_ID']].sum())
    for pos in range(0, len(legend)):
        total_page[legend.loc[pos, 'NGS_ID'] + '%'] = round((total_page[legend.loc[pos, 'NGS_ID']]/NGSID_ASVs[pos])*100, 3)

    total_page['AVG'] = round(total_page.iloc[:,len(total_page.axes)-len(legend):].mean(axis=1), 2)     #colonna 'AVG' 
    """
    ricalcolo la percentuale di abbondanza di ogni ASV per Sample
    può accadere che per ASV che in 'Total' superassero la soglia del 0.01% in un Sample 
    non la superino a seguito del trimming e ricalcolo e in 'Total_order' vengono segnati come 0
    """
    total_page['>0.01%'] = (total_page.loc[:,ngs_id_perc] > 0.01).sum(axis=1)   #colonna 'ASV > 0.01%'
    total_page['>0.5%'] = (total_page.loc[:,ngs_id_perc] > 0.5).sum(axis=1)     #colonna 'ASV > 0.5%'
    return total_page, ngs_id

def samples_table(main_table:pd.DataFrame, legend:pd.DataFrame, Sample_aboundance_path:str):
    sample_table = pd.DataFrame()
    samples_dict = {}
    for pos in range(0, len(legend)):
        sample_table[legend.loc[pos, 'Samples']] = np.NaN
        samples_dict[legend.loc[pos, 'Samples']] = int(legend.loc[pos, 'Order']+1)
    samples = sample_table.columns.tolist()
    column_pos = main_table.columns.get_loc(legend.loc[0, 'NGS_ID'] + '%')
    for pos in range(0, len(samples)):
        sample_table[samples[pos]] = round(main_table.iloc[:,column_pos:column_pos+samples_dict[samples[pos]]].mean(axis=1), 3)
        column_pos += samples_dict[samples[pos]]
    
    
    DF_sample_aboundance_table = sample_table.copy()
    main_table = pd.concat([main_table, DF_sample_aboundance_table], axis=1)

    io.WriteFiles.writeToCsv(Sample_aboundance_path, main_table)
    
    DF_percentages01 = pd.DataFrame(sample_table)
    DF_percentages01 = DF_percentages01> 0.01 # convert bool table into 0/1 values
    DF_percentages01 = DF_percentages01.astype(int)
    DF_percentages01.rename(columns={x: x + ' >1%' for x in DF_percentages01.columns}, inplace=True)
    DF_percentages01['>1%'] = DF_percentages01.sum(axis=1)
    sample_table.replace(0, np.NaN, inplace=True)
    sample_table = pd.concat([sample_table, DF_percentages01], axis=1)
    return sample_table, samples

def log_table(sample_table:pd.DataFrame, samples:list):
    #calcola il log2 tra i campioni, produce 2 tabelle: una con i valori e una con i valori assoluti dei log2
    log_table = pd.DataFrame()
    for pos in range(0, len(samples), 2):
        if pos+1 == len(samples):
            break
        log_table[samples[pos+1] +' Vs '+ samples[pos]] = round(np.log2(sample_table[samples[pos+1]]/sample_table[samples[pos]]), 2)
    
    log_table['GR -/- TRAT'+' Vs '+'WT CTR'] = round(np.log2(sample_table['GR -/- TRAT']/sample_table['WT CTR']), 2)
    log_table['GR -/- TRAT'+' Vs '+'GR +/+ TRAT'] = round(np.log2(sample_table['GR -/- TRAT']/sample_table['GR +/+ TRAT']), 2)
    log_table['GR +/+ TRAT'+' Vs '+'WT TRAT'] = round(np.log2(sample_table['GR +/+ TRAT']/sample_table['WT TRAT']), 2)
    log_table['GR -/- TRAT'+' Vs '+'GR +/+ TRAT'] = round(np.log2(sample_table['GR -/- TRAT']/sample_table['GR +/+ TRAT']), 2)

    log_table_columns = log_table.columns.tolist()
    abs_log_table = log_table.abs()**2
    abs_log_table.rename(columns={x: 'abs('+x+')' for x in abs_log_table.columns}, inplace=True)
    return_table = pd.concat([log_table,abs_log_table], axis=1)
    return return_table,log_table_columns

def stat_page_update(main_table:pd.DataFrame, legend:pd.DataFrame, stat_path:str):
    #prende i dati di Stat e li aggiorna con i nuovi dati
    stat_page = io.ReadFiles.readToDataframe(stat_path, 'Stat')
    for pos in range(0, len(legend)):
        stat_page.loc[pos, 'Divesity ASVs >0.5%'] = main_table[legend.loc[pos,'Samples']].sum()
    return stat_page    
    
def refined_page(final_table:pd.DataFrame, sample_id:list, log_id:list):
    Header = ('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
    refined_table = final_table[sample_id]
    refined_table = pd.concat([refined_table, final_table['>1%']], axis=1)
    refined_table = pd.concat([refined_table, final_table[log_id]], axis=1)
    
    
    refined_table['IDs'] = final_table['Genus']
    
    for pos in reversed(range(len(Header))):
        refined_table['IDs'] = refined_table['IDs'].replace(0, np.NaN)
        refined_table['IDs'] = refined_table['IDs'].fillna(value=final_table[Header[pos]])
    refined_table['ASVs'] = final_table['ASVs']
        
    refined_table = refined_table[refined_table['>1%'] > 0]
    
    return refined_table

def main(Xlsx_output:str, Legend_path:str, TotalCSV_path:str, Sample_aboundance_path:str):
    total_page = io.ReadFiles.readToDataframe(TotalCSV_path, None)
    legend = io.ReadFiles.readToDataframe(Legend_path, None)

    # ricalcolo la percentuale di ASV per ogni NGS_ID
    main_table, ngs_id = data_table(total_page, legend)
    sample_table, sample_id = samples_table(main_table, legend, Sample_aboundance_path)
    log2_table, log_id = log_table(sample_table, sample_id)

    DF_final_table = pd.concat([main_table, sample_table], axis=1)
    DF_final_table = pd.concat([DF_final_table, log2_table], axis=1)
            
    DF_final_table.fillna(0, inplace=True)
    stat_page = stat_page_update(DF_final_table, legend, Xlsx_output)
    io.WriteFiles.writeToXlsx(Xlsx_output, 'Stat', stat_page, 0, 0) #aggiorna la pagina Stat
    
    DF_refined = refined_page(DF_final_table, sample_id, log_id)
    
    with pd.ExcelWriter(
        Xlsx_output, engine='openpyxl', mode='a', if_sheet_exists='replace'
    )as writer:
            #genera pagina Total_order con i dati di DF_final_table colorati
            DF_final_table.style.background_gradient(subset=ngs_id, axis=None, cmap='Greens').\
                background_gradient(subset=sample_id, axis=None, cmap='Greens').\
                    background_gradient(subset=log_id, axis=None, cmap=colormap, vmin=-1, vmax=1).\
                        to_excel(writer, sheet_name='Total_order')
       
            #genera pagina Refined con i dati filtrati di DF_final_table (Total_order)            
            DF_refined.style.background_gradient(subset=sample_id, axis=None, cmap='Greens').\
                background_gradient(subset=log_id, axis=None, cmap=colormap, vmin=-1, vmax=1).\
                    to_excel(writer, sheet_name='Refined')
    print("Total_order element saved in:\t\t", Xlsx_output)
    print("Refined element saved in:\t\t", Xlsx_output)
    
    return True
