import pandas as pd
import InputOutput_class as io
import matplotlib.pyplot as plt

"""
# info
utilizzando i file 'quality_info.csv' 'alpha_diversity.csv' 'NGS_IDtoSample.csv' e 'group_sample_size.csv'
genera la seconda pagina del file excel, Stat page contenente:
1.una tabella con le colonne di 'quality_info' e 'alpha_diversity'
2.una tabella con la media di ASVs e Shannon per campione
2.un grafico con la media di ASVs per campione

# paths
    QualityInfo_path = '/Users/lucas/16S_amplicon_analysis/phyloseq/quality_info.csv'
    AlphaDiversity_path = '/Users/lucas/16S_amplicon_analysis/phyloseq/alpha_diversity.csv'
    Xlsx_output = "/Users/lucas/16S_amplicon_analysis/OUT_16S_Amplicon_Analysis.xlsx"
"""

def createQualityInfoTable(q_info:pd.DataFrame, ngsid_sample:pd.DataFrame):  # crea un dataframe con le colonne di quality_info
    Header = ('Position', 'Samples', 'Replicates', 'NGS ID', 'Raw Reads', 'Lenght', 'Non-chimeric', 'perc of chimeric', 'assigned to ASVs', 'perc of assigned to ASVs')
    q_tb = pd.DataFrame(columns=Header)
    q_tb['Position'] = ngsid_sample['Position']  # numera i campioni
    q_tb['Samples'] = ngsid_sample['Samples']  # copia la colonna Samples di NGS_IDtoSample
    q_tb['NGS ID'] = q_info['NGS_ID']  # copia la colonna NGS_ID di quality_info
    q_tb['Replicates'] = ngsid_sample['Replicates']  # copia la colonna replicates di NGS_IDtoSample
    ##
    q_tb['Raw Reads'] = q_info['input'] 
    q_tb.loc[len(q_info), 'Raw Reads'] = round(q_info['input'].mean())  # calcola la media di Raw Reads
    ##
    q_tb['Lenght'] = q_info['length_REV']  # copia la colonna lenght_REV di quality_info
    ##
    q_tb['Non-chimeric'] = q_info['nonchim']  # copia la colonna nonchim di quality_info
    q_tb.loc[len(q_info), 'Non-chimeric'] = round(q_info['nonchim'].mean())  # calcola la media di non-chimeric
    ##
    q_tb['perc of chimeric'] = round(q_info['nonchim']/q_info['input']*100)
    q_tb.loc[len(q_info), 'perc of chimeric'] = round(q_tb['perc of chimeric'].mean())  # calcola la media di perc of chimeric
    
    return q_tb

def createAlphaDiversityTable(alpha_div:pd.DataFrame):  # crea un dataframe con le colonne di alpha_diversity
    Header = ('ASVs', 'Shannon', 'ASVs >0.01%', 'ASVs >0.5%', 'Divesity ASVs >0.5%')
    al_tb = pd.DataFrame(columns=Header)
    al_tb['ASVs'] = alpha_div['Observed']
    al_tb['Shannon'] = alpha_div['Shannon']
    return al_tb

def createAlphaTable(al_df:pd.DataFrame, legend:pd.DataFrame, Graphs_folder_path:pd.DataFrame):  # crea un dataframe di riassunto delle informazioni di ASVs e Shannon relativi ai campioni
    Header = ('Sample AVG', 'ASVs', 'Shannon')
    avg_tb = pd.DataFrame(columns=Header)
    pointer = 0
    sample_size = legendToDict(legend)
    sample_id = list(sample_size.keys())
    avg_tb['Sample AVG'] = sample_id
    for val in range(0,len(sample_id)):
        num = sample_size[sample_id[val]]
        avg_tb.loc[val,'ASVs'] = round(al_df.loc[pointer:pointer+num-1, 'ASVs'].mean())  # calcola la media di ASVs
        avg_tb.loc[val,'Shannon'] = round(al_df.loc[pointer:pointer+num-1, 'Shannon'].mean(),1)  # calcola la media di Shannon + arrotonda a 1 decimale
        pointer += num
    
    plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
    plt.bar(avg_tb['Sample AVG'], avg_tb['ASVs'], width=0.5, color='#145DA0', alpha=1)
    plt.xlabel('Samples')
    plt.ylabel('ASVs')
    plt.title('Avgerage ASVs per Sample')
    graph_path = Graphs_folder_path + '/Average_ASVs.png'
    plt.savefig(graph_path, dpi=400,  bbox_inches='tight')
    plt.clf()
    return avg_tb, graph_path

def legendToDict(legend:pd.DataFrame):  # crea un dizionario con le colonne di legend
    samples_dict = {}
    for pos in range(0, len(legend)):
        samples_dict[legend.loc[pos, 'Samples']] = int(legend.loc[pos, 'Order']+1)
        
    return samples_dict

# main ###
def main(Xlsx_output:str, Graphs_folder_path:str, QualityInfo_path:str, AlphaDiversity_path:str):
    
    # lettura dei file
    legend_xlxs = io.ReadFiles.readToDataframe(Xlsx_output,'Legend')
    quality_info = io.ReadFiles.readToDataframe(QualityInfo_path, None)  # legge il file quality_info.csv
    alpha_diversity = io.ReadFiles.readToDataframe(AlphaDiversity_path, None)  # legge il file alpha_diversity.csv

    # creazione dei dataframe
    DF_q_info= createQualityInfoTable(quality_info, legend_xlxs)  # crea il dataframe di quality_info
    alTable = createAlphaDiversityTable(alpha_diversity)  # crea il dataframe di alpha_diversity
    table = pd.concat([DF_q_info, alTable], axis=1)  # unisce i due dataframe
    alphaTable, Graph_path = createAlphaTable(alTable, legend_xlxs, Graphs_folder_path)  # crea il dataframe di riassunto di ASVs e Shannon dei campioni per formare l'istogramma

    # scrittura del file excel
    
    table_out = io.WriteFiles.writeToXlsx(Xlsx_output, 'Stat', table, 0, 0)
    x_coordinate = len(table.columns) + 3
    alphaTable_out = io.WriteFiles.writeToXlsx(Xlsx_output, 'Stat', alphaTable, 0, x_coordinate)
    cell_pos = 'W3' # posizione del grafico
    graph_out = io.WriteFiles.writeGraphs(Graph_path, Xlsx_output, 'Stat', cell_pos)
    if table_out and alphaTable_out and graph_out == True:
        return True
    else:
        return False

