import pandas as pd
import matplotlib.pyplot as plt
import InputOutput_class as io

"""
# info
utilizzando i file 'beta_diversity_values.csv', 'beta_diversity_vectors.csv', 'Legend.csv' e 'group_sample_size.csv'
genera la terza pagina del file excel, Beta page contenente:
1.  una tabella (VECTOR_table) composta alla prima riga dei valori di beta_diversity_values.csv (in particolare la colonna 'relative_eig') 
    convertiti in percentuale e rinominati in 'Bray-Curtis', alla seconda riga i nomi delle colonne di beta_diversity_vectors.csv 
2.  una tabella (PCoA_2D_table) con le prime 2 colonne di VECTOR_table per essere utilizzata in un grafico
3.  un grafico composto dai valori numeri di PCoA_2D_table, colorati in base al Sample specificato tramite 'Legend.csv'

#input directories
    Group_size_path = '/Users/lucas/16S_amplicon_analysis/group_sample_size.csv'
    Beta_VECTOR_path= '/Users/lucas/16S_amplicon_analysis/phyloseq/beta_diversity_vectors.csv'
    BetaDiversity_path = '/Users/lucas/16S_amplicon_analysis/phyloseq/beta_diversity_values.csv'
#output directories
    Xlsx_path = '/Users/lucas/16S_amplicon_analysis/OUT_16S_Amplicon_Analysis.xlsx'
    Graph_folder_path = '/Users/lucas/16S_amplicon_analysis/graphs'
"""
def graphScatter(Graph_folder_path:str ,PCoA_2D_table:pd.DataFrame, DF_Legend:pd.DataFrame):
    # inizializzo il grafico
    plt.axhline(0,color='black', linewidth=0.8)
    plt.axvline(0,color='black', linewidth=0.8)
    plt.grid(color='gray', axis='y', linewidth=0.4)
    plt.xlabel('PCo1')
    plt.ylabel('PCo2')
    plt.title('PCoA 2D')

    # creo lo scatter plot
    num = 0
    for pos in range(2,len(PCoA_2D_table)):
        # if else servono ad attribuire un colore diverso per ogni campione e a creare una legenda
        if DF_Legend.loc[pos-2,'Order'] == 0:
            num +=1
            dot_color = 'C'+str(num) # colori diversi per ogni campione
            legend_name = DF_Legend['Samples'][pos-2] # nomi dei campioni
        else: 
            legend_name = None
        plt.scatter(PCoA_2D_table.loc[pos:pos+1,1], PCoA_2D_table.loc[pos:pos+1,2], alpha=0.8, color = dot_color, label = legend_name)


    plt.legend(loc='right', bbox_to_anchor=(1.35, 0.5),fancybox=True)
    graph_path = Graph_folder_path + '/PCoA_2D.png'
    plt.savefig(graph_path, dpi=550, bbox_inches='tight')
    plt.clf()
    return graph_path

def main(Xlsx_output, Graph_folder_path, BetaDiversity_VALUES_path, BetaDiversity_VECTOR_path):
    DF_VALUES_beta_diversity = io.ReadFiles.readToDataframe(BetaDiversity_VALUES_path, None)
    DF_VECTOR_beta_diversity = io.ReadFiles.readToDataframe(BetaDiversity_VECTOR_path, None)
    DF_Legend = io.ReadFiles.readToDataframe(Xlsx_output, 'Legend')

    """
    Genera la prima colonna della tabella contenente i valori di beta_diversity_vectors.csv
    questa prima colonna è composta dai valori di beta_diversity_values.csv
    convertiti in percentuale e sottoposti a trimming (rimozione dei valori negativi)
    successivamente rinominata in 'Bray-Curtis'
    """

    relative_eig = pd.DataFrame(DF_VALUES_beta_diversity['Relative_eig']) # seleziona la colonna 'Relative_eig' di beta_diversity_values.csv
    relative_eig = round(relative_eig*100,2)
    relative_eig = relative_eig[relative_eig['Relative_eig'] > 0] # rimuove i valori negativi 'trimming' dei valori
    relative_eig.rename(columns={'Relative_eig':'0' }, inplace=True)
    relative_eig = relative_eig.transpose()
    relative_eig.insert(0, '0', 'Bray-Curtis')
    relative_eig.columns = range(0,len(relative_eig.columns))
    # relative_eig rappresenta la prima riga futura VECTOR_table, contenenti le % di beta_diversity_values.csv

    DF_VECTOR_beta_diversity.insert(1, 'Samples', DF_Legend['Samples'])  # aggiunta colonna 'Samples' a DF_VECTOR_beta_diversity
    for x in range(2, len(DF_VECTOR_beta_diversity.columns)):
        DF_VECTOR_beta_diversity.rename(columns={DF_VECTOR_beta_diversity.columns[x]: 'PCo' + str(x-1)}, inplace=True)

    DF_VECTOR_beta_diversity =DF_VECTOR_beta_diversity.drop(DF_VECTOR_beta_diversity.columns[0], axis=1) 
    #rimuovo la colonna degli NGS ID siccome ho aggiunto la colonna 'Samples'

    DF_header = pd.DataFrame(data=DF_VECTOR_beta_diversity.columns)
    DF_header = DF_header.transpose()
    VECTOR_table = pd.concat([relative_eig, DF_header ], axis=0) # creo la tabella VECTOR_table
    #prendo gli indici di colonna di DF_VECTOR_beta_diversity e li aggiungo come seconda riga della tabella VECTOR_table

    beta_diversity = pd.DataFrame(data=DF_VECTOR_beta_diversity)
    beta_diversity.columns = range(0,len(relative_eig.columns)) # rinomino le colonne di beta_diversity, che derivano da DF_VECTOR_beta_diversity in numeri

    VECTOR_table = pd.concat([VECTOR_table, beta_diversity], axis=0)
    VECTOR_table.index = range(0,len(VECTOR_table)) # rinomino gli indici(righe) di VECTOR_table

    PCoA_2D_table = pd.DataFrame(VECTOR_table.iloc[:,0:3]) # seleziono le prime 3 colonne di VECTOR_table e le salvo in DF_PCoA_2D
    graph_path = graphScatter(Graph_folder_path, PCoA_2D_table, DF_Legend) # creo il grafico PCoA_2D
    
    VECTOR_table_out = io.WriteFiles.writeToXlsx(Xlsx_output, 'Beta', VECTOR_table, 1, 1)
    y_coordinate = len(VECTOR_table)+3
    PCoA_2D_table_out = io.WriteFiles.writeToXlsx(Xlsx_output, 'Beta', PCoA_2D_table, y_coordinate, 1)
    Graph_out = io.WriteFiles.writeGraphs(graph_path, Xlsx_output, 'Beta', 'V2')
    
    if VECTOR_table_out and PCoA_2D_table_out and Graph_out == True:
        return True
    else:
        return False