import pandas as pd
import InputOutput_class as io


"""
#info
Il programma genera una tabella contente in ordine gli 'NGS_ID' associati ai rispettivi 'Samples'
Il ruolo del file generato è di essere di riferimento per i programmi successivi per mantenere 
l'associazione per le analisi successive, risparmiando tempo e ridondanza di chiamate e complessità
genera un secondo file chiamato 'group_sample_size.csv' che contiene il numero di replicati(size)
per ogni campione(Samples) per semplificare analisi dei campioni

PRODUCE 3 OUTPUT:
- Legend.xlsx
- Legend.csv

#paths
    #input
        NGS_IDtoSample_path = '/Users/lucas/16S_amplicon_analysis/phyloseq/NGS_IDtoSample.csv'
        QualityInfo_path = '/Users/lucas/16S_amplicon_analysis/phyloseq/quality_info.csv'
    #output
        legend_xlsx_path = "/Users/lucas/16S_amplicon_analysis/OUT_16S_Amplicon_Analysis.xlsx/Legend"
        legend_csv_path = '/Users/lucas/16S_amplicon_analysis/Legend.csv'
"""

## main ###
def main(Xlsx_path, Legend_csv_path, QualityInfo_path, NGS_IDtoSample_path):
   
    # read files
    Header = ('Position', 'NGS_ID', 'Samples','Replicates','Order') # columns of the Legend dataframe
    DF_NGS_IDtoSample = io.ReadFiles.readToDataframe(NGS_IDtoSample_path,'Legend')
    DF_q_info = io.ReadFiles.readToDataframe(QualityInfo_path,None)
    
    # create legend dataframe
    DF_Legend = pd.DataFrame(columns=Header)
    DF_Legend['Position'] = range(1, len(DF_q_info)+1)
    DF_Legend['NGS_ID'] = DF_q_info['NGS_ID']

    """
    Ho creato questo ciclo for nested per associare i valori NGS_ID al rispettivo campione (Sample)
    Tenendo conto della presenza di elementi nulli all'interno del file NGS_IDtoSample.csv
    Praticamente esploro ogni singola cella passando da X a Y (coordinate cartesiane) 
    e copio il valore se non è nullo nella colonna Samples
    """ 
    
    pos = 0
    for X in range(len(DF_NGS_IDtoSample.columns)):
        for Y in range(len(DF_NGS_IDtoSample[DF_NGS_IDtoSample.columns[0]])):
            value = DF_NGS_IDtoSample[DF_NGS_IDtoSample.columns[X]][Y]
            if pd.notnull(value):
                DF_Legend.loc[pos, 'Samples'] = DF_NGS_IDtoSample.columns[X]  # copia i nomi dei campioni
                pos += 1
            else:
                break    
    """
    Questo ciclo for serve per associare un valore alfabetico ai replicati di un campione max 26 replicati
    e per calcolare il numero di replicati per ogni campione
    
    riempe la colonna 'Replicates' e 'Order' del dataframe 'DF_Legend'
    crea un dataframe 'group_sample_size' con le colonne 'Samples' e 'size' contenente il numero di replicati
    """
    word =''
    alphabet = 'abcdefghijklmnopqrstuvwxyz' # max 26 replicates
    num = 0
    for sample in DF_Legend['Samples']:
        if sample != word:
            for i in range(0,DF_NGS_IDtoSample[sample].count()):
                pos = DF_Legend.loc[DF_Legend['Samples'] == sample].index[i]
                DF_Legend.loc[pos, 'Replicates'] = alphabet[i]
                DF_Legend.loc[pos, 'Order'] = i
            word = sample
            num += 1
        """
        group_sample_size è una lista che contiene il numero di replicati per ogni campione
        con associazione Samples,size: WT CTR,3
        """
        

    legend_xlsx_out = io.WriteFiles.writeToXlsx(Xlsx_path, 'Legend', DF_Legend, 0, 0) # scrive la leggenda in Legend.xlsx
    legend_csv_out = io.WriteFiles.writeToCsv(Legend_csv_path, DF_Legend) # scrive la leggenda in Legend.csv
    
    if legend_xlsx_out and legend_csv_out == True: # verifico se i file sono stati scritti correttamente
        return True
    else:
        return False
