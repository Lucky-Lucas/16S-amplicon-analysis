import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import InputOutput_class as io

def legendToDict(legend:pd.DataFrame):  # crea un dizionario con le colonne di legend
    samples_dict = {}
    for pos in range(0, len(legend)):
        samples_dict[legend.loc[pos, 'Samples']] = int(legend.loc[pos, 'Order']+1)
        
    return samples_dict


DF_sample_aboundance = io.ReadFiles.readToDataframe('/Users/lucas/16S_amplicon_analysis/sample_aboundance_table.csv', None)
DF_legend = io.ReadFiles.readToDataframe('/Users/lucas/16S_amplicon_analysis/Legend.csv', None)
DF_legend_dict = legendToDict(DF_legend)
Header = ('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')


bar_graph, bar_ax = plt.subplots()
plt.title('Sample Aboundance Table, ASVs > 1%')
plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right', fontsize=12)

bottom = {keys: 0 for keys in DF_legend_dict.keys()}
filtered_samples = pd.DataFrame()
average_samples = pd.DataFrame()

for key in DF_legend_dict.keys():
    filtered_samples = DF_sample_aboundance[DF_sample_aboundance.loc[:, key] > 0.01]
    filtered_samples.sort_values(by=key, ascending=False, inplace=True)
    sample_sum = filtered_samples[key].sum()
    for row in filtered_samples.index:
        perc_sample = DF_sample_aboundance.loc[row,key]
        bar_label = DF_sample_aboundance.loc[row,'ASVs'] + '  ' + round(perc_sample*100, 1).astype(str) + '%'
        for pos in reversed(range(len(Header))):
            if Header[pos] is not np.NaN:
                name = Header[pos]
        if perc_sample > 0.03:
            plt.text(key, bottom[key]+ perc_sample/2, bar_label, fontsize=11, fontweight='book', color='black')
            bar_ax.bar(key, perc_sample, bottom=bottom[key], align='edge', width=0.6)
        else:
            bar_ax.bar(key, perc_sample, bottom=bottom[key], align='edge', label=bar_label, width=0.6)

        bottom[key] += perc_sample
        
    bar_ax.bar(key, 1 - bottom[key], bottom=bottom[key], align='edge', label=bar_label, width=0.6)
    plt.gca().set_prop_cycle(None)
    plt.text(key, bottom[key], 'Others', fontsize=12, fontweight='book', color='black')

bar_graph.set_size_inches(20, 10.5)
bar_graph.savefig('samplefigure', dpi=550, bbox_inches='tight')
bar_graph.legend(loc='right', bbox_to_anchor=(1.35, 0.5),fancybox=True)
plt.show()
plt.clf()
