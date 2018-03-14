import pandas as pd
import matplotlib.pyplot as plt
 
# Create a dataframe
gtex = pd.read_csv('../logs/gtex_v7_hallmark_results.log', sep='\t')
tcga = pd.read_csv('../logs/panTCGA_hallmark_results.log', sep='\t')

df = pd.DataFrame({'group':list(gtex['Num']), 'gtex':gtex['Average'] , 'panTCGA':tcga['Average'] })
 
# Reorder it following the values of the first value:
ordered_df = df.sort_values(by='gtex')
my_range=range(1,len(df.index)+1)
 
# The vertical plot is made using the hline function
# I load the seaborn library only to benefit the nice looking feature
import seaborn as sns
plt.vlines(x=my_range, ymin=ordered_df['gtex'], ymax=ordered_df['panTCGA'], color='grey', alpha=0.4)
plt.scatter(my_range, ordered_df['gtex'], color='skyblue', alpha=1, label='gtex')
plt.scatter(my_range, ordered_df['panTCGA'], color='red', alpha=0.4 , label='panTCGA')
plt.legend()
 
# Add title and axis names
plt.xticks(my_range, ordered_df['group'], rotation=70, fontsize=8, ha='right')
plt.title("Comparison of hallmark accuracies between GTEX and panTCGA", loc='Center')
plt.xlabel('Set')
plt.ylabel('Accuracies')

# plt.tight_layout()
plt.gcf().subplots_adjust(bottom=0.30)
plt.show()