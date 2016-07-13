import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def genbank2pandas(gb_file,feature=None,name_func=None):
    table = pd.read_table(gb_file,header=None)
    table.columns = ['genome', 'source','feature','start','end','score','strand','frame','attr']
    table = table.sort_values(['strand','start'])
    if feature != None:
        table = table[table.feature.isin(feature)]
    return table

def gff2pandas(gff_filename,normalize=True):
    # Import data into pandas table
    raw_data = pd.read_table(gff_filename,header=None)
    raw_data.columns=['genome', 'source','feature','start','end','score','strand','frame','attr']
    raw_data = raw_data[['start','score','strand']]

    # Isolate plus and minus strands
    plus_data = raw_data[raw_data.strand=='+']
    plus_data.index = plus_data.start.values
    plus_series = plus_data.score
    minus_data = raw_data[raw_data.strand=='-']
    minus_data.index = minus_data.start.values
    minus_series = minus_data.score
    
    # Merge together strands into single dataframe
    df = pd.DataFrame(index = range(raw_data.start.max()))
    df = pd.concat([df,plus_series,minus_series],axis=1).fillna(0)
    df.columns = ['score_plus','score_minus']
    if normalize:
        return df*np.true_divide(10**7,df.values.sum())
    else:
        return df

def plot_region(start,end,tracks,annotation=[],figsize=(None,None),labels=None,num_ticks=20):
    # Get total number of tracks
    num_tracks = len(tracks) + int(len(annotation)!=0)
    
    # Calculate figure size
    fig_width,fig_height = figsize
    if fig_height == None:
        fig_height = 1.5*len(tracks) + 0.75*int(len(annotation)!=0)
    if fig_width == None:
        fig_width = 10
    
    # Initialize subplots
    fig,axes = plt.subplots(num_tracks,1,squeeze=False,figsize=(fig_width,fig_height),
                            subplot_kw={'xlim':(start,end),
                                        'xticks':[],
                                        'xticklabels':[],
                                        },
                            gridspec_kw = {'height_ratios':[2]*len(tracks)+[1]}
                            )
    
    
    # Include annotation arrows on bottom plot if included
    if len(annotation)>0:
        # Trim plot height
        axes[-1][0].set_ylim((-0.5,1))
        # Remove y-ticks
        axes[-1][0].set_yticks([])
        
        # Trim annotation dataframe
        trimmed_data = annotation[(annotation['end'] >= start) & (annotation['start'] <= end)]
        for i,row in trimmed_data.iterrows():
            # Flip arrows if on negative strand
            if row.strand =='+':
                gene_start = row.start
                gene_end = row.end
            else:
                gene_start = row.end
                gene_end = row.start
            
            # Add arrows and labels for genes
            axes[-1][0].add_patch(mpatches.Arrow(gene_start,0,gene_end-gene_start,0,fc='lightgray',ec='k'))
            axes[-1][0].text((gene_start+gene_end)/2,0.5,row.name,ha='center')

    # Plot tracks
    for i,track in enumerate(tracks):
        # Add baseline
        axes[i][0].plot(range(start,end),[0]*(end-start),linestyle='--', color='k')
        # Add data
        axes[i][0].plot(track['score_plus'][start:end],'b')
        axes[i][0].plot(-track['score_minus'][start:end],'g')
        # Add labels
        if labels != None:
            axes[i][0].set_ylabel(labels[i],fontsize=12)

    fig.tight_layout()

    # Add x-ticks to bottom axes
    axes[-1][0].set_xticks(range(start,end,(end-start)/num_ticks))
    axes[-1][0].set_xticklabels(range(start,end,(end-start)/num_ticks))
    axes[-1][0].xaxis.set_ticks_position('bottom')
    return axes

def plot_genes(gene_list,tracks,annotation=[],figsize=(None,None),labels=None,num_ticks=20,ranged=False):
    rows = annotation.loc[gene_list]
    start = min(rows.start)
    end = max(rows.end)
    if ranged:
        annot_rows = annotation
    else:
        annot_rows = rows
    return plot_region(start,end,tracks,annotation=annot_rows,figsize=figsize,labels=labels,num_ticks=num_ticks)
