#!/usr/bin/env python
"""
# Author: Xianbin Meng
# File Name: SCMS.py
# Description:
"""
import numpy as np
import pandas as pd
import argparse
import seaborn as sns
import os
import time
import matplotlib.pyplot as plt
from sklearn.preprocessing import MaxAbsScaler

from lib.utils import _process, _filter3, _knn, _plot_embedding, _diff

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Single Cell Mass Spectrometry Analysis')
    parser.add_argument('--file', '-f', type=str, default='')
    parser.add_argument('--metafile', type=str, default=None, help='Pre-processed file for visualization. Default:None')
    parser.add_argument('--group_names', '-g', nargs='+', default=None, help='Group names. Default:None')
    parser.add_argument('--ppm_threshold_peak', type=int, default=10, help='Peak error threshold for peak selection in the same cell. Default:10')
    parser.add_argument('--ppm_threshold_cell', type=int, default=20, help='Peak error threshold for combining different cells. Default:20')
    parser.add_argument('--decrease', default=True, help='Peak selection mode. Default:True')
    parser.add_argument('--peak', default=True, help='Whether select peak for each cell. Default:True')
    parser.add_argument('--filter', type=float, default=0.5, help='Peaks that appear in less than 100*filter cells will be filtered. Default:0.5')
    parser.add_argument('--knn', default=True, help='Whether use knn for missing value imputation. Default:True')
    parser.add_argument('--n_neighbors', type=int, default=5, help='KNN algorithm parameter. Default:5')
    parser.add_argument('--method', '-m', type=str, default='all', help='Dimensional reduction methods, including PLS,PCA,UMAP,tSNE. Default:all')
    parser.add_argument('--q_value', type=float, default=0.05, help='Adjusted p_value threshhold for extracting differential expressed signal among different groups. Default:0.05')
    parser.add_argument('--log2fold', type=float, default=0.5, help='Log-foldchange threshhold for extracting differential expressed signal among different groups. Default:0.5')
    parser.add_argument('--seed', type=int, default=2020, help='Random seed. Default:2020')
    parser.add_argument('--outdir', '-o', type=str, default='./', help='Output path. Default:./')
    args = parser.parse_args()

    np.random.seed(args.seed) 
    
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    
    timestart = time.time()
    
    if args.metafile is None:
        if os.path.isfile(args.file):
            df = pd.read_csv(args.file)
        else:
            raise ValueError("{} file does not exist!")
        print('Preprocessing...')
        meta, m_z = _process(df,
                             groups=args.group_names,
                             peak=args.peak,
                             ppm_threshold_peak=args.ppm_threshold_peak,
                             ppm_threshold_cell=args.ppm_threshold_cell,
                             decrease=args.decrease)

        meta.to_csv(outdir+'/signal.csv')
        m_z.to_csv(outdir+'/m_z.csv')

        plt.figure(figsize=(6,6))
        counts = pd.DataFrame({'counts':(meta>0).sum(axis=1)})
        ax = sns.violinplot(y='counts',data=counts, orient='h')
        ax.tick_params(axis='x',bottom=True, top=False, labelbottom=True, labeltop=False, labelsize=12, length=3, pad=3,labelrotation=0)
        ax.tick_params(axis='y',bottom=True, top=False, labelbottom=True, labeltop=False, labelsize=12, length=3, pad=3,labelrotation=0)
        ax.set_ylabel('m/z signal counts', fontsize=14, labelpad=10, va='center')
        plt.savefig(outdir+"/violinplot.pdf")

        meta = _filter3(meta, filter_peak=args.filter)
        m_z = m_z.loc[meta.columns,:]
        print('After filtering, {} m/z signals are retained.'.format(meta.shape[0]))
        meta.to_csv(outdir+'/signal_filtered.csv')
        m_z.T.to_csv(outdir+'/m_z_filtered.csv')
    
    else:
        print('Using pre-processed file for visualization')
        if '.csv' in args.metafile:
            meta = pd.read_csv(args.metafile, index_col=0)
        elif '.txt' in args.metafile:
            meta = pd.read_csv(args.metafile, sep='\t', index_col=0)
    
    if args.knn:
        print('Performing knn for imputation...')
        meta = _knn(meta,n_neighbors=args.n_neighbors)
        meta.to_csv(outdir+'/signal_filtered_knn.csv')
    
    labels = [item.rsplit('-',1)[0] for item in meta.index]

    if args.method=='all':
        print('Performing visualization...')
        for method in ['PLS', 'PCA', 'UMAP', 'tSNE']:
            emb = _plot_embedding(meta, method=method, labels=labels, save=outdir+"/{}.pdf".format(method), return_emb=True)
            emb = pd.DataFrame(emb, index=meta.index, columns=["{}_{}".format(method,i) for i in np.arange(1, emb.shape[1]+1)])
            emb.to_csv(outdir+"/{}.csv".format(method))
    else:
        print('Performing {} for visualization...'.format(args.method))
        emb = _plot_embedding(meta, method=args.method, labels=labels, save=outdir+"/{}.pdf".format(args.method), return_emb=True)
        emb = pd.DataFrame(emb, index=meta.index, columns=["{}_1".format(args.method),"{}_2".format(args.method)])
        emb.to_csv(outdir+"/{}.csv".format(args.method))
    
    print('Performing differential expression analysis...')
    de = _diff(meta, labels, q_value=args.q_value, log2fold=args.log2fold)
    de.to_csv(outdir+'/de.csv')
    
    meta = pd.DataFrame(MaxAbsScaler().fit_transform(meta), index=meta.index, columns=meta.columns)
    sns.clustermap(meta.loc[:,de.index].T,
                   col_cluster=False,cmap='RdBu_r',
                   figsize=(meta.shape[0]/15, de.shape[0]/8),
                   dendrogram_ratio=0.25)
    plt.savefig(outdir+"/heatmap.pdf")
    
    timeend = time.time()
    print('Running time: {}'.format(timeend-timestart))
    print('Done')
    