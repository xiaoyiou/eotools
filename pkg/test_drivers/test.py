"""
This file serves as the testing driver for pre-processing unit of data
"""
from data import mdata, analysis as ana
from eotools.dtools import rawdata
import cPickle as pickle

ppfix = '/Users/xiaoyiou/projects/pyprojs/eotools/'
prefix = ppfix + 'tracks/'
gprefix = ppfix + 'genelst/'
n = 4
mod_files=[
 'GSM701923_H3K4me2.bed',
 'GSM701924_H3K4me3.bed',
 'GSM701925_H3K9Ac.bed',
 'GSM701926_H3K9me2.bed',
 'GSM701927_H3K18Ac.bed',
 'GSM701928_H3K27me1.bed',
 'GSM701929_H3K27me3.bed',
 'GSM701930_H3K36me2.bed',
 'GSM701931_H3K36me3.bed',
 'GSM701932_H3.bed']
names = [
    'H3K4me2', 'H3K4me3', 'H3K9ac','H3K9me2',
    'H3K18ac', 'H3K27me1','H3K27me3','H3K36me2','H3K36me3','H3'
]
chrlen_path = ppfix + 'genomes/chrlen.tab'
genome_path = ppfix + 'genomes/TAIR8_GFF3_genes.gff'

exp = rawdata.Experiment([prefix + x for x in mod_files[:n]], names[:n])
exp.preprocess(chrlen_path, genome_path)
exp.get_binary()

# Testing exporting and reloading
exp.export('./test.p')
ret = pickle.load(open('./test.p', 'rb'))
young = mdata.Mdata(ret.binary, ret.names)
young.addLabels(ret.genes, 'all')
col_names = ['defense', 'develop', 'stress', 'stimulus', 'flowerN', 'salt', 'heat']
glst_paths = ['defense.glst','development.glst','stress.glst'\
              ,'stimulus.glst','floweringN.glst','salt.glst','heat.glst']
inds = [0,1,1,1,0,0,0]
dfLst = [ana.selectDataGL(ret.binary, gprefix+glst_paths[i], inds[i])
         for i in range(len(glst_paths))]

gLsts = [x.index.tolist() for x in dfLst]
for i in range(len(col_names)):
    young.addLabels(gLsts[i], col_names[i])

young.findPatterns(2)
young.findPScores(mode='ar',ratios=young.ratios)
