"""
This module is the library for conducting analysis and visulizations
All these methods are still hard coded for arabidopsis not flexible
This line is to test the git
"""
import csv
import numpy as np
import math
import pandas as pd
from yxtools import dictInsert
from scipy.spatial.distance import cosine
from scipy.stats import poisson as pois
from scipy.stats import binom
import pylab as pl
import re
import matplotlib.pyplot as plt
from pymining import itemmining


def findCommonGenes(paths):
    
    geneLists=[]
    for path in paths:
        glst=[]
        with open(path,'rb') as file:
            for line in file:
                tokens = line.rstrip().split('\t')
                if tokens[1].startswith('AT'):
                    glst.append(tokens[1])
        file.close()
        geneLists.append(set(glst))

    return geneLists

def getMods(path , countType=np.float):
    result= pd.read_csv(path,sep='\t',names=['chrom','left','right','gene','modid','chromb','start','end','count'],dtype={'left':np.int32,'right':np.int32,'count':countType,'start':np.int32,'end':np.int32})
    result['modid']-=1
    return result
    

def findMax(df):
    """
    This will find the max of each modification on each gene
    and reuturn as a df
    """

    return pd.DataFrame({'count':df.groupby(['gene','chrom','modid'])['count'].max()}).reset_index()
                        

def findBinary(df,baseline,names,alpha):
    df['appear']=df[df['count']> alpha*baseline[(names[df['modid']-1],df['chrom'])]]
    return df

def visMatrix(data, names,cmin=0,cmax=1):
    """
    This is used to visualize the correlation matrix
    """

    column_labels = names
    row_labels = names
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(data,label='big')

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(row_labels, minor=False,rotation=30)
    ax.set_yticklabels(column_labels, minor=False)
    
    

    fig.colorbar(heatmap,ticks=np.arange(cmin,cmax,(cmax-cmin)/float(10)))
    plt.show()

def findCos(df):
    """
    the data frame input here is the binary format of the 
    dataset
    """

    N=df.shape[1]-1
    data = np.zeros((N,N),dtype=float)

    for i in range(N):
        for j in range(N):
            data[i,j]=1-cosine(df[i],df[j])

    return data

def findExc(df):
    
    L=float(df.shape[0])
    N=df.shape[1]-1
    data = np.zeros((N,N),dtype=float)

    for i in range(N):
        for j in range(N):
            data[i,j]=df[((df[i]==1)&(df[j]==0))\
                         |((df[j]==1)&(df[i]==0))][i].count()/L

    return data

def xfindFS(df,sup):
   return __findFS(df, sup)


def __findFS(df,sup):
    """
    This function is used to find the "frequent" itemsets
    for the database
    """
    if type(df) == pd.core.frame.DataFrame:
        L=df.shape[0]
        tss=[]
        for i in range(L):
            tss.append(df.columns[df.iloc[i]>0].tolist())
        relim=itemmining.get_relim_input(tss)
    else:
        relim = itemmining.get_relim_input(df)
    report = itemmining.relim(relim , min_support =sup )    
    return report 

def __visFS(report,names,ilen_thresh,maxReport,filter,N):
    len_sorted = sorted(report.items(),key = \
                        lambda s: len(s[0]))[::-1]
    sup_sorted = sorted(report.items(), key=\
                        lambda s: s[1]/N)[::-1]

    print "+++++++++++++++++++++++++++++++++++++++++++++++++++"
    print "Patterns sorted by pattern length"
    i=0
    for (fset,count) in len_sorted:
        mods=[names[T] for T in fset]

        if filter !=None:
            if len(set(mods).intersection(set(filter)))==0:
                continue



        i+=1
        if maxReport>0 and i==maxReport:
            break
        
        if len(mods)<ilen_thresh:
            continue 
        print mods,'\t',count
    print "++++++++++++++++++++++++++++++++++++++++++++++++"
    print "Patterns sorted by frequency"
        
    i=0
    for (fset, count)in sup_sorted:
        mods=[names[T] for T in fset]

        if filter !=None:
            if len(set(mods).intersection(set(filter)))==0:
                continue
                

        i+=1
        if maxReport>0 and i==maxReport:
            break
        if len(mods)<ilen_thresh:
            continue 

        print mods,'\t',"{0:.2f}".format(count/N)

def findFS(df,names,ilen_thresh=2,min_support=2,maxReport=50,filter=None,ratio=False):

    """
    Visualizer of the frequent itemsets report
    """
    
    
    report=__findFS(df,min_support)

    # Now we want to recover names and sort the result
    if ratio:
        __visFS(report,names,ilen_thresh,maxReport,filter,float(df.shape[0]))
    else:
        __visFS(report,names,ilen_thresh,maxReport,filter,1)

    
def selectData(df,go,action,goid):
    
    """
    select the partial data based on the action and goid

    """
    
    temp=df.copy()
    temp['gene']=df.index
    data = pd.merge(temp,go[(go.action==action)&\
                            (go.goid==goid)],on='gene')
    data=data[temp.columns]

    print "partial data has "+\
        str(data.shape[0])+' data points' 
    
    return data.set_index('gene')



def selectDataGL(df, path,ind=1,sep='\t',start=0):
    
    """
    This function selects the data based on the glst file 
    
    """

    # read the file and create the pandas df
    # because of male formatting of the file, we cannot directly load it

    glst=[]
    i =0
    with open(path,'rb') as file:
        for line in file:
            i+=1
            if i<=start:
                continue
            tokens = line.rstrip().split(sep)
            if tokens[ind].startswith('AT'):
                glst.append(tokens[ind])
        file.close()

    return df.ix[df.index.intersection(glst)]
    
        
  

    print " partial data has " +\
        str( data.shape[0]) + ' data points'

    return data.set_index('gene')
    
def compareFS(df1,df2,names,ilen_thresh=1,min_support=2,maxReport=50,filter=None,norm=True, absolute=True,disjoint=False):
    
    """
    compare the partial data with all genes
    """

    if disjoint:
        df1 = df1.ix[df1.index.difference(df2.index)]
        df2 = df2.ix[df2.index.difference(df1.index)]


    reportBG = __findFS(df1,min_support)
    
    reportPD = __findFS(df2,min_support)
    
    # generating the difference report

    data = dict()
    
    for key in reportPD:
        if not norm:
            if key in reportBG:
                
                data[key]=reportPD[key]-reportBG[key]
            else:
                data[key]=reportPD[key]

        if norm:
            if key in reportBG:
                data[key]= reportPD[key]/float(df2.shape[0])-reportBG[key]/float(df1.shape[0])
            else:
                data[key]= reportPD[key]/float(df2.shape[0])

        if absolute:
            for key in data:
                data[key]=abs(data[key])
    __visFS(data,names,ilen_thresh,maxReport,filter)



def compareCos(df1,df2,disjoint=True):

    if disjoint:
        df1=df1.ix[df1.index.difference(df2.index)]
        df2=df2.ix[df2.index.difference(df1.index)]
        
    
    data1=findCos(df1)
    data2=findCos(df2)
    
    return data2-data1
    


def getEnrich(df,genes,N):
    """
    Get the modification enrichments for a particular gene or a set of genes
    """

    glst = df['gene'].unique()
    result = dict()
    
    for gene in genes:
        if gene in glst:
            result[gene]=dict()
            for i in range(N):
                result[gene][i]=df[(df['gene']==gene)&\
                                   (df['modid']==i)]['count'].tolist()
                

    return result

def getEnrich2(df,grps):

    result = dict()

    for key in grps.keys():
        result[key]=df.ix[grps[key],'count'].tolist()
    return result
            


def getData(df):
    """
    change the weird pandas format into dictionary format data
    
    """
    gp = df.groupby(['gene','modid'])

    return gp.groups

    
def visEnrich(data,gene,left,right,mod_names,N,modLst):
    
    data = getEnrich(data,[gene],len(mod_names))

    x = np.linspace(left,right,N)
    for i in modLst:
        pl.plot(x,data[gene][i],label=mod_names[i])
        
    pl.legend(loc='upper right')        
    plt.show()


def visEnrich2(data,gene,left,right,mod_names,N,modLst,baseline=None):
    # build a rectangle in axes coords
    

    
    data = getEnrich(data,[gene],len(mod_names))

    x = np.linspace(left,right,N)

    f,ax=plt.subplots(len(modLst),sharex=True, sharey=True,figsize=(16,9))

    colors = plt.get_cmap('jet')(np.linspace(0, 1.0, len(modLst)))
    ii=0
    
    for i in modLst:
        ax[ii].bar(x,data[gene][i],width=(right-left)/N,color=colors[ii],align='center')
       
        ax[ii].text(0.1, 0.9,mod_names[i], ha='center', va='center', transform=ax[ii].transAxes,weight='bold')
       
        if baseline!= None:
            ax[ii].plot([left,right],[baseline[i],baseline[i]])
        ax[ii].autoscale(tight=True)
        ii+=1 
        
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.xlabel('Chromoson Position Relative to Gene', fontsize=18)
    plt.show()

    

def calcThresh(data,N,power = 4):
    """
    This function is used to calculate the thresholds of different modifications
    N is the number of modifications

    """
    
    temp = [0.0]*N
    count = [0.0]*N
    distro = []

    pValue = math.pow(10,-power)
    
    for (PL,mod_num) in data:
        temp[mod_num]+=sum(data[PL,mod_num])
        count[mod_num]+=len(data[PL,mod_num])


    for i in range(N):
        distro.append(pois(temp[i]/count[i]))

    return [x.isf(pValue) for x in distro]


def calcThreshAll(data,power=4):
    """
    This function returns the baseline considering entire genome 
    instead of gene regions only
    """

    # The mean value for each modification
    temp = data.groupby(['modid']).mean()['count'].tolist()
    pValue = math.pow(10,-power)
    distro = [pois(x) for x in temp]

    return [x.isf(pValue) for x in distro]

    
def maxBin(data, baseline):
    """
    This is for binarization of the data,
    compatible with old data format and baseline for
    """

    if type(baseline)==list:
        baseline = pd.DataFrame({'modid':range(len(baseline)),'mean':baseline})
    
    
    maxes=findMax(data)

    
    binary=pd.merge(maxes,baseline)
    binary['appear']=0
    binary.ix[binary['count']>binary['mean'],'appear']=1
    binary=binary.pivot(index='gene',columns='modid',values='appear')
    return binary


def __findFSGenes(binary):
    
    """
    This function find the genes for each pattern (frequent Itemset)
    and return the result in the dictionary
    """

    N = binary.shape[1]
    return  binary.groupby(range(N)).groups


def findFSGenes(report,binary):
    N = binary.shape[0]
    eSets = __findFSGenes(binary)
    result = dict()
    
    for key in report:
        temp = [0]*N
        for i in key:
            temp[i]=1
        result[key] = []
        for k in eSets:
            if all([x-y>=0 for x,y in zip(k,temp)]):
    
                result[key]+=(eSets[k].tolist())

    return result


def compareFSS(dfLst,min_support,mod_names,base=0):
    """
    This function determines get the results by ratios
    The last item in the list is the gene names lists related to 
    the pattern
    ===================================
    This method returns the two outputs: 1. the result of frequent itemset, 2. the names asso. with 
    patterns
    """
    
    Ns=[]
    for df in dfLst:
        Ns.append(float(df.shape[0]))

    report = __findFS(dfLst[base],min_support)
    result = dict()
    for key in report:
        result[key]=[report[key]/Ns[base]]
    
    for i in range(1,len(dfLst)):
        report = __findFS(dfLst[i],min_support)
        
        for key in result:
            if key in report:
                result[key].append(report[key]/Ns[i])
            else:
                result[key].append(0)

    gNames = findFSGenes(report,dfLst[base])


    return (result,gNames)


def visCompareFSS(result, mod_names,to_file=None):

    if to_file == None:
        for key in result:
            line = result[key]
            print [mod_names[T] for T in key],'\t','\t'.join(map(str,line))
    
    else:
        with open(to_file, 'wb') as f:
            writer = csv.writer(f,delimiter='\t')
            for key in result:
                line = [[mod_names[T] for T in key]]
                line += result[key]
                writer.writerow(line)
            f.close()
        
    
def createFssDF(result,mod_names,col_names,base,lens):
    """
    This funciton will transform the frequent itemsets
    compairison result (including ratio change and P values)

    Parameters:
    result : the dict baesd result of different lists of genes
    mod_names: the modification names
    col_names: the column names of the result
    base: the reference column
    lens: the lengths of the different lists
    """

    N = len(col_names)
    col_inds = range(N)
    data = dict()
    patterns = []
    for i in col_inds:
        j =0
        for key in result:
            if i ==0:
                patterns.append(key)
            if j ==0:
                data[col_names[i]]=[result[key][i]]
            else:
                data[col_names[i]].append(result[key][i])
            j+=1

    df = pd.DataFrame.from_dict(data).set_index([patterns])
    refCol = df.iloc[:,base]
    diff_data = dict()
    

    for i in col_inds:
        if i == base:
            continue
        diff_data[col_names[i]] = (df.iloc[:,i]\
            -df.iloc[:,base])/df.iloc[:,base]

    df2 = pd.DataFrame.from_dict(diff_data).set_index([patterns])

    P_data =dict()
    for i in col_inds:
        if i == base:
            continue
        k = df.iloc[:,i]*lens[i]
        k = k.astype(int)
        P_data[col_names[i]] = binom.cdf\
            (k,lens[i],df.iloc[:,base])

    df3 = pd.DataFrame.from_dict(P_data).set_index([patterns])
    
    return (df,df2,df3)

def findPatterns(dfP,lowT=0.05,highT=0.05,\
                 ratios=None,threshP=0.05,threshN=0.05):
    """
    dfP is the rPvalue data frame
    patterns is the int->pattern book
    lowT and highT is the pValue thresholds for left/right tail
    ratios: default to be none. If not, we are also considering 
    the absolute ratio in addition to P values
    
    tresh: the ratio threshold
    """
    
    # The positively related patterns 
    pPtns = dict()
    nPtns = dict()
    # The negatively related patterns

    col_names = dfP.columns.tolist()
    
    for name in col_names:

        tempN = dfP.loc[dfP[name]<lowT]
        if ratios is not None:
            tempN=tempN.ix[ratios.loc[ratios['all']>\
        threshN].index.intersection(tempN.index)]
        nlst = tempN.sort_values(by=[name]).index.tolist()

        tempP = dfP.loc[1-dfP[name]<highT]

        if ratios is not None:
            tempP= tempP.ix[ratios[ratios[name]\
                >threshP].index.intersection(tempP.index)]
        plst = tempP.sort_values(by=[name],ascending=False).index.tolist()
        pPtns[name] = plst
        nPtns[name] = nlst

    return (pPtns,nPtns)


def findPatternGlst(ptns,dfLst,mod_names,col_names):

    N=len(mod_names)
    
    result = dict()
    
    for i in range(1,len(col_names)):
        key = col_names[i]
        result[key]=dict()
        eSets = __findFSGenes(dfLst[i])
        for p in ptns[key]:
            temp = [0]*N
            for i in p:
                temp[i]=1
            result[key][p]=[]
            for k in eSets:
                if all([x-y>=0 for x,y in zip(k,temp)]):
                    result[key][p]+=(eSets[k])
    return result
                    

def summarizePtns(pPtns,nPtns):
    """
    summarizePtns will summarize the patterns into sccinct rules 
    """

    # the edges are going from super set to subset
    ppatterns = dict()
    # the edges are going from subset to supber set
    npatterns = dict()

    for key in pPtns:
        sbNet = dict()
        spNet = dict()
        A = pPtns[key]
        B = nPtns[key]
        C = list(set(A+B))
        tempP=[]
        tempN=[]
        for i in C:
            for j in C:
                if i==j:
                    continue
                if i.issubset(j):
                    dictInsert(sbNet,i,j)
                elif i.issuperset(j):
                    dictInsert(spNet,i,j)
        for x in A:
            if x not in sbNet:
                tempP.append(x)
                continue
            plst = sbNet[x]
            print "#####",x,plst
            if all([y in B for y in plst]):
                tempP.append(x)

        for x in B:
            if x not in spNet:
                tempN.append(x)
                continue
            plst = spNet[x]
            if all([y in A for y in plst]):
                tempN.append(x)
        
        ppatterns[key]=tempP
        npatterns[key]=tempN

    return (ppatterns,npatterns)
                
    
def readGeneInfo(path):
    res = {}
    # each entry record contains [direction,type]
    IDreg = re.compile('ID=(.*?);')
    NoteReg = re.compile('Note=(.*?);')
    with open(path,'rb') as f:
        for line in f:
            tokens = line.rstrip().split('\t')
            if tokens[2] == 'gene':
                ID = IDreg.findall(tokens[-1])[0]
                note = NoteReg.findall(tokens[-1])[0]
                res[ID]=[True if tokens[6] == '+' else False,\
                         note]
    return res

def cleanDir(gdata, geneInfo):
    for gene, mod in gdata:
        if gene in geneInfo:
            if geneInfo[gene][0] == False:
                gdata[(gene,mod)] =gdata[(gene,mod)][::-1]

    


    
            
    
