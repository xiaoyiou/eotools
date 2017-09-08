# This module has the class for wrapping the
# histone modification data
import analysis as ana
import pandas as pd
import numpy as np
from scipy.stats import binom, norm


class Mdata:
    """
    Mdata is the encapsule for the histone modificaiton
    data where the only binary data is stored for now.
    we also add the labels
    """
    def __init__(self,binary, mod_names):
        self.binary=binary.copy()
        self.mods= mod_names
        self.N = binary.shape[0]
        self.M = len(mod_names)
        self.labels = pd.DataFrame(index=binary.index)
        self.ratios = None
        self.rPvalue = None
        self.classes=[]
        
    def addLabels(self,labels,lName):
        """
        labels is a Series,list or dataframe
        of one column
        lName is the label name :flowering, defense
        """
        self.labels[lName] = 0
        self.labels.loc[labels,lName]=1
        self.classes.append(lName)
        
        
    def calcAllRatios(self):
        """
        After calling addLabels and findPatterns
        we calculate the ratios,rpvalue,and etc
        """
        patterns =[]
        self.lens = self.labels.sum()
        df = pd.DataFrame(index=self.pgMap.columns)
        ndf = pd.DataFrame(index=self.pgMap.columns)
        # ndf is the ratios in negative data points
        for col in self.labels.columns:
            for ind in df.index:
                df.loc[ind,col]=len(\
            self.binary[(self.pgMap[ind]==1) & \
            (self.labels[col]==1)])/float(self.lens[col])
                if col=='all':
                    continue
                ndf.loc[ind,col]=len(\
            self.binary[(self.pgMap[ind]==1) & \
            (self.labels[col]==0)])/float(self.N-self.lens[col])
        self.ratios=df
        self.nratios=ndf
        
        
        
    def calcrPvalue(self,base='all',mode='norm'):
        
        if self.ratios == None:
            self.calcAllRatios()

        df = pd.DataFrame(index=self.ratios.index)
        if mode =='norm':
            df = self.ratios[self.nratios.columns]-\
                 self.nratios
            
            bg = self.ratios.loc[:,base]
            for col in df.columns:
                df.loc[:,col]/=np.sqrt(bg*(1-bg)*\
                    (float(1)/self.lens[col]+\
                    float(1)/(self.N-self.lens[col])))
            self.diff=df.copy()
            for col in df.columns:
                df.loc[:,col]=norm.cdf(df.loc[:,col])
                
            
        elif mode=='binom':
            for col in self.labels.columns:
                print col
                temp = self.ratios.loc[:,col]\
                       *self.lens[col].astype(int)
                df[col] = binom.cdf(temp,\
                self.lens[col],self.ratios[base])

        self.rPvalue = df
    
        
        
    
    def findPatterns(self,min_support):
        """
        Find the patterns with the min_support
        """
        
        
        report = ana.xfindFS(self.binary,min_support)
        temp=ana.findFSGenes(report,self.binary)
        self.pgMap = pd.DataFrame(index=self.binary.index)

        for p in temp:
            self.pgMap[p]=0
            self.pgMap.loc[temp[p],p]=1


    def findPScores(self,mode='ar',ratios=False):
        """
        """
        if self.rPvalue==None:
            self.calcrPvalue()

        r =self.ratios.loc[:,self.ratios.columns[1:]]
        if mode == 'ar':
            if not ratios:
                self.pScore= self.rPvalue.applymap(lambda x: (np.log(x) - np.log(1-x))if x < 1 else 44 )
            else:
                self.pScore=np.log(self.rPvalue).multiply(self.ratios['all'],axis='index')-np.log(1-self.rPvalue)*r
#            self.pScore=self.pScore.apply(lambda row:row*len(row.name),axis=1)
        elif mode =='bs':
            # Use binary streaming classificaiton method
            if not ratios:
                self.pScore= np.log(self.rPvalue)
            else:
                self.pScore=-np.log(1-self.rPvalue).r

        elif mode =='lr':
            self.pScore = np.log(1-self.rPvalue)/np.log(self.rPvalue)
            
    def predictALL(self,binary,mode='norm'):
        """
        reevaluate all genes
        """
        result =None
        m = self.pgMap
        p = self.pScore
        r = self.ratios.iloc[:,range(1,self.ratios.shape[1])]
        if mode =='sum':
            result = binary.apply(lambda row:\
        p.ix[m.columns[m.ix[row.name].nonzero()]].sum() if len(m.ix[row.name].nonzero()[0])>0 else 0,axis=1)
        elif mode == 'maxmin':
            
            result = binary.apply(lambda row:\
        p.ix[m.columns[m.ix[row.name].nonzero()]].max()+p.ix[m.columns[m.ix[row.name].nonzero()]].min()if len(m.ix[row.name].nonzero()[0])>0 else 0,axis=1)

        elif mode == 'norm':
            result = binary.apply(lambda row:\
    p.ix[m.columns[m.ix[row.name].nonzero()]].sum()/r.ix[m.columns[m.ix[row.name].nonzero()]].sum() if len(m.ix[row.name].nonzero()[0])>0 else 0,axis=1)
        elif mode =='max':
            result = binary.apply(lambda row:\
        p.ix[m.columns[m.ix[row.name].nonzero()]].max()if len(m.ix[row.name].nonzero()[0])>0 else 0,axis=1)

        elif mode == 'ratio':
            result = binary.apply(lambda row:\
            r.ix[m.columns[m.ix[row.name].nonzero()]].sum() if len(m.ix[row.name].nonzero()[0])>0 else 0,axis=1)


        return result
          
    
