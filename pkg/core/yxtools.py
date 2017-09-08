#! /usr/bin/python
# Filename: yxtools.py
# This module finally appears from nowhere to deal with the nasty tab files


def dictInsert(dct,key,data):
    """
    This function takes the one dict, one key and one data as argument
    it checks whether the key exists in dict:
    if not,a new list containing the data will be registered to the key in dict
    otherwise, the data will be appended to the list
    Since in python dct is mutable object, so no need to return
    """
    if not key in dct.keys():
        dct[key]=[data]
    else:
        assert len(dct[key])>0
        dct[key].append(data)
        

# getList is the function used to get one-column list of strings
# numbers, IDs.
def getList(path):
    result=[];
    file=open(path,'rb');
    for line in file:
        line=line.rstrip();
        result.append(line)
        
    file.close()
    return result

# getNames is the function returns the mapping of names/ids.
# it returns a tuple of two dictionaries that map the relations between names
# The input file should have at least two columns, and the mapping column 
# should be specified in the inds arguments, the default value is [0,1]
# WARNING: if the input file doesn't gaurantee a one-to-one mapping.
# The returned dictionaries are not reliable, either.
def getNames(path,inds=[0,1],delim=None,capital=False):
    dictA=dict()
    dictB=dict()
    file=open(path,'rb');
    for line in file:
        line = line.rstrip();
        ind1=inds[0];
        ind2=inds[1];
        if delim==None:
            tokens=line.split();
        else:
            tokens=line.split(delim);
        
        x=tokens[ind1];
        y=tokens[ind2];

        if capital:
            x=x.upper()
            y=y.upper()
        dictA[x]=y;
        dictB[y]=x;
    file.close()
    return (dictA,dictB);

# this function returns the selected columns as dictionary

def getTabData(path,key=0,delim=None,capital=False):
 
    file=open(path,'rb');
    result=dict()
    
    for line in file:
        line=line.rstrip();
        if delim==None:
            tokens=line.split();
        else:
            tokens=line.split(delim);
        
        id=tokens[key]
        if capital:
            id=id.upper();
        result[id]=tokens[key+1:len(tokens)];
    file.close()
    return result


def flatenDict(data):
    '''
    flatenDict transforms a nested dictionary (dict of dicts)
    to a single level dictionary with tuples of two keys 
    {(key1,key2):value}
    '''
    result=dict()
    
    for key in data.keys():
        for key2 in data[key].keys():
            result[(key,key2)]=data[key][key2]

    return result

def flatenDictList(data):
    '''
    flatenDictList transforms dict of lists into a 
    a dict of 
    {(key1,key2):value}
    '''
    result=dict()
    for key in data.keys():
        for elem in data[key]:
            result[(key,elem)]=1

    return result
            



