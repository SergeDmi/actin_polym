#!/usr/bin/env python3
####### PACKAGES
import math
from numpy import *
import pandas as pd
import sys,os

# @TODO : make array from string : make recurent

def make_array_from_str(value):
    # Makes an array from a string
    i=value.find('[')
    j=value.find(']')
    value=value[i+1:j]
    elements=value.split(',')
    return array([float(val) for val in elements])

def compute_dist_from_row(row):
    # Computes distance between two filaments f0 and f1
    pos1=make_array_from_str(row['f0.position'])
    pos2=make_array_from_str(row['f1.position'])
    ori1=make_array_from_str(row['f0.orientation'])
    ori2=make_array_from_str(row['f1.orientation'])
    return distance(pos1,ori1,pos2,ori2)

def compute_volume_from_row(row):
    # Computes volume from box sizes
    box=make_array_from_str(row['box'])
    return box[0]*box[1]*box[2]

def compute_ratio_from_row(row):
    # Computes ratio of reaction events
    return row['f0']/row['f1']

def distance(pos1,ori1,pos2,ori2):
    # Actually computes distance between two filaments
    # Normalizes the distance in units of monomer size
    m1=pos1+get_offset_from_orientation(ori1)
    m2=pos2+get_offset_from_orientation(ori2)-m1
    return sqrt(m2[0]**2.0+m2[1]**2.0)/sqrt(ori1[0]**2.0+ori1[1]**2.0)

def get_offset_from_orientation(ori):
    # offset between filament position and filament center (depends on orientation)
    bori=abs(ori)
    if sum(bori)==2:
        if bori[0]==2:
            return array([1.5,-0.5])
        else:
            return array([0.5,-1.5])
    else:
        if bori[0]==2.0:
            return array([0.5,0])
        else:
            return array([0,-0.5])


if __name__=="__main__":
    """
        summary_diffusion.py : arranges the results from several simulations

        USAGE :
        python summary_diffusion.py folder [FOLDERS] [name=NAME] [-sort KEY [KEYS]]

        EXAMPLES :
        python summary_diffusion.py run* name=res.csv

        Makes a summary of the res.csv files in folders run*

        python summary_diffusion.py run* name=res.csv sort=distance sort=volume

        Makes a summary of the res.csv files in folders run*
        Outputs results sorted by distance in by_distances.csv
        Outputs results sorted by volume in by_volumes.csv

    """
    args=sys.argv[1:]
    options=[]
    folders=[]
    to_sort=[]

    for arg in args:
        if arg.find('=')>0:
            options.append(arg)
        elif arg.startswith('-'):
            options.append(arg)
        else:
            folders.append(arg)

    # Name of the global output file
    res_name='res.csv'
    for opt in options:
        if opt.startswith('name='):
            res_name=opt[5:]
        if opt.startswith('sort='):
            to_sort.append(opt[5:])

    # Checking what's in the results to prepare pandas datagrame
    fname=os.path.join(folders[0],res_name)
    ex_data=pd.read_csv(fname,squeeze=True, index_col=0)
    keys=ex_data.columns.to_list()
    datas=pd.DataFrame(columns=keys)

    # Reading all folders and assembling the dataframe
    nf=len(folders)
    for folder in folders:
        fname=os.path.join(folder,res_name)
        ex_data=pd.read_csv(fname,squeeze=True,index_col=0)
        ex_data.index=[os.path.basename(folder)]
        datas=datas.append(ex_data,)


    # Specific measurements derived from dataframe values
    datas['distance']=datas.apply(lambda row: compute_dist_from_row(row), axis=1)
    datas['volume']=datas.apply(lambda row: compute_volume_from_row(row), axis=1)
    datas['ratio']=datas.apply(lambda row: compute_ratio_from_row(row), axis=1)

    # Exporting gathered  file
    export_file='export.csv'
    datas.to_csv(export_file)

    # Now we sort data with keys defined by user ! Exciting !
    for key in to_sort:
        choices=sort(unique(datas[key].to_numpy()))
        keys=[key,'fil1','fil2','ratio','std_n1','std_n2','std_ratio','n']
        summary=pd.DataFrame(columns=keys)

        for i,choice in enumerate(choices):
            select=datas[datas[key]==choice]
            n1=select['f0'].to_numpy()
            n2=select['f1'].to_numpy()
            ratio=select['ratio'].to_numpy()
            data={ key : choice , 'fil1' : mean(n1) , 'fil2' : mean(n2) , 'ratio' : mean(ratio) , 'std_n1' : std(n1) , 'std_n2' : std(n2) , 'std_ratio' : std(ratio), 'n':len(n1) }
            summary.loc[i]=data


        summary_file='by_%ss.csv' %key
        summary.to_csv(summary_file)
