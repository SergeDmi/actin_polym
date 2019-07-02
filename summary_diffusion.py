####### PACKAGES
import math
from numpy import *
import pandas as pd
import sys,os
#from statistical_tools import multilinear_regression as linfit


def make_array_from_str(value):
    i=value.find('[')
    j=value.find(']')
    value=value[i+1:j]
    elements=value.split(',')
    return array([float(val) for val in elements])

def compute_dist_from_row(row):
	pos1=make_array_from_str(row['f1.position'])
	pos2=make_array_from_str(row['f2.position'])
	ori1=make_array_from_str(row['f1.orientation'])
	ori2=make_array_from_str(row['f2.orientation'])

	return distance(pos1,ori1,pos2,ori2)


def compute_ratio_from_row(row):
	#pos1=make_array_from_str(row['f1.position'])
	#pos2=make_array_from_str(row['f2.position'])
	#ori1=make_array_from_str(row['f1.orientation'])
	#ori1=make_array_from_str(row['f2.orientation'])

	return row['f1']/row['f2']

def  distance(pos1,ori1,pos2,ori2):
	m1=ori1/2.0
	m2=(2*(pos2-pos1)+ori2)/2.0-m1
	return sqrt(m2[0]**2.0+m2[1]**2.0)



if __name__=="__main__":
	args=sys.argv[1:]
	options=[]
	folders=[]

	for arg in args:
		if arg.find('=')>0:
			options.append(arg)
		else:
			folders.append(arg)

	res_name='res.csv'
	for opt in options:
		if opt.startswith('name='):
			res_name=opt[5:]

	fname=os.path.join(folders[0],res_name)

	nf=len(folders)
	ex_data=pd.read_csv(fname,squeeze=True, index_col=0)
	keys=ex_data.columns.to_list()
	datas=pd.DataFrame(columns=keys)

	for folder in folders:
		fname=os.path.join(folder,res_name)
		ex_data=pd.read_csv(fname,squeeze=True,index_col=0)
		ex_data.index=[os.path.basename(folder)]
		#index=[os.path.basename(folder)]
		datas=datas.append(ex_data,)
		#fname.name=os.path.basename(folder)

	datas['distance']=datas.apply(lambda row: compute_dist_from_row(row), axis=1)
	datas['ratio']=datas.apply(lambda row: compute_ratio_from_row(row), axis=1)
	# Now we create stuff
	dists=sort(unique(datas['distance'].to_numpy()))
	nd=len(dists)
	keys=['distance','fil1','fil2','ratio','std_n1','std_n2','std_ratio','n']
	by_dist=pd.DataFrame(columns=keys)

	for i,dist in enumerate(dists):
		select=datas[datas['distance']==dist]
		n1=select['f1'].to_numpy()
		n2=select['f2'].to_numpy()
		ratio=select['ratio'].to_numpy()
		data={ 'distance' : dist , 'fil1' : mean(n1) , 'fil2' : mean(n2) , 'ratio' : mean(ratio) , 'std_n1' : std(n1) , 'std_n2' : std(n2) , 'std_ratio' : std(ratio), 'n':len(n1) }
		by_dist.loc[i]=data

	# Exporting gathered  file
	export_file='export.csv'
	datas.to_csv(export_file)
	dist_file='by_distances.csv'
	by_dist.to_csv(dist_file)
