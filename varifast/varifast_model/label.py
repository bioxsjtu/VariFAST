import argparse
import pandas as pd
def get_label(k_all,gd):       #get label of train data
	#gd = pd.read_csv(good_mut)
	k_all['label'] = 0
	m = list(k_all['Start'])
	m = [int(i) for i in m]
	k_all['Start'] = m
	m = list(gd['Start'])
	m = [int(i) for i in m]
	gd['Start'] = m
	k_all.sort_values(by=['Start','sample','Chr','Ref','Alt'],inplace=True)
	gd.sort_values(by=['Start','sample','Chr','Ref','Alt'], inplace=True)
	
	gd.index = range(len(gd.index))
	k_all.index = range(len(k_all.index))
	c1 = []
	c2 = []
	for i in range(len(gd.index)):
		m = gd.ix[i,'Chr']+str(gd.ix[i,'Start'])+gd.ix[i,'Ref']+gd.ix[i,'Alt']+gd.ix[i,'sample']
		c1.append(m)
	for j in range(len(k_all.index)):
		m = k_all.ix[j,'Chr'] + str(k_all.ix[j,'Start']) + k_all.ix[j,'Ref'] + k_all.ix[j,'Alt'] + k_all.ix[j,'sample']
		c2.append(m)
	k_all.index = c2
	gd.index = c1
	c3 = list(set(c1)&set(c2))
	k_all.ix[c3,'label'] = 1
	#k_all.to_csv('all_mut_feature4.csv',index=False)
	return k_all

def label():
	parser = argparse.ArgumentParser(description='script for igv')
	parser.add_argument('--fvar',required=True,help ='the file name of true variant list')
	parser.add_argument('--file',required=True ,help = 'file name of all variant')
	parser.add_argument('-o',default='labeled_result.csv',type=str,help='output file name of labeled file')
	args = parser.parse_args()
	gd = pd.read_csv(args.fvar,low_memory=False)
	k_all = pd.read_csv(args.file)
	k_all = get_label(k_all, gd)
	k_all.to_csv(args.o,index=False)
	return
if __name__=='__main__':
	label()
