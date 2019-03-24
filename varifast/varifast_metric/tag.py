import argparse
import pandas as pd

def get_tag(feature,vaf=0.2,threshold_h=0.9,threshold_e=0.9,threshold_ni=0.9,threshold_nd=0.9,threshold_dir=0.9,threshold_sse=0.9,threshold_si=0.9,threshold_hdr = 1,threshold_mm=1,threshold_lm=0.2,threshold_mv=0.2,paired=False,rate_tn=0.1):
	result = pd.DataFrame()
	feature.index = range(len(feature.index))
	tag_all =[]
	score_all = []
	l1 = ['h','e','dir','sse','r','ri','mm','ao','si']
	l2 = ['mv','hdr']
	l3 = ['lm','lcr','vafr','ni','nd']
	for i in range(len(feature.index)):
		score = 0
		tag = ''
		if feature.ix[i,'vafr']>=1-vaf:
			tag = tag+'LVF:'
		if feature.ix[i,'lcr']>0:
			tag = tag+'LCT:'
		else:
			if feature.ix[i,'lm']>=threshold_lm and feature.ix[i,'lm_num']>=2:
				tag = tag + 'LM:'
			if feature.ix[i,'mm']>=threshold_mm: #and feature.ix[i,'MM_num']>=10:
				tag = tag + 'MM:'
			if feature.ix[i,'hdr']>=threshold_hdr:
				tag += 'HDR:'
			if feature.ix[i,'h']>=threshold_h:
				tag += 'H:'
			if feature.ix[i,'e']>=threshold_e:
				tag += 'E:'
			if feature.ix[i,'ni']>=threshold_ni:
				tag += 'NI:'
			if feature.ix[i,'nd']>=threshold_nd:
				tag += 'ND:'
			if feature.ix[i,'dir']>=threshold_dir:
				tag += 'D:'
			if feature.ix[i,'ao']==1:
				tag += 'AO:'
			if feature.ix[i, 'sse']>=threshold_sse:
				tag += 'SSE:'
			if feature.ix[i,'mv']>=threshold_mv and feature.ix[i, 'mv_num']>=2:
				tag += 'MV:'
			if feature.ix[i,'r']==1:
				tag += 'RR:'
			if feature.ix[i,'ri']==1:
				tag +='RI:'
			if feature.ix[i,'si']>=threshold_si:
				tag += 'SI:'
		if paired:
			if feature.ix[i,'lncr']==1:
				tag +='NCN:'
			if feature.ix[i,'lncr']>0 and feature.ix[i,'lncr']<1:
				tag +='LCN:'
			if feature.ix[i,'nvaf']!='NAN':
				if float(feature.ix[i,'nvaf'])>=vaf*rate_tn and feature.ix[i,'tn_num']>=2:
					tag +='TN:'
		tag_all.append(tag)
		for j in l1:
			score += feature.ix[i,j]
		for j in l2:
			score += feature.ix[i,j]*2.0
		for j in l3:
			score += feature.ix[i,j]*3.0
		if paired:
			if feature.ix[i,'nvaf']!='NAN':
				score += float(feature.ix[i,'nvaf'])*2.0
			score += feature.ix[i,'lncr']*3.0
		score_all.append(score)
	result['index'] = feature['index']
	result['Chr'] = feature['Chr']
	result['Start'] = feature['Start']
	result['Alt'] = feature['Alt']
	result['Ref'] = feature['Ref']
	result['sample'] = feature['sample']
	result['score'] = score_all
	result['tag'] = tag_all
	#feature['score'] = score_all
	#feature['tag'] = tag_all
	return feature

def tag():
	parser = argparse.ArgumentParser(description='script for tag')
	parser.add_argument('--file',required=True ,type=str,help = 'feature file name')
	parser.add_argument('--vaf',required=True,type=float,help ='at least vaf of variant')
	parser.add_argument('--paired',action='store_true',help='somatic')
	parser.add_argument('--h',default=0.9,type=float,help='threshold for H')
	parser.add_argument('--e',default=0.9,type=float,help='threshold for E')
	parser.add_argument('--ni',default=0.9,type=float,help='threshold for NI')
	parser.add_argument('--nd',default=0.9,type=float,help='threshold for ND')
	parser.add_argument('--dir',default=0.9,type=float,help='threshold for D')
	parser.add_argument('--sse',default=0.9,type=float,help='threshold for SSE')
	parser.add_argument('--mm',default=1,type=float,help='threshold for MM')
	parser.add_argument('--lm',default=0.2,type=float,help='threshold for LM')
	parser.add_argument('--mv',default=0.2,type=float,help='threshold for MV')
	parser.add_argument('--hdr',default=1,type=float,help='threshold for HDR')
	parser.add_argument('--si',default=1,type=float,help='threshold for SI')
	parser.add_argument('--tn',default=0.1,type=float,help='threshold for TN')
	parser.add_argument('-o',default='tag.csv',type=str,help='output file name')
	args = parser.parse_args()
	feature = pd.read_csv(args.file)
	tag = get_tag(feature,vaf=args.vaf,paired=args.paired,threshold_h=args.h,threshold_e=args.e,threshold_ni=args.ni,threshold_nd=args.nd,threshold_dir=args.dir,threshold_sse=args.sse,threshold_si=args.si,threshold_hdr=args.hdr,threshold_mm=args.mm, threshold_lm=args.lm, threshold_mv=args.mv, rate_tn = args.tn)
	tag.to_csv(args.o,index=False, sep='\t')
	print('tag get and all finished')

if __name__=='__main__':
	tag()
