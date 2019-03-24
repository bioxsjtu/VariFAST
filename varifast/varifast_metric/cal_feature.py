import argparse
import pandas as pd
import numpy as np
def make_feature(k_all,train=False,paired=False,hdr_thr=0.9,depth_thr=20,normal_depth_thr=20):  #threshold for at least ratio to sum
	feature = pd.DataFrame()
	k_all = k_all.ix[(k_all['depth']!=0)&(k_all['mut_depth']!=0),:]
	k_all.index = range(len(k_all.index))
	ind = k_all['index']
	ind = ['VAR'+str(i) for i in ind]
	feature['index'] = ind
	feature['Chr'] = k_all['Chr']
	feature['Start'] = k_all['Start']
	feature['Alt'] = k_all['Alt']
	feature['Ref'] = k_all['Ref']
	feature['sample'] = k_all['sample']
	feature['vafr'] = 1 - k_all['mut_depth']*1.0/k_all['depth']   #MUT DEPTH RATIO
	def f1(x,thr):
		return 1-(min(x,thr)*1.0/thr)
	feature['lcr'] = k_all.apply(lambda row:f1(row['depth'],depth_thr),axis=1)
	feature['lm'] = k_all['LM']*1.0/k_all['mut_depth'] 
	feature['lm_num'] = k_all['LM']
	feature['mm'] = k_all['low_d_mut']*100             #*1.0/k_all['mut_depth_cancer'] ###low mut ratio
	#feature['MM_num'] = k_all['low_d_mut']*k_all['num_all_nb']
	t = k_all['HDR']
	def get_hdr(t,hdr_thr):
		p = []
		for i in t:
			try:
				if np.isnan(i):
					p.append([])
			except:
				p.append(i.split(':'))
		t_max = []
		tt =[]
		for i in p:
			m_a = 0     #sum ratio
			m2 = 0      #max ratio
			for j in i:
				m = float(j.split(',')[0])
				if m>=hdr_thr:
					m_a += m
				if m>m2:
					m2 = m
			if m2<hdr_thr:
				m_a = m2
			t_max.append(m2)
			tt.append(m_a)
		return t_max, tt
	t_max, tt = get_hdr(t,hdr_thr = hdr_thr)
	feature['hdr'] = tt
	feature['h'] = k_all['head_mut']*1.0/k_all['mut_depth']   #near head
	feature['e'] = k_all['end_mut']*1.0/k_all['mut_depth']     #near end
	feature['ni'] = k_all['near_insert']*1.0/k_all['mut_depth'] #near insert
	feature['nd'] = k_all['near_del']*1.0/k_all['mut_depth']       # near del
	feature['dir'] = k_all['D']*1.0/k_all['mut_depth']              #directional ratio
	feature['ao'] = k_all['AO']
	feature['sse'] = k_all['SSE']*1.0/k_all['mut_depth']         #same head and end ratio
	feature['mv'] = k_all['MV']*1.0/k_all['mut_depth']           #not same as mut
	feature['mv_num'] = k_all['MV']
	feature['r'] = k_all['R']
	feature['ri'] = k_all['RI']
	feature['si'] = k_all['SI']*1.0/k_all['mut_paired_depth']
	if train:
		feature['label'] = k_all['label']
	if paired:
		def f2(a,b):
			if b==0:
				return 'NAN'
			else:
				return a*1.0/b
		#feature['depth_normal'] = k_all['depth_normal']
		feature['tn_num'] = k_all['mut_depth_normal']
		feature['nvaf'] = k_all.apply(lambda row:f2(row['mut_depth_normal'], row['depth_normal']), axis=1)
		feature['lncr'] = k_all.apply(lambda row:f1(row['depth_normal'],normal_depth_thr),axis=1)
	return feature
def cal_feature():
	parser = argparse.ArgumentParser(description='script for igv')
	parser.add_argument('--file',required=True ,help = 'file name')
	#parser.add_argument('--label','-l',action='store_true',help='get label')
	parser.add_argument('-d',default=20,type=float,help='threshold of depth')
	parser.add_argument('--paired',action='store_true',help='somatic')
	parser.add_argument('-nd',default=20,type=float,help='threshold of normal depth')
	parser.add_argument('--ht',default=0.9,type=float,help='threshold of hdrr to sum together')
	parser.add_argument('-o',default='feature.csv',type=str,help='output file name of feature')
	args = parser.parse_args()
	k_all = pd.read_csv(args.file)
	feature = make_feature(k_all,paired=args.paired,hdr_thr= args.ht, depth_thr = args.d, normal_depth_thr = args.nd)
	feature.to_csv(args.o,index=False, sep='\t')
	print('feature calculation finished')
	return
if __name__=='__main__':
	cal_feature()
