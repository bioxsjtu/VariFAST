import re
import pandas as pd
import os
import time
import ray
import numpy as np
def make_spare_list(f,spare_list=[],end_list=[],duplicated= True, chr_b='$'):   #get spare list for read
	r = f.readline()
	count = 0
	chr_all = set()
	while r!='':
		if r.startswith('@'):
			r = f.readline()
			#row_count += 1
			continue
			
		
		read = r.split('\t')
		chr = read[2]
		if not chr.startswith('chr'):
			chr = 'chr' + chr

		if chr_b=='$' and count==0:
			chr_bef = chr
		elif count==0:
			chr_bef = chr_b
		start = int(read[3])
		if start==0:
			r = f.readline()
			continue
		type = read[5]
		num = re.findall('\d+',type)
		nn = re.findall('[^\d]+',type)
		nn = nn[::-1]
		num = num[::-1]
		if 'N' in nn:
			print(nn)
		if len(num)==0:
			r = f.readline()
			continue
		num = [int(kk) for kk in num]
		end = start
		r = read[9]
		name = read[0]
		rnn = ''
		flag = int(read[1])
		r_start = 0     #get read position
		r_start_ins = 0   #get insert position
		lm = int(read[4])
		quality = read[10]

		######get insert
		label_insert = False
		insert = []
		while len(nn)!=0:
			ty = nn.pop()
			ny = num.pop()
			if ty=='S' and len(nn)!=0:
				r_start = ny
			if ty=='M':
				end = end + ny
				rnn = rnn + r[(r_start):(ny+r_start)]
				r_start = ny + r_start
				r_start_ins = r_start_ins + ny
			if ty=='I':
				#insert.append((r_start_ins+start-1,r[(r_start):(ny+r_start)],'I',type,start))
				insert.append((r_start_ins+start,r[(r_start):(ny+r_start)],'I',type,start))
				r_start = r_start + ny
				label_insert = True
			if ty=='D':
				insert.append((r_start_ins+start,'^'*ny,'D',type,start))
				rnn = rnn + '^'*ny
				end = end + ny
				r_start_ins = r_start_ins + ny
				label_insert = True
		
		#length = sum(num)
		#end = start + length -1
		end = end -1
		if duplicated:
			if ((flag>=1024)and(flag<2048))or(flag>=3072):
				r = f.readline()
				continue
		ff = bin(flag)
		strand = '+'
		if len(ff)>=7:
			if ff[-5]=='1':
				strand = '-'
		#rn.add(rnn)
		if label_insert:
			width = [start,end,rnn,chr,type,lm,strand,name,insert]
		else:
			width = [start,end,rnn,chr,type,lm,strand,name]
		line = 0
		if chr!=chr_bef:
			return spare_list, [end], [[width]], chr, chr_bef
		else:
			if (count==0) and (chr_b=='$'):
				spare_list.append([width])
				end_list.append(end)
			while line<len(end_list):
				if start<=end_list[line]:
					line = line +1
				else:
					break
			if line == len(end_list):
				spare_list.append([width])
				end_list.append(end)
			else:
				spare_list[line].append(width)
				end_list[line] = end
		count = count +1
		r = f.readline()
	return spare_list, [],[], '&', chr


def get_nearest_read(d_start, mut_f):   #for a position get the nearest read index
	l = len(d_start)
	end = l
	start = 0
	ind = l//2
	while True:
		if d_start[ind]>mut_f:
			if ind-1<0:
				t= 'NAN'
				break
			if d_start[ind-1]<=mut_f:
				t = ind -1 
				break
			else:
				end = ind-1
				ind = (start+ind-1)//2 
		elif d_start[ind]<mut_f:
			if ind+1>=l:
				t = l-1
				break
			if(d_start[ind+1]>=mut_f):
				t = ind
				break
			else:
				start = ind+1
				ind = (ind+1+end)//2
		else:
			t = ind
			break
	return t

def quicksearch_nb(d,mut,forward, behind, mut_nb, type): #get all nb numbers at a threshold,and get the read at mut 
	def get_del_num(r):
		l = re.findall('\^+',r)
		l = [len(j) for j in l]
		return sum(l)
	mut_f = mut-forward
	if mut_f<0:
		mut_f = 0
	if type=='D':
		mut_b = mut+behind + len(mut_nb)-1
	else:
		mut_b = mut + behind
	d_start = [i[0] for i in d]
	t = get_nearest_read(d_start,mut_f)
	t2 = get_nearest_read(d_start, mut)
	read = []
	c = 0
	read_mut = []
	if t2=='NAN':
		read_mut = 'NAN'
	else:
		if d[t2][1]>=mut:
			if len(d[t2])==9:
				read_mut.append([d[t2][0],d[t2][2],d[t2][4],d[t2][8],t2])
			else:
				read_mut.append([d[t2][0],d[t2][2],d[t2][4],t2])
	if t=='NAN':
		t = 0
		if d[t][0]>mut_b:
			read = 'NAN'
			all_nb = 0
			return all_nb, read, read_mut
		else:
			if d[t][1]<=mut_b:
				r = d[t][2][:(d[t][1]-d[t][0]+1)]
				del_num = get_del_num(r)
				all_nb = d[t][1] - d[t][0]+1 - del_num
				if len(d[t])==9:
					read.append([d[t][0],d[t][2][:(d[t][1]-d[t][0]+1)],d[t][4],d[t][8],t2])
				else:
					read.append([d[t][0],d[t][2][:(d[t][1]-d[t][0]+1)],d[t][4],t2])
			else:
				r = d[t][2][:(mut_b-d[t][0]+1)]
				del_num = get_del_num(r)
				all_nb = mut_b - d[t][0] +1 - del_num
				if len(d[t])==9:
					read.append([d[t][0],d[t][2][:(mut_b-d[t][0]+1)],d[t][4],d[t][8],t2])
				else:
					read.append([d[t][0],d[t][2][:(mut_b-d[t][0]+1)],d[t][4],t2])
	else:
		if d[t][1]<mut_f:
			all_nb=0
		elif d[t][1]<=mut_b:
			r = d[t][2][(mut_f-d[t][0]):(d[t][1]-d[t][0]+1)]
			del_num = get_del_num(r)
			all_nb = d[t][1]-mut_f+1 - del_num
			if len(d[t])==9:
				read.append([mut_f,d[t][2][(mut_f-d[t][0]):(d[t][1]-d[t][0]+1)],d[t][4],d[t][8],t2])
			else:
				read.append([mut_f,d[t][2][(mut_f-d[t][0]):(d[t][1]-d[t][0]+1)],d[t][4],t2])
		else:
			r = d[t][2][(mut_f-d[t][0]):(mut_b-d[t][0]+1)]
			del_num = get_del_num(r)
			all_nb = mut_b - mut_f+1 - del_num
			if len(d[t])==9:
				read.append([mut_f,d[t][2][(mut_f-d[t][0]):(mut_b-d[t][0]+1)],d[t][4],d[t][8],t2])
			else:
				read.append([mut_f,d[t][2][(mut_f-d[t][0]):(mut_b-d[t][0]+1)],d[t][4],t2])
	t = t+1
	if t<len(d):
		while d[t][0]<=mut_b:
			if d[t][1]<=mut_b:
				all_nb = all_nb + d[t][1] - d[t][0] +1
				if len(d[t])==9:
					read.append([d[t][0],d[t][2][:(d[t][1]-d[t][0]+1)],d[t][4],d[t][8],t2])
				else:
					read.append([d[t][0],d[t][2][:(d[t][1]-d[t][0]+1)],d[t][4],t2])
			else:
				all_nb = all_nb + mut_b - d[t][0] +1
				if len(d[t])==9:
					read.append([d[t][0],d[t][2][:(mut_b-d[t][0]+1)],d[t][4],d[t][8],t2])
				else:
					read.append([d[t][0],d[t][2][:(mut_b-d[t][0]+1)],d[t][4],t2])
			t = t+1
			if t>=len(d):
				break
	return all_nb,read,read_mut     ###read_mut the list of read conatins mut site
	
def get_mut_list_from_fa(fas,read,mut_site,mut_nb,type): ##get mut list from read list
	mut_list = []
	for r in read:
		pos = r[0]
		if len(r)==5:
			insert = r[3]
			for i in insert:
				if (i[0]<r[0])or(i[0]>(r[0]+len(r[1])-1)):
					continue
				mut_list.append((i[0],i[1],i[2]))
		for i in range(len(r[1])):
			if (fas[pos+i-1]!=r[1][i])and(r[1][i]!='^'):
				if (pos+i-1<r[0])or(pos+i-1>(r[0]+len(r[1])-1)):
					continue
				mut_list.append((pos+i,r[1][i],'M'))
	if type=='I':
		mut_site_t = mut_site + 1
	else:
		mut_site_t = mut_site
	try:
		mut_list.remove((mut_site_t,mut_nb,type))
	except:
		pass
	return mut_list
def make_fas(fa_path):        #get reference sequencing
	fa_all = dict()
	f = open(fa_path)
	r = f.readline()
	count = 0
	fas = []
	while r!='':
		if r.startswith('>'):
			r = r.strip()
			chr = r.split(' ')[0][1:]
			if not chr.startswith('chr'):
				chr = 'chr'+chr
			if count!=0:
				fas = ''.join(fas)
				fas = fas.upper()
				fa_all[chr_bef] = ray.put(fas)
				count =count +1
				chr_bef = chr
				fas = []
			else:
				chr_bef = chr
				count = count +1
			r = f.readline()
			continue
		fas.append(r.strip())
		r = f.readline()
	fas = ''.join(fas)
	fas = fas.upper()
	fa_all[chr_bef] = ray.put(fas)
	return fa_all

####mut_nb the mut nb 
def find_repeat(mut_site,mut_list,d,time,l): #if repeat and minimun repeat unit
	start = mut_site - d
	end = mut_site + d
	label = False
	label_repeat = False
	i = start
	site = set()
	while i< end:
		mut_nb = mut_list[i:i+l]
		count =1
		label = True
		while label:
			if i+l*count+1>len(mut_list):
				label = False
				break
			if mut_list[(i+l*count):(i+l*count +l)]==mut_nb:
				count+=1
				if count>=time:
					label_repeat = True
					site.add(mut_nb)
					break
			else:
				label = False
				break
		i=i+1
	return label_repeat,site

	
def make_feature_n(k_all,hdr_thr=0.9,depth_thr=20,normal_depth_thr=20):  #threshold for at least ratio to sum
	feature = pd.DataFrame()
	k_all = k_all.ix[(k_all['depth']!=0)&(k_all['mut_depth']!=0),:]
	k_all.index = range(len(k_all.index))
	feature['index'] = k_all['index']
	#feature['depth'] = k_all['depth']
	#feature['mut_depth'] = k_all['mut_depth']
	feature['vafr'] = 1 - k_all['mut_depth']*1.0/k_all['depth']   #MUT DEPTH RATIO
	def f1(x,thr):
		return 1-(min(x,thr)*1.0/thr)
	feature['lcr'] = k_all.apply(lambda row:f1(row['depth'],depth_thr),axis=1)
	feature['lm'] = k_all['LM']*1.0/k_all['mut_depth'] 
	feature['lm_num'] = k_all['LM']
	feature['mm'] = k_all['low_d_mut']*100             #*1.0/k_all['mut_depth_cancer'] ###low mut ratio
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
			if i==['']:
				m_a = 0
				m2 = 0
			else:
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
	return feature

def get_score(feature,vaf=0.16):
	feature.index = range(len(feature.index))
	score_all = []
	l1 = ['h','e','dir','sse','r','ri','mm','ao','si']
	l2 = ['mv','hdr']
	l3 = ['lm','lcr','vafr','ni','nd']
	for i in range(len(feature.index)):
		score = 0
		for j in l1:
			score += feature.ix[i,j]
		for j in l2:
			score += feature.ix[i,j]*2.0
		for j in l3:
			score += feature.ix[i,j]*3.0
		score_all.append(score)
	feature['score'] = score_all
	return feature
	
@ray.remote
def cal_depth(sampleid,ind,chr, d, mut_site,mut_nb,ref,fas_all,type,thr,forward,behind,head_num,insert_num,nb_quality,repeat_time,distance_repeat,insert_thr,insert_dis, normal_germline, seq_quality):
	#####hdr_thr the threshold near max hdrr will be saved
	#print('start:'+str(time.time()))
	#print("process:%s"%(os.getpid()))
	ttime = time.time()
	#d = ray.get(d_all[chr])
	fas = ray.get(fas_all[chr])
	ref_mv = fas[mut_site-1]
	depth = 0        #all depth
	mut_depth = 0    #mut depth
	count = 0
	mut_d_p = 0   #plus strand
	mut_d_r = 0  #reverse strand
	label_ao = 0
	site_repeat = []
	def mut_nb_in_site(mut_nb,site): ###if insert nb the same as repeat nb
		for s in site:
			s1 = '^'+s+'+'+'$'
			s2 = '^' + mut_nb+'+$'
			if re.search(s1,mut_nb) or re.search(s2,s):
				return True
		return False
###########near repeat region?
	for i in range(len(repeat_time)):
		label_r,site_r = find_repeat(mut_site,fas,distance_repeat,repeat_time[i],i+1)
		site_repeat.append(site_r)
		if i==0:
			label_repeat = label_r
			if type=='I':
				label_repeat_insert = mut_nb_in_site(mut_nb,site_r)
			else:
				label_repeat_insert = False
		else:
			label_repeat = label_repeat or label_r
			if type=='I':
				label_insert = mut_nb_in_site(mut_nb,site_r)
				label_repeat_insert = label_repeat_insert or label_insert

	if type=='I':
		if len(mut_nb)>=insert_thr:
			if re.search(mut_nb,fas[(mut_site-insert_dis):(mut_site+insert_dis)]):
				label_ao = 1
	if label_repeat_insert==True:
		label_repeat_insert = 1
	else:
		label_repeat_insert = 0
	if label_repeat==True:
		label_repeat = 1
	else:
		label_repeat = 0
	label_repeat = [label_repeat,label_repeat_insert,label_ao]
	#site_repeat = [site_dn,site_mn,site_tr]
	#depth_normal = 0
	mut_sse = 0  #mut read has same start and end
	mut_sse_d = dict()
	#mut_lm = 0   #low varinat mapping quality
	mut_sequence_lm = 0  #low sequence mapping quality
	mut_end_depth = 0   #near end
	mut_head_depth = 0   #near head
	mut_si = 0         #si feature
	mut_name_list = []
	def get_site(read,mut_site,head_num,site_num=1):
		start = read[0]
		site = mut_site - start
		if site<head_num:
			label_head = True
		else:
			label_head = False
		if mut_site+site_num-1>(read[1]-head_num):
			label_end = True
		else:
			label_end = False
		if (site+site_num)>(read[1]+1):
			return 'NAN','NAN',False,False
		else:
			return read[2][site:(site+site_num)],label_head,label_end
	mut_read_ind = [] ##the index of mut happen
	mut_read_id_ind = []
	dd_ind = 0
	num_all_nb = 0
	read_all = []
	read_mut_all = []
	#fa = fa_path+chr+'.fa'
	read_mut_ind = []  #index correspond with read
	if type=='I':
		mut_site_t = mut_site + 1
	else:
		mut_site_t = mut_site
#############cal mut depth 
	for dd in d:
		d_start = [i[0] for i in dd]
		t = get_nearest_read(d_start, mut_site_t)
		nb,read,read_mut = quicksearch_nb(dd,mut_site,forward,behind, mut_nb,type)
		#mut_site_insert_size_abnormal_all += mut_site_insert_size_abnormal
		#mut_arround_insert_size_abnormal_all += mut_arround_insert_size_abnormal
		if read_mut!='NAN' and len(read_mut)!=0:
			read_mut_ind.append([dd_ind,read_mut[0][-1]])
			read_mut_all.append(read_mut)
		if read!='NAN' and len(read)!=0:
			read_all.append(read)
		num_all_nb = num_all_nb + nb
		if t=='NAN':
			dd_ind = dd_ind +1
			continue
		if mut_site_t<=dd[t][1]:
			depth = depth +1
			site, label_head,label_end = get_site(dd[t], mut_site_t,site_num=len(mut_nb),head_num=head_num)
			if site=='NAN':
				dd_ind = dd_ind+1
				continue
			p_sse = (dd[t][0],dd[t][1])
			if (len(mut_nb)==1) and type=='M':
				if site==mut_nb:
					if dd[t][5]<=seq_quality:
						mut_sequence_lm +=1
					if p_sse in mut_sse_d:
						mut_sse_d[p_sse] +=1
					else:
						mut_sse_d[p_sse] = 1
					mut_depth = mut_depth + 1
					mut_read_ind.append(dd_ind)
					mut_read_id_ind.append((dd_ind,t))
					if label_head:
						mut_head_depth +=1
					if label_end:
						mut_end_depth +=1
					if dd[t][6]=='-':
						mut_d_r +=1
					elif dd[t][6]=='+':
						mut_d_p +=1
					mut_name_list.append(dd[t][7])
			if type=='D':
				if site==mut_nb:
					if p_sse in mut_sse_d:
						mut_sse_d[p_sse] +=1
					else:
						mut_sse_d[p_sse] = 1
					mut_depth += 1
					mut_read_ind.append(dd_ind)
					mut_read_id_ind.append((dd_ind,t))
					if dd[t][5]<=seq_quality:
						mut_sequence_lm += 1
					if label_head:
						mut_head_depth +=1
					if label_end:
						mut_end_depth+=1
					if dd[t][6]=='-':
						mut_d_r +=1
					elif dd[t][6]=='+':
						mut_d_p +=1
					mut_name_list.append(dd[t][7])
			if type=='I':
				if len(dd[t])==9:
					insert = dd[t][8]
					for i in range(len(insert)):
						if ((insert[i][0]==mut_site_t) and (mut_nb==insert[i][1])):
							if p_sse in mut_sse_d:
								mut_sse_d[p_sse] +=1
							else:
								mut_sse_d[p_sse] = 1
							mut_depth = mut_depth+1
							mut_read_ind.append(dd_ind)
							mut_read_id_ind.append((dd_ind,t))
							if dd[t][5]<=seq_quality:
								mut_sequence_lm+=1
							if label_head:
								mut_head_depth +=1
							if label_end:
								mut_end_depth +=1
							if dd[t][6]=='-':
								mut_d_r +=1
							elif dd[t][6]=='+':
								mut_d_p +=1
							mut_name_list.append(dd[t][7])
		dd_ind = dd_ind +1
	for i in mut_sse_d:
		if mut_sse_d[i]>mut_sse:
			mut_sse = mut_sse_d[i]
#######cal SIO
	mut_paired_depth = len(set(mut_name_list))
	mut_name_list.sort()
	j = 0
	while j<(len(mut_name_list)-1):
		if mut_name_list[j]==mut_name_list[j+1]:
			mut_si +=1
			j = j+2
		else:
			j =j +1
#########get mut list
	count = 0
	if len(read_all)==0:
		mut_list = []
	for read in read_all:
		if count==0:
			mut_list = get_mut_list_from_fa(fas,read,mut_site,mut_nb, type)
			count = count+1
		else:
			m = get_mut_list_from_fa(fas,read,mut_site,mut_nb,type)
			mut_list = list((set(m))|(set(mut_list)))
###########get the read which has mut as mut site
	count2 = 0
	if len(read_mut_all)==0:
		mut_list_one_read = []
	for read2 in read_mut_all:
		if count2==0:
			mut_list_one_read = get_mut_list_from_fa(fas, read2,mut_site,mut_nb,type)
			count2+=1
		else:
			m = get_mut_list_from_fa(fas, read2,mut_site,mut_nb, type)
			mut_list_one_read = list((set(m)|(set(mut_list_one_read))))
#######cal sub mut depth
	#d_s = [d[i] for i in mut_read_ind]
	d_s = d
	#read_mut_ind = [i for i in read_mut_ind if i[0] in mut_read_ind]
	#d_s2 = [[d[i[0]][i[1]]] for i in read_mut_ind]
	def cal_sub_mut(d_s,chr, ref,mut_list,mut_site,insert_num,mut_nb,type,normal_germline):
		all_mut = 0       #all mut number(except mut site)
		#l_max = 0         #use for  not same read hdr max
		l_insert_max = 0   ##near insert
		l_del_max = 0      ##near del
		low_depth_mut = 0  #low frequency mutant depth
		#cor_cite = ''      
		hdr_site = []
		hdr_cor_cite = []  #use for same read max
		mv = dict()           #mut_mv  number of mut at same site but different nb or insert
		mut_mv_type = []
		mv_max = 0
		#min_hdr = 0 #########min num fro sub mut
		#num_hdr = 0    #number of hdr 
		if type!='I':
			scope = len(mut_nb)
		else:
			scope = 1 
		for m in mut_list:
			m_sub_ind = []
			d_ind = 0
			sub_mut = 0
			mut_mv = 0
			m_near_insert_ind = []
			m_near_del_ind = []
			label_mv = False
			if m[2]=='I':
				ref2 = '-'
			else:
				ref2 = ref
			mut_n = str(ref2) + str(chr) + str(m[1]) + str(m[0])
			if normal_germline!='NULL':
				if mut_n in normal_germline:
					label_mv2 = False
				else:
					label_mv2 = True
			else:
				label_mv2 = True
			if m[2]=='I':
				if m[0]>=(mut_site+1) and m[0]<=(mut_site+scope):
					m2 = list(m)
					m2.append(ref2)
					m2.append(chr)
					mut_mv_type.append(m2)
					label_mv = True
			else:
				if m[0]>=(mut_site) and m[0]<=(mut_site+scope-1):
					m2 = list(m)
					m2.append(ref2)
					m2.append(chr)
					mut_mv_type.append(m2)
					label_mv = True
			label_mv = label_mv2&label_mv
			for d2 in d_s:
				d_start2 = [i2[0] for i2 in d2]
				t2 = get_nearest_read(d_start2,m[0])
				if t2=='NAN':
					d_ind = d_ind +1
					continue
				if m[0]<=d2[t2][1]:
					r = d2[t2]
					lm = len(m[1])
					site,_,_ = get_site(r, m[0], head_num=0,site_num=lm)
					if site=='NAN':
						d_ind = d_ind + 1
						continue
					if m[2]=='I':
						if len(r)==9:
							ins = r[8]
							for i in ins:
								if i[2]=='I':
									if (i[0]==m[0]) and (i[1]==m[1]):
										if label_mv:
											mut_mv +=1
										m_sub_ind.append((d_ind,t2))
										all_mut = all_mut+1
										sub_mut = sub_mut + 1
										if(mut_site>=(m[0]-insert_num)) and (mut_site<=(m[0]+insert_num)):
											m_near_insert_ind.append(d_ind)
					elif m[2]=='D':
						if len(r)==9:
							ins = r[8]
							for i in ins:
								if i[2]=='D':
									if (i[0]==m[0]) and (i[1]==m[1]):
										if label_mv:
											mut_mv+=1
										m_sub_ind.append((d_ind,t2))
										all_mut = all_mut +1
										sub_mut = sub_mut+1
										if (mut_site>=(m[0]-insert_num)) and (mut_site<=(m[0]+len(m[1])-1+insert_num)):
											m_near_del_ind.append(d_ind)
					else:
						if site==m[1]:
							if label_mv:
								mut_mv +=1
							m_sub_ind.append((d_ind,t2))
							all_mut = all_mut +1
							sub_mut = sub_mut+1
				d_ind = d_ind +1
			if m[2]=='I':
				if m[0]>=(mut_site+1) and m[0]<=(mut_site+scope):
					if m[0]-1 not in mv:
						mv[m[0]-1] = mut_mv
					else:
						mv[m[0]-1] += mut_mv
			else:
				if m[0]>=(mut_site) and m[0]<=(mut_site+scope-1):
					if m[0] not in mv:
						mv[m[0]] = mut_mv
					else:
						mv[m[0]] += mut_mv
			
			if sub_mut<=max(thr*depth,1):
				low_depth_mut += sub_mut
			ll = len((set(mut_read_id_ind))&(set(m_sub_ind)))
			l_sub = len(m_sub_ind)
			if (l_sub + len(mut_read_ind))==0:
				h_r =0
			else:
				h_r = (ll*2.0)/(l_sub+len(mut_read_ind))
			hdr_site.append([str(h_r),str(l_sub),str(len(mut_read_ind)),str(ll)])
			hdr_cor_cite.append(m)
			ll_del = len((set(mut_read_ind))&(set(m_near_del_ind)))
			ll_insert = len((set(mut_read_ind))&(set(m_near_insert_ind)))
			if l_insert_max<ll_insert:         #near insert
				l_insert_max = ll_insert
			if l_del_max<ll_del:               #near del
				l_del_max = ll_del
		for i in mv:
			if mv[i]>mv_max:
				mv_max = mv[i]
		#return mut_hdr,hdr_cor_cite,num_hdr,low_depth_mut,all_mut,l_insert_max,l_del_max,mut_mv,mut_mv_type
		return hdr_site,hdr_cor_cite,low_depth_mut,all_mut,l_insert_max,l_del_max,mv_max,mut_mv_type
	_,_,low_d_mut,all_m,_,_,mut_mv,mut_mv_type = cal_sub_mut(d_s,chr, ref_mv, mut_list,mut_site,insert_num=insert_num, mut_nb=mut_nb, type=type, normal_germline= normal_germline)
	mut_hdr,hdr_cor_cite,_,_,ll_insert,ll_del,_,_ = cal_sub_mut(d_s, chr, ref_mv, mut_list_one_read,mut_site,insert_num=insert_num, mut_nb=mut_nb, type=type, normal_germline= normal_germline)
	mut_D = max(mut_d_p, mut_d_r)
	if num_all_nb==0:
		low_d_mut = 1000
	else:
		low_d_mut = low_d_mut*1.0/num_all_nb
	mut_hdr = [','.join(j) for j in mut_hdr]
	mut_hdr = ':'.join(mut_hdr)
	return [ttime,ind,sampleid,chr,mut_site,ref,mut_nb,depth,mut_depth,mut_paired_depth,low_d_mut,all_m, hdr_cor_cite,mut_hdr, mut_head_depth,mut_end_depth,ll_del, ll_insert,mut_mv,mut_mv_type,mut_D,label_repeat[2], mut_sse,site_repeat,label_repeat[0],label_repeat[1], mut_si,mut_sequence_lm,num_all_nb]

def cal_all(sam_file, mut_info,fas_all,duplicated, sampleid,thread_num=2,thr = 0.1,forward=10,behind=10,head_num=1,insert_num=1,nb_quality=30,repeat_time=[10,5,5,3,3],distance_repeat=10,insert_thr=10,insert_dis=20, normal_germline='NULL',seq_quality=20):
	mut_info['index'] = range(len(mut_info.index))
	mm = mut_info.ix[mut_info['sample']==sampleid,]
	mm.index = mm['index']
	def sep_chr_result(mm, sampleid, chr, dat, fas_all, thr, forward, behind, head_num, insert_num, nb_quality, repeat_time, distance_repeat, insert_thr, insert_dis, normal_germline, seq_quality):
		k = pd.DataFrame(columns=('time','index','sample','Chr','Start','Ref','Alt','depth','mut_depth','mut_paired_depth','low_d_mut','all_mut','HDR_cor_cite','HDR','head_mut','end_mut','near_del','near_insert','MV','mv_type','D','AO','SSE','repeat_type','R','RI','SI','LM','num_all_nb'))
		pro = []
		for i in list(mm.index):
			index = int(i)
			pos = int(mm.ix[i,'Start'])
			alt = mm.ix[i,'Alt']
			ref = mm.ix[i,'Ref']
			if ref=='-':
				type = 'I'
			if alt=='-':
				alt = '^'*(len(ref))
				type = 'D'
			elif ((len(alt)==1) and (alt!='-')and(ref!='-') and (len(ref)==1)):
				type = 'M'
			#depth,mut_depth,low_d_mut,all_mut,mut_hdr, hdr_cor_cite,maa,msa,near_insert,near_del,del_abnormal_insert_depth,del_abnormal_insert_mut_depth,head_mut,end_mut,mut_mv,mut_mv_type,mut_sse,label_repeat,site_repeat,mut_D, mut_si, mut_paired_depth, mut_lm, num_all_nb  = cal_depth(d,pos,alt,fas,type,thr=thr,forward=forward,behind=behind,head_num=head_num, insert_num=insert_num, nb_quality=nb_quality, repeat_time=repeat_time, distance_repeat=distance_repeat, insert_thr= insert_thr, insert_dis=insert_dis, hdr_thr=hdr_thr)
			pro.append(cal_depth.remote(sampleid,index,chr,dat,pos,alt,ref,fas_all,type,thr,forward,behind,head_num, insert_num, nb_quality, repeat_time, distance_repeat, insert_thr, insert_dis, normal_germline, seq_quality))
			#pro.append(cal_depth(sampleid,index,chr,d_all,pos,alt,ref,fas_all,type,thr,forward,behind,head_num, insert_num, nb_quality, repeat_time, distance_repeat, insert_thr, insert_dis, hdr_thr, normal_germline))
		r = ray.get(pro)
		count = 0
		for i in r:
			k.loc[count] =i
			count = count +1
		return k
	pro = []
	f = open(sam_file,'r+')
	dat, end_list, spare_list, chr_next, chr= make_spare_list(f, spare_list=[], end_list = [],chr_b='$', duplicated=duplicated)
	count_chr = 0
	while chr_next!='&':
		mm_chr = mm.ix[mm['Chr']==chr,]
		if len(mm_chr.index)==0:
			dat, end_list, spare_list, chr_next, chr = make_spare_list(f, spare_list = spare_list, end_list = end_list, chr_b = chr_next, duplicated=duplicated)
			continue
		if chr not in fas_all.keys():
			dat, end_list, spare_list, chr_next, chr = make_spare_list(f, spare_list = spare_list, end_list = end_list, chr_b = chr_next, duplicated=duplicated)
			continue
		if count_chr==0:
			k_all = sep_chr_result(mm_chr, sampleid, chr, dat, fas_all, thr, forward, behind, head_num, insert_num, nb_quality, repeat_time, distance_repeat, insert_thr, insert_dis, normal_germline, seq_quality)
			count_chr += 1
		else:
			k = sep_chr_result(mm_chr, sampleid, chr, dat, fas_all, thr, forward, behind, head_num, insert_num, nb_quality, repeat_time, distance_repeat, insert_thr, insert_dis, normal_germline, seq_quality)
			k_all = pd.concat([k_all,k])
			count_chr += 1
		dat, end_list, spare_list, chr_next, chr= make_spare_list(f, spare_list = spare_list, end_list = end_list, chr_b = chr_next, duplicated = duplicated)
	mm_chr = mm.ix[mm['Chr']==chr,]
	if (len(mm_chr.index)!=0) and (chr in fas_all.keys()):
		k = sep_chr_result(mm_chr, sampleid, chr, dat, fas_all, thr, forward, behind, head_num, insert_num, nb_quality, repeat_time, distance_repeat, insert_thr, insert_dis, normal_germline, seq_quality)
		k_all = pd.concat([k_all,k])
	f.close()
	print(sampleid+':end:'+str(time.time()))
	#print('all num' + str(len(k_all.index)))
	return k_all

def multi_cal(param,d,vaf,thread_num,fas_all,thr,forward,behind,duplicated,nb_quality,head_num,insert_num, repeat_time, distance_repeat, insert_thr, insert_dis,score_thr,seq_quality, ht):
	file = param[0]
	sam_file_cancer=param[1]
	sam_file_normal = param[2]
	sampleid = param[3]
	print(sampleid+':start:'+str(time.time()))
	mut_info = pd.read_csv(file)
	if sam_file_normal!='NULL':
		print(sampleid+':dataframe made:'+str(time.time()))
		m_normal = cal_all(sam_file_normal ,mut_info,fas_all, duplicated, sampleid,thread_num=thread_num, nb_quality=nb_quality,thr=thr,forward=forward,behind=behind,head_num=head_num, insert_num=insert_num, repeat_time = repeat_time,distance_repeat=distance_repeat,insert_thr=insert_thr, insert_dis=insert_dis,seq_quality=seq_quality)
		mv = list(m_normal['mv_type'])
		mv_d = pd.DataFrame()
		ref = []
		alt = []
		chr = []
		pos = []
		for i in mv:
			for j in i:
				if j[2]=='D':
					alt.append('-')
				else:
					alt.append(j[1])
				pos.append(j[0])
				chr.append(j[4])
				ref.append(j[3])
		mv_d['Chr'] = chr
		mv_d['Start'] = pos
		mv_d['Alt'] = alt
		mv_d['Ref'] = ref
		mv_d['sample'] = sampleid
		#print(mv_d)
		normal_ger = set()
		if len(mv_d.index)!=0:
			k_normal = cal_all(sam_file_normal, mv_d,fas_all,duplicated, sampleid,thread_num=thread_num, nb_quality=nb_quality,thr=thr,forward=forward,behind=behind,head_num=head_num, insert_num=insert_num, repeat_time = repeat_time,distance_repeat=distance_repeat,insert_thr=insert_thr, insert_dis=insert_dis,seq_quality=seq_quality)
			#k_normal.to_csv('test.csv')
			feature_n = make_feature_n(k_normal, depth_thr=d, hdr_thr=ht)
			tag_n = get_score(feature_n,vaf = vaf)
			#tag_n.to_csv('test_tag.csv')
			tag_n = tag_n.ix[tag_n['score']<=score_thr,:]
			k_normal.index = k_normal['index']
			kk = k_normal.ix[list(tag_n['index']),:]
			for i in kk.index:
				ss = str(kk.ix[i,'Ref']) + str(kk.ix[i,'Chr']) + str(kk.ix[i,'Alt']) + str(kk.ix[i,'Start'])
				normal_ger.add(ss)
		#d_all = make_spare_list(sam_file_cancer,duplicated=duplicated)
		k_cancer = cal_all(sam_file_cancer, mut_info,fas_all, duplicated, sampleid,thread_num=thread_num,nb_quality=nb_quality,thr=thr,forward=forward,behind=behind,head_num=head_num, insert_num=insert_num, repeat_time = repeat_time,distance_repeat=distance_repeat,insert_thr=insert_thr, insert_dis=insert_dis,normal_germline = normal_ger, seq_quality=seq_quality)
		k_normal = m_normal.ix[:,['index','depth','mut_depth']]
		k_normal.columns = ['index','depth_normal','mut_depth_normal']
		k = pd.merge(k_normal,k_cancer,on='index')
		k.to_csv('./result/' + sampleid+'.csv',index=False)
	else:
		#d_all = make_spare_list(sam_file_cancer,duplicated=duplicated)
		print(sampleid+':dataframe made:'+str(time.time()))
		k_cancer = cal_all(sam_file_cancer, mut_info,fas_all, duplicated, sampleid,thread_num=thread_num, nb_quality=nb_quality,thr=thr,forward=forward,behind=behind,head_num=head_num, insert_num=insert_num, repeat_time = repeat_time,distance_repeat=distance_repeat,insert_thr=insert_thr, insert_dis=insert_dis, seq_quality=seq_quality)
		k_cancer.to_csv('./result/'+sampleid+'.csv', index=False)
def get_all_sam(file, sam_file_name,fa_path,d=20, vaf = 0.2, thr=0.1,forward=10,behind=10,duplicated=True,nb_quality=20,thread_num=1,head_num=1,insert_num=1,repeat_time = [10,5,5,3,3],distance_repeat=20,insert_thr=10,insert_dis=20, score_thr=4, seq_quality=20, ht=0.9):        ###r1t at least r1t number nb as repeat region   ####insert_thr at least insert_thr number of nb regard as ao     
	ray.init(num_cpus=thread_num)
	print('start')
	sam_file_list = pd.read_csv(sam_file_name)
	sam_case = list(sam_file_list['case'])
	size = []
	for i in sam_case:
		size.append(os.path.getsize(i))
	sam_file_list['size'] = size
	sam_file_list.sort_values(by='size',inplace=True,ascending=False)
	sam_file_list.index = range(len(sam_file_list.index))
	print(sam_file_list)
	count = 0
	fas_all = make_fas(fa_path)
	try:
		os.system('mkdir result')
	except:
		pass
	for i in range(len(sam_file_list.index)):
		sampleid = sam_file_list.ix[i,'sample']
		sam_file_cancer = sam_file_list.ix[i,'case']
		if len(sam_file_list.columns)==3:
			sam_file_normal='NULL'
		else:
			sam_file_normal = sam_file_list.ix[i,'control']
		param = [file,sam_file_cancer,sam_file_normal,sampleid]
		multi_cal(param,d, vaf,thread_num,fas_all,thr,forward,behind,duplicated,nb_quality,head_num, insert_num,repeat_time,distance_repeat,insert_thr,insert_dis, score_thr, seq_quality, ht)
	ray.shutdown()
	f = os.listdir('./result/')
	for i in f:
		if count==0:
			k_all = pd.read_csv('./result/' + i)
			count +=1
		else:
			k = pd.read_csv('./result/' + i)
			k_all = pd.concat([k_all,k],ignore_index=True)
	print('finished')
	print(time.time())
	return k_all

if __name__=='__main__':
	k_all = get_all_sam('mut_info2.csv','sam_file.csv','/share/home/hangzhang/igv/fa/',thr=0.1,duplicated=True,forward=10,behind=10,thread_num=1,insert_percentile=0.005)
	k_all.to_csv('all_no_label_test.csv')
	gd = pd.read_csv('mut_label.csv')
	k_all = get_label(k_all,gd)
	k_all.to_csv('all_mut_label.csv')
	feature = make_feature(k_all,paired=True)
	feature.to_csv('feature.csv')
	tag = get_tag(feature,vaf=0.16,depth=20,paired=True)
	tag.to_csv('tag.csv')
	#filter(k_all,paired=True)
