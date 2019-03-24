from varifast_metric import *
import argparse
def main():
	parser = argparse.ArgumentParser(description='script for igv')
	parser.add_argument('-p',default=2,type=int,help='number of additional threads to use')
	parser.add_argument('--vcf',required=True,type=str,help='file name of vcf file,should be csv')
	parser.add_argument('--sam',required=True,type=str,help='file name of sam file,should be csv')
	parser.add_argument('--fasta',required=True ,type=str,help = 'fasta file name')
	parser.add_argument('--vaf',required=True,type=float,help ='at least vaf of variant')
	parser.add_argument('-d',required=True,type=int,help ='at least sequencing depth')
	parser.add_argument('-nd',default=1,type=int,help ='at least sequencing depth for normal')
	parser.add_argument('--duplicated',action='store_true',help ='if remove duplicated read')
	parser.add_argument('--score',default=4,type=float,help='threshold of score for normal germline')
	parser.add_argument('--paired',action='store_true',help='somatic')
	parser.add_argument('--quality',default=10,type=float,help='sequence low mapping quality')
	
	####model paramter
	parser.add_argument('--model', action='store_true', help='use xgboost model')
	parser.add_argument('--pmodel', default='./model/A_trio_germline.model',type=str, help='path of xgboost model')
	parser.add_argument('--train', action='store_true', help ="use user's data to train an xgb model")
	parser.add_argument('--fvar',default='NULL',help ='if train, the file name of golden variants')
	
	####tools paramter
	parser.add_argument('--lt',default=0.05,type=float,help='low frequnt threshold')
	parser.add_argument('--f',default=20,type=int,help='forward number of nb')
	parser.add_argument('--repeat',default=[10,5,5,3,3],nargs='+',help='as least repeat time')
	parser.add_argument('--rd',default=20,type=int,help='forward number of nb to determine repeat')
	parser.add_argument('--hn',default=10,type=int,help='num of nb to head or end')
	parser.add_argument('--ins',default=1,type=int,help='num of nb to insert or del')
	parser.add_argument('--ao',default=10,type=int,help='num of insert nb for ao')
	parser.add_argument('--aod',default=20,type=int,help='distance for ao')
	
	###metric paramter
	parser.add_argument('--ht',default=0.9,type=float,help='threshold of hdrr to sum together')
	
	######tag paramter
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
	
	parser.add_argument('-o',default='Mid-Result.csv',type=str,help='output file name of result')
	parser.add_argument('--so',default='Score.csv',type=str,help='output file name of score')
	parser.add_argument('--mo',default='Metric.csv',type=str,help='output file name of metric')
	args = parser.parse_args()
	r = args.repeat
	r = [int(i) for i in r]
	k_all = get_all_sam(args.vcf,args.sam,args.fasta,vaf = args.vaf, d = args.d,thr=args.lt,duplicated=args.duplicated,forward=args.f,behind=args.f,thread_num=args.p,head_num=args.hn, insert_num=args.ins,repeat_time=r,distance_repeat=args.rd,insert_thr=args.ao,insert_dis=args.aod,score_thr = args.score, seq_quality=args.quality)
	k_all.to_csv(args.o,index=False)
	print('basic finished')
	feature = make_feature(k_all,paired=args.paired,hdr_thr= args.ht,depth_thr=args.d, normal_depth_thr = args.nd)
	feature.to_csv(args.mo,index=False, sep='\t')
	print('feature calculation')
	tag = get_tag(feature,vaf=args.vaf,paired=args.paired,threshold_h=args.h,threshold_e=args.e,threshold_ni=args.ni,threshold_nd=args.nd,threshold_dir=args.dir,threshold_sse=args.sse,threshold_si=args.si,threshold_hdr=args.hdr,threshold_mm=args.mm, threshold_lm=args.lm, threshold_mv=args.mv, rate_tn = args.tn)
	if args.train:
		from varifast_model import *
		assert args.fvar!='NULL'
		feature_labeled = get_label(feature, args.fvar)
		model = XGBOOST(feature_labeled)
		model.train()
	if args.model:
		from varifast_model import *
		r = predict(args.pmodel, feature)
		tag = pd.merge(tag, r, on='index')
	tag.to_csv(args.so,index=False, sep='\t')
	print('tag get and all finished')
	return
if __name__=='__main__':
	main()
