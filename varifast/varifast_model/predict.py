import pandas as pd
import xgboost as xgb
import numpy as np
def predict(model_path, data):
	k = ['lm_num', 'mv_num', 'index', 'Chr', 'Start', 'Alt', 'Ref', 'sample']
	r = data.ix[:,'index']
	n = set(data.columns) - set(k)
	data = data.ix[:,n]
	model = xgb.Booster({'nthread': 4})
	model.load_model(model_path)
	p = model.predict(data)
	r['xgb_proba'] = p
	p = np.array(p)
	p[p>=0.5] = 1
	p[p<0.5] = 0
	r['xgb_label'] = p
	return r