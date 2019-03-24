import xgboost as xgb
from xgboost import XGBClassifier
import pandas as pd
import numpy as np
class XGBOOST():
	def __init__(self, train_data):
		#k = ['lm_num', 'mv_num', 'index', 'Chr', 'Start', 'Alt', 'Ref', 'sample','label']
		k = list(set(data_train.columns) - set(['label','index','mv_num', 'lm_num']))
		#k = list(set(train_data.columns) - set(k))
		x_train = np.array(train_data.ix[:,k])
		y_train = np.array(train_data['label'])
		x_train,x_valid,y_train,y_valid = train_test_split(x_train,y_train,test_size=0.2,random_state=0)
		self.x_train = x_train
		self.x_valid = x_valid
		self.y_train = y_train
		self.y_valid = y_valid

	def train(self, round=500,early_stopping=20):
		param = {'gamma':[0.1,1,10],'learning_rate':[0.01,0.05,0.1], 'max_depth':[3,6,9],'objective':['binary:logistic']}
		model = XGBClassifier(n_estimators=round,random_state=0,early_stopping_rounds=early_stopping)
		clf = GridSearchCV(model, param,cv=3)
		clf.fit(self.x_train, self.y_train)
		p = clf.best_params_
		max_depth = p['max_depth']
		eta = p['learning_rate']
		gamma = p['gamma']
		d_train = xgb.DMatrix(self.x_train, self.y_train)
		d_valid = xgb.DMatrix(self.x_valid, self.y_valid)
		param = {'max_depth':max_depth,'eta':eta,'objective':'binary:logistic','gamma':gamma}
		watchlist = [(d_train,'train'),(d_valid,'valid')]
		model = xgb.train(param,d_train,num_boost_round=round,evals=watchlist,early_stopping_rounds=early_stopping)
		model.save_model('./model/user.model')

