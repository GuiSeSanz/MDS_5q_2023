
from simiclasso.clus_regression import simicLASSO_op
from simiclasso.weighted_AUC_mat import main_fn
import time
import os
import subprocess as sp
# 

DATA_PATH = '/root/SimiC/Data_MDS'
RES_PATH = '/root/SimiC/Data_MDS'

# lambda1 = 0.001, lamda2 = 0.01, done
# ----> adjusetd R2 = 0.7972
lambda1 = 0.001
lambda2 = 0.01


similarity = True
k_cluster = None
num_TFs = -1
num_target_genes = -1
max_rcd_iter = 100000
df_with_label = False
percent_of_target = 1
cross_val = False

SAMPLE = 'ThreeWayComparison_1000'
# SAMPLE = 'Non5qVsElder_1000'
# SAMPLE = '5qVsElder_1000'
p2assignment = os.path.join(DATA_PATH, 'ThreeWayComparison' +'.clustAssign.txt')
# We need the following files:
p2df = os.path.join(DATA_PATH, SAMPLE+'.DF.pickle')
p2tf = os.path.join(DATA_PATH, SAMPLE+'.TF.pickle')
assert os.path.exists(p2df), 'File not found: ' + p2df
assert os.path.getsize(p2df) > 0, 'File is empty: ' + p2df
assert os.path.exists(p2assignment), 'File not found: ' + p2assignment
assert os.path.getsize(p2assignment) > 0, 'File is empty: ' + p2assignment
assert os.path.exists(p2tf), 'File not found: ' + p2tf
assert os.path.getsize(p2tf) > 0, 'File is empty: ' + p2tf

lambda1= 0.01
lambda2_list = [0.01]

if cross_val:
	print('CROSS VALIDATION')
	lambda1='CrossVal'
	lambda2='CrossVal'
	p2saved_file = os.path.join(RES_PATH, SAMPLE + '_L1'+str(lambda1)+'_L2'+str(lambda2)+'_Ws.pickle') #Weigths
	simicLASSO_op(p2df, p2assignment, similarity, p2tf, p2saved_file,  k_cluster, num_TFs, num_target_genes, 
				max_rcd_iter = max_rcd_iter, df_with_label = df_with_label, cross_val=cross_val)
else:
	for lambda2 in lambda2_list:
		print(lambda2)
		p2saved_file = os.path.join(RES_PATH, SAMPLE + '_L1'+str(lambda1)+'_L2'+str(lambda2)+'_Ws.pickle') #Weigths
		ts_simic = time.time()

		# # SET LAMBDA1 AND LAMBDA2
		print('\n\nSET LAMBDA1 AND LAMBDA2\n\n')
		print('\nCalculating weights\n')
		simicLASSO_op(p2df, p2assignment, similarity, p2tf, p2saved_file,  k_cluster, num_TFs, num_target_genes, 
			max_rcd_iter = max_rcd_iter, df_with_label = df_with_label,
			lambda1=lambda1, lambda2 = lambda2)
		te_simic = time.time()
		t_simic = te_simic - ts_simic 
		p2saved_file = os.path.join(RES_PATH, SAMPLE + '_L1'+str(lambda1)+'_L2' +str(lambda2)+'_Ws.pickle') #Weigths

		# filter the weigths
		sp.check_call('Rscript /root/SimiC/Filter_weigths.R ' + p2saved_file, shell=True) 
		p2saved_file = os.path.join(RES_PATH, SAMPLE + '_L1'+str(lambda1)+'_L2'+str(lambda2)+'_Ws_filtered_BIC.pickle') #Weigths

		time_pass = lambda x: '{}h{}min'.format(x // 3600, x// 60 - x//3600)
		print('simic uses {}'.format(time_pass(t_simic)))
		print(p2saved_file)
		p2AUC = os.path.join(RES_PATH, SAMPLE + '_L1'+str(lambda1)+'_L2'+str(lambda2)+'_AUCs.pickle') #Auc
		ts_auc = time.time()
		print('\nCalculating AUCs\n')
		main_fn(p2df, p2saved_file, p2AUC, percent_of_target = percent_of_target)

		te_auc = time.time()
		t_auc = te_auc - ts_auc
		t_total = te_auc - ts_simic
		print('simic uses {}'.format(time_pass(t_simic)))
		print('auc uses {}'.format(time_pass(t_auc)))
		print('total uses {}'.format(time_pass(t_total)))

print('Yay!')


