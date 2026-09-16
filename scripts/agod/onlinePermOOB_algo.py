import os
import argparse
import numpy as np
import pandas as pd
import shap
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import train_test_split
#from model_registry_class import ModelRegistry
from river import drift
from bayesian_changepoint_detection.online_changepoint_detection import (
    online_changepoint_detection,
    constant_hazard,
    StudentT)
from bocpd.bayesian_models import online_changepoint_detection as bocpd_online
from bocpd.online_likelihoods import StudentT as bocpd_studentT
from frouros.detectors.concept_drift import DDM, STEPD, HDDMA, ECDDWT
from river.anomaly import HalfSpaceTrees
from bocpd.bayesian_models import online_changepoint_detection as bocpd_online
from bocpd.online_likelihoods import StudentT as BocpdStudentT
from bocpd.hazard_functions import constant_hazard as bocpd_constant_hazard
from frouros.detectors.concept_drift import DDM, STEPD, HDDMA, ECDDWT
from skmultiflow.drift_detection import ADWIN as SKM_ADWIN
from skmultiflow.drift_detection import KSWIN as SKM_KSWIN
#from fitness_anomaly_detection import FitnessGaussian
from bocpd.bayesian_models import online_changepoint_detection as bocpd_online
from frouros.detectors.concept_drift import DDM, STEPD, HDDMA, ECDDWT
import numpy as np
import scipy
import rrcf
from model_registry_class import *
model_factory = ModelRegistry(
ntree = 150, ridge_alpha = 0.25,
nthread = 1, maxit = 500,
max_depth = 5, gamma = 0.25,
eta = 0.15, mlp_hidden_size = 4,
mlp_decay = 1e-4, mlp_max_iter = 500
)
MODEL_REGISTRY = model_factory.as_r_style_dict()
def onlinePermOOB_wholedf(df, 
    model_registry = MODEL_REGISTRY,
    model_m = 'rf_regression',
    ref_batch_size = 2000, batch_size = 50,
    burnin = 5, seed = 2026, alpha = 0.05,
    frouros_burnin = 500, frouros_quantile = 90):
    results = {}
    model_registry = MODEL_REGISTRY
    model_m = 'rf_regressor'
    ref_batch_size = 2000
    batch_size = 1
    burin = 0
    p = df.shape[1]
    df_ref = df[:ref_batch_size, :]
    df_ref_X = np.asarray(df_ref[:ref_batch_size,:(p-2)]).astype(float)
    df_ref_Y = np.asarray(df_ref[:ref_batch_size, (p-1)]).astype(float)
    model = MODEL_REGISTRY[model_m]
    model_fit = model['fit'](
        df_ref_X, df_ref_Y, seed = seed
    )
    y_hat = model['predict'](
        model_fit, df_ref_X
    )
    mse_ref = np.mean((df_ref_Y - y_hat) ** 2)
    df_trail = df[(ref_batch_size+1):, :]
    MSE_list = []
    #Conduct the following inference procedure to get the next predictive value:
    for i in range(0, len(df_trail), batch_size):
        if i > (len(df_trail) - batch_size + 1):
            df_batch = df_trail[i:,:]
        else:
            df_batch = df_trail[i:(i + batch_size), :]
        df_batch_X = np.asarray(df_batch[:, :(p-2)]).astype(float)
        df_batch_Y = np.asarray(df_batch[:, p-1]).astype(float)
        ###Then fit on the new batch, record the new MSE:
        y_hat_batch = model['predict'](
            model_fit, df_batch_X
        )
        MSE_list.append(np.mean((df_batch_Y - y_hat_batch)**2))
    ###Update the empirical p-value:
    pval_list = empirical_pval(MSE_list, burnin = burnin)
    rej_addis = method_addis(pval_list)
    rej_saffron = method_saffron(pval_list)
    rej_fixed = method_fix_alpha(pval_list, alpha = alpha)
    #extract the first rejection point:
    scale = max(np.std(MSE_list[:100]), 1e-10) if len(MSE_list) > 100 else 1e-6
    ph_update = PageHinkley(delta=0.005 * scale, threshold=5.0 * scale, min_instances=30)
    detetection_result = [ph_update.update(x) for x in MSE_list]
    results['PH_1'] = first_k_consecutive_rej(detetection_result, 1)
    results['PH_2'] = first_k_consecutive_rej(detetection_result, 2)
    results['PH_3'] = first_k_consecutive_rej(detetection_result, 3)
    #Bayesian Change-Point Detection(from the standard implementation)
    mu0 = float(np.mean(MSE_list[:min(50, len(MSE_list))]))
    obs = StudentT(alpha = 1.0, beta = 1.0, kappa = 1.0, mu = mu0)
    hazard = lambda r: constant_hazard(100, r)
    _, maxes = online_changepoint_detection(MSE_list, hazard, obs)
    det = np.zeros(len(MSE_list), dtype = bool)
    for i in range(50, len(maxes)):
        if maxes[i] < 3 or (i > 0 and maxes[i] < maxes[i - 1] - 5):
            det[i - 1] = True
    results['BOCPD_SUM'] = np.sum(det)
    results['BOCPD_bcpd_2'] = first_k_consecutive_rej(det, 2)
    results['BOCPD_bcpd_1'] = first_k_consecutive_rej(det, 1)
    results['BOCPD_bcpd_3'] = first_k_consecutive_rej(det, 3)
    model = EWMA(r = 0.1, burnin = 25)
    model.process(MSE_list)
    change_point_list = model.changepoints
    results['EWMA_SUM'] = len(change_point_list)
    results['EWMA_3'] = first_k_consecutive_rej_ind(change_point_list, 3)  
    results['EWMA_2'] = first_k_consecutive_rej_ind(change_point_list, 2)
    results['EWMA_1'] = first_k_consecutive_rej_ind(change_point_list, 1)
    #TwoSample:
    model = TwoSample(statistic = 'Lepage', threshold = 2.5)
    model.process(MSE_list)
    change_point_list = model.changepoints
    results['TS_SUM'] = np.sum(det)
    results['TS_3'] = first_k_consecutive_rej_ind(change_point_list, 3)
    results['TS_2'] = first_k_consecutive_rej_ind(change_point_list, 2)
    results['TS_1'] = first_k_consecutive_rej_ind(change_point_list, 1)
    #ADWIN:
    adwin = drift.ADWIN(delta = 0.02, clock = 1, grace_period = 30)
    det = []
    for MSE in MSE_list:
        adwin.update(float(MSE))
        det.append(adwin.drift_detected)
    results['ADWIN_SUM'] = np.sum(det)
    results['ADWIN_3'] = first_k_consecutive_rej(det, 3)
    results['ADWIN_2'] = first_k_consecutive_rej(det, 2)
    results['ADWIN_1'] = first_k_consecutive_rej(det, 1)
    #HDMMA:
    detector = SKM_ADWIN(delta = 0.05)
    det = []
    for MSE in MSE_list:
        detector.add_element(float(MSE))
        det.append(detector.detected_change())
    results['HDMMA_SUM'] = np.sum(det)
    results["HDMMA_3"] = first_k_consecutive_rej(det, 3)
    results["HDMMA_2"] = first_k_consecutive_rej(det, 2)
    results["HDMMA_1"] = first_k_consecutive_rej(det, 1)
    #ECDDWT: 
    detector = SKM_KSWIN(window_size = 50, stat_size = 20)
    det = []
    for i, MSE in enumerate(MSE_list):
        detector.add_element(float(MSE))
        det.append(detector.detected_change())
        print(i)
    results['ECDDWT_SUM'] = np.sum(det)
    results["ECDDWT_3"] = first_k_consecutive_rej(det, 3)
    results["ECDDWT_2"] = first_k_consecutive_rej(det, 2)
    results["ECDDWT_1"] = first_k_consecutive_rej(det, 1) 
    #HalfSpace Tree based methods, adapted from the original implementation
    init_len = min(30, len(MSE_list))
    mi, ma = float(np.min(MSE_list[:init_len])), float(np.max(MSE_list[:init_len]))
    if ma - mi < 1e-12:
        ma = mi + 1.0
    scaled = [np.clip((t - mi) / (ma - mi), 0.0, 1.0) for t in MSE_list]
    hst = HalfSpaceTrees(n_trees=150, height=4, window_size=200)
    scores = []
    for i, v in enumerate(scaled):
        x = {"r": float(v)}
        scores.append(hst.score_one(x))
        hst.learn_one(x)
        print(i)
    threshold = np.percentile(scores[:min(200, len(scores))], 95) if len(scores) >= 20 else (np.max(scores) if scores else 0.0)
    det = [s > threshold for s in scores]
    results['HalfSpaceTrees_SUM'] = np.sum(det)
    results["HalfSpaceTrees_1"] = first_k_consecutive_rej(det, 1)
    results["HalfSpaceTrees_2"] = first_k_consecutive_rej(det, 2)
    results["HalfSpaceTrees_3"] = first_k_consecutive_rej(det, 3)
    #Martingale based methods:
    pvals = empirical_pval(MSE_list)
    rej_martingale = method_martingale(pvals, alpha = alpha)
    rej_martingale2 = method_martingale2(pvals, alpha = alpha)
    results['martingale_SUM'] = np.sum(rej_martingale)
    results['martingale2_SUM'] = np.sum(rej_martingale2)    
    results['martingale_1'] = first_k_consecutive_rej(rej_martingale, 1)
    results['martingale_2'] = first_k_consecutive_rej(rej_martingale, 2)
    results['martingale_3'] = first_k_consecutive_rej(rej_martingale, 3)
    results['martingale2_1'] = first_k_consecutive_rej(rej_martingale2, 1)
    results['martingale2_2'] = first_k_consecutive_rej(rej_martingale2, 2)
    results['martingale2_3'] = first_k_consecutive_rej(rej_martingale2, 3)
    rej_fix = (np.asarray(pvals) <= alpha)
    results['fix_SUM'] = np.sum(rej_fix)  
    results['fix_1'] = first_k_consecutive_rej(rej_fix, 1)
    results['fix_2'] = first_k_consecutive_rej(rej_fix, 2)
    results['fix_3'] = first_k_consecutive_rej(rej_fix, 3)
    rej_addis = method_addis(pvals, alpha = alpha)
    rej_saffron = method_saffron(pvals, alpha = alpha)
    results['addis_SUM'] = np.sum(rej_fix)  
    results['addis_1'] = first_k_consecutive_rej(rej_addis, 1)
    results['addis_2'] = first_k_consecutive_rej(rej_addis, 2)
    results['addis_3'] = first_k_consecutive_rej(rej_addis, 3)
    results['saffron_SUM'] = np.sum(rej_fix)  
    results['saffron_1'] = first_k_consecutive_rej(rej_saffron, 1)
    results['saffron_2'] = first_k_consecutive_rej(rej_saffron, 2)
    results['saffron_3'] = first_k_consecutive_rej(rej_saffron, 3)
    #rrcf building trees incrementally here:
    #Streaming Anomaly Detection Here:
    num_trees = 30
    shingle_size = 5
    tree_size = 128
    for _ in range(num_trees):
        forest.append(tree)
    points = rrcf.shingle(df[:])
    for idx, row in df.iterrows():
        row.update({})
    #Extract the feature importance of the online anomaly detection procedure.
    #frouros packages, assume stable?
    #What is the burn-in period?
    threshold = np.percentile(MSE_list[:frouros_burnin], frouros_quantile)
    #Require the Error Stream:
    MSE_stream = (MSE_list > threshold).astype(int)
    #Then update them:
    for detector_name, Detector in [('DDM', DDM), ('STEPD', STEPD),
        ('HDDMA', HDDMA), ('ECDDWT', ECDDWT)]:
        detector = Detector()
        detection_result = []
        for mse in MSE_stream:
            detector.update(value = mse)
            detection_result.append(detector.drift)
        results[f"{detector_name}_SUM"] = np.sum(detection_result)
        results[f"{detector_name}_1"] = first_k_consecutive_rej(detection_result, 1)
        results[f"{detector_name}_2"] = first_k_consecutive_rej(detection_result, 2)
        results[f"{detector_name}_3"] = first_k_consecutive_rej(detection_result, 3)        
    return results

