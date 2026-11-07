#Post_Hoc Subset Localization Software:
from LOCO_R_risk import R_risk  
from MMD_LOCO import MMD_Loco
from benchmark_config import MODEL_REGISTRY

#Inputs is the dataframes as well as the 

#Ind_fea:    List/Indices of the features to be mapped to the split
#Subset_Loc: The granularity 
#Statistic:  The list of statistics for visualize -- include the 
#B:          The number of resampling procedures 
#Split:      A list of splitting mechanisms for how to split each of the features,
#            it can either be the quantile cut or the 
#Output:     The MMD 
STATISTIC_LIST = {
    'mean': np.mean(),
    'PSI': psi_func(),
    
}
def post_hoc_localization(df, Subset_Loc, Statistic, B):
    df_subset1 = df
    
    