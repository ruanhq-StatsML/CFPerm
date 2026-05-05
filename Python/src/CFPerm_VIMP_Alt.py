#CFperm alternative variable importance:

from __future__ import annotations
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple, Literal
import numpy as np
import pandas as pd
import numpy as np
import pandas as pd
from typing import Dict, Optional, Sequence, Literal, Tuple, Any
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from econml.dml import CausalForestDML


'''
Utility Functions:
- _as_1d(x): return the 1d array (p, ) -> (p)
- _safe_normalize(v): return the normalized version of vector v: v/np.sum(v)
- _fit_riesz_representer(X_train, M_train): fit the outcome regression for the Riesz Representer(M~X)
'''
def _as_1d(X):
    return np.asarray(X).reshape(-1)


class CFPerm_Models:
    


'''
Variable Importances for the Causal Forest with the impurity gain
'''

def _fit

























