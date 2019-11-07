import numpy as np
pred_values = ["TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP", "TP" , "FN", "FN", "FN", "FN", "FN", "FN", "FN"]
len(pred_values)
27
np.random.seed(116)
boot_sensi = []

for _ in range(10000):
    bootsample = np.random.choice(pred_values, 20, replace=True)
    se = list(bootsample).count("TP") * 1.0/20
    boot_sensi.append(se)

np.percentile(boot_sensi, 2.5)
# 0.55
np.percentile(boot_sensi, 97.5)
# 0.9


## ref https://cloud.tencent.com/developer/ask/213904
## ref https://www.jianshu.com/p/f7dc2da8dd13

import statsmodels.api as sm
from statsmodels import stats
stats.proportion.proportion_confint(27, 32, alpha=0.05, method="beta")
## (0.6721212343179582, 0.9472494355978659)
stats.proportion.proportion_confint(27, 32, alpha=0.05, method="binom_test")
## (0.6742326654458024, 0.9363472682486428)

## ref https://stackoverflow.com/questions/51794473/calculating-confidence-interval-for-a-proportion-in-one-sample

import pandas as pd
import scipy
import scipy.stats as st
data = pd.DataFrame({
     "exp1":[34, 41, 39]
    ,"exp2":[45, 51, 52]
    ,"exp3":[29, 31, 35]}).T

data.loc[:,"row_mean"] = data.mean(axis=1)
data.loc[:,"row_std"] = data.std(axis=1)
data
#        0   1   2   row_mean   row_std
# exp1  34  41  39  38.000000  2.943920
# exp2  45  51  52  49.333333  3.091206
# exp3  29  31  35  31.666667  2.494438
tscore = st.t.ppf(1-0.025, data.shape[0]-1)
mean_of_means = data.row_mean.mean()
std_of_means = data.row_mean.std()

lower_bound = mean_of_means - (tscore*std_of_means/(data.shape[0]**0.5))
upper_bound = mean_of_means + (tscore*std_of_means/(data.shape[0]**0.5))
lower_bound
## 17.432439139464606
upper_bound
## 61.90089419386874
data.shape[0]
## 3
