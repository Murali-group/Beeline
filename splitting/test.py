# import importlib
# import pandas as pd
# graph_split = importlib.import_module("graph-split")

import pandas as pd
from split_script import *
df = pd.DataFrame({'source':[1,2,3,3,4,4,3,5], 'target':[2,3, 1, 1,2,4,1,1]})
df.sort_values(by=['target'], inplace=True)

print(df)
train, test = split_train_test(df, 'edge', 0.4, seed=None)

print('train: \n', df.iloc[train])
print('test: \n', df.iloc[test])
new_samples = generate_negative_samples(df)
print(new_samples)