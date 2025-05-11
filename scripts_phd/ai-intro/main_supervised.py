#%%
# Import Python libraries for data manipuation and visualization
import pandas as pd
import numpy as np
import matplotlib.pyplot as pyplot

# Import the Python machine learning libraries we need
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score

# Import some convenience functions.  This can be found on the course github
import sys
sys.path.append('data_supervised')
from functions import *
#%%

df=pd.read_csv('data_supervised/world_data_really_tiny.csv')
df.head()
df.shape
df.describe()
# %%
