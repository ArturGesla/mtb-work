#%% Core libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Sklearn processing
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import train_test_split

# Sklearn regression algorithms
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import KNeighborsRegressor
from sklearn.tree import DecisionTreeRegressor

# Sklearn regression model evaluation function
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import r2_score

from sklearn.svm import SVR

# Convenience functions.  This can be found on the course github
import sys
sys.path.append('data_supervised')
from functions import *

#%%
dataset=pd.read_csv("data_supervised/world_data.csv")
dataset=dataset.drop(["murder","urbanpopulation","unemployment"],axis=1)
means=dataset.iloc[:,1:].mean().to_dict()
for m in means:
    dataset[m]=dataset[m].fillna(value=means[m])

#%%
print(plt.style.available)
plt.style.use("fast")
dataset.hist(figsize=(20,20))
# scattermatrix(dataset)
plt.show()

#%%

dataset.columns
y=dataset.lifeexp
X=dataset[['happiness', 'income', 'sanitation', 'water',
       'literacy', 'inequality', 'energy', 'childmortality', 'fertility',
       'hiv', 'foodsupply', 'population']]

       # Rescale the data
scaler = MinMaxScaler(feature_range=(0,1))
rescaledX = scaler.fit_transform(X)

# Convert X back to a Pandas DataFrame, for convenience
X = pd.DataFrame(rescaledX, index=X.index, columns=X.columns)
test_size = 0.33
seed = 1
X_train, X_test, Y_train, Y_test = train_test_split(X, y, test_size=test_size, random_state=seed)

models = [LinearRegression(), KNeighborsRegressor(), SVR()]
for model in models:
    model.fit(X_train, Y_train)
    predictions = model.predict(X_train)
    print(type(model).__name__, mean_absolute_error(Y_train, predictions))

#%%
for model in models:
    predictions = model.predict(X_test)
    print(type(model).__name__, mean_absolute_error(Y_test, predictions))

#%% quiz
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
df=pd.read_csv("data_supervised/kc_house_data.csv")
df.drop(columns = ['id', 'zipcode', 'lat', 'long', 'date'], inplace = True)


df.columns
y=df.price
# X=df[['grade']]
X=df[['grade','sqft_living']]

       # Rescale the data
scaler = MinMaxScaler(feature_range=(0,1))
rescaledX = scaler.fit_transform(X)
X = pd.DataFrame(rescaledX, index=X.index, columns=X.columns)
m=LinearRegression()
m.fit(X,y)
print(m.coef_)
ypred=m.predict(X)
print(r2_score(y,ypred))

#%%
from sklearn.metrics import mean_squared_error, mean_absolute_error

from sklearn.model_selection import train_test_split

import numpy as np


def rmse(y_pred, y_true):

    return(np.sqrt( mean_squared_error(y_pred, y_true)))

predictors = ['bedrooms', 'bathrooms', 'sqft_living', 'sqft_lot','floors', 'waterfront', 'view', 'condition', 'grade', 'sqft_above', 'sqft_living15','sqft_lot15',]


# store the score for each  random_state

scores = []

# random state values

for seed in range(0,40, 1):

    # split the dataset train (90%), test (10%)

    X_train, X_test, y_train, y_test = train_test_split(

        df[predictors],

        df.price,

        test_size=0.10,

        random_state=seed)

    # fit a linear regression model

    model = LinearRegression().fit(X_train,y_train)

    # predicted values for the test set

    yhat = model.predict(X_test)

    # record the scores

    scores.append({

        'seed': seed,

        'rmse': np.sqrt(mean_squared_error(yhat, y_test)),

    })


# transform the list of dictionnaries into a dataframe for ease of use and sort the dataframe

scores = pd.DataFrame(scores).sort_values(by = 'rmse').reset_index()    
df['basement'] = 0 # by default no basement 

df.loc[ df.sqft_basement > 0 , 'basement'] = 1 # except when the surface is not zero


df['renovated'] = 0 # by default no renovation 

df.loc[ df.yr_renovated > 0 , 'renovated'] = 1 # except when the year of renovation is not zero

predictors = ['bedrooms', 'bathrooms', 'floors']

y=df.price
# X=df[['grade']]
# X=df[predictors] 
# X=df[predictors+ ['sqft_basement', 'yr_renovated']]
X=df[predictors+  ['sqft_basement', 'yr_renovated', 'basement', 'renovated']]

       # Rescale the data
scaler = MinMaxScaler(feature_range=(0,1))
rescaledX = scaler.fit_transform(X)
X = pd.DataFrame(rescaledX, index=X.index, columns=X.columns)
m=LinearRegression()
m.fit(X,y)
print(m.coef_)
ypred=m.predict(X)
print(r2_score(y,ypred))

#%%
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import  StandardScaler, MinMaxScaler
predictors = ['sqft_living', 'sqft_lot','sqft_above', 'sqft_living15','sqft_lot15']

X_train, X_test, y_train, y_test = train_test_split(

        df[predictors],

        df.price,

        test_size=0.33,

        random_state=8)
model = KNeighborsRegressor().fit(X_train,y_train)
