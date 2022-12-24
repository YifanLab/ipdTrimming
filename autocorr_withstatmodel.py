import pandas as pd
from matplotlib import pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf

data = pd.read_table('./ecog2andamt1_auto.xls')
print(data)
ecog2 = data[['Distance', 'EcoG2']].set_index(['Distance'])
amt1 = data[['Distance', 'AMT1']].set_index(['Distance'])
print(data)
plot_acf(ecog2, lags=999)
#plot_acf(amt1, lags=500)
plt.ylim(-0.1,0.1)
plt.show()