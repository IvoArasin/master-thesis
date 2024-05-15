import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.dates as mdates
import datetime as dt
import matplotlib.ticker as ticker
from fontTools.ttLib import TTFont


yield_data = pd.read_excel("/Users/ivoarasin/Desktop/Master/Semester Four/thesis/master_thesis_code_R/Finished DNS model files/zero_rates08til23.xlsx")
lower = 25
upper = 145
yield_data.dates = yield_data.iloc[lower:upper,0]
yield_data_dates = [dt.datetime.strptime(str(d), '%Y-%m-%d %H:%M:%S').date().strftime('%Y') for d in yield_data.dates]

yield_data = yield_data.iloc[:,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,22]]
yield_data = yield_data.iloc[lower:upper,1:]
#yield_data.columns = np.array([1/365, 7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15])
yield_data.columns = np.array([7/365, 14/365, 1/12, 2/12, 3/12, 6/12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15])
#yield_data = yield_data/100
#yield_data = yield_data.iloc[1:152,:]

range_ = range(len(yield_data))
X = np.array([yield_data.columns for i in range_])
#Y = np.array([np.array([a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a]) for a in range_])
Y = np.array([np.array([a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a]) for a in range_])
Z = np.array([ [None]*len(yield_data.columns) for a in range_])
for x in range_:
    Z[x,:] = pd.DataFrame(yield_data).iloc[x,:]

# Plot a basic wireframe.
import matplotlib.font_manager as font_manager
font_path = "/Users/ivoarasin/Downloads/latin-modern-roman/lmroman5-regular.otf"
font_manager.fontManager.addfont(font_path)
prop = font_manager.FontProperties(fname=font_path)

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = prop.get_name()
plt.rcParams["figure.figsize"] = (10,10)
plt.rcParams['font.family'] = 'Latin Modern Roman'
#print(plt.rcParams["font.sans-serif"])
#print(plt.rcParams["font.monospace"])

fig = plt.figure()
ax = fig.add_subplot(projection='3d', proj_type="ortho")
plt.gca().invert_yaxis()

ax.plot_surface(X, Y, Z, rstride=1, cstride=1, linewidth=0.5, edgecolor="black", cmap="spring")
ax.yaxis.set_ticks(range_)
ax.yaxis.set_ticklabels(yield_data_dates)
plt.locator_params(axis="y", nbins=7)
ax.view_init(elev=20, azim=-135, roll=0)
ax.set_ylabel('Date', fontsize=13, rotation=-22, labelpad=8)
ax.set_xlabel('Maturity in Years', fontsize=13)
ax.set_zlabel('Yield in %', fontsize=13)

plt.savefig("Yield_CalmYears_only19.pdf", transparent=True, pad_inches=0)
plt.show()


