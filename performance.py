import pandas as pd
import matplotlib.pyplot as plt
import model
import numpy

ordernumbers = [7,8,10,17,20,22,23,24]
performanceAREA = []
performanceGAP = []
performanceTIME= []
ordernum =[]
result = []
for a in ordernumbers:
    result = model.Optimize(a,1,1)
    ordernum.append(a)
    performanceGAP.append(result[1])
    performanceTIME.append(round(result[2],2))
    performanceAREA.append(int(result[0]/1000000))
"""datapath = "Data/input_data.ods"
datasheet = pd.ExcelFile(datapath, engine="odf")
orders = []
a= 1
while a <=8:
    orders.append(Order(datasheet,a))
    a+=1
widths = []
for itemtype in orders:
    widths.append(itemtype.Width)
plt.boxplot(widths)
#plt.show()

lengths = []
for itemtype in orders:
    lengths.append(itemtype.Length)
plt.boxplot(lengths)
#plt.show()
demands = []
for itemtype in orders:
    demands.append(itemtype.Demand)
plt.boxplot(demands)
#plt.show() """
fig, ax = plt.subplots()
fig.patch.set_visible(False)
ax.axis('off')
ax.axis('tight')
performance = [ordernum,performanceAREA,performanceGAP,performanceTIME]
performance = numpy.array(performance).transpose()
df = pd.DataFrame(data=performance,columns=['ORDER','AREA','GAP','TIME'])
ax.table(cellText=df.values, colLabels=df.columns, loc='center')

fig.tight_layout()
plt.show()
