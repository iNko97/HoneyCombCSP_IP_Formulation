import pandas as pd
import matplotlib.pyplot as plt
import model
import numpy

#Available: [7,8,10,17,20,22,23,24]
ordernumbers = [10]
performanceAREA = []
performanceGAP = []
performanceTIME= []
ordernum =[]
result = []
settings = []
for scenario in [1,2,3]:
    for nsmax in [1,2]:
        if scenario==3 and nsmax == 1:
            continue
        for a in ordernumbers:
            stringa = "Ordine:" + str(a) + "Scenario:" + str(scenario) + " NsMax:" + str(nsmax)
            print(stringa)
            settings.append(stringa)
            result = model.optimise(a,scenario,nsmax)
            ordernum.append(a)
            performanceGAP.append(round(result[1]*100,2))
            performanceTIME.append(round(result[2],2))
            performanceAREA.append(round(result[0]/1000000))
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
performance = [settings,ordernum,performanceAREA,performanceGAP,performanceTIME]
performance = numpy.array(performance).transpose()
df = pd.DataFrame(data=performance,columns=['SETTINGS','ORDER','AREA','GAP','TIME'])
table =ax.table(cellText=df.values, colLabels=df.columns, loc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)


fig.tight_layout()
plt.show()
