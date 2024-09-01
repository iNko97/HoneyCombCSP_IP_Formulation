import pandas as pd
import matplotlib.pyplot as plt
import model
import numpy

class Order:
    def __init__(self, datasheet, number):
        self.sheet = pd.read_excel(datasheet, number)
        self.Width = self.sheet['Width']
        self.Length = self.sheet['Length']
        self.Demand = self.sheet['Demand']
ordernumbers = [7,8,10,17,20,22,23,24]
performanceGAP = []
performanceTIME= []
result = [0.5,100]
for a in range(8):
    performanceGAP.append(result[0])
    performanceTIME.append(result[1])
datapath = "Data/input_data.ods"
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
#plt.show()
fig, ax = plt.subplots()

# hide axes
fig.patch.set_visible(False)
ax.axis('off')
ax.axis('tight')
performance = [performanceGAP,performanceTIME]
performance = numpy.array(performance).transpose()
df = pd.DataFrame(data=performance,columns=['GAP','TIME'])
#df.rename({"GAP","TIME"})
ax.table(cellText=df.values, colLabels=df.columns, loc='center')

fig.tight_layout()
plt.show()
