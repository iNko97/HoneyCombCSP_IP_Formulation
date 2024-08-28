import pandas as pd
import matplotlib.pyplot as plt

class Order:
    def __init__(self, datasheet, number):
        self.sheet = pd.read_excel(datasheet, number)
        self.Width = self.sheet['Width']
        self.Length = self.sheet['Length']
        self.Demand = self.sheet['Demand']

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
plt.show()

lengths = []
for itemtype in orders:
    lengths.append(itemtype.Length)
plt.boxplot(lengths)
plt.show()
demands = []
for itemtype in orders:
    demands.append(itemtype.Demand)
plt.boxplot(demands)
plt.show()