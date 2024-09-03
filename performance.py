import pandas as pd
import matplotlib.pyplot as plt
import model
import numpy
import os

def addToCSV(model,ordernum,stringa):
    output_path='Data/result.csv'
    performance = {
        "Settings": stringa,
        "Ordine": ordernum,
        "Area": model.ObjBoundC/1000000,
        "Gap": model.MIPGap,
        "Time": model.Runtime}
    print(performance)
    df = pd.DataFrame([performance])
    df.to_csv(path_or_buf=output_path, mode='a',index=False, header=None)

#Available: [7,8,10,17,20,22,23,24]
ordernumbers = [8]
performanceAREA = []
performanceGAP = []
performanceTIME= []
result = []
settings = []
for scenario in [1]:
    for nsmax in [1,2]:
        if scenario==3 and nsmax == 1:
            continue
        for a in ordernumbers:
            stringa = "Ordine:" + str(a) + "Scenario:" + str(scenario) + " NsMax:" + str(nsmax)
            print(stringa)
            result = model.optimise(a,scenario,nsmax)
            addToCSV(result,a,stringa)
            
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



