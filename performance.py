import pandas as pd
import matplotlib.pyplot as plt
import model
import numpy
import os
def getPercentage(n1,n2):
    diff = n1-n2
    sum = n1+n2
    avrg = sum/2
    percentage = (diff/avrg)*100
    percentage = round(percentage,1)
    return str(percentage) + "%"
def addToCSV(model,ordernum,stringa,gapAc):
    output_path='Data/result.csv'
    GapPerc = getPercentage(model.ObjVal, gapAc)
    if(model.runTime <1800):
        print("Here")
        GapPerc = str(model.MIPGAP) +"%"
    performance = {
        "Settings": stringa,
        "Ordine": ordernum,
        "Area": model.ObjVal/1000000,
        "Gap": GapPerc,
        "Time": round(model.Runtime,1)}
    print(performance)
    df = pd.DataFrame([performance])
    df.to_csv(path_or_buf=output_path, mode='a',index=False, header=None)
    print(performance)

#Available: [7,8,10,17,20,22,23,24]
ordernumbers = [7,8,10,17]
for scenario in [1,2,3]:
    for nsmax in [1,2]:
        if scenario==3 and nsmax == 1:
            continue
        for a in ordernumbers:
            stringa ="Scenario:" + str(scenario) + " NsMax:" + str(nsmax)
            print(stringa)
            result, Ac_lower = model.optimise(a,scenario,nsmax)
            addToCSV(result,a,stringa,Ac_lower)
            
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



