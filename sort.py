import sys
import csv
import pandas as pd

csvData = pd.read_csv("data/result.csv")
csvData.sort_values(['Ordine'],axis=0,inplace=True)
print(csvData)