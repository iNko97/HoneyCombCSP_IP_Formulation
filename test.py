import pandas as pd
import order

datapath = "Data/input_data.ods"
datasheet = pd.ExcelFile(datapath, engine="odf")
ordine = order.Order(datasheet)
print(ordine.I)