import pandas as pd


class code_dict:
    def __init__(self, code_df = pd.read_csv("code_lookup.csv"), 
                 entrycol = "sgg", keycol = "district_name"):
        self.df = code_df
        self.lookupdict = dict(zip(self.df[keycol], self.df[entrycol]))
    def lookup(self, name):
        return self.lookupdict[name]


              
