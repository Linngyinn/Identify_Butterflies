from xbbg import blp
import pandas as pd
from datetime import date, timedelta
import numpy as np
from itertools import combinations
from xbbg.core import conn,process
import os

def _init_instrument_request(req,ticker): #Fill the request data
    req.set('ticker',ticker)
    req.set('partialMatch',False)
    req.set('maxResults',1000)

def _process_instruments(msg): #Process the response
    for elt in msg.asElement().getElement('results').values():
        yield elt.getElementAsString('parseky')

def allGovts(ticker): #Return all govts with the given ticker, matured or not
    request = process.create_request(service = '//blp/instruments',request='govtListRequest')
    _init_instrument_request(request,ticker)
    conn.send_request(request=request)
    return process.rec_events(func=_process_instruments)

def liveGovts(ticker): #Just return 'live' bonds, ordered by maturity

    dfAll = blp.bdp([g for g in allGovts(ticker)],['id_isin','maturity'])

    if cty == "ACGB":
        exclusion = [] #["BK467634 Corp", "BK467634 Corp"]
        dfAll = dfAll[dfAll['maturity'] >= today + timedelta(days = 365*2)]
    elif cty == "T":
        exclusion = list(blp.bds("TUH3 COMB Comdty", "FUT_DLVRBLE_BNDS_ISINS")["isin_of_deliverable_bonds"]) + list(blp.bds("FVH3 COMB Comdty", "FUT_DLVRBLE_BNDS_ISINS")["isin_of_deliverable_bonds"]) + list(blp.bds("UXYH3 COMB Comdty", "FUT_DLVRBLE_BNDS_ISINS")["isin_of_deliverable_bonds"]) + list(blp.bds("TYH3 COMB Comdty", "FUT_DLVRBLE_BNDS_ISINS")["isin_of_deliverable_bonds"]) + list(blp.bds("WNH3 COMB Comdty", "FUT_DLVRBLE_BNDS_ISINS")["isin_of_deliverable_bonds"])
        exclusion = [a[:-5] for a in exclusion]
        dfAll = dfAll[dfAll['maturity'] >= today + timedelta(days = 365*3)]
    elif cty == "JGB":
        exclusion = list(blp.bds("JBH3 COMB Comdty", "FUT_DLVRBLE_BNDS_ISINS")["isin_of_deliverable_bonds"])
        exclusion = [a[:-7] for a in exclusion]
        dfAll = dfAll[dfAll['maturity'] >= today + timedelta(days = 365*3)]
    
    dfAll = dfAll[~dfAll['id_isin'].isin(exclusion)]

    #filter inflation-linked and amt outstanding
    Inf_Ind = [blp.bdp(dfAll.index[i],"INFLATION_LINKED_INDICATOR") for i in range(len(dfAll.index))]
    AMT_OUT = [blp.bdp(dfAll.index[i],"amt_outstanding") for i in range(len(dfAll.index))]
    GB_df = pd.merge(dfAll, pd.concat(Inf_Ind, axis = 0), left_index=True, right_index=True)
    GB_df = pd.merge(GB_df, pd.concat(AMT_OUT, axis = 0), left_index=True, right_index=True)
    GB_df = GB_df[(GB_df["inflation_linked_indicator"] == "N") & (GB_df["amt_outstanding"] > 0)]

    return GB_df.sort_values('maturity')

def beta(x, y):
    covariance = np.cov(y,x) 
    return covariance[0,1]/covariance[1,1]

class Flys:
    
    def __init__(self, GB_df, start_date, History_Window, Repo, Rev_Repo, cty):
        
        self.GB_df = GB_df
        self.start_date = start_date
        self.History_Window = History_Window
        self.Repo = Repo
        self.Rev_Repo = Rev_Repo
        self.cty = cty
        
        self.get_yields(self.GB_df)
        self.GB_data(self.GB_df, self.GB_prices)
        fly_combo_base = self.create_flys(self.GB_df)
        self.run_stats(fly_combo_base)
        
    def get_yields(self, GB_df):
        
        today_yld = [blp.bdp(name, "yas_bond_yld") for name in GB_df.index]
        today_yld = pd.concat(today_yld, axis = 0).rename(columns={"yas_bond_yld": today}).T
        past_yld = [blp.bdh(name, "yld_ytm_mid", self.start_date, yesterday, Fill='P') for name in GB_df.index]
        past_yld_clean = []
        for df in past_yld:
            if not df.empty:
                df.columns = df.columns.droplevel(1)
                past_yld_clean.append(df)
            
        past_yld_clean = pd.concat(past_yld_clean, axis = 1)
        
        GB_prices = pd.concat([today_yld, past_yld_clean], axis = 0)
        #GB_prices = today_yld.merge(past_yld_clean, how="outer", on=list(today_yld.columns), left_index=True, right_index=True)
        GB_prices = GB_prices.sort_index(ascending=0)
        
        def non_match_elements(list_a, list_b):
            non_match = []
            for i in list_a:
                if i not in list_b:
                    non_match.append(i)
            return non_match
        
        non_match = non_match_elements(list(GB_df.index), list(past_yld_clean.columns))

        GB_df.drop(index = non_match,axis=1, inplace = True)

        self.GB_prices = GB_prices
        self.GB_df = GB_df
        
    def GB_data(self, GB_df, GB_prices):
        
        GB_df["datapoints"] = list(GB_prices.count(axis=0))
        GB_df["duration"] = [blp.bdp(col,"DUR_MID")["dur_mid"][0] for col in GB_prices]
        GB_df["time_maturity"] = list(((GB_df["maturity"] - today)/365)/np.timedelta64(1, 'D'))
        GB_df.drop(['inflation_linked_indicator', 'amt_outstanding'], axis=1, inplace = True)
        
        #drop if datapoints < History_Window
        GB_prices.drop(list(GB_df[GB_df['datapoints'] < self.History_Window].index), axis=1, inplace = True)
        GB_df.drop(GB_df[GB_df['datapoints'] < self.History_Window].index, inplace = True)
        
        GB_df["yield"] = list(GB_prices[GB_prices.index == today].values[0])
        GB_df["long_carry"] = (GB_df["yield"] - GB_df["duration"])/self.Repo/12*3
        GB_df["short_carry"] = (GB_df["yield"] - GB_df["duration"])/self.Rev_Repo/12*3
        
        self.GB_df = GB_df
        
    def create_flys(self, GB_df):
        
        lst = list(GB_df.index)
        comb = combinations(lst, 3)
        combi = [', '.join(i) for i in comb]
        combine = [i.split(", ") for i in combi]
        fly_combo = pd.DataFrame(combine, columns=["Bond 1", "Bond 2", "Bond 3"])
        
        if self.cty == "T" or self.cty == "JGB":
            fly_combo = pd.merge(fly_combo, GB_df["time_maturity"].to_frame().reset_index().rename(columns = {"index": "Bond 1"}), on = "Bond 1").rename(columns = {"time_maturity": "time_maturity 1"})
            fly_combo = pd.merge(fly_combo, GB_df["time_maturity"].to_frame().reset_index().rename(columns = {"index": "Bond 2"}), on = "Bond 2").rename(columns = {"time_maturity": "time_maturity 2"})
            fly_combo = pd.merge(fly_combo, GB_df["time_maturity"].to_frame().reset_index().rename(columns = {"index": "Bond 3"}), on = "Bond 3").rename(columns = {"time_maturity": "time_maturity 3"})
            fly_combo["steepener"] = (fly_combo["time_maturity 2"] - fly_combo["time_maturity 1"]).round(0)
            fly_combo["flattener"] = (fly_combo["time_maturity 3"] - fly_combo["time_maturity 2"]).round(0)
            fly_combo = fly_combo.loc[fly_combo["steepener"] == fly_combo["flattener"], :][["Bond 1", "Bond 2", "Bond 3"]]
    
        for num in fly_combo.columns[::-1]:
            fly_combo = pd.merge(GB_df["yield"].reset_index().rename(columns = {"index": "Bond " + str(num.split(" ")[-1])}), fly_combo, on = "Bond " + str(num.split(" ")[-1]), how = "right").rename(columns = {"yield": "Yield " + str(num.split(" ")[-1])})
            fly_combo = pd.merge(GB_df["long_carry"].reset_index().rename(columns = {"index": "Bond " + str(num.split(" ")[-1])}), fly_combo, on = "Bond " + str(num.split(" ")[-1]), how = "right").rename(columns = {"long_carry": "Long_Carry " + str(num.split(" ")[-1])})
            fly_combo = pd.merge(GB_df["short_carry"].reset_index().rename(columns = {"index": "Bond " + str(num.split(" ")[-1])}), fly_combo, on = "Bond " + str(num.split(" ")[-1]), how = "right").rename(columns = {"short_carry": "Short_Carry " + str(num.split(" ")[-1])})
            
        return fly_combo
            
    def run_stats(self, fly_combo):
        
        average, Stdev, Min, Max = [], [], [], []
        ten_perc, ninety_perc = [], []
        corr_b2, corr_b3_b1 = [], []
        beta_b2, beta_b3_b1, beta_b2_b1, beta_b3_b2 = [], [], [], []
        corr_1d, corr_5d, corr_25d = [], [], []
        for i in fly_combo.index:
            bonds = fly_combo[["Bond 1", "Bond 2", "Bond 3"]].iloc[i, :].values
            relevant_prices = self.GB_prices[bonds].iloc[:self.History_Window, :]
            series = 2*relevant_prices[relevant_prices.columns[1]] - relevant_prices[relevant_prices.columns[0]] - relevant_prices[relevant_prices.columns[2]]
            average.append(series.mean())
            Stdev.append(series.std())
            Min.append(series.min())
            Max.append(series.max())
            ten_perc.append(np.percentile(series, 10))
            ninety_perc.append(np.percentile(series, 90))
            corr_b2.append(series.corr(relevant_prices[relevant_prices.columns[1]]))
            corr_b3_b1.append(series.corr(relevant_prices[relevant_prices.columns[2]] - relevant_prices[relevant_prices.columns[0]]))
            beta_b2.append(beta(series, relevant_prices[relevant_prices.columns[1]]))
            beta_b3_b1.append(beta(series, relevant_prices[relevant_prices.columns[2]] - relevant_prices[relevant_prices.columns[0]]))
            beta_b2_b1.append(beta(series, relevant_prices[relevant_prices.columns[1]] - relevant_prices[relevant_prices.columns[0]]))
            beta_b3_b2.append(beta(series, relevant_prices[relevant_prices.columns[2]] - relevant_prices[relevant_prices.columns[1]]))
            
            relevant_prices_1d = self.GB_prices[bonds].iloc[1:1+self.History_Window, :]
            series1d = 2*relevant_prices_1d[relevant_prices_1d.columns[1]] - relevant_prices_1d[relevant_prices_1d.columns[0]] - relevant_prices_1d[relevant_prices_1d.columns[2]]
            corr_1d.append(series.corr(series1d))
            relevant_prices_5d = self.GB_prices[bonds].iloc[5:5+self.History_Window, :]
            series5d = 2*relevant_prices_5d[relevant_prices_5d.columns[1]] - relevant_prices_5d[relevant_prices_5d.columns[0]] - relevant_prices_5d[relevant_prices_5d.columns[2]]
            corr_5d.append(series.corr(series5d))
            relevant_prices_25d = self.GB_prices[bonds].iloc[25:25+self.History_Window, :]
            series25d = 2*relevant_prices_25d[relevant_prices_25d.columns[1]] - relevant_prices_25d[relevant_prices_25d.columns[0]] - relevant_prices_25d[relevant_prices_25d.columns[2]]
            corr_25d.append(series.corr(series25d))
            
        fly_combo["Fly Live"] = -100*(fly_combo["Yield 1"]+fly_combo["Yield 3"]-2*fly_combo["Yield 2"])
        fly_combo["Average"] = 100*np.array(average)
        fly_combo["Stdev"] = 100*np.array(Stdev)
        fly_combo["ZScore"] = (fly_combo["Fly Live"] - fly_combo["Average"])/fly_combo["Stdev"]
        fly_combo["Min"] = 100*np.array(Min)
        fly_combo["Max"] = 100*np.array(Max)
        fly_combo["10th Per"] = 100*np.array(ten_perc)
        fly_combo["90th Per"] = 100*np.array(ninety_perc)
        fly_combo["Corr with B2"] = corr_b2
        fly_combo["Corr with B3-B1"] = corr_b3_b1
        fly_combo["Beta of B2"] = beta_b2
        fly_combo["Beta of B3-B1"] = beta_b3_b1
        fly_combo["Beta of B2-B1"] = beta_b2_b1
        fly_combo["Beta of B3-B2"] = beta_b3_b2
        fly_combo["Serial Correl (1D)"] = corr_1d
        fly_combo["Serial Correl (5D)"] = corr_5d
        fly_combo["Serial Correl (25D)"] = corr_25d
        fly_combo["Long Belly Carry"] = 100*(2*fly_combo["Long_Carry 2"]-fly_combo["Long_Carry 1"]-fly_combo["Long_Carry 3"])
        fly_combo["Short Belly Carry"] = -100*(2*fly_combo["Short_Carry 2"]-fly_combo["Short_Carry 1"]-fly_combo["Short_Carry 3"])
        
        fly_combo.drop(['Short_Carry 1', 'Long_Carry 1', 'Yield 1', 'Short_Carry 2', 'Long_Carry 2', 'Yield 2', 'Short_Carry 3', 'Long_Carry 3', 'Yield 3'], axis=1, inplace = True)
        
        self.fly_combo = fly_combo
    
    #Receive Trades
    def filter_receive_trades(self):
        
        fly_combo = self.fly_combo
        
        fly_combo = fly_combo[(fly_combo["ZScore"] >= 1) & (fly_combo["Long Belly Carry"]>= -1) & (fly_combo["Stdev"]>=5)]
        fly_combo = fly_combo[fly_combo["Beta of B2-B1"].abs() <= 1]
        fly_combo = fly_combo[fly_combo["Beta of B3-B2"].abs() <= 1]
        fly_combo = fly_combo[fly_combo["Beta of B2"].abs() <= 1]
        #fly_combo = fly_combo[fly_combo["Serial Correl (25D)"].abs() <= 0.6]
        return fly_combo
    
    #Pay Trades
    def filter_pay_trades(self):
        
        fly_combo = self.fly_combo
        
        fly_combo = fly_combo[(fly_combo["ZScore"] <= -1) & (fly_combo["Short Belly Carry"]>= -1) & (fly_combo["Stdev"]>=5)]
        fly_combo = fly_combo[fly_combo["Beta of B2-B1"].abs() <= 1]
        fly_combo = fly_combo[fly_combo["Beta of B3-B2"].abs() <= 1]
        fly_combo = fly_combo[fly_combo["Beta of B2"].abs() <= 1]
        #fly_combo = fly_combo[fly_combo["Serial Correl (25D)"].abs() <= 0.6]
        return fly_combo
        
        
    
if __name__ == "__main__":
    
    today = date.today()
    yesterday = today - timedelta(days = 1)
    
    '''get history of data'''
    start_date = "1/12/2016"
    
    cty = 'ACGB' #ACGB, T, JGB
    
    GB_df = liveGovts(cty) 

    History_Window = 252*2
    Repo = 3.16 #3.16/3.8/-0.05
    Rev_Repo = 3 #3/3.7/-0.15
    
    obj = Flys(GB_df, start_date, History_Window, Repo, Rev_Repo, cty)
    unfiltered = obj.fly_combo
    receive_trades = obj.filter_receive_trades()
    pay_trades = obj.filter_pay_trades()
    
    if os.path.exists(r"Z:\Flys_Ideas\Output/" + str(today)):
            
        with pd.ExcelWriter(r"Z:\Flys_Ideas\Output/" + str(today) + "/Fly_" + str(cty) + "_" +  str(today) + ".xlsx") as writer:  
            receive_trades.to_excel(writer, sheet_name='Receive Flys')
            pay_trades.to_excel(writer, sheet_name='Pay Flys')
            unfiltered.to_excel(writer, sheet_name='Unfiltered')
                
    else:
        os.mkdir(r"Z:\Flys_Ideas\Output/" + str(today))

        with pd.ExcelWriter(r"Z:\Flys_Ideas\Output/" + str(today) + "/Fly_" + str(cty) + "_" + str(today) + ".xlsx") as writer:  
            receive_trades.to_excel(writer, sheet_name='Receive Flys')
            pay_trades.to_excel(writer, sheet_name='Pay Flys')
            unfiltered.to_excel(writer, sheet_name='Unfiltered')
    
