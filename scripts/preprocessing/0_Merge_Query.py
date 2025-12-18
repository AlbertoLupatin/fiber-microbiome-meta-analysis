import pandas as pd
import sys

query1_file = sys.argv[1]
query2_file = sys.argv[2]

query1 = pd.read_csv(query1_file)
query2 = pd.read_csv(query2_file)

all_query = pd.merge(query1, query2, on = "PMID", how = "outer")

print(all_query.columns)

all_query["First Author"] = all_query["First Author_x"].fillna(all_query["First Author_y"]).astype(str).str.replace(" ", "", regex=False) + "_" + all_query["Publication Year"].astype(str).str.replace(".0", "")
all_query["First Author"] = all_query["First Author"].str.replace("_nan", "")
all_query["DOI"] = all_query["DOI_y"].fillna(all_query["DOI_x"]).astype(str)
all_query["DOI"] = all_query["DOI"].apply(lambda x: x if x.startswith("http") else "https://doi.org/" + x)
all_query["Title"] = all_query["Title_y"].fillna(all_query["Title_x"])

columns_to_select = ["PMID", "First Author", "Title", "N", "Disease", "Intervention", "Total Fibers [g]", "Fiber Type", 
                     "Duration [weeks]", "Arms", "Crossover", "Country", "Region", "Data Availability", "Metadata Availability", "Filter", "Notes", "DOI"]
all_query = all_query[columns_to_select]

all_query.fillna("", inplace = True)

all_query.to_csv("../../Data/Query/Query_Final.csv")

print("Querys merged\nSaving and Exiting\n")
