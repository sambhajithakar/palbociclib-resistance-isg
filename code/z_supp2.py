import pandas as pd, json, openpyxl
X = "/mnt/user-data/outputs/Palbociclib_Resistance_Study/Supplementary_Tables.xlsx"
O = "/home/claude/an/out"
s12 = pd.read_csv(f"{O}/feline_10x_scores.csv")
s12 = s12.rename(columns={s12.columns[0]: "biopsy"})
R = json.load(open(f"{O}/feline_10x_results.json"))
rows = []
for k, v in R.items():
    p = k.split("|")
    d = dict(analysis=p[0], arm=p[1], score=p[2], comparison=p[3] if len(p) > 3 else "")
    d.update(v); rows.append(d)
s13 = pd.DataFrame(rows)
with pd.ExcelWriter(X, engine="openpyxl", mode="a", if_sheet_exists="replace") as w:
    s12.to_excel(w, sheet_name="S12 FELINE biopsy scores", index=False)
    s13.to_excel(w, sheet_name="S13 FELINE tests", index=False)
print(openpyxl.load_workbook(X).sheetnames)
print(s12.shape, s13.shape)
