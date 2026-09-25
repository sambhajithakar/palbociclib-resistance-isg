"""Append the qPCR supplementary tables to the workbook."""
import pandas as pd, openpyxl
X = "/mnt/user-data/outputs/Palbociclib_Resistance_Study/Supplementary_Tables.xlsx"
Q = "/mnt/user-data/outputs/Palbociclib_Resistance_Study/analysis/qpcr"
s10 = pd.read_csv(f"{Q}/qpcr_primers_final.csv")
s11 = pd.read_csv(f"{Q}/qpcr_insilico_pcr.csv")
s11 = s11.rename(columns={"n_F_sites": "forward primer annealing sites",
                          "n_R_sites": "reverse primer annealing sites",
                          "n_amplicons": "transcripts amplified",
                          "accessions": "RefSeq accessions amplified"})
with pd.ExcelWriter(X, engine="openpyxl", mode="a", if_sheet_exists="replace") as w:
    s10.to_excel(w, sheet_name="S10 qPCR primers", index=False)
    s11.to_excel(w, sheet_name="S11 in-silico PCR", index=False)
print(openpyxl.load_workbook(X).sheetnames)
