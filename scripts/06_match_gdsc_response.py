import os
import pandas as pd

ccle_path = "data/sample_metadata.csv"
gdsc1_path = "data/GDSC1_fitted_dose_response.xlsx"
gdsc2_path = "data/GDSC2_fitted_dose_response.xlsx"
outdir = "shiny_data"


def normalize_name(name):
    if pd.isna(name):
        return ""
    return str(name).upper().replace("-", "").replace(" ", "").replace("_", "").replace(".", "")


def match_gdsc(gdsc_path, ccle_lookup):
    gdsc = pd.read_excel(gdsc_path)
    gdsc["name_norm"] = gdsc["CELL_LINE_NAME"].apply(normalize_name)
    matched = gdsc.merge(ccle_lookup, on="name_norm", how="inner")

    keep_cols = [
        "CellLine", "DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET",
        "PATHWAY_NAME", "LN_IC50", "AUC", "MIN_CONC", "MAX_CONC", "Z_SCORE"
    ]
    if "RMSE" in matched.columns:
        keep_cols.append("RMSE")
    if "CANCER_TYPE" in matched.columns:
        keep_cols.append("CANCER_TYPE")
    if "SANGER_MODEL_ID" in matched.columns:
        keep_cols.append("SANGER_MODEL_ID")

    return matched[[c for c in keep_cols if c in matched.columns]].copy()


ccle = pd.read_csv(ccle_path)
ccle["name_norm"] = ccle["CellLine"].apply(normalize_name)
ccle_lookup = ccle[["CellLine", "name_norm"]].drop_duplicates("name_norm")

gdsc1_matched = match_gdsc(gdsc1_path, ccle_lookup)
gdsc2_matched = match_gdsc(gdsc2_path, ccle_lookup)

gdsc1_matched.to_csv(os.path.join(outdir, "gdsc1_matched.csv"), index=False)
gdsc2_matched.to_csv(os.path.join(outdir, "gdsc2_matched.csv"), index=False)

print(f"GDSC1: {gdsc1_matched['CellLine'].nunique()} cell lines, {gdsc1_matched['DRUG_NAME'].nunique()} drugs")
print(f"GDSC2: {gdsc2_matched['CellLine'].nunique()} cell lines, {gdsc2_matched['DRUG_NAME'].nunique()} drugs")
