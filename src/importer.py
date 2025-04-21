# /Users/yuewang/Documents/uniprot_importer/importer.py

import pandas as pd

def process_go_annotations(input_path: str, output_path: str):
    df = pd.read_csv("go_annotations.csv")

    columns_to_drop = ["DB_Object_Symbol", "Aspect", "DB_Object_Name", 
                       "Synonym", "DB_Object_Type", "Taxon", 
                       "Annotation_Extension", "Gene_Product_Form_ID"]
    
    df = df.drop(columns=columns_to_drop)

    df["subject"] = df["DB"] + ":" + df["DB_Object_ID"]
    df = df.drop(columns=["DB", "DB_Object_ID"])
    cols = ["subject"] + [col for col in df.columns if col != "subject"]
    df = df[cols]

    df = df.rename(columns={
        "Qualifier": "predicate", 
        "GO_ID": "object",
        "DB_Reference": "publications", 
        "With_From": "supporting_objects",
        "Date": "annotation_date",
        "Assigned_By": "primary_knowledge_source"
    })

    df["annotation_date"] = pd.to_datetime(df["annotation_date"], format='%Y%m%d').dt.strftime('%Y-%m-%d')
    df["annotation_date"] = df["annotation_date"].astype(str)
   

    df["aggregator"] = "UniProt"
    df["protocol_id"] = "NULL"
    df["supporting_objects"] = df["supporting_objects"].fillna("NULL")

    allowed_predicates = ["enables", "contributes_to", "acts_upstream_of_or_within", "involved_in",
                          "acts_upstream_of", "acts_upstream_of_positive_effect",
                          "acts_upstream_of_negative_effect", "acts_upstream_of_or_within_negative_effect",
                          "acts_upstream_of_or_within_positive_effect", "located_in", "part_of", 
                          "is_active_in", "colocalizes_with"]

    df["negated"] = ~df["predicate"].isin(allowed_predicates)
    df["predicate"] = df["predicate"].str.replace(r"^NOT\|", "", regex=True)

    # Evidence mapping
    url = "https://raw.githubusercontent.com/evidenceontology/evidenceontology/master/gaf-eco-mapping.txt"
    mapping_columns = ["Evidence_Code", "DB_Reference", "evidence_type"]
    df2 = pd.read_csv(url, sep="\t", comment="#", header=None, names=mapping_columns)

    df["publications"] = df["publications"].str.strip().str.upper()
    df2["DB_Reference"] = df2["DB_Reference"].str.strip().str.upper()
    df2["Evidence_Code"] = df2["Evidence_Code"].str.strip().str.upper()

    merged = df.merge(
        df2,
        how="left",
        left_on=["Evidence_Code", "publications"],
        right_on=["Evidence_Code", "DB_Reference"]
    )
    
    unmatched_mask = merged["evidence_type"].isna()
    fallback_df = df2[df2["DB_Reference"].str.upper() == "DEFAULT"]

    fallback_match = merged.loc[unmatched_mask, ["Evidence_Code"]].merge(
        fallback_df,
        how="left",
        on="Evidence_Code"
    )

    merged.loc[unmatched_mask, "evidence_type"] = fallback_match["evidence_type"].values
    merged = merged.drop(columns=["DB_Reference"])

    merged.to_csv(output_path, index=False)

if __name__ == "__main__":
    process_go_annotations(
        input_path="go_annotations.csv",  
        output_path="new_annotations1.csv"
    )
#     df = pd.read_csv("new_annotations1.csv")
#     print(df.head())
#     print(df.dtypes)  
