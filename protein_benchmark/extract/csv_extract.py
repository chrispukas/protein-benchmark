import os
import polars as pf

def create_length_index(link_to_database, save_dir, seq_columns_remap: dict[str: str]):
    if os.path.exists(save_dir):
        print(f"Index {save_dir} already exists. Skipping creation.")
        return

    print(f"Reading database from {link_to_database}")
    df = pf.read_csv(link_to_database)

    print(df.columns)
    print(df.head())

    rows_to_export = ['row_index']

    for col, len_name in seq_columns_remap.items():
        rows_to_export.append(len_name)

        df[len_name] = df[col].str.len()
    
    df['row_index'] = df.index
    df[rows_to_export].to_parquet(save_dir)
    print(df.columns)

    print(f"Saved length index to {save_dir}")

def reformat_db(db, seq_columns: list[str]):
   return db.groupby(seq_columns)['row_index'].agg(list).reset_index().set_index(seq_columns)