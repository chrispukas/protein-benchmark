import os
import datetime as dt

import incito_pipeline.util.query.rcsbapi as rcs
import incito_pipeline.util.data_utility as du



def get_pdb_dirs_by_date(pdb_dir: str, cutoff_date: dt.date = dt.datetime(2023, 1, 6)):
    print(f"Filtering dir: {pdb_dir} with cutoffdate {str(cutoff_date)}")
    files = du.get_file_names(pdb_dir, ".pdb")

    pdb_dirs = []

    for f in files:
        file_path = os.path.join(pdb_dir, f)

        pdb_id = du.get_pdb_id(f)
        add_date = rcs.rcsbapi_query_date(pdb_id) if cutoff_date else None
        if add_date is None or add_date < cutoff_date.strftime('%Y-%m-%d'):
            continue

        print(f"Appending {pdb_id}, date: {add_date}, cutoff_date: {cutoff_date}")

        pdb_dirs.append(file_path)

    return pdb_dirs