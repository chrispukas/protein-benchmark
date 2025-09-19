

import requests


def rcsbapi_query_date(pdb_id):
    base_url = "https://data.rcsb.org/rest/v1/core/entry/"
    url = f"{base_url}{pdb_id.upper()}"
    response = requests.get(url)
    print(f"Querying {url}")

    if response.status_code == 200:
        data = response.json()
        accession_info = data.get("rcsb_accession_info", {})
        date = accession_info.get("deposit_date")
        return date
    else:
        print(f"Failed to fetch {pdb_id}: {response.status_code}")
        return None


def rcsbapi_query_dates(pdb_ids):
    deposit_dates = []

    for pdb_id in pdb_ids:
        date = rcsbapi_query_date(pdb_id)
        deposit_dates.append(date)

    return deposit_dates