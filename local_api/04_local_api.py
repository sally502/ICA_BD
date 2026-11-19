
import requests
import csv

# ==========================================
# 1. Configuration
# ==========================================

BASE_URL = "http://bioinfmsc6.bio.ed.ac.uk:8080"
SUMMARY_URL = BASE_URL + "/summary"

GENES = [
    "BRCA1",
    "BRCA2",
    "TP53",
    "EGFR",
    "MYC",
    "APOE",
    "CFTR",
    "ACE2",
    "TNF",
    "VEGFA"
]

OUTPUT_FILE = "local_api_data.tsv"


# ==========================================
# 2. Fetch data from Local API
# ==========================================

def fetch_local_data(genes):

    params = {
        "identifiers": ",".join(genes),
        "sleeptime": 1
    }

    print("[INFO] Fetching Local API data...")

    try:
        response = requests.get(
            SUMMARY_URL,
            params=params,
            timeout=300
        )

        response.raise_for_status()

        data = response.json()

        if isinstance(data, dict) and "error" in data:
            print("[ERROR] API returned:", data["error"])
            return None

        if not isinstance(data, list):
            print("[ERROR] Unexpected data format.")
            return None

        print(f"[INFO] Retrieved {len(data)} records.")

        return data

    except requests.RequestException as error:
        print("[ERROR] Request failed:", error)
        return None

    except ValueError:
        print("[ERROR] Response is not valid JSON.")
        return None


# ==========================================
# 3. Save data as TSV
# ==========================================

def save_tsv(data):

    if not data:
        print("[ERROR] No data to save.")
        return

    fieldnames = [
        "ncbi_gene_id",
        "official_symbol",
        "description",
        "status",
        "current_gene_id",
        "chromosome",
        "map_location",
        "organism_scientific",
        "organism_common",
        "other_aliases",
        "other_designations"
    ]

    with open(
        OUTPUT_FILE,
        "w",
        newline="",
        encoding="utf-8"
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore"
        )

        writer.writeheader()

        for record in data:
            writer.writerow(record)

    print(f"[SUCCESS] TSV saved to {OUTPUT_FILE}")


# ==========================================
# 4. Main
# ==========================================

def main():

    data = fetch_local_data(GENES)

    if data is None:
        return

    save_tsv(data)

    returned_genes = {
        record.get("official_symbol")
        for record in data
    }

    missing_genes = set(GENES) - returned_genes

    if missing_genes:
        print("[WARNING] Missing genes:",
              ", ".join(sorted(missing_genes)))
    else:
        print("[SUCCESS] All 10 genes retrieved!")


if __name__ == "__main__":
    main()
