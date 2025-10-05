import requests
import os
import re
import pandas as pd
import time
from requests.adapters import HTTPAdapter, Retry
import gzip
# ================= CONFIG =================
INPUT_CSV = "phytochemicals_turmeric.csv"
OUTPUT_GZ = "pubchem_results.txt.gz"
TMP_FILE = "names_only.txt"  # temp file for upload
PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi"
TIMEOUT = 30  # seconds for requests
RETRIES = 3
SLEEP_BETWEEN_RETRIES = 5
# =========================================

def prepare_temp_file(csv_file, tmp_file):
    """Extract two columns and write to temporary file for PubChem."""
    df = pd.read_csv(csv_file)
    required_cols = ["Phytochemical_name"]
    for col in required_cols:
        if col not in df.columns:
            raise SystemExit(f"❌ Column '{col}' not found in CSV")

    # Write tab-separated two-column file
    df_subset = df[required_cols].dropna()  
    df_subset.to_csv(tmp_file, sep="\t", index=False, header=False)
    print(f"✅ Temporary upload file '{tmp_file}' created with {len(df_subset)} rows.")

def submit_job(tmp_file):
    """Submit job to PubChem with retries."""
    session = requests.Session()
    retries = Retry(total=RETRIES, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))

    data = {
        "inputtype": "synofiltered",
        "inputdsn": "",
        "idinput": "file",
        "operatortype": "samecid",
        "outputtype": "smiles",
        "outputdsn": "",
        "method": "file-pair",
        "compression": "gzip",
        "submitjob": "Submit Job",
        "xmlfile": "",
    }

    with open(tmp_file, "rb") as f:
        files = {"idfile": (os.path.basename(tmp_file), f, "application/octet-stream")}
        print("Submitting job with file upload...")
        resp = session.post(PUBCHEM_URL, data=data, files=files, allow_redirects=True, timeout=TIMEOUT)

    if resp.status_code != 200:
        print(resp.text[:1000])
        raise SystemExit("❌ Submission failed")

    return resp, session

def extract_download_link(resp):
    """Extract download link or reqid page."""
    m = re.search(r"https://pubchem\.ncbi\.nlm\.nih\.gov/rest/download/[^\s\"']+", resp.text)
    if m:
        return m.group(0)
    elif "idexchange.cgi?reqid=" in resp.text:
        print("✅ Job submitted, PubChem is processing. Check the page:")
        print(resp.url)
        return None
    else:
        raise Exception("Could not find download link or reqid in response:\n" + resp.text[:800])

def download_results(download_url, session, output_file):
    """Download results from PubChem."""
    print("Downloading results...")
    with session.get(download_url, stream=True, timeout=TIMEOUT) as r:
        r.raise_for_status()
        with open(output_file, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
    print(f"✅ Done. Results saved to '{output_file}'")

def main():
    prepare_temp_file(INPUT_CSV, TMP_FILE)
    resp, session = submit_job(TMP_FILE)
    download_url = extract_download_link(resp)
    if download_url:
        download_results(download_url, session, OUTPUT_GZ)
    # Cleanup temp file
    if os.path.exists(TMP_FILE):
        os.remove(TMP_FILE)
        print(f"🗑️ Temporary file '{TMP_FILE}' removed.")

if __name__ == "__main__":
    start_time = time.time()
    try:
        main()
    except Exception as e:
        print("❌ Error:", e)
    finally:
        elapsed = time.time() - start_time
        print(f"⏱️ Total elapsed time: {elapsed:.2f} seconds")

import gzip
import pandas as pd
import os

# --- After the download logic ---

# Load original CSV
df_original = pd.read_csv(INPUT_CSV)

# Path to the extracted txt inside the gz
extracted_txt = "pubchem_results.txt"  # temporary extracted file

# Extract the gz file
with gzip.open(OUTPUT_GZ, "rt", encoding="utf-8") as f_in, open(extracted_txt, "w", encoding="utf-8") as f_out:
    f_out.write(f_in.read())

# Load the SMILES results
# Assuming tab-separated: Phytochemical_name\tSMILES
df_smiles = pd.read_csv(extracted_txt, sep="\t", header=None, names=["Phytochemical_name", "SMILES"])

# Strip whitespace to match correctly
df_smiles["Phytochemical_name"] = df_smiles["Phytochemical_name"].str.strip()

# Merge with original CSV on Phytochemical_name
df_smiles_grouped = (
    df_smiles.groupby("Phytochemical_name")["SMILES"]
    .apply(lambda x: "|".join(str(s) for s in x if pd.notna(s)))
    .reset_index()
)
df_merged = df_original.merge(df_smiles_grouped, on="Phytochemical_name", how="left")

# Save merged file
merged_file = "phytochemicals_with_smiles.csv"
df_merged.to_csv(merged_file, index=False)

print(f"✅ Merged file with SMILES saved to {merged_file}")

# Optional: clean up temporary extracted txt
os.remove(extracted_txt)
