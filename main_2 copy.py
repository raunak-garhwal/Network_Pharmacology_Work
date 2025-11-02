"""
PubChem SMILES Fetcher
Fetch SMILES structures from PubChem for compounds in a CSV file
"""

import requests
import os
import re
import pandas as pd
import time
from requests.adapters import HTTPAdapter, Retry
import gzip

PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi"
TIMEOUT = 30
RETRIES = 3


def fetch_smiles(input_csv):
    """
    Fetch SMILES from PubChem and merge with input CSV
    Args:
        input_csv (str): Path to input CSV with 'Phytochemical_name' column
    Returns:
        str: Path to output file with SMILES
    """
    TMP_FILE = "names_only.txt"
    OUTPUT_GZ = "pubchem_results.txt.gz"
    EXTRACTED_TXT = "pubchem_results.txt"
    
    # Setup session
    session = requests.Session()
    retries = Retry(total=RETRIES, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    
    try:
        # Step 1: Prepare temp file
        print(f"Loading {input_csv}")
        df_original = pd.read_csv(input_csv)
        
        if "Phytochemical_name" not in df_original.columns:
            raise ValueError("Column 'Phytochemical_name' not found in CSV")
        
        df_subset = df_original[["Phytochemical_name"]].dropna()
        df_subset.to_csv(TMP_FILE, sep="\t", index=False, header=False)
        print(f"Prepared {len(df_subset)} compounds for upload")
        
        # Step 2: Submit job to PubChem
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
        
        with open(TMP_FILE, "rb") as f:
            files = {"idfile": (os.path.basename(TMP_FILE), f, "application/octet-stream")}
            print("Submitting job to PubChem...")
            resp = session.post(PUBCHEM_URL, data=data, files=files, allow_redirects=True, timeout=TIMEOUT)
        
        if resp.status_code != 200:
            raise Exception(f"Submission failed with status {resp.status_code}")
        
        # Step 3: Extract download link
        m = re.search(r"https://pubchem\.ncbi\.nlm\.nih\.gov/rest/download/[^\s\"']+", resp.text)
        if m:
            download_url = m.group(0)
        elif "idexchange.cgi?reqid=" in resp.text:
            print("Job submitted, PubChem is processing")
            print(f"Check status at: {resp.url}")
            return None
        else:
            raise Exception("Could not find download link in response")
        
        # Step 4: Download results
        print("Downloading results...")
        with session.get(download_url, stream=True, timeout=TIMEOUT) as r:
            r.raise_for_status()
            with open(OUTPUT_GZ, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
        print("Download complete")
        
        # Step 5: Extract and process SMILES
        print("Processing SMILES data...")
        with gzip.open(OUTPUT_GZ, "rt", encoding="utf-8") as f_in:
            with open(EXTRACTED_TXT, "w", encoding="utf-8") as f_out:
                f_out.write(f_in.read())
        
        df_smiles = pd.read_csv(EXTRACTED_TXT, sep="\t", header=None, 
                               names=["Phytochemical_name", "SMILES"])
        
        df_smiles["Phytochemical_name"] = df_smiles["Phytochemical_name"].str.strip()
        
        df_smiles_grouped = (
            df_smiles.groupby("Phytochemical_name")["SMILES"]
            .apply(lambda x: "|".join(str(s) for s in x if pd.notna(s)))
            .reset_index()
        )
        
        df_merged = df_original.merge(df_smiles_grouped, on="Phytochemical_name", how="left")
        
        # Step 6: Save output
        output_file = input_csv.replace(".csv", "_with_smiles.csv")
        df_merged.to_csv(output_file, index=False)
        print(f"Output saved to: {output_file}")
        
        return output_file
        
    finally:
        # Cleanup temp files
        for temp_file in [TMP_FILE, OUTPUT_GZ, EXTRACTED_TXT]:
            if os.path.exists(temp_file):
                os.remove(temp_file)


if __name__ == "__main__":
    input_csv = "phytochemicals_tulsi.csv"
    
    start_time = time.time()
    try:
        result = fetch_smiles(input_csv)
        if result:
            print(f"Success: {result}")
        else:
            print("Job submitted, check PubChem later")
    except Exception as e:
        print(f"Error: {e}")
    finally:
        elapsed = time.time() - start_time
        print(f"Total time: {elapsed:.2f} seconds")