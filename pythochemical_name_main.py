"""
PubChem SMILES Fetcher - Optimized Production Version
"""

import requests
import pandas as pd
import time
import os
import re
import gzip
from requests.adapters import HTTPAdapter, Retry

PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi"


class PubChemSMILES:
    def __init__(self):
        self.session = requests.Session()
        self.session.mount("https://", HTTPAdapter(
            max_retries=Retry(total=3, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
        ))

    def fetch_smiles(self, compound_names):
        """Fetch SMILES for compound names"""
        if isinstance(compound_names, str):
            compound_names = [compound_names]
        
        tmp_file = "temp_compounds.txt"
        gz_file = "temp_smiles.txt.gz"
        
        try:
            # Create input file
            with open(tmp_file, "w", encoding="utf-8") as f:
                f.write("\n".join(compound_names))
            
            # Submit job
            with open(tmp_file, "rb") as f:
                resp = self.session.post(
                    PUBCHEM_URL,
                    data={
                        "inputtype": "synofiltered",
                        "idinput": "file",
                        "operatortype": "samecid",
                        "outputtype": "smiles",
                        "method": "file-pair",
                        "compression": "gzip",
                        "submitjob": "Submit Job",
                        "inputdsn": "",
                        "outputdsn": "",
                        "xmlfile": "",
                    },
                    files={"idfile": (os.path.basename(tmp_file), f, "application/octet-stream")},
                    allow_redirects=True,
                    timeout=30
                )
            
            if resp.status_code != 200:
                return None
            
            # Extract download link
            m = re.search(r"https://pubchem\.ncbi\.nlm\.nih\.gov/rest/download/[^\s\"']+", resp.text)
            if not m:
                return None
            
            # Download results
            with self.session.get(m.group(0), stream=True, timeout=30) as r:
                r.raise_for_status()
                with open(gz_file, "wb") as f:
                    for chunk in r.iter_content(8192):
                        if chunk:
                            f.write(chunk)
            
            # Extract SMILES
            results = []
            with gzip.open(gz_file, "rt", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t", 1)
                    if len(parts) == 2:
                        results.append({"Phytochemical_name": parts[0], "SMILES": parts[1]})
            
            return results
            
        except:
            return None
        finally:
            for f in [tmp_file, gz_file]:
                if os.path.exists(f):
                    os.remove(f)


def get_smiles(compound_names, max_retries=3):
    """
    Get SMILES for compound name(s) with retry logic
    
    Args:
        compound_names: str or list of compound names
        max_retries: number of retry attempts (default: 3)
    
    Returns:
        tuple: (DataFrame, elapsed_time) or (None, 0)
    """
    fetcher = PubChemSMILES()
    start = time.time()
    
    for attempt in range(max_retries):
        results = fetcher.fetch_smiles(compound_names)
        if results:
            elapsed = time.time() - start
            return pd.DataFrame(results), elapsed
        
        if attempt < max_retries - 1:
            time.sleep(2 ** attempt)  # Exponential backoff: 1s, 2s, 4s
    
    return None, 0


def main():
    compound = input("Enter compound name(s) [comma-separated]: ").strip()
    
    if not compound:
        print("❌ Empty input")
        return
    
    compounds = [c.strip() for c in compound.split(",")] if "," in compound else [compound]
    
    df, elapsed = get_smiles(compounds)
    
    if df is not None:
        filename = "compound_smiles.csv"
        df.to_csv(filename, index=False)
        print(f"✅ {filename} | {len(df)} compounds | {elapsed:.1f}s")
    else:
        print("❌ Failed - PubChem may be processing, retry in a moment")


if __name__ == "__main__":
    main()