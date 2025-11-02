"""
IMPPAT Phytochemical Scraper with PubChem SMILES Integration
"""

import requests
from bs4 import BeautifulSoup
import time
import re
import pandas as pd
from requests.adapters import HTTPAdapter, Retry
import gzip
import os
from io import StringIO

PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi"
TIMEOUT = 30
RETRIES = 3
MAX_FETCH_RETRIES = 3
RETRY_DELAY = 5


class IMPPATScraper:
    def __init__(self):
        self.session = requests.Session()
        retries = Retry(total=3, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
        self.session.mount("https://", HTTPAdapter(max_retries=retries))
        self.headers = {
            "Content-Type": "application/x-www-form-urlencoded; charset=UTF-8",
            "X-Requested-With": "XMLHttpRequest",
            "Referer": "https://cb.imsc.res.in/imppat/basicsearch/phytochemical",
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
        }
        self.url = "https://cb.imsc.res.in/imppat/ajax"

    def fetch_phytochemicals(self, plant_name: str):
        payload = {
            "combine": plant_name,
            "pid": "",
            "pname": "",
            "phytochemical": "",
            "action": "phytochemical"
        }

        try:
            r = self.session.post(self.url, data=payload, headers=self.headers, timeout=25)
            r.raise_for_status()
            data = r.json()
        except:
            return []

        html = data.get("output", "")
        if not html:
            return []

        soup = BeautifulSoup(html, "lxml")
        table = soup.find("table", {"id": "table_id"})
        if not table or not table.find("tbody"):
            return []

        return [[td.get_text(strip=True) for td in tr.find_all("td")] 
                for tr in table.find("tbody").find_all("tr") if tr.find_all("td")]


def fetch_smiles(data, plant_name):
    TMP_FILE = "temp_names.txt"
    OUTPUT_GZ = "temp_results.txt.gz"
    
    session = requests.Session()
    retries = Retry(total=RETRIES, backoff_factor=1, status_forcelist=[500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    
    try:
        headers = ["Indian_medicinal_plant", "Plant_part", "IMPPAT_Phytochemical_identifier", 
                   "Phytochemical_name", "References"]
        
        df_original = pd.DataFrame(data, columns=headers)
        df_subset = df_original[["Phytochemical_name"]].dropna()
        df_subset.to_csv(TMP_FILE, sep="\t", index=False, header=False)
        
        data_payload = {
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
            files = {"idfile": (TMP_FILE, f, "application/octet-stream")}
            resp = session.post(PUBCHEM_URL, data=data_payload, files=files, 
                              allow_redirects=True, timeout=TIMEOUT)
        
        if resp.status_code != 200:
            return None
        
        m = re.search(r"https://pubchem\.ncbi\.nlm\.nih\.gov/rest/download/[^\s\"']+", resp.text)
        if not m:
            return None
        
        download_url = m.group(0)
        
        with session.get(download_url, stream=True, timeout=TIMEOUT) as r:
            r.raise_for_status()
            with open(OUTPUT_GZ, "wb") as f:
                for chunk in r.iter_content(chunk_size=16384):
                    if chunk:
                        f.write(chunk)
        
        with gzip.open(OUTPUT_GZ, "rt", encoding="utf-8") as f:
            smiles_text = f.read()
        
        df_smiles = pd.read_csv(StringIO(smiles_text), sep="\t", header=None,
                               names=["Phytochemical_name", "SMILES"])
        
        df_smiles["Phytochemical_name"] = df_smiles["Phytochemical_name"].str.strip()
        
        df_smiles_grouped = (
            df_smiles.groupby("Phytochemical_name")["SMILES"]
            .apply(lambda x: "|".join(str(s) for s in x if pd.notna(s)))
            .reset_index()
        )
        
        df_merged = df_original.merge(df_smiles_grouped, on="Phytochemical_name", how="left")
        
        safe_name = "".join(c for c in plant_name if c.isalnum() or c in (" ", "-", "_")).replace(" ", "_")
        output_file = f"phytochemicals_{safe_name}_with_smiles.csv"
        df_merged.to_csv(output_file, index=False)
        
        return output_file
        
    except:
        return None
        
    finally:
        for temp_file in [TMP_FILE, OUTPUT_GZ]:
            if os.path.exists(temp_file):
                try:
                    os.remove(temp_file)
                except:
                    pass


def main():
    scraper = IMPPATScraper()
    plant = input("Enter plant name: ").strip()

    if not plant:
        print("Error: Empty plant name")
        return

    print(f"Searching IMPPAT for: {plant}")
    start = time.time()

    data = scraper.fetch_phytochemicals(plant)
    if not data:
        print("No data found")
        return
    
    print(f"Found {len(data)} compounds ({time.time() - start:.2f}s)")
    print("Fetching SMILES from PubChem...")
    
    result = None
    for attempt in range(1, MAX_FETCH_RETRIES + 1):
        if attempt > 1:
            print(f"Retry {attempt}/{MAX_FETCH_RETRIES} (waiting {RETRY_DELAY}s...)")
            time.sleep(RETRY_DELAY)
        
        result = fetch_smiles(data, plant)
        if result:
            print(f"Success! Output: {result} ({time.time() - start:.2f}s total)")
            break
    
    if not result:
        print("Failed after all retries. Try again later.")


if __name__ == "__main__":
    main()