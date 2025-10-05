#!/usr/bin/env python3
"""
idexchange_automate.py

Automate PubChem Identifier Exchange (idexchange.cgi) to convert names -> SMILES in bulk.

Usage:
    python idexchange_automate.py input.csv

Assumptions:
- Input CSV has a column named 'Phytochemical_name' (adjust COLUMN_NAME if different).
- Script uses Selenium to interact with the idexchange web UI and downloads results.
"""

import os
import sys
import time
import gzip
import shutil
import tempfile
import logging
from pathlib import Path
from typing import List

import pandas as pd
import requests
from tqdm import tqdm

# Selenium imports
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException, NoSuchElementException, ElementClickInterceptedException
from webdriver_manager.chrome import ChromeDriverManager
# Add this import at top of file

from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager

def start_driver(headless: bool = False):
    """Start Chrome webdriver via webdriver-manager (works with modern selenium)."""
    options = webdriver.ChromeOptions()
    if headless:
        # modern selenium supports "--headless=new"; fallback if needed
        try:
            options.add_argument("--headless=new")
        except Exception:
            options.add_argument("--headless")
        options.add_argument("--disable-gpu")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument(f"user-agent={USER_AGENT}")
    options.add_experimental_option("excludeSwitches", ["enable-logging"])

    # IMPORTANT: create Service object and pass it with service=...
    service = Service(ChromeDriverManager().install())
    # driver = webdriver.Chrome(service=service, options=options)
    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=options)
    driver.set_page_load_timeout(60)
    return driver


# ----------------- CONFIG -----------------
IDEXCHANGE_URL = "https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi"
COLUMN_NAME = "Phytochemical_name"    # change if your CSV uses a different column
CHUNK_SIZE = 2000                      # number of names per idexchange job (tune if needed)
WAIT_TIMEOUT = 300                     # seconds to wait for job completion/download link
DOWNLOAD_RETRY_DELAY = 3               # seconds between checking for result link
OUTPUT_SUFFIX = "_with_smiles_idex.csv"
FAILURES_SUFFIX = "_failures_idex.csv"
USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) AutoPubChemIdExchange/1.0"
# ------------------------------------------

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("idexchange_automate")

def prepare_text_file(names: List[str], path: str):
    """Write one name per line (Synonyms mode expects newline-separated)."""
    with open(path, "w", encoding="utf-8") as f:
        for n in names:
            # strip tabs/newlines inside names, preserve commas if needed
            line = str(n).strip().replace("\r", " ").replace("\n", " ")
            if line:
                f.write(line + "\n")

def start_driver(headless: bool = False):
    """Start Chrome webdriver via webdriver-manager."""
    options = webdriver.ChromeOptions()
    if headless:
        options.add_argument("--headless=new")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument(f'user-agent={USER_AGENT}')
    # optional: disable images / css for speed (not necessary)
    # driver = webdriver.Chrome(ChromeDriverManager().install(), options=options)
    service = Service(ChromeDriverManager().install())
    driver = webdriver.Chrome(service=service, options=options)
    driver.set_page_load_timeout(60)
    return driver

def find_select_by_option_text(driver, option_text_substr):
    """
    Search all <select> elements and return the first Select where any option text contains option_text_substr (case-insensitive).
    """
    selects = driver.find_elements(By.TAG_NAME, "select")
    for sel in selects:
        try:
            options = sel.find_elements(By.TAG_NAME, "option")
            for opt in options:
                if option_text_substr.lower() in opt.text.lower():
                    return Select(sel)
        except Exception:
            continue
    return None

def set_select_option_by_text(select_obj: Select, option_text_substr: str):
    """Choose the option whose visible text contains option_text_substr."""
    for opt in select_obj.options:
        if option_text_substr.lower() in opt.text.lower():
            select_obj.select_by_visible_text(opt.text)
            return True
    return False

def submit_job_and_get_result_link(driver, upload_filepath: str) -> str:
    """
    Interact with idexchange page:
    - set Input format = Synonyms (if available)
    - set Output type = SMILES
    - set Output method = Two column file showing input-output correspondence
    - upload file, submit job, wait for result URL and return it
    """
    driver.get(IDEXCHANGE_URL)
    wait = WebDriverWait(driver, 30)

    # Wait for the file input to appear (page is JS-driven)
    try:
        # wait until any file input present
        file_input = wait.until(EC.presence_of_element_located((By.XPATH, "//input[@type='file']")))
    except TimeoutException:
        raise RuntimeError("Timeout waiting for idexchange page to load the upload control.")

    # Set Input ID type -> Synonyms
    sel_input = find_select_by_option_text(driver, "Synonyms")
    if sel_input:
        set_select_option_by_text(sel_input, "Synonyms")
        logger.info("Set input type -> Synonyms")
    else:
        logger.warning("Could not find input-type select; leaving default")

    # Set Output type -> SMILES
    sel_output_type = find_select_by_option_text(driver, "SMILES")
    if sel_output_type:
        set_select_option_by_text(sel_output_type, "SMILES")
        logger.info("Set output type -> SMILES")
    else:
        logger.warning("Could not find output-type select; leaving default")

    # Set Output method -> two column mapping (if available)
    sel_outmethod = find_select_by_option_text(driver, "Two column")
    if sel_outmethod:
        set_select_option_by_text(sel_outmethod, "Two column")
        logger.info("Set output method -> Two column file")
    else:
        logger.warning("Could not find output-method select; leaving default")

    time.sleep(0.3)  # small pause

    # Upload file
    try:
        file_input = driver.find_element(By.XPATH, "//input[@type='file']")
        file_input.send_keys(os.path.abspath(upload_filepath))
        logger.info("Uploaded file: %s", upload_filepath)
    except Exception as e:
        raise RuntimeError(f"Failed to upload file: {e}")

    time.sleep(0.5)

    # Click Submit Job (button text may vary—try different selectors)
    submit_btn = None
    try:
        # try button element containing 'Submit Job'
        submit_btn = driver.find_element(By.XPATH, "//button[contains(., 'Submit Job') or contains(., 'Submit')]")
    except Exception:
        # fallback: input[type=submit]
        try:
            submit_btn = driver.find_element(By.XPATH, "//input[@type='submit' and (contains(@value, 'Submit') or contains(@value,'submit'))]")
        except Exception:
            submit_btn = None

    if not submit_btn:
        raise RuntimeError("Cannot find 'Submit Job' button on idexchange page.")

    try:
        submit_btn.click()
    except ElementClickInterceptedException:
        driver.execute_script("arguments[0].click();", submit_btn)

    logger.info("Job submitted, now waiting for result link...")

    # After submitting, the page becomes a self-refreshing waiting page.
    # Poll the page until an anchor with a downloadable result appears (.gz or .txt)
    result_href = None
    start = time.time()
    while True:
        # look for anchors with .gz, .bz2, .txt in href
        anchors = driver.find_elements(By.XPATH, "//a[contains(@href, '.gz') or contains(@href, '.bz2') or contains(@href, '.txt')]")
        if anchors:
            # choose first that is likely the result (heuristic)
            for a in anchors:
                href = a.get_attribute("href") or ""
                if href and ("download" in href or href.endswith(".gz") or href.endswith(".txt") or "idexchange" in href):
                    result_href = href
                    break
            if result_href:
                break

        # Some pages display the result URL as plain text; try to find pre or textarea
        try:
            pre = driver.find_element(By.TAG_NAME, "pre")
            txt = pre.text.strip()
            if txt and (txt.startswith("http") or "\nhttp" in txt):
                # extract first http...
                for token in txt.split():
                    if token.startswith("http"):
                        result_href = token
                        break
            if result_href:
                break
        except Exception:
            pass

        if time.time() - start > WAIT_TIMEOUT:
            raise RuntimeError("Timeout waiting for idexchange result link (increase WAIT_TIMEOUT).")
        time.sleep(DOWNLOAD_RETRY_DELAY)

    logger.info("Found result URL: %s", result_href)
    return result_href

def download_result(url: str, out_path: str, session=None):
    """Download result file via requests (faster & robust than browser click)"""
    s = session or requests.Session()
    s.headers.update({"User-Agent": USER_AGENT})
    r = s.get(url, stream=True, timeout=60)
    r.raise_for_status()
    with open(out_path, "wb") as fh:
        for chunk in r.iter_content(1024*32):
            if chunk:
                fh.write(chunk)
    logger.info("Downloaded result to %s", out_path)
    return out_path

def extract_if_gz(src_path: str, dst_path: str):
    """If gz, extract to dst_path, else move"""
    if src_path.endswith(".gz"):
        with gzip.open(src_path, "rb") as f_in:
            with open(dst_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        shutil.copy(src_path, dst_path)

def parse_two_column_result(txt_path: str) -> pd.DataFrame:
    """
    Parse two-column TSV style file. Format usually:
      input_value<TAB>output_value
    or same separated by whitespace.
    """
    # try tab-delimited first
    try:
        df = pd.read_csv(txt_path, sep="\t", header=None, names=["input", "output"], dtype=str, engine="python")
        # If only one column, retry whitespace split
        if df.shape[1] == 1 or df['output'].isnull().all():
            df = pd.read_csv(txt_path, sep=r"\s+", header=None, names=["input", "output"], dtype=str, engine="python")
    except Exception:
        # fallback: read lines and split first whitespace
        lines = []
        with open(txt_path, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    input_token = parts[0]
                    output_token = parts[-1]
                    lines.append((input_token, output_token))
        df = pd.DataFrame(lines, columns=["input", "output"])
    # normalize
    df['input'] = df['input'].astype(str).str.strip()
    df['output'] = df['output'].astype(str).str.strip().replace("", pd.NA)
    return df

def run_idexchange_for_csv(input_csv: str, output_csv: str = None, headless: bool = False):
    df = pd.read_csv(input_csv, dtype=str)
    if COLUMN_NAME not in df.columns:
        raise SystemExit(f"Column '{COLUMN_NAME}' not found in {input_csv}. Columns: {list(df.columns)}")

    names_series = df[COLUMN_NAME].fillna("").astype(str)
    unique_names = list(dict.fromkeys([n.strip() for n in names_series if n.strip()]))
    logger.info("Unique names to process: %d", len(unique_names))

    # chunking
    chunks = [unique_names[i:i+CHUNK_SIZE] for i in range(0, len(unique_names), CHUNK_SIZE)]
    logger.info("Divided into %d chunk(s) of up to %d names", len(chunks), CHUNK_SIZE)

    driver = start_driver(headless=headless)
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})

    merged_map = {}  # name -> smiles (or None)
    try:
        for idx, chunk in enumerate(chunks, start=1):
            logger.info("Processing chunk %d/%d (%d names)", idx, len(chunks), len(chunk))
            with tempfile.TemporaryDirectory() as tmpdir:
                txt_file = os.path.join(tmpdir, f"chunk_{idx}.txt")
                prepare_text_file(chunk, txt_file)

                # submit job and get result link
                try:
                    result_href = submit_job_and_get_result_link(driver, txt_file)
                except Exception as e:
                    logger.error("Failed to submit job or obtain result link for chunk %d: %s", idx, e)
                    # mark names as failed in this chunk
                    for n in chunk:
                        merged_map[n] = None
                    continue

                # download result
                local_dl = os.path.join(tmpdir, f"result_{idx}.gz" if result_href.endswith(".gz") else f"result_{idx}.txt")
                try:
                    download_result(result_href, local_dl, session=session)
                except Exception as e:
                    logger.error("Failed to download result for chunk %d: %s", idx, e)
                    for n in chunk:
                        merged_map[n] = None
                    continue

                # extract if needed
                extract_path = os.path.join(tmpdir, f"result_{idx}.txt")
                try:
                    extract_if_gz(local_dl, extract_path)
                except Exception as e:
                    logger.warning("Could not extract gz, trying to read raw file: %s", e)
                    extract_path = local_dl

                # parse two-column output into df
                try:
                    out_df = parse_two_column_result(extract_path)
                except Exception as e:
                    logger.error("Failed to parse result for chunk %d: %s", idx, e)
                    for n in chunk:
                        merged_map[n] = None
                    continue

                # merge results into merged_map
                for _, row in out_df.iterrows():
                    inp = str(row['input']).strip()
                    outv = None if pd.isna(row['output']) or row['output'] == '' else str(row['output']).strip()
                    merged_map[inp] = outv

                # any names in chunk not present in out_df -> None
                for n in chunk:
                    if n not in merged_map:
                        merged_map[n] = None

                logger.info("Chunk %d done: retrieved %d mappings", idx, out_df['output'].notna().sum())

    finally:
        try:
            driver.quit()
        except Exception:
            pass

    # Map back to original df
    df['SMILES_idexchange'] = df[COLUMN_NAME].map(lambda x: merged_map.get(str(x).strip(), None))
    # Save output
    if output_csv is None:
        output_csv = os.path.splitext(input_csv)[0] + OUTPUT_SUFFIX
    df.to_csv(output_csv, index=False)
    # failures
    failures = df[df['SMILES_idexchange'].isna() | (df['SMILES_idexchange'].astype(str).str.len() == 0)]
    failures_csv = os.path.splitext(input_csv)[0] + FAILURES_SUFFIX
    failures.to_csv(failures_csv, index=False)
    logger.info("Finished. Output saved to %s. Failures saved to %s", output_csv, failures_csv)
    return df

if __name__ == "__main__":
    # if len(sys.argv) < 2:
    #     print("Usage: python idexchange_automate.py <input.csv>")
    #     sys.exit(1)
    input_csv = "phytochemicals_haldi.csv"
    # optional: set headless True to run without UI
    df_out = run_idexchange_for_csv(input_csv, headless=False)
    print("Done. Rows with SMILES:", df_out['SMILES_idexchange'].notna().sum())
