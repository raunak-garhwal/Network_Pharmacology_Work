import requests
from bs4 import BeautifulSoup
import csv
import time
from requests.adapters import HTTPAdapter, Retry

class IMPPATScraper:
    def __init__(self):
        # Use a persistent session for speed (keeps TCP connection alive)
        self.session = requests.Session()

        # Retries for robustness (network hiccups, 500 errors, etc.)
        retries = Retry(total=3, backoff_factor=1,
                        status_forcelist=[500, 502, 503, 504])
        self.session.mount("https://", HTTPAdapter(max_retries=retries))

        # Common headers (mimic browser)
        self.headers = {
            "Content-Type": "application/x-www-form-urlencoded; charset=UTF-8",
            "X-Requested-With": "XMLHttpRequest",
            "Referer": "https://cb.imsc.res.in/imppat/basicsearch/phytochemical",
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Safari/537.36"
        }

        self.url = "https://cb.imsc.res.in/imppat/ajax"

    def fetch_phytochemicals(self, plant_name: str):
        """Fetch phytochemicals for a single plant using the backend AJAX endpoint"""
        payload = {
            "combine": plant_name,
            "pid": "",
            "pname": "",
            "phytochemical": "",
            "action": "phytochemical"
        }

        try:
            r = self.session.post(self.url, data=payload, headers=self.headers, timeout=10)
            r.raise_for_status()
        except requests.RequestException as e:
            print(f"❌ Network error for {plant_name}: {e}")
            return []

        try:
            data = r.json()
        except ValueError:
            print(f"❌ Invalid JSON for {plant_name}")
            return []

        html = data.get("output", "")
        if not html:
            return []

        soup = BeautifulSoup(html, "lxml")
        table = soup.find("table", {"id": "table_id"})
        if not table or not table.find("tbody"):
            return []

        results = []
        for tr in table.find("tbody").find_all("tr"):
            cols = [td.get_text(strip=True) for td in tr.find_all("td")]
            if cols:
                results.append(cols)

        return results

    @staticmethod
    def save_to_csv(data, plant_name: str):
        """Save results to a CSV file"""
        if not data:
            return None

        safe_name = "".join(c for c in plant_name if c.isalnum() or c in (" ", "-", "_")).replace(" ", "_")
        filename = f"phytochemicals_{safe_name}.csv"

        headers = [
            "Indian_medicinal_plant",
            "Plant_part",
            "IMPPAT_Phytochemical_identifier",
            "Phytochemical_name",
            "References"
        ]
        if len(data[0]) != len(headers):
            headers = [f"Column_{i+1}" for i in range(len(data[0]))]

        try:
            with open(filename, "w", newline="", encoding="utf-8", buffering=16384) as f:
                writer = csv.writer(f)
                writer.writerow(headers)
                writer.writerows(data)
            return filename
        except Exception as e:
            print(f"❌ CSV save error: {e}")
            return None


def main():
    scraper = IMPPATScraper()
    plant = input("Enter plant name: ").strip()

    if not plant:
        print("❌ Error: Empty plant name")
        return

    print(f"🔍 Searching IMPPAT for: {plant} ...")
    start_time = time.time()

    data = scraper.fetch_phytochemicals(plant)

    if data:
        file = scraper.save_to_csv(data, plant)
        elapsed = time.time() - start_time
        print(f"✅ {len(data)} compounds saved to {file} (took {elapsed:.2f}s)")
    else:
        print("❌ No phytochemical data found")

if __name__ == "__main__":
    main()