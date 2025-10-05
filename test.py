import requests
from bs4 import BeautifulSoup
import csv

def fetch_phytochemicals_direct(plant_name: str):
    url = "https://cb.imsc.res.in/imppat/ajax"

    # Form payload - match exactly what you saw in DevTools
    payload = {
        "combine": plant_name,      # search term
        "pid": "",
        "pname": "",
        "phytochemical": "",
        "action": "phytochemical"   # ← sometimes it must match exactly
    }

    headers = {
        "Content-Type": "application/x-www-form-urlencoded; charset=UTF-8",
        "X-Requested-With": "XMLHttpRequest",   # required for many AJAX backends
        "Referer": "https://cb.imsc.res.in/imppat/basicsearch/phytochemical",
        "User-Agent": "Mozilla/5.0"
    }

    response = requests.post(url, data=payload, headers=headers)
    response.raise_for_status()

    data = response.json()

    html = data.get("output", "")
    if not html:
        return []

    soup = BeautifulSoup(html, "lxml")
    table = soup.find("table", {"id": "table_id"})
    if not table:
        return []

    results = []
    for tr in table.find("tbody").find_all("tr"):
        cols = [td.get_text(strip=True) for td in tr.find_all("td")]
        if cols:
            results.append(cols)

    return results


def save_to_csv(data, plant_name: str):
    if not data:
        return None

    safe_name = "".join(c for c in plant_name if c.isalnum() or c in (" ", "-", "_")).replace(" ", "_")
    filename = f"phytochemicals_{safe_name}.csv"

    headers = ["Indian_medicinal_plant", "Plant_part", "IMPPAT_Phytochemical_identifier", "Phytochemical_name", "References"]
    if len(data[0]) != len(headers):
        headers = [f"Column_{i+1}" for i in range(len(data[0]))]

    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(data)

    return filename


if __name__ == "__main__":
    plant = input("Enter plant name: ").strip()
    data = fetch_phytochemicals_direct(plant)
    if data:
        file = save_to_csv(data, plant)
        print(f"✅ {len(data)} compounds saved to {file}")
    else:
        print("❌ No data found (check spelling or payload fields)")
