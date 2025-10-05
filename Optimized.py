from playwright.sync_api import sync_playwright
from bs4 import BeautifulSoup
import csv
import time

def fetch_phytochemicals(plant_name: str, max_pages: int = 200):
    results = []
    
    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True,
            args=[
                "--no-sandbox", "--disable-dev-shm-usage", "--disable-gpu",
                "--disable-web-security", "--disable-features=VizDisplayCompositor",
                "--disable-extensions", "--disable-plugins", "--no-first-run",
                "--no-default-browser-check", "--disable-logging",
                "--disable-background-timer-throttling",
                "--disable-renderer-backgrounding",
                "--disable-blink-features=AutomationControlled",
                "--disable-ipc-flooding-protection"
            ]
        )
        
        context = browser.new_context(
            viewport={"width": 1280, "height": 720},
            user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
        )
        page = context.new_page()
        
        # Block unnecessary resources
        page.route("**/*", lambda route: (
            route.abort() if route.request.resource_type in [
                "image", "stylesheet", "font", "media", "websocket",
                "manifest", "other"
            ] or any(x in route.request.url for x in ["favicon", "analytics", "google"])
            else route.continue_()
        ))
        
        try:
            # Load IMPPAT
            page.goto("https://cb.imsc.res.in/imppat/basicsearch/phytochemical", timeout=15000)
            page.fill("input[name=combine]", plant_name)
            page.click("input[type=submit][value=Search]")
            
            # Wait for table to load
            page.wait_for_selector("table#table_id.dataTable", timeout=15000)
            page.wait_for_function("""
                () => {
                    const info = document.querySelector('.dataTables_info');
                    return info && info.textContent.includes('entries');
                }
            """, timeout=10000)
            
            page_num, previous_count, no_data_streak = 1, 0, 0
            
            while page_num <= max_pages:
                try:
                    # Fast inner_html extraction
                    html = page.locator("table#table_id").inner_html()
                except Exception:
                    break

                try:
                    soup = BeautifulSoup(html, "lxml")
                except Exception:
                    soup = BeautifulSoup(html, "html.parser")

                tbody = soup.find("tbody")
                if not tbody:
                    break
                rows = tbody.find_all("tr")

                page_data = []
                for row in rows:
                    tds = row.find_all("td")
                    if len(tds) >= 4:
                        cols = [td.get_text(strip=True) for td in tds]
                        if cols[0] and cols[3]:
                            page_data.append(cols)

                if not page_data:
                    no_data_streak += 1
                    if no_data_streak >= 3:  # consecutive empty pages
                        break
                else:
                    no_data_streak = 0
                    results.extend(page_data)

                # Prevent infinite loop
                if len(results) == previous_count:
                    break
                previous_count = len(results)

                # Next page check
                has_next = page.evaluate("""
                    () => {
                        const btn = document.querySelector('.dataTables_paginate .paginate_button.next');
                        return btn && !btn.classList.contains('disabled');
                    }
                """)

                if not has_next:
                    break

                # Go next
                page.evaluate("document.querySelector('.dataTables_paginate .paginate_button.next').click()")
                page.wait_for_function("""
                    () => {
                        const tbody = document.querySelector('table#table_id tbody');
                        return tbody && tbody.children.length > 0;
                    }
                """, timeout=8000)

                page_num += 1
                time.sleep(0.3)  # tiny delay for stability

        except Exception as e:
            print(f"Scraping error: {e}")
        finally:
            context.close()
            browser.close()

    return results


def save_to_csv(data, plant_name: str):
    """Save results to CSV with safe headers"""
    if not data:
        return None
    
    safe_name = "".join(c for c in plant_name if c.isalnum() or c in (" ", "-", "_")).strip().replace(" ", "_")
    filename = f"phytochemicals_{safe_name}.csv"

    # Default IMPPAT headers
    headers = ["Indian_medicinal_plant", "Plant_part", "IMPPAT_Phytochemical_identifier", "Phytochemical_name", "References"]

    if len(data[0]) != len(headers):
        headers = [f"Column_{i+1}" for i in range(len(data[0]))]

    try:
        with open(filename, "w", newline="", encoding="utf-8", buffering=16384) as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(data)
        return filename
    except Exception as e:
        print(f"File save error: {e}")
        return None


def main():
    print("=== IMPPAT Phytochemical Scraper ===")
    plant_name = input("Enter plant name: ").strip()

    if not plant_name:
        print("❌ Error: Empty plant name")
        return

    print(f"🔍 Searching IMPPAT for: {plant_name} ...")

    data = fetch_phytochemicals(plant_name)

    if data:
        filename = save_to_csv(data, plant_name)
        if filename:
            print(f"✅ Success: {len(data)} compounds saved to {filename}")
        else:
            print("❌ Could not save data")
    else:
        print("❌ No phytochemical data found")


if __name__ == "__main__":
    main()
