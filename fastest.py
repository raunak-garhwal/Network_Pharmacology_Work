from playwright.sync_api import sync_playwright
from bs4 import BeautifulSoup
import csv

def fetch_phytochemicals(plant_name: str):
    results = []

    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True,
            args=[
                '--no-sandbox', '--disable-dev-shm-usage', '--disable-gpu', 
                '--disable-web-security', '--disable-features=VizDisplayCompositor',
                '--disable-extensions', '--disable-plugins', '--no-first-run', 
                '--no-default-browser-check', '--disable-logging', '--disable-background-timer-throttling',
                '--disable-backgrounding-occluded-windows', '--disable-renderer-backgrounding',
                '--disable-blink-features=AutomationControlled', '--disable-ipc-flooding-protection'
            ]
        )

        context = browser.new_context(
            viewport={'width': 1280, 'height': 720},
            user_agent='Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        )
        page = context.new_page()
        
        # Block all unnecessary resources for maximum speed
        page.route('**/*', lambda route: (
            route.abort() if route.request.resource_type in
            ['image', 'stylesheet', 'font', 'media', 'websocket', 'manifest', 'other'] 
            or 'favicon' in route.request.url
            or 'analytics' in route.request.url
            or 'google' in route.request.url
            else route.continue_()
        ))

        try:
            # Navigate directly to search page
            page.goto("https://cb.imsc.res.in/imppat/basicsearch/phytochemical", timeout=8000)
            
            # Fill and submit form immediately
            page.fill("input[name=combine]", plant_name)
            page.click("input[type=submit][value=Search]")
            
            # Wait for DataTables to initialize - IMPPAT specific
            page.wait_for_selector("table#table_id.dataTable", timeout=8000)
            
            # IMPPAT-specific: Wait for initial data load
            page.wait_for_function("""
                () => {
                    const info = document.querySelector('.dataTables_info');
                    return info && info.textContent.includes('entries');
                }
            """, timeout=5000)

            page_count = 0
            while True:
                page_count += 1
                
                # Get table HTML using faster method
                html = page.locator("table#table_id").inner_html()
                soup = BeautifulSoup(html, 'html.parser')
                
                # IMPPAT uses tbody structure
                tbody = soup.find("tbody")
                if not tbody:
                    break
                
                rows = tbody.find_all("tr")
                if not rows:
                    break

                # Process rows with minimal validation
                page_data = []
                for row in rows:
                    tds = row.find_all("td")
                    if tds and len(tds) >= 4:  # IMPPAT has at least 4-5 columns
                        cols = [td.get_text(strip=True) for td in tds]
                        if cols[0] and cols[3]:  # Plant name and phytochemical name should exist
                            page_data.append(cols)

                if not page_data:
                    break

                results.extend(page_data)

                # IMPPAT-specific pagination check
                next_enabled = page.evaluate("""
                    () => {
                        const nextBtn = document.querySelector('.dataTables_paginate .paginate_button.next');
                        return nextBtn && !nextBtn.classList.contains('disabled');
                    }
                """)
                
                if not next_enabled:
                    break

                # Click next with immediate execution
                page.evaluate("document.querySelector('.dataTables_paginate .paginate_button.next').click()")
                
                # Minimal wait for IMPPAT DataTables refresh
                page.wait_for_function("""
                    () => {
                        const tbody = document.querySelector('table#table_id tbody');
                        return tbody && tbody.children.length > 0;
                    }
                """, timeout=6000)

        except Exception as e:
            if "timeout" in str(e).lower():
                print(f"Request timeout: Plant '{plant_name}' may not exist in IMPPAT database.")
            else:
                print(f"Scraping error: {e}")
        finally:
            try:
                context.close()
                browser.close()
            except:
                pass

    return results

def save_to_csv(data, plant_name: str):
    if not data:
        print(f"No phytochemical data found for: {plant_name}")
        return None
    
    # Quick filename generation
    safe_name = "".join(c for c in plant_name if c.isalnum() or c in (' ', '-', '_')).strip().replace(' ', '_')
    filename = f"phytochemicals_{safe_name}.csv"
    
    # IMPPAT standard headers
    headers = ["Indian_medicinal_plant", "Plant_part", "IMPPAT_Phytochemical_identifier", "Phytochemical_name", "References"]
    
    # Auto-adjust for different column counts
    if data and len(data[0]) != 5:
        headers = [f"Col_{i+1}" for i in range(len(data[0]))]
    
    try:
        with open(filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(data)
        return filename
    except Exception as e:
        print(f"File save error: {e}")
        return None

def main():
    print("IMPPAT Phytochemical Scraper")
    print("-" * 28)
    
    plant_name = input("Plant name: ").strip()
    
    if not plant_name:
        print("Error: Empty plant name")
        return
    
    print(f"Searching IMPPAT for: {plant_name}")
    
    try:
        data = fetch_phytochemicals(plant_name)
        
        if not data:
            print("No compounds found. Possible reasons:")
            print("- Plant not in IMPPAT database")
            print("- Incorrect spelling")
            print("- No phytochemical data available")
            return
        
        filename = save_to_csv(data, plant_name)
        
        if filename:
            print(f"Success: {len(data)} compounds saved to {filename}")
        else:
            print("Error: Could not save data")
            
    except KeyboardInterrupt:
        print("\nCancelled by user")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()