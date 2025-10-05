from playwright.sync_api import sync_playwright
from bs4 import BeautifulSoup
import csv

def fetch_phytochemicals(plant_name):
    results = []
    
    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True,
            args=[
                '--no-sandbox', 
                '--disable-dev-shm-usage', 
                '--disable-gpu',
                '--disable-web-security',
                '--disable-features=VizDisplayCompositor',
                '--disable-extensions',
                '--disable-plugins',
                '--disable-images',
                '--disable-javascript',
                '--no-first-run',
                '--no-default-browser-check',
                '--disable-logging',
                '--silent'
            ]
        )
        
        context = browser.new_context(
            viewport={'width': 1280, 'height': 720},
            user_agent='Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        )
        page = context.new_page()
        
        # Block all unnecessary resources
        page.route('**/*', lambda route: (
            route.abort() if route.request.resource_type in ['image', 'stylesheet', 'font', 'media'] 
            else route.continue_()
        ))
        
        try:
            page.goto("https://cb.imsc.res.in/imppat/basicsearch/phytochemical", timeout=15000)
            page.fill("input[name=combine]", plant_name)
            page.evaluate("document.querySelector('input[type=submit][value=Search]').click()")
            page.wait_for_selector("table#table_id", timeout=15000)
            
            page_num = 1
            previous_count = 0
            max_pages = 200
            no_data_streak = 0
            
            while page_num <= max_pages:
                try:
                    html = page.evaluate("document.querySelector('table#table_id').innerHTML")
                except:
                    break
                
                try:
                    soup = BeautifulSoup(html, 'lxml')
                except:
                    soup = BeautifulSoup(html, 'html.parser')
                
                rows = soup.find_all("tr")[1:]
                
                if not rows:
                    break
                
                page_data = [[td.get_text(strip=True) for td in row.find_all("td")] 
                           for row in rows if len(row.find_all("td")) > 1]
                
                if not page_data:
                    no_data_streak += 1
                    if no_data_streak >= 3:
                        break
                else:
                    no_data_streak = 0
                    results.extend(page_data)
                
                if len(results) == previous_count:
                    break
                previous_count = len(results)
                
                try:
                    has_next = page.evaluate("""
                        () => {
                            const btn = document.querySelector('.dataTables_paginate .paginate_button.next:not(.disabled)');
                            return btn && !btn.classList.contains('disabled');
                        }
                    """)
                    
                    if has_next:
                        page.evaluate("document.querySelector('.dataTables_paginate .paginate_button.next').click()")
                        page.wait_for_selector("table#table_id", timeout=6000)
                        page_num += 1
                    else:
                        break
                except:
                    break
                    
        except:
            pass
        finally:
            context.close()
            browser.close()
    
    return results

def save_to_csv(data, plant_name):
    """Auto-save data to CSV"""
    if not data:
        return None
    
    safe_name = "".join(c for c in plant_name if c.isalnum() or c in (' ', '-', '_')).strip().replace(' ', '_')
    filename = f"phytochemicals_{safe_name}.csv"
    
    headers = [
        "Indian_medicinal_plant", "Plant_part", "IMPPAT_Phytochemical_identifier", 
        "Phytochemical_name", "References"
    ]
    
    if data and len(data[0]) != len(headers):
        headers = [f"Column_{i+1}" for i in range(len(data[0]))]
    
    try:
        with open(filename, 'w', newline='', encoding='utf-8', buffering=16384) as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(headers)
            writer.writerows(data)
        return filename
    except:
        return None

def main():
    """Main function"""
    while True:
        plant_name = input("Plant name: ").strip()
        if plant_name:
            break
    
    data = fetch_phytochemicals(plant_name)
    
    if data:
        filename = save_to_csv(data, plant_name)
        if filename:
            print(f"✓ {len(data)} phytochemicals → {filename}")
        else:
            print("❌ Save failed")
    else:
        print("❌ No data found")

if __name__ == "__main__":
    main()