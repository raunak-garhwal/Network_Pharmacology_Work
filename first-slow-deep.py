from playwright.sync_api import sync_playwright
from bs4 import BeautifulSoup
import time
import csv

def fetch_phytochemicals(plant_name):
    results = []
    
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        
        try:
            # Navigate to the search page
            print(f"Searching for phytochemicals in: {plant_name}")
            page.goto("https://cb.imsc.res.in/imppat/basicsearch/phytochemical", timeout=60000)
            
            # Fill the search form
            page.fill("input[name=combine]", plant_name)
            page.click("input[type=submit][value=Search]")
            
            # Wait for results table to load
            page.wait_for_selector("table#table_id", timeout=60000)
            
            page_num = 1
            previous_results_count = 0
            max_pages = 50  # Safety limit to prevent infinite loops
            no_new_data_count = 0  # Track consecutive pages with no new data
            
            while page_num <= max_pages:
                print(f"Processing page {page_num}...")
                
                # Get the HTML content of the table
                try:
                    html = page.inner_html("table#table_id")
                    soup = BeautifulSoup(html, "html.parser")
                except Exception as e:
                    print(f"Error getting table content: {e}")
                    break
                
                # Extract data from all rows (skip header row)
                rows = soup.find_all("tr")[1:]
                if not rows:
                    print("No data rows found on this page.")
                    break
                
                current_page_data = []
                for row in rows:
                    cells = [td.get_text(strip=True) for td in row.find_all("td")]
                    if cells and len(cells) > 1:  # Only add non-empty rows with actual data
                        current_page_data.append(cells)
                
                # Check if we got new data
                if not current_page_data:
                    no_new_data_count += 1
                    print(f"No new data on page {page_num}")
                    if no_new_data_count >= 2:  # Stop if 2 consecutive pages have no data
                        print("No new data found on consecutive pages. Stopping.")
                        break
                else:
                    no_new_data_count = 0  # Reset counter if we found data
                    results.extend(current_page_data)
                    print(f"Extracted {len(current_page_data)} compounds from page {page_num}")
                
                # Check if we're getting duplicate data (same as previous page)
                if len(results) == previous_results_count:
                    print("No new results added. Likely reached end of data.")
                    break
                
                previous_results_count = len(results)
                
                # Look for pagination controls
                pagination_exists = False
                next_btn = None
                
                # Try to find the next button with various selectors
                selectors = [
                    ".dataTables_paginate .paginate_button.next:not(.disabled)",
                    ".pagination .next:not(.disabled)",
                    "a.next:not(.disabled)",
                    ".paginate_button.next:not(.disabled)"
                ]
                
                for selector in selectors:
                    try:
                        next_btn = page.query_selector(selector)
                        if next_btn and next_btn.is_visible():
                            class_attr = next_btn.get_attribute("class") or ""
                            if "disabled" not in class_attr.lower():
                                pagination_exists = True
                                break
                    except:
                        continue
                
                # If no specific next button, look for pagination links
                if not pagination_exists:
                    try:
                        pagination_links = page.query_selector_all(".dataTables_paginate a, .pagination a")
                        for link in pagination_links:
                            if link.is_visible():
                                link_text = link.inner_text().strip().lower()
                                class_attr = link.get_attribute("class") or ""
                                if link_text in ['next', '>', '»', str(page_num + 1)] and "disabled" not in class_attr.lower():
                                    next_btn = link
                                    pagination_exists = True
                                    break
                    except:
                        pass
                
                # Try to click next button
                if pagination_exists and next_btn:
                    try:
                        print(f"Clicking next button to go to page {page_num + 1}")
                        next_btn.click()
                        
                        # Wait for page to load
                        time.sleep(3)
                        
                        # Wait for table to update
                        try:
                            page.wait_for_selector("table#table_id", timeout=15000)
                            # Additional wait for data to load
                            time.sleep(2)
                        except:
                            print("Timeout waiting for next page to load")
                            break
                            
                        page_num += 1
                        
                    except Exception as e:
                        print(f"Error navigating to next page: {e}")
                        break
                else:
                    print("No more pages available or next button not found.")
                    break
                    
        except Exception as e:
            print(f"Error during scraping: {e}")
        finally:
            browser.close()
    
    return results

def save_to_csv(data, plant_name, filename=None):
    """Save the scraped data to a CSV file"""
    if not data:
        print("No data to save.")
        return
    
    # Generate filename if not provided
    if not filename:
        safe_plant_name = "".join(c for c in plant_name if c.isalnum() or c in (' ', '-', '_')).strip()
        safe_plant_name = safe_plant_name.replace(' ', '_')
        filename = f"phytochemicals_{safe_plant_name}.csv"
    
    # Create headers based on typical IMPPAT structure
    # You may need to adjust these based on the actual table structure
    headers = [
        "Phytochemical_ID", 
        "Phytochemical_Name", 
        "Molecular_Formula", 
        "Molecular_Weight", 
        "Plant_Name",
        "Plant_Part",
        "Activity"
    ]
    
    # If we have data, use the first row to determine the number of columns
    if data and len(data[0]) != len(headers):
        # Adjust headers based on actual data structure
        headers = [f"Column_{i+1}" for i in range(len(data[0]))]
    
    try:
        with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.writer(csvfile)
            
            # Write headers
            writer.writerow(headers)
            
            # Write data
            writer.writerows(data)
        
        print(f"Data saved to {filename}")
        print(f"Total records saved: {len(data)}")
        
    except Exception as e:
        print(f"Error saving to CSV: {e}")

def main():
    """Main function to handle user input and orchestrate the scraping"""
    print("Phytochemical Database Scraper")
    print("=" * 40)
    
    # Get plant name from user
    while True:
        plant_name = input("Enter the plant name (e.g., 'Curcuma longa'): ").strip()
        if plant_name:
            break
        print("Please enter a valid plant name.")
    
    # Optional: Ask for custom filename
    custom_filename = input("Enter custom filename (press Enter for auto-generated): ").strip()
    if not custom_filename:
        custom_filename = None
    elif not custom_filename.endswith('.csv'):
        custom_filename += '.csv'
    
    print(f"\nStarting search for: {plant_name}")
    print("This may take several minutes depending on the number of results...")
    
    # Fetch the data
    try:
        data = fetch_phytochemicals(plant_name)
        
        if data:
            print(f"\nSuccessfully scraped {len(data)} phytochemical records!")
            
            # Display first few records as preview
            print("\nPreview of first 3 records:")
            for i, record in enumerate(data[:3]):
                print(f"Record {i+1}: {record}")
            
            # Save to CSV
            save_to_csv(data, plant_name, custom_filename)
            
        else:
            print("No data found for the specified plant.")
            
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()