"""
GeneCards Production Scraper
Optimized for reliability and completeness
"""

try:
    import cloudscraper
except ImportError:
    print("ERROR: Install cloudscraper with: pip install cloudscraper")
    exit(1)

from bs4 import BeautifulSoup
import pandas as pd
import time

BASE_URL = "https://www.genecards.org/Search/Keyword"
PAGE_SIZES = [10000, 5000, 2000, 1000, 500, 100]
REQUEST_DELAY = 1.0


class GeneCardsScraper:
    def __init__(self):
        self.scraper = cloudscraper.create_scraper(
            browser={'browser': 'chrome', 'platform': 'windows', 'mobile': False}
        )
        self.page_size = None

    def fetch(self, query, page):
        """Fetch single page"""
        try:
            r = self.scraper.get(
                BASE_URL,
                params={"queryString": query, "pageSize": self.page_size, "startPage": page},
                timeout=60
            )
            if r.status_code == 200 and "searchResults" in r.text:
                return r.text
        except:
            pass
        return None

    def parse(self, html):
        """Parse genes from HTML"""
        soup = BeautifulSoup(html, "lxml")
        table = soup.find("table", {"id": "searchResults"})
        if not table or not table.find("tbody"):
            return []
        
        results = []
        for tr in table.find("tbody").find_all("tr"):
            cols = tr.find_all("td")
            if len(cols) >= 9:
                try:
                    results.append([
                        cols[0].get_text(strip=True),
                        (cols[2].find("a") or cols[2]).get_text(strip=True),
                        cols[3].get_text(strip=True),
                        cols[4].get_text(strip=True),
                        cols[5].get_text(strip=True),
                        cols[6].get_text(strip=True),
                        cols[7].get_text(strip=True),
                        cols[8].get_text(strip=True)
                    ])
                except:
                    continue
        return results

    def get_total(self, html):
        """Extract total results"""
        soup = BeautifulSoup(html, "lxml")
        for link in soup.find_all("a", {"data-ga-total-results": True}):
            total = link.get("data-ga-total-results")
            if total:
                return int(total)
        return None

    def find_page_size(self, query):
        """Find largest working page size"""
        for size in PAGE_SIZES:
            try:
                r = self.scraper.get(
                    BASE_URL,
                    params={"queryString": query, "pageSize": size, "startPage": 0},
                    timeout=60
                )
                if r.status_code == 200 and "searchResults" in r.text:
                    soup = BeautifulSoup(r.text, "lxml")
                    if soup.find("table", {"id": "searchResults"}):
                        self.page_size = size
                        return r.text
            except:
                continue
        return None

    def scrape(self, query):
        """Scrape all pages sequentially"""
        html = self.find_page_size(query)
        if not html:
            return []
        
        all_results = self.parse(html)
        total = self.get_total(html)
        
        if not total:
            return all_results
        
        pages = (total // self.page_size) + (1 if total % self.page_size else 0)
        
        for page in range(1, pages):
            time.sleep(REQUEST_DELAY)
            html = self.fetch(query, page)
            if html:
                all_results.extend(self.parse(html))
        
        return all_results

    def save(self, data, query):
        """Save to CSV"""
        if not data:
            return None
        
        filename = f"genecards_{''.join(c for c in query if c.isalnum() or c in ' -_').replace(' ', '_')}.csv"
        
        df = pd.DataFrame(data, columns=[
            "Index", "Symbol", "Description", "Category", 
            "UniProt_ID", "GIFtS", "GC_id", "Score"
        ])
        
        original = len(df)
        df = df.drop_duplicates(subset=["Symbol", "GC_id"], keep="first")
        df["_score"] = pd.to_numeric(df["Score"], errors="coerce")
        df = df.sort_values("_score", ascending=False).drop("_score", axis=1)
        df.to_csv(filename, index=False)
        
        return filename, len(df), original - len(df)


def scrape_disease(disease):
    """
    Scrape disease genes from GeneCards
    
    Args:
        disease: Disease name to search
    
    Returns:
        tuple: (filename, gene_count, duplicates_removed, time_elapsed)
    """
    scraper = GeneCardsScraper()
    start = time.time()
    
    data = scraper.scrape(disease)
    
    if not data:
        return None, 0, 0, 0
    
    result = scraper.save(data, disease)
    elapsed = time.time() - start
    
    if result:
        filename, count, dupes = result
        return filename, count, dupes, elapsed
    
    return None, 0, 0, elapsed


def main():
    disease = input("Enter disease name: ").strip()
    if not disease:
        print("Error: Empty query")
        return
    
    print(f"Scraping GeneCards for: {disease}")
    
    filename, count, dupes, elapsed = scrape_disease(disease)
    
    if filename:
        print(f"\n✅ Success: {filename}")
        print(f"Genes: {count:,} | Time: {elapsed:.1f}s")
        if dupes:
            print(f"Duplicates removed: {dupes:,}")
    else:
        print("\n❌ Failed: No data retrieved")


if __name__ == "__main__":
    main()