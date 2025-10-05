import pandas as pd
import requests
import time
import json
from urllib.parse import quote, unquote
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from collections import defaultdict
import asyncio
import aiohttp
import nest_asyncio
from typing import Optional, Tuple, Dict, List
import logging

# Enable nested event loops for Jupyter compatibility
try:
    nest_asyncio.apply()
except:
    pass

class OptimizedSMILESRetriever:
    def __init__(self, max_workers=50, use_async=True):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
        self.max_workers = max_workers
        self.use_async = use_async
        self.lock = threading.Lock()
        self.processed_count = 0
        self.success_count = 0
        self.cache = {}  # In-memory cache for processed compounds
        
        # Setup logging
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
        # Setup session with optimized connection pooling
        self.session = requests.Session()
        
        retry_strategy = Retry(
            total=3,
            backoff_factor=0.5,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET"]
        )
        
        adapter = HTTPAdapter(
            pool_connections=30,
            pool_maxsize=30,
            max_retries=retry_strategy
        )
        
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)
        self.session.headers.update({
            'User-Agent': 'Optimized SMILES Retriever v2.0',
            'Accept': 'application/json',
            'Connection': 'keep-alive'
        })
        
    def clean_compound_name(self, name: str) -> Optional[str]:
        """Enhanced compound name cleaning with better validation"""
        if pd.isna(name) or not name:
            return None
            
        name = str(name).strip()
        if not name or len(name) < 2:
            return None
            
        # Enhanced cleaning
        # Remove extra whitespace
        name = re.sub(r'\s+', ' ', name)
        
        # Remove common problematic characters but preserve chemical notation
        name = re.sub(r'[^\w\s\-\+\(\)\[\],\.\'\"\:\;]', '', name)
        
        # Handle common naming variations
        name = name.replace('_', '-')
        name = re.sub(r'\s*\([^)]*\)$', '', name)  # Remove trailing parentheses with descriptions
        
        return name.strip() if name.strip() else None
    
    def get_multiple_search_terms(self, compound_name: str) -> List[str]:
        """Generate multiple search variants for better matching"""
        cleaned = self.clean_compound_name(compound_name)
        if not cleaned:
            return []
            
        search_terms = [cleaned]
        
        # Add variations
        if '-' in cleaned:
            search_terms.append(cleaned.replace('-', ' '))
            search_terms.append(cleaned.replace('-', ''))
            
        if ' ' in cleaned:
            search_terms.append(cleaned.replace(' ', '-'))
            search_terms.append(cleaned.replace(' ', ''))
            
        # Remove duplicates while preserving order
        seen = set()
        unique_terms = []
        for term in search_terms:
            if term not in seen:
                seen.add(term)
                unique_terms.append(term)
                
        return unique_terms[:3]  # Limit to top 3 variants
    
    async def async_get_smiles(self, session: aiohttp.ClientSession, compound_name: str) -> Tuple[str, Optional[str], str]:
        """Async SMILES retrieval with multiple search strategies"""
        # Check cache first
        if compound_name in self.cache:
            return compound_name, self.cache[compound_name], 'cached'
            
        search_terms = self.get_multiple_search_terms(compound_name)
        if not search_terms:
            return compound_name, None, 'invalid_name'
        
        # Try multiple search strategies
        for search_term in search_terms:
            try:
                # Strategy 1: Name search
                encoded_name = quote(search_term)
                url = f"{self.base_url}/name/{encoded_name}/property/CanonicalSMILES/json"
                
                async with session.get(url, timeout=aiohttp.ClientTimeout(total=10)) as response:
                    if response.status == 200:
                        data = await response.json()
                        smiles = self._extract_smiles_from_response(data)
                        if smiles:
                            self.cache[compound_name] = smiles
                            return compound_name, smiles, 'success_name'
                            
                # Strategy 2: Try synonym search if name search fails
                url = f"{self.base_url}/name/{encoded_name}/synonyms/json"
                async with session.get(url, timeout=aiohttp.ClientTimeout(total=8)) as response:
                    if response.status == 200:
                        data = await response.json()
                        # If we find synonyms, try the first few
                        synonyms = data.get('InformationList', {}).get('Information', [])
                        if synonyms and 'Synonym' in synonyms[0]:
                            first_synonym = synonyms[0]['Synonym'][0]
                            # Recursively search with synonym
                            encoded_syn = quote(first_synonym)
                            smiles_url = f"{self.base_url}/name/{encoded_syn}/property/CanonicalSMILES/json"
                            async with session.get(smiles_url, timeout=aiohttp.ClientTimeout(total=8)) as smiles_response:
                                if smiles_response.status == 200:
                                    smiles_data = await smiles_response.json()
                                    smiles = self._extract_smiles_from_response(smiles_data)
                                    if smiles:
                                        self.cache[compound_name] = smiles
                                        return compound_name, smiles, 'success_synonym'
                                        
            except asyncio.TimeoutError:
                continue
            except Exception as e:
                self.logger.debug(f"Error processing {search_term}: {str(e)}")
                continue
        
        return compound_name, None, 'not_found'
    
    def sync_get_smiles(self, compound_name: str) -> Tuple[str, Optional[str], str]:
        """Synchronous SMILES retrieval with fallback strategies"""
        # Check cache first
        if compound_name in self.cache:
            return compound_name, self.cache[compound_name], 'cached'
            
        search_terms = self.get_multiple_search_terms(compound_name)
        if not search_terms:
            return compound_name, None, 'invalid_name'
        
        for search_term in search_terms:
            try:
                # Primary search
                encoded_name = quote(search_term)
                url = f"{self.base_url}/name/{encoded_name}/property/CanonicalSMILES/json"
                
                response = self.session.get(url, timeout=10)
                
                if response.status_code == 200:
                    data = response.json()
                    smiles = self._extract_smiles_from_response(data)
                    if smiles:
                        self.cache[compound_name] = smiles
                        with self.lock:
                            self.success_count += 1
                        return compound_name, smiles, 'success'
                
                # Fallback: Try CID-based search
                cid_url = f"{self.base_url}/name/{encoded_name}/cids/json"
                cid_response = self.session.get(cid_url, timeout=8)
                
                if cid_response.status_code == 200:
                    cid_data = cid_response.json()
                    cids = cid_data.get('IdentifierList', {}).get('CID', [])
                    if cids:
                        cid = cids[0]  # Use first CID
                        smiles_url = f"{self.base_url}/cid/{cid}/property/CanonicalSMILES/json"
                        smiles_response = self.session.get(smiles_url, timeout=8)
                        
                        if smiles_response.status_code == 200:
                            smiles_data = smiles_response.json()
                            smiles = self._extract_smiles_from_response(smiles_data)
                            if smiles:
                                self.cache[compound_name] = smiles
                                with self.lock:
                                    self.success_count += 1
                                return compound_name, smiles, 'success_cid'
                
            except Exception as e:
                self.logger.debug(f"Error processing {search_term}: {str(e)}")
                continue
        
        return compound_name, None, 'not_found'
    
    def _extract_smiles_from_response(self, data: dict) -> Optional[str]:
        """Robust SMILES extraction from PubChem response"""
        try:
            # Handle different response structures
            if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                properties = data['PropertyTable']['Properties']
                if isinstance(properties, list) and properties:
                    smiles = properties[0].get('CanonicalSMILES')
                    if smiles and isinstance(smiles, str) and len(smiles) > 3:
                        return smiles.strip()
            
            # Alternative structure
            if 'Properties' in data and isinstance(data['Properties'], list):
                for prop in data['Properties']:
                    if 'CanonicalSMILES' in prop:
                        smiles = prop['CanonicalSMILES']
                        if smiles and isinstance(smiles, str) and len(smiles) > 3:
                            return smiles.strip()
            
        except Exception as e:
            self.logger.debug(f"Error extracting SMILES: {str(e)}")
        
        return None
    
    def update_progress(self, compound_name: str):
        """Thread-safe progress tracking with enhanced stats"""
        with self.lock:
            self.processed_count += 1
            if self.processed_count % 25 == 0:  # More frequent updates
                success_rate = (self.success_count / self.processed_count) * 100
                cache_hits = len([v for v in self.cache.values() if v is not None])
                print(f"⚡ {self.processed_count} processed | ✓ {self.success_count} SMILES | 📋 {cache_hits} cached | {success_rate:.1f}% success")
    
    async def async_process_compounds(self, compounds: List[str]) -> Dict[str, Optional[str]]:
        """Async processing of compounds"""
        results = {}
        
        # Setup async session with connection limits
        connector = aiohttp.TCPConnector(
            limit=50,
            limit_per_host=25,
            ttl_dns_cache=300,
            use_dns_cache=True
        )
        
        timeout = aiohttp.ClientTimeout(total=15)
        
        async with aiohttp.ClientSession(
            connector=connector,
            timeout=timeout,
            headers={'User-Agent': 'Optimized Async SMILES Retriever v2.0'}
        ) as session:
            
            # Create semaphore to limit concurrent requests
            semaphore = asyncio.Semaphore(self.max_workers)
            
            async def process_with_semaphore(compound):
                async with semaphore:
                    result = await self.async_get_smiles(session, compound)
                    self.update_progress(compound)
                    return result
            
            # Process compounds concurrently
            tasks = [process_with_semaphore(compound) for compound in compounds]
            
            for coro in asyncio.as_completed(tasks):
                name, smiles, status = await coro
                results[name] = smiles
                
                if status.startswith('success'):
                    with self.lock:
                        self.success_count += 1
        
        return results
    
    def sync_process_compounds(self, compounds: List[str]) -> Dict[str, Optional[str]]:
        """Synchronous processing with threading"""
        results = {}
        
        def process_compound(compound):
            result = self.sync_get_smiles(compound)
            self.update_progress(compound)
            return result
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_compound = {
                executor.submit(process_compound, compound): compound 
                for compound in compounds
            }
            
            for future in as_completed(future_to_compound):
                compound = future_to_compound[future]
                try:
                    name, smiles, status = future.result()
                    results[name] = smiles
                except Exception as e:
                    self.logger.error(f"Error processing {compound}: {str(e)}")
                    results[compound] = None
        
        return results
    
    def process_dataset(self, csv_file_path: str, compound_column: str = 'Phytochemical_name', 
                       output_file_path: Optional[str] = None) -> pd.DataFrame:
        """Main processing function with enhanced error handling"""
        
        print(f"🚀 OPTIMIZED SMILES RETRIEVER v2.0")
        print(f"📁 Processing: {csv_file_path}")
        print("="*70)
        
        # Read CSV with error handling
        try:
            start_time = time.time()
            print("📖 Reading CSV...")
            df = pd.read_csv(csv_file_path, encoding='utf-8')
            read_time = time.time() - start_time
            print(f"✅ CSV loaded in {read_time:.2f}s | Rows: {len(df):,}")
        except Exception as e:
            print(f"❌ Error reading CSV: {str(e)}")
            return None
        
        # Validate column
        if compound_column not in df.columns:
            print(f"❌ Column '{compound_column}' not found!")
            print(f"Available columns: {list(df.columns)}")
            return None
        
        # Analyze data
        print("🔍 Analyzing dataset...")
        valid_compounds = df[compound_column].dropna()
        unique_compounds = valid_compounds.unique()
        unique_count = len(unique_compounds)
        
        print(f"📊 Dataset Analysis:")
        print(f"   • Total rows: {len(df):,}")
        print(f"   • Valid compound names: {len(valid_compounds):,}")
        print(f"   • Unique compounds: {unique_count:,}")
        
        if unique_count == 0:
            print("❌ No valid compounds found!")
            return df
        
        # Processing mode selection
        processing_mode = "Async" if self.use_async else "Threading"
        estimated_time = unique_count / (self.max_workers * 1.5)
        print(f"⚙️  Processing mode: {processing_mode}")
        print(f"🔧 Workers: {self.max_workers}")
        print(f"⏱️  Estimated time: {estimated_time/60:.1f} minutes")
        print()
        
        # Reset counters
        self.processed_count = 0
        self.success_count = 0
        
        # Process compounds
        print("🚀 Starting SMILES retrieval...")
        processing_start = time.time()
        
        try:
            if self.use_async:
                # Use asyncio for better performance
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)
                smiles_results = loop.run_until_complete(
                    self.async_process_compounds(unique_compounds.tolist())
                )
                loop.close()
            else:
                smiles_results = self.sync_process_compounds(unique_compounds.tolist())
        except Exception as e:
            print(f"❌ Error during processing: {str(e)}")
            print("🔄 Falling back to synchronous processing...")
            smiles_results = self.sync_process_compounds(unique_compounds.tolist())
        
        processing_time = time.time() - processing_start
        
        # Add results to dataframe
        print("\n📝 Adding SMILES to dataset...")
        df['CanonicalSMILES'] = df[compound_column].map(smiles_results)
        
        # Generate output filename
        if output_file_path is None:
            base_name = csv_file_path.replace('.csv', '')
            output_file_path = f"{base_name}_with_smiles_optimized.csv"
        
        # Save results
        print("💾 Saving results...")
        save_start = time.time()
        df.to_csv(output_file_path, index=False)
        save_time = time.time() - save_start
        
        # Enhanced statistics
        total_time = time.time() - start_time
        compounds_with_smiles = df['CanonicalSMILES'].notna().sum()
        success_percentage = (self.success_count / unique_count) * 100
        processing_rate = unique_count / processing_time
        cache_efficiency = len(self.cache) / unique_count * 100
        
        print("\n" + "="*70)
        print("🎉 PROCESSING COMPLETE!")
        print("="*70)
        print(f"📊 RESULTS:")
        print(f"   • Total rows processed: {len(df):,}")
        print(f"   • Unique compounds: {unique_count:,}")
        print(f"   • SMILES retrieved: {self.success_count:,} ({success_percentage:.1f}%)")
        print(f"   • Rows with SMILES: {compounds_with_smiles:,}")
        print(f"   • Cache hits: {len(self.cache):,}")
        print()
        print(f"⚡ PERFORMANCE:")
        print(f"   • Processing rate: {processing_rate:.1f} compounds/second")
        print(f"   • Processing time: {processing_time:.1f}s ({processing_time/60:.1f} min)")
        print(f"   • Total time: {total_time:.1f}s ({total_time/60:.1f} min)")
        print(f"   • Cache efficiency: {cache_efficiency:.1f}%")
        print()
        print(f"💾 OUTPUT: {output_file_path}")
        
        # Create summary report
        self.create_detailed_summary(df, compounds_with_smiles, unique_count)
        
        return df
    
    def create_detailed_summary(self, df: pd.DataFrame, compounds_with_smiles: int, unique_count: int):
        """Create detailed summary with analysis"""
        summary_path = "smiles_retrieval_summary.txt"
        
        # Analyze SMILES
        smiles_data = df[df['CanonicalSMILES'].notna()]['CanonicalSMILES']
        
        summary = f"""
OPTIMIZED SMILES RETRIEVAL SUMMARY
==================================

Dataset Statistics:
• Total rows: {len(df):,}
• Rows with SMILES: {compounds_with_smiles:,} ({100*compounds_with_smiles/len(df):.1f}%)
• Unique compounds processed: {unique_count:,}
• Success rate: {100*self.success_count/unique_count:.1f}%
• Cache efficiency: {100*len(self.cache)/unique_count:.1f}%

SMILES Analysis:
• Unique SMILES found: {smiles_data.nunique():,}
• Average SMILES length: {smiles_data.str.len().mean():.1f} characters
• SMILES length range: {smiles_data.str.len().min()}-{smiles_data.str.len().max()} characters

Top 10 Most Common Compounds:
{df['Phytochemical_name'].value_counts().head(10).to_string()}

Failed Retrievals (sample):
{df[df['CanonicalSMILES'].isna()]['Phytochemical_name'].value_counts().head(10).to_string()}

Processing Statistics:
• Cache hits: {len(self.cache):,}
• Network requests made: ~{self.processed_count:,}
• Average processing time per compound: {(time.time()/self.processed_count if self.processed_count > 0 else 0)*1000:.1f}ms
        """
        
        with open(summary_path, 'w') as f:
            f.write(summary.strip())
        
        print(f"📋 Detailed summary saved: {summary_path}")

# Convenient wrapper functions
def process_file_optimized(file_path: str, compound_column: str = 'Phytochemical_name', 
                          output_path: Optional[str] = None, max_workers: int = 50, 
                          use_async: bool = True) -> pd.DataFrame:
    """
    Optimized file processing with enhanced error handling
    
    Args:
        file_path: Path to CSV file
        compound_column: Column containing compound names
        output_path: Output file path (optional)
        max_workers: Number of parallel workers
        use_async: Use async processing for better performance
    """
    retriever = OptimizedSMILESRetriever(max_workers=max_workers, use_async=use_async)
    return retriever.process_dataset(file_path, compound_column, output_path)

# Example usage
if __name__ == "__main__":
    # Example with optimized settings
    file_path = "phytochemicals_haldi.csv"
    
    print("🚀 Running Optimized SMILES Retriever...")
    df = process_file_optimized(
        file_path=file_path,
        compound_column="Phytochemical_name",
        output_path="phytochemicals_with_smiles_optimized.csv",
        max_workers=50,  # Increased workers
        use_async=True   # Enable async processing
    )
    
    if df is not None:
        print(f"\n✅ Success! {df['CanonicalSMILES'].notna().sum():,} compounds now have SMILES!")
    else:
        print("❌ Processing failed. Please check the error messages above.")