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
from typing import Optional, Tuple, Dict, List, Set
import logging
from difflib import SequenceMatcher
import warnings
warnings.filterwarnings('ignore')

# Enable nested event loops for Jupyter compatibility
try:
    nest_asyncio.apply()
except:
    pass

class ComprehensiveSMILESRetriever:
    def __init__(self, max_workers=50, use_async=True, enable_fuzzy_search=True):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
        self.max_workers = max_workers
        self.use_async = use_async
        self.enable_fuzzy_search = enable_fuzzy_search
        self.lock = threading.Lock()
        self.processed_count = 0
        self.success_count = 0
        self.cache = {}  # In-memory cache
        self.failed_compounds = set()  # Track permanently failed compounds
        
        # Setup detailed logging
        self.setup_logging()
        
        # Setup optimized session
        self.setup_session()
        
        # Search strategy tracking
        self.strategy_stats = defaultdict(int)
        
    def setup_logging(self):
        """Setup detailed logging for debugging"""
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)
        
        # Create file handler for detailed logs
        try:
            file_handler = logging.FileHandler('smiles_retrieval.log')
            file_handler.setLevel(logging.DEBUG)
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)
        except:
            pass  # File logging is optional
        
    def setup_session(self):
        """Setup optimized requests session"""
        self.session = requests.Session()
        
        retry_strategy = Retry(
            total=4,
            backoff_factor=0.5,
            status_forcelist=[429, 500, 502, 503, 504, 520, 521, 522, 523, 524],
            allowed_methods=["GET", "POST"]
        )
        
        adapter = HTTPAdapter(
            pool_connections=50,
            pool_maxsize=50,
            max_retries=retry_strategy
        )
        
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)
        self.session.headers.update({
            'User-Agent': 'Comprehensive SMILES Retriever v3.0 - Academic Research',
            'Accept': 'application/json',
            'Connection': 'keep-alive'
        })
    
    def clean_compound_name(self, name: str) -> List[str]:
        """Generate multiple cleaned variants of compound name"""
        if pd.isna(name) or not name:
            return []
            
        name = str(name).strip()
        if not name or len(name) < 2:
            return []
        
        variants = []
        
        # Original cleaned name
        cleaned = re.sub(r'\s+', ' ', name)
        cleaned = cleaned.strip()
        if cleaned:
            variants.append(cleaned)
        
        # Remove parenthetical information
        no_parens = re.sub(r'\s*\([^)]*\)', '', cleaned).strip()
        if no_parens and no_parens != cleaned:
            variants.append(no_parens)
            
        # Remove bracketed information
        no_brackets = re.sub(r'\s*\[[^\]]*\]', '', cleaned).strip()
        if no_brackets and no_brackets != cleaned:
            variants.append(no_brackets)
        
        # Handle common separators
        for sep in ['-', '_', '/', ',', ';']:
            if sep in cleaned:
                # Replace with space
                spaced = cleaned.replace(sep, ' ')
                spaced = re.sub(r'\s+', ' ', spaced).strip()
                if spaced not in variants:
                    variants.append(spaced)
                
                # Remove separator
                no_sep = cleaned.replace(sep, '')
                no_sep = re.sub(r'\s+', ' ', no_sep).strip()
                if no_sep not in variants:
                    variants.append(no_sep)
        
        # Handle Greek letters and special characters
        greek_replacements = {
            'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'δ': 'delta',
            'ε': 'epsilon', 'ζ': 'zeta', 'η': 'eta', 'θ': 'theta',
            'κ': 'kappa', 'λ': 'lambda', 'μ': 'mu', 'ν': 'nu',
            'ξ': 'xi', 'π': 'pi', 'ρ': 'rho', 'σ': 'sigma',
            'τ': 'tau', 'υ': 'upsilon', 'φ': 'phi', 'χ': 'chi',
            'ψ': 'psi', 'ω': 'omega'
        }
        
        greek_variant = cleaned
        for greek, english in greek_replacements.items():
            if greek in greek_variant:
                greek_variant = greek_variant.replace(greek, english)
        
        if greek_variant != cleaned and greek_variant not in variants:
            variants.append(greek_variant)
        
        # Remove duplicates while preserving order
        unique_variants = []
        seen = set()
        for variant in variants:
            if variant and variant.lower() not in seen:
                seen.add(variant.lower())
                unique_variants.append(variant)
        
        return unique_variants[:5]  # Limit to top 5 variants
    
    def extract_smiles_comprehensive(self, data: dict) -> Optional[str]:
        """Comprehensive SMILES extraction from various response formats"""
        try:
            # Strategy 1: Standard PropertyTable format
            if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                properties = data['PropertyTable']['Properties']
                if isinstance(properties, list) and properties:
                    prop = properties[0]
                    # Try multiple SMILES property names
                    for smiles_key in ['CanonicalSMILES', 'IsomericSMILES', 'SMILES']:
                        if smiles_key in prop:
                            smiles = prop[smiles_key]
                            if self.validate_smiles(smiles):
                                return smiles.strip()
            
            # Strategy 2: Direct Properties array
            if 'Properties' in data and isinstance(data['Properties'], list):
                for prop in data['Properties']:
                    for smiles_key in ['CanonicalSMILES', 'IsomericSMILES', 'SMILES']:
                        if smiles_key in prop:
                            smiles = prop[smiles_key]
                            if self.validate_smiles(smiles):
                                return smiles.strip()
            
            # Strategy 3: Nested structure search
            if isinstance(data, dict):
                for key, value in data.items():
                    if isinstance(value, dict):
                        smiles = self.extract_smiles_comprehensive(value)
                        if smiles:
                            return smiles
                    elif isinstance(value, list):
                        for item in value:
                            if isinstance(item, dict):
                                smiles = self.extract_smiles_comprehensive(item)
                                if smiles:
                                    return smiles
            
        except Exception as e:
            self.logger.debug(f"Error extracting SMILES: {str(e)}")
        
        return None
    
    def validate_smiles(self, smiles: str) -> bool:
        """Validate SMILES string"""
        if not smiles or not isinstance(smiles, str):
            return False
        
        smiles = smiles.strip()
        
        # Basic validation checks
        if len(smiles) < 3:
            return False
        
        # Check for common SMILES characters
        valid_chars = set('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789()[]@+-=#/\\.')
        if not all(c in valid_chars for c in smiles):
            return False
        
        # Check balanced parentheses and brackets
        parens = smiles.count('(') - smiles.count(')')
        brackets = smiles.count('[') - smiles.count(']')
        
        if parens != 0 or brackets != 0:
            return False
        
        return True
    
    async def search_strategy_1_name_direct(self, session: aiohttp.ClientSession, compound_name: str) -> Tuple[Optional[str], str]:
        """Strategy 1: Direct name-to-SMILES search"""
        variants = self.clean_compound_name(compound_name)
        
        for variant in variants:
            try:
                encoded_name = quote(variant)
                url = f"{self.base_url}/name/{encoded_name}/property/CanonicalSMILES/json"
                
                async with session.get(url, timeout=aiohttp.ClientTimeout(total=12)) as response:
                    if response.status == 200:
                        data = await response.json()
                        smiles = self.extract_smiles_comprehensive(data)
                        if smiles:
                            self.strategy_stats['name_direct'] += 1
                            return smiles, f'name_direct:{variant}'
                            
            except Exception as e:
                self.logger.debug(f"Strategy 1 failed for {variant}: {str(e)}")
                continue
        
        return None, 'strategy_1_failed'
    
    async def search_strategy_2_cid_based(self, session: aiohttp.ClientSession, compound_name: str) -> Tuple[Optional[str], str]:
        """Strategy 2: Name-to-CID-to-SMILES search"""
        variants = self.clean_compound_name(compound_name)
        
        for variant in variants:
            try:
                encoded_name = quote(variant)
                cid_url = f"{self.base_url}/name/{encoded_name}/cids/json"
                
                async with session.get(cid_url, timeout=aiohttp.ClientTimeout(total=10)) as response:
                    if response.status == 200:
                        data = await response.json()
                        cids = data.get('IdentifierList', {}).get('CID', [])
                        
                        if cids:
                            # Try multiple CIDs if available
                            for cid in cids[:3]:  # Try up to 3 CIDs
                                smiles_url = f"{self.base_url}/cid/{cid}/property/CanonicalSMILES/json"
                                async with session.get(smiles_url, timeout=aiohttp.ClientTimeout(total=10)) as smiles_response:
                                    if smiles_response.status == 200:
                                        smiles_data = await smiles_response.json()
                                        smiles = self.extract_smiles_comprehensive(smiles_data)
                                        if smiles:
                                            self.strategy_stats['cid_based'] += 1
                                            return smiles, f'cid_based:{cid}'
                                            
            except Exception as e:
                self.logger.debug(f"Strategy 2 failed for {variant}: {str(e)}")
                continue
        
        return None, 'strategy_2_failed'
    
    async def search_strategy_3_synonym_based(self, session: aiohttp.ClientSession, compound_name: str) -> Tuple[Optional[str], str]:
        """Strategy 3: Synonym-based search"""
        variants = self.clean_compound_name(compound_name)
        
        for variant in variants:
            try:
                encoded_name = quote(variant)
                synonym_url = f"{self.base_url}/name/{encoded_name}/synonyms/json"
                
                async with session.get(synonym_url, timeout=aiohttp.ClientTimeout(total=10)) as response:
                    if response.status == 200:
                        data = await response.json()
                        info_list = data.get('InformationList', {}).get('Information', [])
                        
                        if info_list and 'Synonym' in info_list[0]:
                            synonyms = info_list[0]['Synonym'][:5]  # Try first 5 synonyms
                            
                            for synonym in synonyms:
                                try:
                                    encoded_synonym = quote(synonym)
                                    smiles_url = f"{self.base_url}/name/{encoded_synonym}/property/CanonicalSMILES/json"
                                    
                                    async with session.get(smiles_url, timeout=aiohttp.ClientTimeout(total=8)) as smiles_response:
                                        if smiles_response.status == 200:
                                            smiles_data = await smiles_response.json()
                                            smiles = self.extract_smiles_comprehensive(smiles_data)
                                            if smiles:
                                                self.strategy_stats['synonym_based'] += 1
                                                return smiles, f'synonym_based:{synonym}'
                                except:
                                    continue
                                    
            except Exception as e:
                self.logger.debug(f"Strategy 3 failed for {variant}: {str(e)}")
                continue
        
        return None, 'strategy_3_failed'
    
    async def search_strategy_4_fuzzy_search(self, session: aiohttp.ClientSession, compound_name: str) -> Tuple[Optional[str], str]:
        """Strategy 4: Fuzzy/similarity search using PubChem's search API"""
        if not self.enable_fuzzy_search:
            return None, 'fuzzy_disabled'
            
        variants = self.clean_compound_name(compound_name)
        
        for variant in variants:
            try:
                # Use PubChem's search endpoint for fuzzy matching
                search_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/cids/json"
                encoded_name = quote(variant + "*")  # Wildcard search
                
                url = search_url.format(encoded_name)
                
                async with session.get(url, timeout=aiohttp.ClientTimeout(total=15)) as response:
                    if response.status == 200:
                        data = await response.json()
                        cids = data.get('IdentifierList', {}).get('CID', [])
                        
                        if cids:
                            # Try the first few CIDs from fuzzy search
                            for cid in cids[:2]:
                                smiles_url = f"{self.base_url}/cid/{cid}/property/CanonicalSMILES/json"
                                async with session.get(smiles_url, timeout=aiohttp.ClientTimeout(total=10)) as smiles_response:
                                    if smiles_response.status == 200:
                                        smiles_data = await smiles_response.json()
                                        smiles = self.extract_smiles_comprehensive(smiles_data)
                                        if smiles:
                                            self.strategy_stats['fuzzy_search'] += 1
                                            return smiles, f'fuzzy_search:{cid}'
                                            
            except Exception as e:
                self.logger.debug(f"Strategy 4 failed for {variant}: {str(e)}")
                continue
        
        return None, 'strategy_4_failed'
    
    async def search_strategy_5_post_request(self, session: aiohttp.ClientSession, compound_name: str) -> Tuple[Optional[str], str]:
        """Strategy 5: POST request for complex names with special characters"""
        variants = self.clean_compound_name(compound_name)
        
        for variant in variants:
            try:
                # Use POST request for names with special characters
                url = f"{self.base_url}/name/property/CanonicalSMILES/json"
                
                # Prepare POST data
                post_data = variant
                
                async with session.post(
                    url,
                    data=post_data,
                    headers={'Content-Type': 'text/plain'},
                    timeout=aiohttp.ClientTimeout(total=15)
                ) as response:
                    if response.status == 200:
                        data = await response.json()
                        smiles = self.extract_smiles_comprehensive(data)
                        if smiles:
                            self.strategy_stats['post_request'] += 1
                            return smiles, f'post_request:{variant}'
                            
            except Exception as e:
                self.logger.debug(f"Strategy 5 failed for {variant}: {str(e)}")
                continue
        
        return None, 'strategy_5_failed'
    
    async def comprehensive_smiles_search(self, session: aiohttp.ClientSession, compound_name: str) -> Tuple[str, Optional[str], str]:
        """Execute all search strategies for maximum success rate"""
        
        # Check cache first
        if compound_name in self.cache:
            return compound_name, self.cache[compound_name], 'cached'
        
        # Skip if permanently failed
        if compound_name in self.failed_compounds:
            return compound_name, None, 'permanently_failed'
        
        # Execute search strategies in order of efficiency
        strategies = [
            self.search_strategy_1_name_direct,
            self.search_strategy_2_cid_based,
            self.search_strategy_3_synonym_based,
            self.search_strategy_5_post_request,
            self.search_strategy_4_fuzzy_search  # Most expensive, try last
        ]
        
        for strategy in strategies:
            try:
                smiles, method = await strategy(session, compound_name)
                if smiles:
                    self.cache[compound_name] = smiles
                    return compound_name, smiles, method
            except Exception as e:
                self.logger.debug(f"Strategy failed for {compound_name}: {str(e)}")
                continue
        
        # Mark as permanently failed to avoid repeated attempts
        self.failed_compounds.add(compound_name)
        return compound_name, None, 'all_strategies_failed'
    
    def sync_comprehensive_search(self, compound_name: str) -> Tuple[str, Optional[str], str]:
        """Synchronous version of comprehensive search"""
        # Check cache first
        if compound_name in self.cache:
            return compound_name, self.cache[compound_name], 'cached'
        
        if compound_name in self.failed_compounds:
            return compound_name, None, 'permanently_failed'
        
        variants = self.clean_compound_name(compound_name)
        
        # Try each variant with different strategies
        for variant in variants:
            # Strategy 1: Direct name search
            try:
                encoded_name = quote(variant)
                url = f"{self.base_url}/name/{encoded_name}/property/CanonicalSMILES/json"
                
                response = self.session.get(url, timeout=12)
                if response.status_code == 200:
                    data = response.json()
                    smiles = self.extract_smiles_comprehensive(data)
                    if smiles:
                        self.cache[compound_name] = smiles
                        self.strategy_stats['sync_name_direct'] += 1
                        return compound_name, smiles, f'sync_name_direct:{variant}'
            except:
                pass
            
            # Strategy 2: CID-based search
            try:
                cid_url = f"{self.base_url}/name/{encoded_name}/cids/json"
                response = self.session.get(cid_url, timeout=10)
                
                if response.status_code == 200:
                    data = response.json()
                    cids = data.get('IdentifierList', {}).get('CID', [])
                    
                    for cid in cids[:2]:
                        smiles_url = f"{self.base_url}/cid/{cid}/property/CanonicalSMILES/json"
                        smiles_response = self.session.get(smiles_url, timeout=10)
                        
                        if smiles_response.status_code == 200:
                            smiles_data = smiles_response.json()
                            smiles = self.extract_smiles_comprehensive(smiles_data)
                            if smiles:
                                self.cache[compound_name] = smiles
                                self.strategy_stats['sync_cid_based'] += 1
                                return compound_name, smiles, f'sync_cid_based:{cid}'
            except:
                pass
        
        self.failed_compounds.add(compound_name)
        return compound_name, None, 'sync_all_failed'
    
    def update_progress(self, compound_name: str, success: bool = False):
        """Enhanced progress tracking"""
        with self.lock:
            self.processed_count += 1
            if success:
                self.success_count += 1
                
            if self.processed_count % 20 == 0:
                success_rate = (self.success_count / self.processed_count) * 100
                cache_hits = len(self.cache)
                failed_count = len(self.failed_compounds)
                
                print(f"⚡ {self.processed_count} processed | ✓ {self.success_count} SMILES ({success_rate:.1f}%) | "
                      f"📋 {cache_hits} cached | ❌ {failed_count} failed")
    
    async def async_process_compounds(self, compounds: List[str]) -> Dict[str, Optional[str]]:
        """Async processing with comprehensive search strategies"""
        results = {}
        
        # Enhanced async session configuration
        connector = aiohttp.TCPConnector(
            limit=60,
            limit_per_host=30,
            ttl_dns_cache=600,
            use_dns_cache=True,
            keepalive_timeout=60,
            enable_cleanup_closed=True
        )
        
        timeout = aiohttp.ClientTimeout(total=20, connect=5)
        
        async with aiohttp.ClientSession(
            connector=connector,
            timeout=timeout,
            headers={
                'User-Agent': 'Comprehensive Async SMILES Retriever v3.0',
                'Accept': 'application/json',
                'Connection': 'keep-alive'
            }
        ) as session:
            
            semaphore = asyncio.Semaphore(self.max_workers)
            
            async def process_with_semaphore(compound):
                async with semaphore:
                    result = await self.comprehensive_smiles_search(session, compound)
                    success = result[1] is not None
                    self.update_progress(compound, success)
                    return result
            
            tasks = [process_with_semaphore(compound) for compound in compounds]
            
            for coro in asyncio.as_completed(tasks):
                name, smiles, status = await coro
                results[name] = smiles
        
        return results
    
    def sync_process_compounds(self, compounds: List[str]) -> Dict[str, Optional[str]]:
        """Enhanced synchronous processing"""
        results = {}
        
        def process_compound(compound):
            result = self.sync_comprehensive_search(compound)
            success = result[1] is not None
            self.update_progress(compound, success)
            return result
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_compound = {
                executor.submit(process_compound, compound): compound 
                for compound in compounds
            }
            
            for future in as_completed(future_to_compound):
                try:
                    name, smiles, status = future.result(timeout=30)
                    results[name] = smiles
                except Exception as e:
                    compound = future_to_compound[future]
                    self.logger.error(f"Error processing {compound}: {str(e)}")
                    results[compound] = None
        
        return results
    
    def process_dataset(self, csv_file_path: str, compound_column: str = 'Phytochemical_name', 
                       output_file_path: Optional[str] = None) -> pd.DataFrame:
        """Main processing function with comprehensive search strategies"""
        
        print(f"🚀 COMPREHENSIVE SMILES RETRIEVER v3.0")
        print(f"📁 Processing: {csv_file_path}")
        print(f"🔍 Search strategies: Name Direct, CID-based, Synonym-based, POST Request, Fuzzy Search")
        print("="*80)
        
        # Read and validate CSV
        try:
            start_time = time.time()
            print("📖 Reading CSV...")
            
            # Try different encodings
            for encoding in ['utf-8', 'latin1', 'cp1252']:
                try:
                    df = pd.read_csv(csv_file_path, encoding=encoding)
                    break
                except UnicodeDecodeError:
                    continue
            else:
                raise Exception("Could not read CSV with any encoding")
                
            read_time = time.time() - start_time
            print(f"✅ CSV loaded in {read_time:.2f}s | Rows: {len(df):,} | Encoding: {encoding}")
            
        except Exception as e:
            print(f"❌ Error reading CSV: {str(e)}")
            return None
        
        if compound_column not in df.columns:
            print(f"❌ Column '{compound_column}' not found!")
            print(f"Available columns: {list(df.columns)}")
            return None
        
        # Enhanced data analysis
        print("🔍 Analyzing dataset...")
        valid_compounds = df[compound_column].dropna()
        unique_compounds = valid_compounds.drop_duplicates()
        unique_count = len(unique_compounds)
        
        # Analyze compound name complexity
        avg_length = valid_compounds.str.len().mean()
        special_chars = valid_compounds.str.contains(r'[^\w\s-]', regex=True).sum()
        
        print(f"📊 Dataset Analysis:")
        print(f"   • Total rows: {len(df):,}")
        print(f"   • Valid compound names: {len(valid_compounds):,}")
        print(f"   • Unique compounds: {unique_count:,}")
        print(f"   • Average name length: {avg_length:.1f} characters")
        print(f"   • Names with special characters: {special_chars:,} ({100*special_chars/len(valid_compounds):.1f}%)")
        
        if unique_count == 0:
            print("❌ No valid compounds found!")
            return df
        
        # Processing configuration
        processing_mode = "Async Multi-Strategy" if self.use_async else "Sync Multi-Strategy"
        estimated_time = unique_count / (self.max_workers * 1.2)  # More conservative estimate
        
        print(f"⚙️  Processing mode: {processing_mode}")
        print(f"🔧 Workers: {self.max_workers}")
        print(f"🎯 Fuzzy search: {'Enabled' if self.enable_fuzzy_search else 'Disabled'}")
        print(f"⏱️  Estimated time: {estimated_time/60:.1f} minutes")
        print()
        
        # Reset counters
        self.processed_count = 0
        self.success_count = 0
        self.cache.clear()
        self.failed_compounds.clear()
        self.strategy_stats.clear()
        
        # Process compounds
        print("🚀 Starting comprehensive SMILES retrieval...")
        print("   Using multiple search strategies for maximum success rate...")
        processing_start = time.time()
        
        try:
            if self.use_async:
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)
                smiles_results = loop.run_until_complete(
                    self.async_process_compounds(unique_compounds.tolist())
                )
                loop.close()
            else:
                smiles_results = self.sync_process_compounds(unique_compounds.tolist())
                
        except Exception as e:
            print(f"⚠️  Error in primary processing: {str(e)}")
            print("🔄 Falling back to synchronous processing...")
            smiles_results = self.sync_process_compounds(unique_compounds.tolist())
        
        processing_time = time.time() - processing_start
        
        # Add results to dataframe
        print(f"\n📝 Mapping SMILES to {len(df):,} rows...")
        df['CanonicalSMILES'] = df[compound_column].map(smiles_results)
        
        # Generate output filename
        if output_file_path is None:
            base_name = csv_file_path.replace('.csv', '')
            output_file_path = f"{base_name}_comprehensive_smiles.csv"
        
        # Save results
        print("💾 Saving enhanced dataset...")
        save_start = time.time()
        df.to_csv(output_file_path, index=False)
        save_time = time.time() - save_start
        
        # Comprehensive statistics
        total_time = time.time() - start_time
        compounds_with_smiles = df['CanonicalSMILES'].notna().sum()
        success_percentage = (self.success_count / unique_count) * 100
        processing_rate = unique_count / processing_time
        
        print("\n" + "="*80)
        print("🎉 COMPREHENSIVE PROCESSING COMPLETE!")
        print("="*80)
        print(f"📊 FINAL RESULTS:")
        print(f"   • Total rows processed: {len(df):,}")
        print(f"   • Unique compounds: {unique_count:,}")
        print(f"   • SMILES retrieved: {self.success_count:,} ({success_percentage:.1f}%)")
        print(f"   • Rows with SMILES: {compounds_with_smiles:,}")
        print(f"   • Cache efficiency: {len(self.cache):,} compounds cached")
        print(f"   • Failed compounds: {len(self.failed_compounds):,}")
        print()
        print(f"🎯 STRATEGY PERFORMANCE:")
        for strategy, count in self.strategy_stats.items():
            percentage = (count / self.success_count * 100) if self.success_count > 0 else 0
            print(f"   • {strategy.replace('_', ' ').title()}: {count:,} ({percentage:.1f}%)")
        print()
        print(f"⚡ PERFORMANCE METRICS:")
        print(f"   • Processing rate: {processing_rate:.1f} compounds/second")
        print(f"   • Processing time: {processing_time:.1f}s ({processing_time/60:.1f} min)")
        print(f"   • Total time: {total_time:.1f}s ({total_time/60:.1f} min)")
        print(f"   • Average time per compound: {processing_time/unique_count:.2f}s")
        print()
        print(f"💾 OUTPUT: {output_file_path}")
        
        # Create comprehensive analysis report
        self.create_comprehensive_report(df, compounds_with_smiles, unique_count, processing_time)
        
        return df
    
    def create_comprehensive_report(self, df: pd.DataFrame, compounds_with_smiles: int, 
                                  unique_count: int, processing_time: float):
        """Create detailed analysis report"""
        report_path = "comprehensive_smiles_report.txt"
        
        # Analyze SMILES data
        smiles_data = df[df['CanonicalSMILES'].notna()]['CanonicalSMILES']
        failed_compounds = df[df['CanonicalSMILES'].isna()]['Phytochemical_name'].dropna()
        
        # SMILES complexity analysis
        if not smiles_data.empty:
            smiles_lengths = smiles_data.str.len()
            avg_length = smiles_lengths.mean()
            complexity_score = smiles_data.str.count(r'[()[\]@#+\-=]').mean()
        else:
            avg_length = 0
            complexity_score = 0
        
        # Generate report
        report = f"""
COMPREHENSIVE SMILES RETRIEVAL REPORT
====================================

PROCESSING SUMMARY:
• Date/Time: {time.strftime('%Y-%m-%d %H:%M:%S')}
• Processing Time: {processing_time/60:.1f} minutes
• Success Rate: {100*self.success_count/unique_count:.1f}%

DATASET STATISTICS:
• Total rows: {len(df):,}
• Rows with SMILES: {compounds_with_smiles:,} ({100*compounds_with_smiles/len(df):.1f}%)
• Unique compounds processed: {unique_count:,}
• Successful retrievals: {self.success_count:,}
• Failed retrievals: {len(self.failed_compounds):,}

SEARCH STRATEGY EFFECTIVENESS:
{chr(10).join([f'• {strategy.replace("_", " ").title()}: {count:,} compounds ({100*count/self.success_count:.1f}%)' 
               for strategy, count in sorted(self.strategy_stats.items(), key=lambda x: x[1], reverse=True)])}

SMILES QUALITY ANALYSIS:
• Unique SMILES found: {smiles_data.nunique() if not smiles_data.empty else 0:,}
• Average SMILES length: {avg_length:.1f} characters
• SMILES complexity score: {complexity_score:.2f}
• Length distribution:
  - Min: {smiles_lengths.min() if not smiles_data.empty else 0} characters
  - Max: {smiles_lengths.max() if not smiles_data.empty else 0} characters
  - Median: {smiles_lengths.median() if not smiles_data.empty else 0:.0f} characters

TOP 10 MOST COMMON COMPOUNDS:
{df['Phytochemical_name'].value_counts().head(10).to_string()}

TOP 10 FAILED COMPOUNDS (for manual review):
{failed_compounds.value_counts().head(10).to_string()}

PERFORMANCE METRICS:
• Processing rate: {unique_count/processing_time:.1f} compounds/second
• Cache efficiency: {100*len(self.cache)/unique_count:.1f}%
• Network requests saved: ~{len(self.cache):,}
• Average response time: {processing_time/unique_count:.3f}s per compound

RECOMMENDATIONS:
1. Failed compounds may need manual curation or alternative databases
2. Consider cross-validation with ChEMBL or ChEBI for failed compounds
3. Review compounds with very short (<10 chars) or very long (>200 chars) SMILES
4. Cache can be saved for future runs to improve performance

TECHNICAL DETAILS:
• Max workers: {self.max_workers}
• Async processing: {self.use_async}
• Fuzzy search: {self.enable_fuzzy_search}
• Total API calls made: ~{self.processed_count * 2:,} (estimated)
        """
        
        try:
            with open(report_path, 'w', encoding='utf-8') as f:
                f.write(report.strip())
            print(f"📋 Comprehensive report saved: {report_path}")
        except Exception as e:
            print(f"⚠️  Could not save report: {str(e)}")
        
        return report

# Enhanced convenience functions
def process_file_comprehensive(file_path: str, compound_column: str = 'Phytochemical_name', 
                             output_path: Optional[str] = None, max_workers: int = 50, 
                             use_async: bool = True, enable_fuzzy_search: bool = True) -> pd.DataFrame:
    """
    Comprehensive file processing with multiple search strategies
    
    Args:
        file_path: Path to CSV file
        compound_column: Column containing compound names
        output_path: Output file path (optional)
        max_workers: Number of parallel workers
        use_async: Use async processing for better performance
        enable_fuzzy_search: Enable fuzzy/similarity search (slower but more comprehensive)
    
    Returns:
        Enhanced DataFrame with SMILES data
    """
    retriever = ComprehensiveSMILESRetriever(
        max_workers=max_workers, 
        use_async=use_async,
        enable_fuzzy_search=enable_fuzzy_search
    )
    return retriever.process_dataset(file_path, compound_column, output_path)

def quick_smiles_lookup(compound_names: List[str], max_workers: int = 15) -> Dict[str, Optional[str]]:
    """
    Quick SMILES lookup for a list of compounds
    
    Args:
        compound_names: List of compound names
        max_workers: Number of parallel workers
    
    Returns:
        Dictionary mapping compound names to SMILES
    """
    retriever = ComprehensiveSMILESRetriever(max_workers=max_workers, use_async=True)
    
    print(f"🔍 Looking up SMILES for {len(compound_names)} compounds...")
    
    try:
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        results = loop.run_until_complete(retriever.async_process_compounds(compound_names))
        loop.close()
    except:
        results = retriever.sync_process_compounds(compound_names)
    
    success_count = sum(1 for v in results.values() if v is not None)
    print(f"✅ Found SMILES for {success_count}/{len(compound_names)} compounds ({100*success_count/len(compound_names):.1f}%)")
    
    return results

def batch_process_files(file_paths: List[str], compound_column: str = 'Phytochemical_name',
                       max_workers: int = 50) -> Dict[str, pd.DataFrame]:
    """
    Process multiple CSV files in batch
    
    Args:
        file_paths: List of CSV file paths
        compound_column: Column containing compound names
        max_workers: Number of parallel workers
    
    Returns:
        Dictionary mapping file paths to processed DataFrames
    """
    results = {}
    
    print(f"🚀 Batch processing {len(file_paths)} files...")
    
    for i, file_path in enumerate(file_paths, 1):
        print(f"\n📁 Processing file {i}/{len(file_paths)}: {file_path}")
        
        try:
            df = process_file_comprehensive(
                file_path=file_path,
                compound_column=compound_column,
                max_workers=max_workers,
                use_async=True,
                enable_fuzzy_search=True
            )
            results[file_path] = df
            print(f"✅ File {i} completed successfully")
            
        except Exception as e:
            print(f"❌ Error processing {file_path}: {str(e)}")
            results[file_path] = None
    
    successful = sum(1 for df in results.values() if df is not None)
    print(f"\n🎉 Batch processing complete: {successful}/{len(file_paths)} files processed successfully")
    
    return results

# Example usage and testing
if __name__ == "__main__":
    # Example 1: Comprehensive processing
    print("🧪 COMPREHENSIVE SMILES RETRIEVER - DEMO")
    print("="*50)
    
    # Test with a sample file
    file_path = "phytochemicals_haldi.csv"  # Replace with your file
    
    print("🚀 Running comprehensive SMILES retrieval...")
    df = process_file_comprehensive(
        file_path=file_path,
        compound_column="Phytochemical_name",
        output_path="phytochemicals_comprehensive_smiles.csv",
        max_workers=10,
        use_async=True,
        enable_fuzzy_search=True
    )
    
    if df is not None:
        success_count = df['CanonicalSMILES'].notna().sum()
        total_count = len(df)
        print(f"\n🎯 FINAL RESULT: {success_count:,} of {total_count:,} compounds now have SMILES!")
        print(f"   Success rate: {100*success_count/total_count:.1f}%")
        
        # Show sample of successful retrievals
        successful_samples = df[df['CanonicalSMILES'].notna()][['Phytochemical_name', 'CanonicalSMILES']].head()
        print(f"\n📋 Sample successful retrievals:")
        for _, row in successful_samples.iterrows():
            name = row['Phytochemical_name']
            smiles = row['CanonicalSMILES']
            print(f"   • {name}: {smiles[:50]}{'...' if len(smiles) > 50 else ''}")
    
    # Example 2: Quick lookup for specific compounds
    print(f"\n🔍 Testing quick lookup...")
    test_compounds = ["curcumin", "quercetin", "resveratrol", "caffeine", "aspirin"]
    quick_results = quick_smiles_lookup(test_compounds, max_workers=10)
    
    print("📋 Quick lookup results:")
    for compound, smiles in quick_results.items():
        status = "✅" if smiles else "❌"
        print(f"   {status} {compound}: {smiles or 'Not found'}")
    
    print(f"\n✅ Demo completed! Check output files for detailed results.")