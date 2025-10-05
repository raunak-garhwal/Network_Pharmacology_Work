import pandas as pd
import requests
import time
import json
from urllib.parse import quote
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

class UltraFastSMILESRetriever:
    def __init__(self, max_workers=15):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
        self.max_workers = max_workers
        self.lock = threading.Lock()
        self.processed_count = 0
        self.success_count = 0
        
        # Setup session with connection pooling and retries
        self.session = requests.Session()
        
        retry_strategy = Retry(
            total=2,
            backoff_factor=0.3,
            status_forcelist=[429, 500, 502, 503, 504],
        )
        
        adapter = HTTPAdapter(
            pool_connections=20,
            pool_maxsize=20,
            max_retries=retry_strategy
        )
        
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)
        self.session.headers.update({'User-Agent': 'Ultra Fast SMILES Retriever'})
        
    def clean_compound_name(self, name):
        """Ultra-fast compound name cleaning"""
        if pd.isna(name) or not name:
            return None
        name = str(name).strip()
        if not name:
            return None
        # Minimal cleaning for speed
        name = re.sub(r'\s+', ' ', name)  # Multiple spaces to single
        return name.strip()
    
    def get_smiles_ultra_fast(self, compound_name):
        """Ultra-optimized single compound SMILES retrieval"""
        cleaned_name = self.clean_compound_name(compound_name)
        if not cleaned_name:
            return compound_name, None, 'invalid'
        
        try:
            encoded_name = quote(cleaned_name)
            url = f"{self.base_url}/name/{encoded_name}/property/CanonicalSMILES/json"
            
            response = self.session.get(url, timeout=8)
            
            if response.status_code == 200:
                data = response.json()
                properties = data.get('PropertyTable', {}).get('Properties', [])
                if properties:
                    smiles = properties[0].get('CanonicalSMILES')
                    if smiles:
                        with self.lock:
                            self.success_count += 1
                        return compound_name, smiles, 'success'
            
            return compound_name, None, 'not_found'
            
        except Exception:
            return compound_name, None, 'error'
    
    def update_progress(self, compound_name):
        """Thread-safe progress tracking"""
        with self.lock:
            self.processed_count += 1
            if self.processed_count % 50 == 0:  # Update every 50 compounds
                success_rate = (self.success_count / self.processed_count) * 100
                print(f"⚡ {self.processed_count} processed | ✓ {self.success_count} SMILES | {success_rate:.1f}% success")
    
    def process_compound_with_progress(self, compound_name):
        """Process compound and update progress"""
        result = self.get_smiles_ultra_fast(compound_name)
        self.update_progress(compound_name)
        return result
    
    def process_large_dataset(self, csv_file_path, compound_column='Phytochemical_name', output_file_path=None, chunk_save=True):
        """Optimized for large datasets with high compound counts"""
        
        print(f"🚀 ULTRA-FAST PROCESSING: {csv_file_path}")
        print("="*60)
        
        # Read CSV efficiently
        start_time = time.time()
        print("📖 Reading CSV...")
        df = pd.read_csv(csv_file_path)
        read_time = time.time() - start_time
        print(f"✅ CSV loaded in {read_time:.2f}s | Rows: {len(df):,}")
        
        if compound_column not in df.columns:
            available_cols = list(df.columns)
            print(f"❌ Column '{compound_column}' not found!")
            print(f"Available columns: {available_cols}")
            return None
        
        # Get unique compounds for processing
        print("🔍 Analyzing unique compounds...")
        unique_compounds = df[compound_column].dropna().unique()
        unique_count = len(unique_compounds)
        print(f"📊 Dataset: {len(df):,} total rows | {unique_count:,} unique compounds")
        
        if unique_count == 0:
            print("❌ No compounds found to process!")
            return df
        
        # Estimate processing time
        estimated_time = unique_count / (self.max_workers * 2)  # Conservative estimate
        print(f"⏱️  Estimated time: {estimated_time/60:.1f} minutes")
        print(f"🔧 Using {self.max_workers} parallel workers")
        print()
        
        # Reset counters
        self.processed_count = 0
        self.success_count = 0
        
        # Process compounds with ultra-high concurrency
        print("🚀 Starting SMILES retrieval...")
        processing_start = time.time()
        
        smiles_results = {}
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            future_to_compound = {
                executor.submit(self.process_compound_with_progress, compound): compound 
                for compound in unique_compounds
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_compound):
                compound = future_to_compound[future]
                try:
                    name, smiles, status = future.result()
                    smiles_results[name] = smiles
                except Exception as e:
                    smiles_results[compound] = None
        
        processing_time = time.time() - processing_start
        
        # Add SMILES column to dataframe
        print("\n📝 Adding SMILES to dataset...")
        df['CanonicalSMILES'] = df[compound_column].map(smiles_results)
        
        # Generate output filename if not provided
        if output_file_path is None:
            base_name = csv_file_path.replace('.csv', '')
            output_file_path = f"{base_name}_with_smiles.csv"
        
        # Save results
        print("💾 Saving results...")
        save_start = time.time()
        df.to_csv(output_file_path, index=False)
        save_time = time.time() - save_start
        
        # Final summary
        total_time = time.time() - start_time
        compounds_with_smiles = df['CanonicalSMILES'].notna().sum()
        success_percentage = (self.success_count / unique_count) * 100
        processing_rate = unique_count / processing_time
        
        print("\n" + "="*60)
        print("🎉 PROCESSING COMPLETE!")
        print("="*60)
        print(f"📊 RESULTS:")
        print(f"   • Total rows processed: {len(df):,}")
        print(f"   • Unique compounds: {unique_count:,}")
        print(f"   • SMILES retrieved: {self.success_count:,} ({success_percentage:.1f}%)")
        print(f"   • Rows with SMILES: {compounds_with_smiles:,}")
        print()
        print(f"⚡ PERFORMANCE:")
        print(f"   • Processing rate: {processing_rate:.1f} compounds/second")
        print(f"   • Processing time: {processing_time:.1f}s ({processing_time/60:.1f} min)")
        print(f"   • Total time: {total_time:.1f}s ({total_time/60:.1f} min)")
        print(f"   • Save time: {save_time:.2f}s")
        print()
        print(f"💾 OUTPUT: {output_file_path}")
        
        return df
    
    def create_smiles_summary(self, df, output_summary_path=None):
        """Create a summary of SMILES retrieval results"""
        if output_summary_path is None:
            output_summary_path = "smiles_summary.txt"
        
        total_rows = len(df)
        rows_with_smiles = df['CanonicalSMILES'].notna().sum()
        unique_compounds = df['Phytochemical_name'].nunique()
        unique_smiles = df['CanonicalSMILES'].dropna().nunique()
        
        summary = f"""
SMILES RETRIEVAL SUMMARY
========================
Dataset Statistics:
• Total rows: {total_rows:,}
• Rows with SMILES: {rows_with_smiles:,} ({100*rows_with_smiles/total_rows:.1f}%)
• Unique compounds: {unique_compounds:,}
• Unique SMILES: {unique_smiles:,}

Top 10 Most Common Compounds:
{df['Phytochemical_name'].value_counts().head(10).to_string()}

SMILES Length Distribution:
{df[df['CanonicalSMILES'].notna()]['CanonicalSMILES'].str.len().describe().to_string()}
        """
        
        with open(output_summary_path, 'w') as f:
            f.write(summary)
        
        print(f"📋 Summary saved: {output_summary_path}")
        return summary

# Ultra-simple usage function
def process_large_file(file_path, output_path=None, max_workers=15):
    """
    Ultra-fast processing for large files with many compounds
    
    Args:
        file_path: Path to your CSV file
        output_path: Output file path (optional)
        max_workers: Number of parallel threads (default: 15 for speed)
    """
    retriever = UltraFastSMILESRetriever(max_workers=max_workers)
    df = retriever.process_large_dataset(file_path, output_file_path=output_path)
    
    if df is not None:
        # Create summary
        retriever.create_smiles_summary(df)
        
        # Quick stats
        smiles_count = df['CanonicalSMILES'].notna().sum()
        total_count = len(df)
        print(f"\n🎯 QUICK STATS: {smiles_count:,} of {total_count:,} rows now have SMILES!")
        
    return df

# Example usage
if __name__ == "__main__":
    # Process your large file
    file_path = "phytochemicals_tulsi.csv"  # Replace with your actual file
    
    # Ultra-fast processing with maximum parallelization
    df = process_large_file(
        file_path=file_path,
        output_path="phytochemicals_with_smiles.csv",
        max_workers=50
    )
    
    print("\n✅ Done! Your file now has SMILES data.")