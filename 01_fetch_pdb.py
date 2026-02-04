import pandas as pd
import requests
import os
import time
import argparse

def get_af_metadata(uniprot_id, session):
    """Queries the AFDB API to get the correct CIF download URL."""
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = session.get(api_url, timeout=15)
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0:
                # Extract the cifUrl from the first prediction result
                return data[0].get('cifUrl')
    except Exception as e:
        print(f"API Error for {uniprot_id}: {e}")
    return None

def download_cifs(csv_file, output_dir):
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return

    # Use exact column name from your CSV
    col_name = 'UniProt_code'
    if col_name not in df.columns:
        print(f"Error: '{col_name}' not found. Available: {list(df.columns)}")
        return
    
    uniprot_ids = df[col_name].dropna().unique()
    os.makedirs(output_dir, exist_ok=True)
    
    session = requests.Session()
    session.headers.update({'User-Agent': 'Mozilla/5.0 (Bioinformatics Pipeline)'})

    stats = {"success": 0, "failed": 0, "skipped": 0}

    print(f"Starting API-verified download for {len(uniprot_ids)} IDs...")

    for up_id in uniprot_ids:
        clean_id = str(up_id).replace('*', '').strip()
        file_path = os.path.join(output_dir, f"{clean_id}.cif")
        
        if os.path.exists(file_path):
            stats["skipped"] += 1
            continue

        # Step 1: Get the correct URL from the API
        cif_url = get_af_metadata(clean_id, session)
        
        if not cif_url:
            print(f"[FAIL] No AFDB entry found for: {clean_id}")
            stats["failed"] += 1
            continue

        # Step 2: Download from the verified URL
        try:
            r = session.get(cif_url, timeout=20)
            if r.status_code == 200:
                with open(file_path, 'wb') as f:
                    f.write(r.content)
                print(f"[OK] Downloaded: {clean_id}")
                stats["success"] += 1
            else:
                print(f"[FAIL] Download failed for {clean_id} (Status {r.status_code})")
                stats["failed"] += 1
        except Exception as e:
            print(f"[ERROR] Connection error for {clean_id}: {e}")
            stats["failed"] += 1
        
        time.sleep(0.1) 

    print(f"\nFinished! Success: {stats['success']} | Failed: {stats['failed']} | Skipped: {stats['skipped']}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("--out", default="plant_cyp_cifs")
    args = parser.parse_args()
    download_cifs(args.input, args.out)