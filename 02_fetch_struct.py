import pandas as pd
import requests
import os
import time
import argparse

def get_af_metadata(uniprot_id, session):
    """Queries the AFDB API for the CIF URL and the average pLDDT score."""
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = session.get(api_url, timeout=15)
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0:
                # The API returns a list; we take the first (usually latest) prediction
                prediction = data[0]
                return {
                    'cifUrl': prediction.get('cifUrl'),
                    'avgPlddt': prediction.get('avgPlddt')
                }
    except Exception as e:
        print(f"API Error for {uniprot_id}: {e}")
    return None

def download_cifs(csv_file, output_dir, plddt_threshold):
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")
        return

    col_name = 'UniProt_code'
    if col_name not in df.columns:
        print(f"Error: '{col_name}' not found. Available: {list(df.columns)}")
        return
    
    uniprot_ids = df[col_name].dropna().unique()
    os.makedirs(output_dir, exist_ok=True)
    
    session = requests.Session()
    session.headers.update({'User-Agent': 'Mozilla/5.0 (Bioinformatics Pipeline)'})

    stats = {"success": 0, "failed": 0, "skipped_low_plddt": 0, "exists": 0}

    print(f"Starting verified download for {len(uniprot_ids)} IDs (Threshold: pLDDT > {plddt_threshold})...")

    for up_id in uniprot_ids:
        clean_id = str(up_id).replace('*', '').strip()
        file_path = os.path.join(output_dir, f"{clean_id}.cif")
        
        if os.path.exists(file_path):
            stats["exists"] += 1
            continue

        # Step 1: Get metadata (URL + pLDDT)
        metadata = get_af_metadata(clean_id, session)
        
        if not metadata or not metadata['cifUrl']:
            print(f"[FAIL] No AFDB entry found: {clean_id}")
            stats["failed"] += 1
            continue

        # Step 2: Check pLDDT Score (The Quality Control Step)
        avg_plddt = metadata['avgPlddt']
        if avg_plddt < plddt_threshold:
            print(f"[LOW CONFIDENCE] {clean_id} skipped (pLDDT: {avg_plddt:.2f})")
            stats["skipped_low_plddt"] += 1
            continue

        # Step 3: Download from the verified URL
        try:
            r = session.get(metadata['cifUrl'], timeout=20)
            if r.status_code == 200:
                with open(file_path, 'wb') as f:
                    f.write(r.content)
                print(f"[OK] {clean_id} downloaded (pLDDT: {avg_plddt:.2f})")
                stats["success"] += 1
            else:
                print(f"[FAIL] Download failed for {clean_id} (Status {r.status_code})")
                stats["failed"] += 1
        except Exception as e:
            print(f"[ERROR] Connection error for {clean_id}: {e}")
            stats["failed"] += 1
        
        time.sleep(0.1) 

    print(f"\nFinal Summary:")
    print(f"Successfully Downloaded: {stats['success']}")
    print(f"Skipped (Low pLDDT):      {stats['skipped_low_plddt']}")
    print(f"Already Existed:         {stats['exists']}")
    print(f"Failed (API/Network):    {stats['failed']}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to kuvek_s2.csv")
    parser.add_argument("--out", default="plant_cyp_cifs")
    parser.add_argument("--plddt", type=float, default=70.0, help="Minimum average pLDDT score (default: 70)")
    
    args = parser.parse_args()
    download_cifs(args.input, args.out, args.plddt)