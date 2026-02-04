import pandas as pd
import requests
import os
import time
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def get_af_metadata(uniprot_id, session):
    """Queries the AFDB API for the CIF URL and the average pLDDT score."""
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = session.get(api_url, timeout=15)
        if response.status_code == 200:
            data = response.json()
            if data and len(data) > 0:
                prediction = data[0]
                return {
                    'cifUrl': prediction.get('cifUrl'),
                    'avgPlddt': prediction.get('avgPlddt')
                }
    except Exception as e:
        print(f"API Error for {uniprot_id}: {e}")
    return None

def download_and_plot(csv_file, output_dir):
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

    plddt_scores = []
    stats = {"success": 0, "failed": 0, "exists": 0}

    print(f"Starting structure retrieval for {len(uniprot_ids)} IDs...")

    for up_id in uniprot_ids:
        clean_id = str(up_id).replace('*', '').strip()
        file_path = os.path.join(output_dir, f"{clean_id}.cif")
        
        # Get metadata to capture pLDDT even if file exists locally
        metadata = get_af_metadata(clean_id, session)
        
        if metadata:
            plddt_scores.append({'ID': clean_id, 'pLDDT': metadata['avgPlddt']})
            
            if os.path.exists(file_path):
                stats["exists"] += 1
            else:
                try:
                    r = session.get(metadata['cifUrl'], timeout=20)
                    if r.status_code == 200:
                        with open(file_path, 'wb') as f:
                            f.write(r.content)
                        stats["success"] += 1
                    else:
                        stats["failed"] += 1
                except:
                    stats["failed"] += 1
        else:
            stats["failed"] += 1
        
        time.sleep(0.1)

    # Convert results to DataFrame for plotting
    results_df = pd.DataFrame(plddt_scores)
    
    if not results_df.empty:
        # Create the distribution plot
        plt.figure(figsize=(10, 6))
        sns.histplot(results_df['pLDDT'], kde=True, color='skyblue', bins=20)
        
        # Add Kuvek et al. threshold markers for visual reference
        plt.axvline(70, color='orange', linestyle='--', label='Confident (70)')
        plt.axvline(90, color='green', linestyle='--', label='Very High Confidence (90)')
        
        plt.title('Distribution of pLDDT Scores: 343 Plant CYPs')
        plt.xlabel('Average pLDDT')
        plt.ylabel('Frequency')
        plt.legend()
        plt.grid(axis='y', alpha=0.3)
        
        plot_path = "plddt_distribution.png"
        plt.savefig(plot_path)
        print(f"\nDistribution plot saved as: {plot_path}")
        plt.show()

    print(f"\nFinal Summary:")
    print(f"Downloaded: {stats['success']} | Existed: {stats['exists']} | Failed: {stats['failed']}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to kuvek_s2.csv")
    parser.add_argument("--out", default="plant_cyp_cifs")
    args = parser.parse_args()
    download_and_plot(args.input, args.out)