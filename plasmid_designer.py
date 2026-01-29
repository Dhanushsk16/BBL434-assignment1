import argparse
import json
import os

# --- DEFAULTS ---
DEFAULT_INPUT = "Input.fa"
DEFAULT_DESIGN = "Design.txt"
FALLBACK_INPUT = "pUC19.fa"
FALLBACK_DESIGN = "Design_pUC19.txt"
MARKER_DB_FILE = "markers.json"
DEFAULT_GENES_FILE = "defaultgenes.json"
OUTPUT_FILE = "Output.fa"
BIOBRICK_SCAR = "TACTAGAG"  # Standard spacer
ORI_KEY_IN_DESIGN = "ori_pMB1"

def read_fasta(filename):
    """Reads FASTA file, returns raw sequence string."""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        return "".join([line.strip() for line in lines if not line.startswith(">")]).upper()
    except FileNotFoundError:
        print(f"Error: The input file '{filename}' was not found.")
        return None

def get_gc_skew_ori(sequence, window_size=600):
    """Finds the index of minimum GC skew and extracts a window around it."""
    skew = 0
    min_skew = 0
    min_idx = 0
    
    for i, base in enumerate(sequence):
        if base == 'G':
            skew += 1
        elif base == 'C':
            skew -= 1
        
        if skew < min_skew:
            min_skew = skew
            min_idx = i

    length = len(sequence)
    start = min_idx - (window_size // 2)
    end = min_idx + (window_size // 2)
    
    if start < 0:
        ori_seq = sequence[start:] + sequence[:end]
    elif end > length:
        ori_seq = sequence[start:] + sequence[:end - length]
    else:
        ori_seq = sequence[start:end]
        
    return ori_seq, min_idx

def load_json_db(filename, db_name="database"):
    """Generic JSON loader for markers and default genes."""
    try:
        with open(filename, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Warning: {db_name} '{filename}' not found.")
        return {}

def parse_design_file(filename):
    keys = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split(',')
                    keys.append(parts[0].strip())
    except FileNotFoundError:
        print(f"Error: Design file '{filename}' not found.")
        return []
    return keys

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plasmid Designer Tool")
    parser.add_argument("--input", default=DEFAULT_INPUT, help="Path to input FASTA")
    parser.add_argument("--design", default=DEFAULT_DESIGN, help="Path to design file")
    args = parser.parse_args()

    # Fallback logic: if default files don't exist, try the backup files
    input_file = args.input
    design_file = args.design

    if not os.path.exists(input_file) and input_file == DEFAULT_INPUT:
        if os.path.exists(FALLBACK_INPUT):
            print(f"'{DEFAULT_INPUT}' not found. Using fallback: '{FALLBACK_INPUT}'")
            input_file = FALLBACK_INPUT
        else:
            print(f"Error: Neither '{DEFAULT_INPUT}' nor '{FALLBACK_INPUT}' found.")
            exit()

    if not os.path.exists(design_file) and design_file == DEFAULT_DESIGN:
        if os.path.exists(FALLBACK_DESIGN):
            print(f"'{DEFAULT_DESIGN}' not found. Using fallback: '{FALLBACK_DESIGN}'")
            design_file = FALLBACK_DESIGN
        else:
            print(f"Error: Neither '{DEFAULT_DESIGN}' nor '{FALLBACK_DESIGN}' found.")
            exit()

    print(f"--- Plasmid Designer Tool ---")

    # 1. Read Input
    host_genome = read_fasta(input_file)
    if not host_genome: exit()

    # 2. Find ORI
    print("1. Identifying Origin of Replication (ORI)...")
    found_ori_seq, ori_idx = get_gc_skew_ori(host_genome)
    print(f"   -> ORI identified at index {ori_idx}. Extracted {len(found_ori_seq)}bp.")

    # 3. Load Data
    markers_db = load_json_db(MARKER_DB_FILE, "Marker database")
    default_genes_db = load_json_db(DEFAULT_GENES_FILE, "Default genes database")
    design_keys = parse_design_file(design_file)
    if not markers_db or not design_keys: exit()

    # 4. SORT PARTS into [ORI] + [Markers] + [MCS] + [Default Genes]
    print("\n2. Assembling Plasmid (Enforcing Order: ORI -> Markers -> MCS -> Default Genes)...")
    
    # We create lists to hold the sequences for each category
    ori_part = ""
    marker_parts = []
    mcs_parts = []
    
    for key in design_keys:
        part_seq = ""
        
        # A. Retrieve Sequence
        if key == ORI_KEY_IN_DESIGN:
            part_seq = found_ori_seq
            print(f"   + Found ORI ({key}) -> Scheduled for Position 1")
            ori_part = part_seq # Store directly
            continue # Done with this key
            
        elif key in markers_db:
            part_seq = markers_db[key]
        else:
            print(f"   [WARNING] '{key}' missing in markers.json. SKIPPING.")
            continue

        # B. Categorize (Marker vs MCS)
        if key.endswith("_site"):
            # It is a Restriction Site -> Goes to MCS (End)
            mcs_parts.append(part_seq)
            print(f"   + {key} -> Scheduled for MCS (Position 3)")
        else:
            # It is likely a Gene/Marker -> Goes to Middle
            marker_parts.append(part_seq)
            print(f"   + {key} -> Scheduled for Markers (Position 2)")

    # 5. CONCATENATE
    # Order: ORI -> SCAR -> MARKERS -> SCAR -> MCS
    
    final_plasmid = ori_part
    
    # Append Markers
    for seq in marker_parts:
        if final_plasmid: final_plasmid += BIOBRICK_SCAR
        final_plasmid += seq
        
    # Append MCS
    for seq in mcs_parts:
        if final_plasmid: final_plasmid += BIOBRICK_SCAR
        final_plasmid += seq

    # Append Default Genes
    if default_genes_db:
        print("\n3. Adding Default Genes...")
        for gene_name, gene_seq in default_genes_db.items():
            if final_plasmid: final_plasmid += BIOBRICK_SCAR
            final_plasmid += gene_seq
            print(f"   + {gene_name} -> Added (Position 4: Default Genes)")

    # 6. Output
    if final_plasmid:
        with open(OUTPUT_FILE, 'w') as f:
            f.write(f">Synthetic_Plasmid_Output\n")
            for i in range(0, len(final_plasmid), 80):
                f.write(final_plasmid[i:i+80] + "\n")
        print(f"\n[SUCCESS] Plasmid saved to: {OUTPUT_FILE}")
        
        if "GAATTC" not in final_plasmid:
            print("[TEST CHECK] Verified: EcoRI site (GAATTC) is ABSENT.")
        else:
            print("[TEST CHECK] Warning: EcoRI site is present.")
    else:
        print("\n[FAILURE] Plasmid generation failed.")