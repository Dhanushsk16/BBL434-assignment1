# Plasmid Designer Tool

This tool generates a plasmid sequence by extracting a replication origin (ORI) from an unknown organism and combining it with user-defined genetic modules.

## Methodology

1.  **ORI Identification**: 
    - Reads the input genome (`Input.fa` or `pUC19.fa`).
    - Uses **GC Skew Analysis** to find the replication start point.
    - Extracts 600bp around the minimum skew as the "newly found ORI".

2.  **Assembly (White-listing Logic)**:
    - Reads the `Design.txt` file.
    - Constructs a **brand new plasmid** by concatenating ONLY the parts listed in the design file.
    - **Concatenation Order**: `[ORI] → [Markers/Genes] → [MCS (Restriction Sites)] → [Default Genes]`
    - **Separator**: Uses the **BioBrick Scar** sequence (`TACTAGAG`) between each part.
    - **Default Genes**: Automatically appends all genes from `defaultgenes.json` at the end.
    - **EcoRI Deletion**: If `EcoRI` is present in the input but NOT in the design file, it is effectively deleted because it is never added to the new plasmid.

3.  **Handling Missing Data**:
    - If a marker requested in the design file is missing from the database (`markers.json`), the tool logs a warning and skips it without crashing.

## Usage

**Default Mode (Uses Input.fa and Design.txt):**
```bash
python plasmid_designer.py
```

**Fallback Mode:**
If `Input.fa` or `Design.txt` are not found, the tool automatically falls back to:
- `pUC19.fa` (for the input genome)
- `Design_pUC19.txt` (for the design file)

**Custom Files:**
```bash
python plasmid_designer.py --input my_genome.fa --design my_design.txt
```