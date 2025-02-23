#!/usr/bin/env python3
import sys
import re
import pysam

def get_csq_header_keys(vcf_header):
    """
    Parse the VCF header to extract the CSQ field format keys.
    """
    csq_line = None
    for record in vcf_header.records:
        if record.key == "INFO" and record.get("ID") == "CSQ":
            csq_line = record.get("Description")
            break
    if not csq_line:
        sys.stderr.write("ERROR: CSQ header not found in VCF.\n")
        sys.exit(1)
    # Extract the format specification from the description.
    m = re.search(r'Format:\s*(.+)', csq_line)
    if not m:
        sys.stderr.write("ERROR: Could not parse CSQ format from header.\n")
        sys.exit(1)
    keys = m.group(1).strip().split("|")
    return keys

def main():
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: python parse_annotations.py <annotated_vcf> <output_all_tsv> <output_protein_tsv>\n")
        sys.exit(1)
    
    annotated_vcf = sys.argv[1]
    output_all_tsv = sys.argv[2]
    output_protein_tsv = sys.argv[3]
    
    try:
        vcf = pysam.VariantFile(annotated_vcf, "r")
    except Exception as e:
        sys.stderr.write(f"ERROR: Could not open VCF file {annotated_vcf}: {e}\n")
        sys.exit(1)
    
    csq_keys = get_csq_header_keys(vcf.header)
    
    try:
        gene_idx = csq_keys.index("SYMBOL")
    except ValueError:
        sys.stderr.write("ERROR: CSQ header does not contain SYMBOL field. Re-run VEP with '--symbol'.\n")
        sys.exit(1)
        
    try:
        consequence_idx = csq_keys.index("Consequence")
    except ValueError:
        sys.stderr.write("ERROR: CSQ header does not contain Consequence field.\n")
        sys.exit(1)
        
    try:
        protein_pos_idx = csq_keys.index("Protein_position")
    except ValueError:
        sys.stderr.write("ERROR: CSQ header does not contain Protein_position field. Re-run VEP with '--protein' or '--hgvs'.\n")
        sys.exit(1)
        
    try:
        amino_acids_idx = csq_keys.index("Amino_acids")
    except ValueError:
        sys.stderr.write("ERROR: CSQ header does not contain Amino_acids field. Re-run VEP with '--protein' or '--hgvs'.\n")
        sys.exit(1)
    
    # Collect summary rows in a list
    rows = []
    
    for record in vcf.fetch():
        pos = record.pos
        ref = record.ref
        alts = record.alts
        if not alts:
            continue
        
        for alt in alts:
            if "CSQ" not in record.info:
                continue
            annotations = record.info["CSQ"]
            selected = None
            for ann in annotations:
                fields = ann.split("|")
                if len(fields) < len(csq_keys):
                    continue
                if fields[gene_idx]:
                    selected = fields
                    break
            if not selected:
                continue
            
            gene = selected[gene_idx]
            consequence = selected[consequence_idx]
            protein_position = selected[protein_pos_idx]
            amino_acids = selected[amino_acids_idx]
            
            # Use the first number if a range is given
            aa_position = protein_position.split("-")[0] if protein_position else "NA"
            if amino_acids and "/" in amino_acids:
                aa_ref, aa_alt = amino_acids.split("/")
            else:
                aa_ref, aa_alt = "NA", "NA"
            
            rows.append([gene, str(pos), ref, alt, consequence, aa_position, aa_ref, aa_alt])
    
    # Sort rows by gene name (first column)
    rows.sort(key=lambda x: x[0])
    
    # Filter for protein-changing annotations:
    # Only include rows where both aa_ref and aa_alt are provided and differ.
    protein_rows = [row for row in rows if row[6] != "NA" and row[7] != "NA" and row[6] != row[7]]
    
    header = "\t".join(["gene", "position", "reference", "alt", "mutation_category", "aa_position", "aa_ref", "aa_alt"])
    
    # Write all annotations to the first output file
    with open(output_all_tsv, "w") as allfh:
        allfh.write(header + "\n")
        for row in rows:
            allfh.write("\t".join(row) + "\n")
    
    # Write protein-changing annotations to the second output file
    with open(output_protein_tsv, "w") as proteinfh:
        proteinfh.write(header + "\n")
        for row in protein_rows:
            proteinfh.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()
