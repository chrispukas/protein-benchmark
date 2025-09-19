import protein_benchmark.io.pdb as pdb



def fasta_from_pdb(pdb_file: str,
                   output_file: str,
                   chain_map: dict[str, str],
                   chain_filter: set[str]):
    records, _ = pdb.reformat_pdb(pdb.extract_lines(pdb_file),
                                        chain_map=chain_map,
                                        chain_filter=chain_filter)
    chains: dict[str, str] = pdb.get_sequences(records)
    with open(output_file, "w") as file_handle:
        for chain_id, sequence in chains.items():
                write_record(chain_id, sequence, file_handle)
    print("Processed sequence for:", pdb_file)


def write_record(chain_id: str, sequence: str, file_handle) -> None:
    file_handle.write(f">{chain_id}|protein|\n")
    file_handle.write(sequence + "\n")