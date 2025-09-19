# Protein Benchmarking Pipeline

General utility benchmarking pipeline for protein prediction ML models (in current configuration, is used for Boltz-2).

To setup, run the following commands:
```
git clone https://github.com/chrispukas/protein-benchmark.git;
cd ./protein-benchmark;
bash setup.sh
```

To run this project, the key notebooks are: ```./scripts/notebooks/main/notebook_input-generation.ipynb``` for input preperation (outputs YAML files for Boltz-2, but can be augmented for other models) and ```./scripts/notebooks/main/notebook_input-generation.ipynb``` for output processing (converting cif to pdb for dockQ scoring etc)



Possible benchmarking configurations include:
| Configuration name                  | Description                                                       |
|-------------------------------------|-------------------------------------------------------------------|
| `epitope_cutoff`                    | The minimum distance of heavy atoms in the antigen and antibody to be considered an eptiope                                 |
| `msa`                               | Path to the MSA (Multisequence Alignment) file, if None, then will use msa server, if 'empty', no msa will be used at all                                |
| `enforce_templates`                 | Enforce fixed steering on the template                              |
| `template_constraint_threshold`     | Max-distance deviation from the template                                  |
| `template_paths`                    | Path to the template pdb files                                 |
| `chain_ids`                         | Chain mapping in the actual sequence e.g [A, H, L]                                |
| `template_ids`                      | Chain mapping in the template e.g [A1, H1, L1]                                 |
| `enforce_pocket`                    | Enforce fixed steering on the pocket (epitope constraints)                                 |
| `pocket_constraint_threshold`       | Max-distance devation from the pocket                                 |
| `pocket_percentage`                 | The percentace of pockets that will be passed into the model inputs                                  |



<br><br>

# Examples:

Running this benchmark on Boltz-2 with the following configurations, will result in the following data:

## Templates, with no epitopes

| Configuration name                  | Setting                                                       |
|-------------------------------------|-------------------------------------------------------------------|
| `epitope_cutoff`                    | 4.50 A                                 |
| `msa`                               | None                                |
| `enforce_templates`                 | True                              |
| `template_constraint_threshold`     | 0.01                                  |
| `template_paths`                    | [<path/to/template(s)>, ...]                                |
| `chain_ids`                         | [A, H, L]                                |
| `template_ids`                      | [A1, H1, L1]                                 |
| `enforce_pocket`                    | False                                |
| `pocket_constraint_threshold`       | n/a                                 |
| `pocket_percentage`                 | n/a                                  |


![Alt text](/res/img/graphs/Screenshot%202025-09-19%20at%2013.49.45.png)


## Restricted Epitopes, with no templates


| Configuration name                  | Setting                                                       |
|-------------------------------------|-------------------------------------------------------------------|
| `epitope_cutoff`                    | 4.50 A                                 |
| `msa`                               | None                                |
| `enforce_templates`                 | False                              |
| `template_constraint_threshold`     | n/a                                  |
| `template_paths`                    | n/a                                |
| `chain_ids`                         | n/a                               |
| `template_ids`                      | n/a                                 |
| `enforce_pocket`                    | True                                |
| `pocket_constraint_threshold`       | 0.01 A                                 |
| `pocket_percentage`                 | 1.00                                  |


![Alt text](/res/img/graphs/Screenshot%202025-09-19%20at%2013.49.50.png)

## Loose Epitopes, with no templates


| Configuration name                  | Setting                                                       |
|-------------------------------------|-------------------------------------------------------------------|
| `epitope_cutoff`                    | 4.50 A                                 |
| `msa`                               | None                                |
| `enforce_templates`                 | False                              |
| `template_constraint_threshold`     | n/a                                  |
| `template_paths`                    | n/a                                |
| `chain_ids`                         | n/a                               |
| `template_ids`                      | n/a                                 |
| `enforce_pocket`                    | True                                |
| `pocket_constraint_threshold`       | 1.00 A                                 |
| `pocket_percentage`                 | 1.00                                  |


![Alt text](/res/img/graphs/Screenshot%202025-09-19%20at%2013.49.55.png)

## Epitopes by percentage, with templates


| Configuration name                  | Setting                                                       |
|-------------------------------------|-------------------------------------------------------------------|
| `epitope_cutoff`                    | 4.50 A                                 |
| `msa`                               | None                                |
| `enforce_templates`                 | True                              |
| `template_constraint_threshold`     | 0.01 A                                  |
| `template_paths`                    | [<path/to/template(s)>, ...]                                |
| `chain_ids`                         | [A, H, L]                                |
| `template_ids`                      | [A1, H1, L1]                                 |
| `enforce_pocket`                    | True                                |
| `pocket_constraint_threshold`       | 1.00 A                              |
| `pocket_percentage`                 | [0, 0.1, 0.25, 0.5, 1.0]                                  |


![Alt text](/res/img/graphs/Screenshot%202025-09-19%20at%2013.50.04.png)