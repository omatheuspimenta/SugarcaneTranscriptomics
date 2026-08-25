# RGA prediction report -- Saccharum officinarum x spontaneum R570

*Generated 2026-08-25 01:23:21 UTC by `rgas_prediction.py` v1.0.0 (config v1.0.0).*

## 1. What this report shows

299,731 proteins were examined and 29,151 (9.73%) carry at least one feature associated with plant immune receptors. These are *candidates* identified from protein domains and topology: they are not experimentally validated resistance genes.

## 2. How the call was made

Six independent annotation tools are harmonised into a single controlled vocabulary of protein features (NB-ARC, TIR, RPW8, CC, LRR, kinase, LysM, transmembrane helix, signal peptide). Protein domains come from InterProScan and are matched by accession only, never by description text, because descriptions change between releases and match unrelated entries. Overlapping hits reported by several signature databases for the same region are merged before anything is counted, so one LRR seen by Pfam, SMART and Gene3D counts once. Transmembrane helices are taken from Phobius and DeepTMHMM, signal peptides from SignalP 6.0 and Phobius, and coiled coils from DeepCoil2 with InterProScan Coils as corroboration. A helix predicted inside the signal peptide is discarded, because signal peptides are routinely mistaken for transmembrane helices. Each protein is then passed through an ordered list of mutually exclusive rules and receives the first class that fits, together with a written justification citing the exact signatures behind the call. Subcellular localisation from DeepLoc 2.0 never decides a class; it only raises or lowers the reported confidence and flags inconsistencies.

## 3. Run metadata

- Command: `code/rgas_prediction.py --input-dir data/rgas --outdir results/rgas/SaccharumR570 --organism-name Saccharum officinarum x spontaneum R570 --workers 6 --log-level INFO`
- Output directory: `results/rgas/SaccharumR570`
- Consensus policies: TM `union`, SP `signalp`, CC `deepcoil`
- Coiled-coil calling: threshold 0.5, minimum length 21 residues, maximum gap 2 residues
- Minimum LRR copies: 1

### Input files

| tool | path | available | size_bytes | n_lines | sha256 |
| --- | --- | --- | --- | --- | --- |
| interproscan | data/rgas/InterProScan/r570_interpro.tsv | True | 545439273 | 3180078 | cffdac56cf60b370d4534c4aac63285a293abe5fb002d30cda7087d87b225a14 |
| phobius | data/rgas/phobius/r570.phobius | True | 14166316 | 299732 | 7e4310b0f57960f716bf4c1ed7e8dfcbd382903bea7ddc025e60a430a6408e1b |
| deeptmhmm | data/rgas/DeepTMHMM/TMRs.gff3 | True | 71427548 | 1706445 | 0588bfc6a34ad03c4a53a3dabae5a59d02a27cd7d0571b04411ae37855835c20 |
| signalp | data/rgas/SignalP6/prediction_results.txt | True | 83538711 | 299733 | 13dbcb3ee54484985ee82af0eb1b17dffa11a16525001dce31ff7d47fb34fa8c |
| deeploc | data/rgas/DeepLoc2/results_20260812-223612.csv | True | 105642186 | 299732 | 83ba3f8add0b001b78d82452872dd8e563ae0779ea7498efb4aedd465293240c |
| deepcoil | data/rgas/DeepCoil | True | 184077592 | 5060 | cf627355ce4d0177060afe6f19bbfc1c057fe14dafd1fa06c5f4110ca1b33738 |

### Evidence channels

| channel | available |
| --- | --- |
| interproscan | True |
| phobius | True |
| deeptmhmm | True |
| signalp | True |
| deeploc | True |
| deepcoil | True |

## 4. Rules applied

| priority | rule_id | family | subclass | requires | requires_one_of | forbids | description |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | CNL | NLR | CNL | NB-ARC;CC;LRR | - | TIR;RPW8 | Coiled-coil NLR: CC + NB-ARC + LRR |
| 2 | TNL | NLR | TNL | NB-ARC;TIR;LRR | - | RPW8 | TIR NLR: TIR + NB-ARC + LRR |
| 3 | RNL | NLR | RNL | NB-ARC;RPW8;LRR | - | - | Helper NLR (ADR1/NRG1 type): RPW8 + NB-ARC + LRR |
| 4 | NL | NLR | NL | NB-ARC;LRR | - | CC;TIR;RPW8 | NB-ARC + LRR without an N-terminal CC/TIR/RPW8 domain |
| 5 | CN | NLR | CN | NB-ARC;CC | - | LRR;TIR;RPW8 | CC + NB-ARC, LRR not detected |
| 6 | TN | NLR | TN | NB-ARC;TIR | - | LRR;RPW8 | TIR + NB-ARC, LRR not detected |
| 7 | RN | NLR | RN | NB-ARC;RPW8 | - | LRR | RPW8 + NB-ARC, LRR not detected |
| 8 | N | NLR | N | NB-ARC | - | CC;TIR;RPW8;LRR | NB-ARC only |
| 9 | TX | NLR-associated | TX | TIR | - | NB-ARC;LRR | TIR-X / TIR-only: TIR without NB-ARC and without LRR |
| 10 | RX | NLR-associated | RX | RPW8 | - | NB-ARC;TIR | RPW8-X: RPW8 without NB-ARC |
| 11 | LRR-RLK | RLK | LRR-RLK | STTK;LRR | (TM OR SP) | NB-ARC;TIR;RPW8 | Kinase + LRR ectodomain + (TM or SP), no NB-ARC |
| 12 | LysM-RLK | RLK | LysM-RLK | STTK;LysM | (TM OR SP) | NB-ARC;TIR;RPW8;LRR | Kinase + LysM ectodomain + (TM or SP), no NB-ARC |
| 13 | other-RLK | RLK | other-RLK | STTK | (TM OR SP) | NB-ARC;TIR;RPW8;LRR;LysM | Kinase + (TM or SP), no recognised ectodomain, no NB-ARC |
| 14 | LRR-RLP | RLP | LRR-RLP | LRR | (TM OR SP) | NB-ARC;TIR;RPW8;STTK | LRR ectodomain + (TM or SP), no kinase, no NB-ARC |
| 15 | LysM-RLP | RLP | LysM-RLP | LysM | (TM OR SP) | NB-ARC;TIR;RPW8;STTK;LRR | LysM ectodomain + (TM or SP), no kinase, no NB-ARC |
| 16 | other-RLP | RLP | other-RLP | - | (TM OR SP) AND (LRR OR LysM) | NB-ARC;TIR;RPW8;STTK;LRR;LysM | Non-LRR/non-LysM ectodomain + (TM or SP), no kinase, no NB-ARC |
| 17 | TM-CC | TM-CC | TM-CC | TM;CC | - | NB-ARC;STTK;LRR;LysM;TIR;RPW8 | Transmembrane + coiled coil, no NB-ARC/kinase/LRR/LysM |
| 18 | Other | Other | Other | - | - | - | Carries at least one core immune feature but fits no rule above |
| 19 | Non-RGA | Non-RGA | NA | - | - | - | No core immune feature detected |

## 5. Counts

### By family

| rga_family | n_proteins | percent_of_proteome |
| --- | --- | --- |
| Non-RGA | 270580 | 90.2743 |
| Other | 11266 | 3.7587 |
| RLK | 8564 | 2.8572 |
| NLR | 4023 | 1.3422 |
| TM-CC | 3960 | 1.3212 |
| RLP | 1318 | 0.4397 |
| NLR-associated | 20 | 0.0067 |

### By subclass

| rga_family | rga_subclass | n_proteins | percent_of_proteome |
| --- | --- | --- | --- |
| NLR | NL | 3038 | 1.0136 |
| NLR | N | 470 | 0.1568 |
| NLR | CNL | 396 | 0.1321 |
| NLR | CN | 79 | 0.0264 |
| NLR | TN | 33 | 0.011 |
| NLR | RNL | 7 | 0.0023 |
| NLR-associated | TX | 20 | 0.0067 |
| Non-RGA | NA | 270580 | 90.2743 |
| Other | Other | 11266 | 3.7587 |
| RLK | other-RLK | 5527 | 1.844 |
| RLK | LRR-RLK | 2992 | 0.9982 |
| RLK | LysM-RLK | 45 | 0.015 |
| RLP | LRR-RLP | 1238 | 0.413 |
| RLP | LysM-RLP | 80 | 0.0267 |
| TM-CC | TM-CC | 3960 | 1.3212 |

### Confidence of RGA calls

| confidence | n_proteins |
| --- | --- |
| high | 23585 |
| medium | 3707 |
| low | 1859 |

### Most frequent domain architectures among RGAs

| domain_architecture | n_proteins |
| --- | --- |
| STTK | 7247 |
| LRR | 3507 |
| SP-TM-STTK | 3025 |
| NB-ARC-LRR | 2602 |
| SP-LRR-TM-STTK | 2166 |
| TM-STTK | 1812 |
| CC-TM | 1287 |
| TM-CC | 1094 |
| LRR-TM-STTK | 708 |
| STTK-TM | 550 |
| TM-CC-TM | 524 |
| SP-LRR | 478 |
| NB-ARC | 451 |
| CC-NB-ARC-LRR | 356 |
| SP-LRR-TM | 341 |
| LRR-TM | 318 |
| CC-STTK | 239 |
| TM-CC-TM-CC-TM-CC-TM-CC-TM | 179 |
| STTK-CC | 170 |
| NB-ARC-LRR-TM | 164 |

## 6. Coiled-coil evidence

Coiled coils are the least reliable feature in every published RGA pipeline. InterProScan's Coils/ncoils module under-detects them, which is why the NLRtracker benchmark (Kourelis et al. 2021) reports CC as the domain most frequently missed. This pipeline therefore uses DeepCoil2 as the primary coiled-coil channel and keeps InterProScan Coils as corroborating evidence. The two tables below show how much the two methods disagree, and how sensitive the NLR subclass counts are to that choice.

### DeepCoil2 versus InterProScan Coils (whole proteome)

| InterProScan Coils | DeepCoil2 | n_proteins |
| --- | --- | --- |
| CC called | CC called | 14165 |
| no CC | CC called | 4579 |
| CC called | no CC | 28135 |
| no CC | no CC | 252852 |

### Subclass counts under each `--cc-policy`

| rga_subclass | deepcoil | union | intersection | coils |
| --- | --- | --- | --- | --- |
| CN | 79 | 230 | 40 | 191 |
| CNL | 396 | 1790 | 244 | 1638 |
| LRR-RLK | 2992 | 2992 | 2992 | 2992 |
| LRR-RLP | 1238 | 1238 | 1238 | 1238 |
| LysM-RLK | 45 | 45 | 45 | 45 |
| LysM-RLP | 80 | 80 | 80 | 80 |
| N | 470 | 319 | 509 | 358 |
| NA | 270580 | 266439 | 272468 | 268327 |
| NL | 3038 | 1644 | 3190 | 1796 |
| Other | 11266 | 11266 | 11266 | 11266 |
| RNL | 7 | 7 | 7 | 7 |
| TM-CC | 3960 | 8101 | 2072 | 6213 |
| TN | 33 | 33 | 33 | 33 |
| TX | 20 | 20 | 20 | 20 |
| other-RLK | 5527 | 5527 | 5527 | 5527 |

### Sensitivity to the segment-calling parameters

| threshold | min_length | n_proteins_with_cc | n_segments |
| --- | --- | --- | --- |
| 0.2 | 14 | 63371 | 177535 |
| 0.2 | 21 | 28319 | 55427 |
| 0.2 | 28 | 16933 | 31658 |
| 0.5 | 14 | 29558 | 58118 |
| 0.5 | 21 | 18744 | 36390 |
| 0.5 | 28 | 13744 | 25126 |

## 7. Identifier reconciliation

| tool | n_ids | n_shared_with_proteome | n_absent_from_tool | n_not_in_proteome |
| --- | --- | --- | --- | --- |
| deepcoil | 299731 | 299731 | 0 | 0 |
| deeploc | 299731 | 299731 | 0 | 0 |
| deeptmhmm | 299731 | 299731 | 0 | 0 |
| interproscan | 284207 | 284207 | 15524 | 0 |
| phobius | 299731 | 299731 | 0 | 0 |
| signalp | 299731 | 299731 | 0 | 0 |

## 8. Warnings

| warning | n_proteins |
| --- | --- |
| 1 TM helix/helices discarded as signal peptide | 4875 |
| CC segment overlaps a predicted TM helix (possible artefact) | 1250 |
| DeepLoc localisation (Cell membrane) is inconsistent with class NL | 890 |
| DeepLoc localisation (Endoplasmic reticulum) is inconsistent with class TM-CC | 887 |
| DeepLoc localisation (Cytoplasm) is inconsistent with class other-RLK | 429 |
| DeepLoc localisation (Lysosome/Vacuole) is inconsistent with class TM-CC | 423 |
| DeepLoc localisation (Cytoplasm) is inconsistent with class TM-CC | 301 |
| DeepLoc localisation (Golgi apparatus) is inconsistent with class TM-CC | 276 |
| DeepLoc localisation (Mitochondrion) is inconsistent with class TM-CC | 268 |
| DeepLoc localisation (Nucleus) is inconsistent with class TM-CC | 253 |
| DeepLoc localisation (Plastid) is inconsistent with class TM-CC | 157 |
| DeepLoc localisation (Cell membrane) is inconsistent with class CNL | 123 |
| 2 TM helix/helices discarded as signal peptide | 113 |
| DeepLoc localisation (Cytoplasm) is inconsistent with class LRR-RLP | 109 |
| DeepLoc localisation (Cell membrane) is inconsistent with class N | 77 |
| DeepLoc localisation (Nucleus) is inconsistent with class LRR-RLP | 64 |
| DeepLoc localisation (Mitochondrion) is inconsistent with class other-RLK | 54 |
| DeepLoc localisation (Nucleus) is inconsistent with class other-RLK | 48 |
| DeepLoc localisation (Plastid) is inconsistent with class other-RLK | 29 |
| DeepLoc localisation (Peroxisome) is inconsistent with class TM-CC | 26 |
| DeepLoc localisation (Endoplasmic reticulum) is inconsistent with class other-RLK | 19 |
| DeepLoc localisation (Endoplasmic reticulum) is inconsistent with class LRR-RLP | 12 |
| DeepLoc localisation (Plastid) is inconsistent with class LRR-RLP | 11 |
| DeepLoc localisation (Cell membrane) is inconsistent with class CN | 9 |
| DeepLoc localisation (Endoplasmic reticulum) is inconsistent with class LysM-RLP | 4 |
| CC lies C-terminal to the NB-ARC domain | 3 |
| DeepLoc localisation (Lysosome/Vacuole) is inconsistent with class other-RLK | 2 |
| DeepLoc localisation (Cytoplasm) is inconsistent with class LRR-RLK | 1 |
| DeepLoc localisation (Extracellular) is inconsistent with class N | 1 |
| DeepLoc localisation (Lysosome/Vacuole) is inconsistent with class LRR-RLP | 1 |
| DeepLoc localisation (Nucleus) is inconsistent with class LRR-RLK | 1 |

## 9. Top 50 RGA candidates

| protein_id | rga_family | rga_subclass | domain_architecture | n_lrr | predicted_localization | confidence |
| --- | --- | --- | --- | --- | --- | --- |
| SoffiXsponR570.10Bg047900.2.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.03Bg092200.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.03Cg107900.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.03Cg107900.2.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.03Fg121100.2.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.03Fg121200.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.03Gg009600.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.03Gg010200.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.03Gg010400.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.03Gg011000.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.03Gg011100.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.06Ag204700.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.09Ag033800.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.09Ag067700.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.09Ag067700.2.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.09Ag067700.3.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.09Ag067700.4.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.09Eg064600.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.09Eg064600.2.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.09Eg064600.3.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.10Bg048200.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.10Bg048200.2.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.10Cg033500.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.10Cg033500.2.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.10Cg033500.3.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.3os1g002100.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.3os1g002600.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.3os1g003200.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.5_9Ag287800.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.6us88g079500.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Nucleus | high |
| SoffiXsponR570.6us88g113400.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.6us88g113900.1.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.6us88g113900.2.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.6us88g113900.4.p | NLR | CNL | CC-NB-ARC-LRR | 3 | Cytoplasm | high |
| SoffiXsponR570.02Bg270000.1.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.02Bg270000.2.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.02Cg297200.1.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.02Cg297200.2.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.02Cg297200.3.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.02Dg296800.4.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.02Eg269400.3.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.02Fg279600.1.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.02Fg279600.2.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.02Gg200300.1.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.05Dg017800.1.p | NLR | CNL | STTK-CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.05Dg017800.2.p | NLR | CNL | STTK-CC-NB-ARC-LRR | 2 | Nucleus | high |
| SoffiXsponR570.05Dg017800.3.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Nucleus | high |
| SoffiXsponR570.06Dg138400.1.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.06Dg138400.2.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |
| SoffiXsponR570.09Ag034200.1.p | NLR | CNL | CC-NB-ARC-LRR | 2 | Cytoplasm | high |

## 10. References

1. Rody HVS, Bombardelli RGH, Creste S, Camargo LEA, Van Sluys M-A, Monteiro-Vitorello CB (2019). Genome survey of resistance gene analogs in sugarcane: genomic features and differential expression of the innate immune system from a smut-resistant genotype. BMC Genomics 20:809. doi:10.1186/s12864-019-6207-y
2. Li P, Quan X, Jia G, Xiao J, Cloutier S, You FM (2016). RGAugury: a pipeline for genome-wide prediction of resistance gene analogs (RGAs) in plants. BMC Genomics 17:852. doi:10.1186/s12864-016-3197-x
3. Sekhwal MK, Li P, Lam I, Wang X, Cloutier S, You FM (2015). Disease resistance gene analogs (RGAs) in plants. Int J Mol Sci 16:19248-19290. doi:10.3390/ijms160819248
4. Kourelis J, Sakai T, Adachi H, Kamoun S (2021). RefPlantNLR is a comprehensive collection of experimentally validated plant disease resistance proteins from the NLR family. PLoS Biology 19(10):e3001124. doi:10.1371/journal.pbio.3001124
5. Smith M, Jones JT, Hein I (2025). Resistify: a novel NLR classifier that reveals Helitron-associated NLR expansion in Solanaceae. Bioinform Biol Insights 19:11779322241308944. doi:10.1177/11779322241308944
6. Shiu S-H, Bleecker AB (2003). Expansion of the receptor-like kinase/Pelle gene family and receptor-like proteins in Arabidopsis. Plant Physiol 132:530-543. doi:10.1104/pp.103.021964
7. Jones JDG, Dangl JL (2006). The plant immune system. Nature 444:323-329. doi:10.1038/nature05286
8. Jones P et al. (2014). InterProScan 5: genome-scale protein function classification. Bioinformatics 30:1236-1240. doi:10.1093/bioinformatics/btu031
9. Blum M et al. (2025). InterPro: the protein sequence classification resource in 2025. Nucleic Acids Res 53:D444-D456. doi:10.1093/nar/gkae1082
10. Paysan-Lafosse T et al. (2025). The Pfam protein families database: embracing AI/ML. Nucleic Acids Res 53:D523-D534. doi:10.1093/nar/gkae997
11. Kall L, Krogh A, Sonnhammer ELL (2004). A combined transmembrane topology and signal peptide prediction method. J Mol Biol 338:1027-1036. doi:10.1016/j.jmb.2004.03.016
12. Hallgren J et al. (2022). DeepTMHMM predicts alpha and beta transmembrane proteins using deep neural networks. bioRxiv. doi:10.1101/2022.04.08.487609
13. Teufel F et al. (2022). SignalP 6.0 predicts all five types of signal peptides using protein language models. Nat Biotechnol 40:1023-1025. doi:10.1038/s41587-021-01156-3
14. Thumuluri V et al. (2022). DeepLoc 2.0: multi-label subcellular localization prediction using protein language models. Nucleic Acids Res 50:W228-W234. doi:10.1093/nar/gkac278
15. Ludwiczak J, Winski A, Szczepaniak K, Alva V, Dunin-Horkawicz S (2019). DeepCoil - a fast and accurate prediction of coiled-coil domains in protein sequences. Bioinformatics 35(16):2790-2795. doi:10.1093/bioinformatics/bty1062
16. Lupas A, Van Dyke M, Stock J (1991). Predicting coiled coils from protein sequences. Science 252:1162-1164. doi:10.1126/science.252.5009.1162
