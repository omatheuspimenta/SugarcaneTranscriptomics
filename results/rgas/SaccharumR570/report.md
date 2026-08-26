# RGA prediction report -- SaccharumR570

*Generated 2026-08-26 19:49:56 UTC by `rgas_prediction.py` v0.0.1 (config v0.0.1).*

## 1. What this report shows

299,731 proteins were examined and 33,296 (11.11%) carry at least one feature associated with plant immune receptors. These are *candidates* identified from protein domains and topology: they are not experimentally validated resistance genes.

## 2. How the call was made

Six independent annotation tools are harmonised into a single controlled vocabulary of protein features (NB-ARC, TIR, RPW8, CC, LRR, kinase, LysM, transmembrane helix, signal peptide). Protein domains come from InterProScan and are matched by accession only, never by description text, because descriptions change between releases and match unrelated entries. Overlapping hits reported by several signature databases for the same region are merged before anything is counted, so one LRR seen by Pfam, SMART and Gene3D counts once. Transmembrane helices are taken from Phobius and DeepTMHMM, signal peptides from SignalP 6.0 and Phobius, and coiled coils from three channels: a domain-level profile HMM plus the DeepCoil2 and InterProScan Coils predictors. A helix predicted inside the signal peptide is discarded, because signal peptides are routinely mistaken for transmembrane helices. Each protein is then passed through an ordered list of mutually exclusive rules and receives the first class that fits, together with a written justification citing the exact signatures behind the call. Subcellular localisation from DeepLoc 2.0 never decides a class; it only raises or lowers the reported confidence and flags inconsistencies.

## 3. Run metadata

### Reproduce this run

The exact command that produced this report, quoted as it was invoked, so it can be pasted back into a shell from the repository root:

```bash
uv run python code/rgas_prediction.py \
    --input-dir data/rgas/ \
    --outdir results/rgas/SaccharumR570/ \
    --organism-name SaccharumR570
```

### Settings

- Output directory: `results/rgas/SaccharumR570`
- Consensus policies: TM `union`, SP `signalp`, CC `union`
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
| deepcoil | data/rgas/DeepCoil | True | 123688840 | 60 | 6733277c760b56236d50741c0e8ac58d101ce92488f38d88b1a9f1f68654b947 |

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
| Non-RGA | 266435 | 88.8914 |
| Other | 16793 | 5.6027 |
| TM-CC | 8105 | 2.7041 |
| NLR | 4023 | 1.3422 |
| RLK | 3037 | 1.0132 |
| RLP | 1318 | 0.4397 |
| NLR-associated | 20 | 0.0067 |

### By subclass

| rga_family | rga_subclass | n_proteins | percent_of_proteome |
| --- | --- | --- | --- |
| NLR | CNL | 2648 | 0.8835 |
| NLR | NL | 786 | 0.2622 |
| NLR | CN | 358 | 0.1194 |
| NLR | N | 191 | 0.0637 |
| NLR | TN | 33 | 0.011 |
| NLR | RNL | 7 | 0.0023 |
| NLR-associated | TX | 20 | 0.0067 |
| Non-RGA | NA | 266435 | 88.8914 |
| Other | Other | 16793 | 5.6027 |
| RLK | LRR-RLK | 2992 | 0.9982 |
| RLK | LysM-RLK | 45 | 0.015 |
| RLP | LRR-RLP | 1238 | 0.413 |
| RLP | LysM-RLP | 80 | 0.0267 |
| TM-CC | TM-CC | 8105 | 2.7041 |

### Confidence of RGA calls

| confidence | n_proteins |
| --- | --- |
| high | 24059 |
| low | 6236 |
| medium | 3001 |


### Confidence by subclass

| rga_subclass | high | medium | low | n_proteins |
| --- | --- | --- | --- | --- |
| Other | 16793 | 0 | 0 | 16793 |
| TM-CC | 420 | 1751 | 5934 | 8105 |
| LRR-RLK | 2990 | 2 | 0 | 2992 |
| CNL | 1617 | 782 | 249 | 2648 |
| LRR-RLP | 1041 | 197 | 0 | 1238 |
| NL | 613 | 173 | 0 | 786 |
| CN | 240 | 65 | 53 | 358 |
| N | 164 | 27 | 0 | 191 |
| LysM-RLP | 76 | 4 | 0 | 80 |
| LysM-RLK | 45 | 0 | 0 | 45 |
| TN | 33 | 0 | 0 | 33 |
| TX | 20 | 0 | 0 | 20 |
| RNL | 7 | 0 | 0 | 7 |

### Most frequent domain architectures among RGAs

| domain_architecture | n_proteins |
| --- | --- |
| STTK | 6381 |
| LRR | 3280 |
| SP-TM-STTK | 2950 |
| TM-CC | 2928 |
| CC-TM | 2584 |
| CC-NB-ARC-LRR | 2230 |
| SP-LRR-TM-STTK | 2154 |
| TM-STTK | 1750 |
| TM-CC-TM | 1291 |
| LRR-TM-STTK | 703 |
| STTK-CC | 692 |
| NB-ARC-LRR | 665 |
| CC-STTK | 524 |
| SP-LRR | 474 |
| STTK-TM | 455 |
| SP-LRR-TM | 338 |
| CC-NB-ARC | 321 |
| LRR-TM | 301 |
| NB-ARC | 188 |
| TM-CC-TM-CC-TM-CC-TM-CC-TM | 179 |

## 6. Coiled-coil evidence

The coiled coil is the least reliable feature in every published RGA pipeline, and it is the one that decides CNL against NL. Three channels are used here, and they are not of equal weight. The leading one is a curated profile HMM for a named domain (the Rx N-terminal domain, PF18052 / IPR041118), which carries the same kind of evidence as the NB-ARC model every NLR call already rests on. The other two, DeepCoil2 and InterProScan Coils, are biophysical propensity predictors: neither publishes a recommended score cut-off, and Simm et al. (2021), benchmarking coiled-coil predictors against the whole PDB, found a 30-fold spread in how many coiled coils they call and agreement with structure close to random. They are kept because they cover proteins no domain model reaches, and a call resting on them alone is graded down rather than hidden. The tables below show how much the channels disagree and how far the subclass counts move with the policy.

### DeepCoil2 versus InterProScan Coils (whole proteome)

| InterProScan Coils | DeepCoil2 | Rx domain | n_proteins |
| --- | --- | --- | --- |
| CC called | CC called | -- | 14165 |
| no CC | CC called | -- | 4579 |
| CC called | no CC | -- | 28135 |
| no CC | no CC | -- | 252852 |
| -- | -- | CC called | 2827 |
| no CC | no CC | CC called | 1056 |

### Subclass counts under each `--cc-policy`

| rga_subclass | rx_domain | deepcoil | coils | union | intersection |
| --- | --- | --- | --- | --- | --- |
| CN | 283 | 79 | 191 | 358 | 30 |
| CNL | 2328 | 396 | 1638 | 2648 | 185 |
| LRR-RLK | 2992 | 2992 | 2992 | 2992 | 2992 |
| LRR-RLP | 1238 | 1238 | 1238 | 1238 | 1238 |
| LysM-RLK | 45 | 45 | 45 | 45 | 45 |
| LysM-RLP | 80 | 80 | 80 | 80 | 80 |
| N | 266 | 470 | 358 | 191 | 519 |
| NA | 274529 | 270580 | 268327 | 266435 | 274540 |
| NL | 1106 | 3038 | 1796 | 786 | 3249 |
| Other | 16793 | 16793 | 16793 | 16793 | 16793 |
| RNL | 7 | 7 | 7 | 7 | 7 |
| TM-CC | 11 | 3960 | 6213 | 8105 | 0 |
| TN | 33 | 33 | 33 | 33 | 33 |
| TX | 20 | 20 | 20 | 20 | 20 |

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
| CC supported only by InterProScan Coils | 4436 |
| DeepLoc localisation (Endoplasmic reticulum) is inconsistent with class TM-CC | 1580 |
| CC segment overlaps a predicted TM helix (possible artefact) | 1252 |
| DeepLoc localisation (Cell membrane) is inconsistent with class CNL | 840 |
| DeepLoc localisation (Cytoplasm) is inconsistent with class TM-CC | 704 |
| DeepLoc localisation (Lysosome/Vacuole) is inconsistent with class TM-CC | 676 |
| DeepLoc localisation (Mitochondrion) is inconsistent with class TM-CC | 655 |
| DeepLoc localisation (Nucleus) is inconsistent with class TM-CC | 639 |
| DeepLoc localisation (Plastid) is inconsistent with class TM-CC | 535 |
| DeepLoc localisation (Golgi apparatus) is inconsistent with class TM-CC | 509 |
| DeepLoc localisation (Cell membrane) is inconsistent with class NL | 173 |
| 2 TM helix/helices discarded as signal peptide | 113 |
| DeepLoc localisation (Cytoplasm) is inconsistent with class LRR-RLP | 109 |
| CC lies C-terminal to the NB-ARC domain | 105 |
| DeepLoc localisation (Nucleus) is inconsistent with class LRR-RLP | 64 |
| DeepLoc localisation (Cell membrane) is inconsistent with class CN | 59 |
| DeepLoc localisation (Peroxisome) is inconsistent with class TM-CC | 38 |
| DeepLoc localisation (Cell membrane) is inconsistent with class N | 27 |
| DeepLoc localisation (Endoplasmic reticulum) is inconsistent with class LRR-RLP | 12 |
| DeepLoc localisation (Plastid) is inconsistent with class LRR-RLP | 11 |
| DeepLoc localisation (Endoplasmic reticulum) is inconsistent with class LysM-RLP | 4 |
| DeepLoc localisation (Cytoplasm) is inconsistent with class LRR-RLK | 1 |
| DeepLoc localisation (Extracellular) is inconsistent with class CN | 1 |
| DeepLoc localisation (Lysosome/Vacuole) is inconsistent with class LRR-RLP | 1 |
| DeepLoc localisation (Nucleus) is inconsistent with class LRR-RLK | 1 |

## 9. Top 50 RGA candidates

| protein_id | rga_family | rga_subclass | domain_architecture | n_lrr | predicted_localization | confidence |
| --- | --- | --- | --- | --- | --- | --- |
| SoffiXsponR570.5_9Ag288500.1.p | NLR | CNL | CC-NB-ARC-LRR | 9 | Nucleus | high |
| SoffiXsponR570.04Cg195600.1.p | NLR | CNL | CC-NB-ARC-LRR | 7 | Cytoplasm | high |
| SoffiXsponR570.9us90g054500.1.p | NLR | CNL | CC-NB-ARC-LRR | 7 | Nucleus | high |
| SoffiXsponR570.02Cg150300.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Nucleus | high |
| SoffiXsponR570.02Cg150300.2.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Nucleus | high |
| SoffiXsponR570.02Cg150300.3.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Nucleus | high |
| SoffiXsponR570.02Dg152200.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Nucleus | high |
| SoffiXsponR570.02Dg152800.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Nucleus | high |
| SoffiXsponR570.02Gg053900.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Nucleus | high |
| SoffiXsponR570.02Gg053900.2.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Nucleus | high |
| SoffiXsponR570.04Ag202900.2.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Cytoplasm | high |
| SoffiXsponR570.04Ag202900.3.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Cytoplasm | high |
| SoffiXsponR570.04Ag202900.4.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Cytoplasm | high |
| SoffiXsponR570.04Dg174500.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Nucleus | high |
| SoffiXsponR570.04Eg130500.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Nucleus | high |
| SoffiXsponR570.07Cg214600.2.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Nucleus | high |
| SoffiXsponR570.07Cg214600.3.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Nucleus | high |
| SoffiXsponR570.08Dg029300.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Nucleus | high |
| SoffiXsponR570.10Ag119400.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Cytoplasm | high |
| SoffiXsponR570.10Bg073500.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Cytoplasm | high |
| SoffiXsponR570.10Cg100200.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Cytoplasm | high |
| SoffiXsponR570.10Dg021500.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Cytoplasm | high |
| SoffiXsponR570.10Dg098400.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Cytoplasm | high |
| SoffiXsponR570.10os1g051800.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Cytoplasm | high |
| SoffiXsponR570.8_10Ag037600.1.p | NLR | CNL | CC-NB-ARC-LRR | 5 | Cytoplasm | high |
| SoffiXsponR570.02Gg056600.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.02Gg056600.2.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.02Gg298300.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.03Ag004300.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.03Ag137500.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Cytoplasm | high |
| SoffiXsponR570.03Bg004500.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.03Bg119800.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Cytoplasm | high |
| SoffiXsponR570.03Cg133900.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Cytoplasm | high |
| SoffiXsponR570.03Dg144500.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Cytoplasm | high |
| SoffiXsponR570.03Eg106300.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.03Fg120900.1.p | NLR | CNL | CC-NB-ARC-TM-LRR-TM-NB-ARC-LRR | 4 | Cytoplasm | high |
| SoffiXsponR570.03Fg151800.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Cytoplasm | high |
| SoffiXsponR570.04Ag252200.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.04Bg011200.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Cytoplasm | high |
| SoffiXsponR570.04Bg011200.2.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.04Cg244400.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.04Cg244400.2.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.04Cg244400.3.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.04Fg132600.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.04Fg132600.2.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.05Dg145200.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Cytoplasm | high |
| SoffiXsponR570.07Bg024300.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Cytoplasm | high |
| SoffiXsponR570.07Bg215900.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.07Cg214600.1.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |
| SoffiXsponR570.07Cg214600.4.p | NLR | CNL | CC-NB-ARC-LRR | 4 | Nucleus | high |

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
17. Simm D, Hatje K, Waack S, Kollmar M (2021). Critical assessment of coiled-coil predictions based on protein structure data. Scientific Reports 11:12439. doi:10.1038/s41598-021-91886-w
