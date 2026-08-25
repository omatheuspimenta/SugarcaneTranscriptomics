import pandas as pd
import regex as re
import os



# Phobius
predictions_phobius = pd.read_csv('../data/rgas/phobius/r570.phobius', sep=r"\s+",
    skiprows=1,
    names=["SEQUENCE_ID", "TM", "SP", "PREDICTION"])

# DeepTMHMM
# Parse domains predictions
def parse_gff(path):
    """Retorna uma lista de dicts, um por sequência/bloco."""
    with open(path, "r") as fh:
        content = fh.read()
 
    blocks = content.strip().split("//")
    records = []
 
    for block in blocks:
        block = block.strip()
        if not block:
            continue
 
        lines = [l for l in block.splitlines() if l.strip()]
 
        seq_id = None
        length = None
        n_tmrs = None
        regions = []
 
        for line in lines:
            if line.startswith("#"):
                m_len = re.search(r"#\s*(\S+)\s+Length:\s*(\d+)", line)
                m_tmr = re.search(r"Number of predicted TMRs:\s*(\d+)", line)
                if m_len:
                    seq_id = m_len.group(1)
                    length = int(m_len.group(2))
                if m_tmr:
                    n_tmrs = int(m_tmr.group(1))
            else:
                # linha de dados: SEQ_ID  regiao  start  end  ...
                parts = line.split("\t")
                if len(parts) >= 4:
                    region = parts[1].strip()
                    start = int(parts[2])
                    end = int(parts[3])
                    regions.append((region, start, end))
 
        if seq_id is None:
            continue
 
        records.append(
            {
                "SEQUENCE_ID": seq_id,
                "Length": length,
                "N_TMRs": n_tmrs,
                "regions": regions,
            }
        )
 
    return records
 
 
def to_wide_dataframe(records):
    """Converte a lista de records (com regions em lista) para um DataFrame wide."""
    rows = []
 
    for rec in records:
        row = {
            "SEQUENCE_ID": rec["SEQUENCE_ID"],
            "Length": rec["Length"],
            "N_TMRs": rec["N_TMRs"],
        }
 
        # agrupa as regiões por tipo, mantendo a ordem em que aparecem
        by_type = {}
        for region, start, end in rec["regions"]:
            by_type.setdefault(region, []).append((start, end))
 
        for region_type, occurrences in by_type.items():
            row[f"has_{region_type}"] = 1
            for i, (start, end) in enumerate(occurrences, start=1):
                row[f"{region_type}_{i}_start"] = start
                row[f"{region_type}_{i}_end"] = end
 
        rows.append(row)
 
    df = pd.DataFrame(rows)
 
    # preenche as colunas dummy que não existiram para uma dada sequência com 0
    has_cols = [c for c in df.columns if c.startswith("has_")]
    df[has_cols] = df[has_cols].fillna(0).astype(int)
 
    return df


tmhmm_records = parse_gff('../data/rgas/DeepTMHMM/TMRs.gff3')
tmhmm_filtered = [r for r in tmhmm_records if (r["N_TMRs"] or 0) >= 1]
predictions_deeptmhmm = to_wide_dataframe(tmhmm_filtered)

# DeepCoil2
"""
Parser para saídas do DeepCoil2 (um arquivo .out por proteína).

Formato de entrada (cada arquivo):
    aa   cc      raw_cc  prob_a  prob_d
    M    0.000   0.000   0.000   0.000
    A    0.000   0.000   0.000   0.000
    ...

Colunas:
    aa      - aminoácido na posição (resíduo a resíduo)
    cc      - propensão de coiled-coil "sharpened" (já pós-processada/
              detecção de picos feita pelo próprio DeepCoil2). Fora de um
              segmento candidato o valor é 0; dentro dele, o valor fica
              constante (plateau) igual à propensão do segmento.
    raw_cc  - propensão bruta, resíduo a resíduo (antes do sharpening)
    prob_a  - probabilidade do resíduo estar na posição "a" do heptad repeat
    prob_d  - probabilidade do resíduo estar na posição "d" do heptad repeat

O que este script faz:
1. Lê todos os .out de um diretório (um por SEQUENCE_ID = nome do arquivo
   sem extensão)
2. Identifica segmentos contíguos onde cc > 0 (segmentos candidatos, já
   definidos pelo próprio DeepCoil2)
3. Aplica dois filtros de qualidade, ajustáveis:
     - CC_THRESHOLD: valor mínimo de cc para considerar o segmento "forte"
       o suficiente (literatura usa 0.2 como permissivo e 0.5 como rígido)
     - MIN_SEGMENT_LEN: comprimento mínimo do segmento em resíduos
       (recomendado >= 14, ou seja, pelo menos 2 heptads)
4. Mantém só as SEQUENCE_ID que têm ao menos 1 domínio após os filtros
5. Salva em formato wide: uma linha por proteína, com dummy has_coiled_coil
   e colunas de posição/estatísticas por domínio (podem existir vários
   domínios coiled-coil na mesma proteína)

Ajuste os parâmetros no topo antes de rodar.
"""
DEEPCOIL_RESULTS_DIR = '../data/rgas/DeepCoil/SofficinarumxspontaneumR570_771_v2.1.protein.part_001'
CC_THRESHOLD = 0.20      # corte mínimo em 'cc' pra considerar o segmento válido
                          # (0.2 = permissivo, 0.5 = rígido - ver literatura)
MIN_SEGMENT_LEN = 14     # comprimento mínimo do domínio em resíduos (~2 heptads)

def parse_out_file(path):
    """Le um arquivo .out do DeepCoil2 e retorna um DataFrame por resíduo."""
    rows = []
    with open(path) as fh:
        header = fh.readline()  # "aa  cc  raw_cc  prob_a  prob_d"
        for i, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            parts = re.split(r"\s+", line)
            if len(parts) < 5:
                continue
            aa, cc, raw_cc, prob_a, prob_d = parts[:5]
            if aa == "*":  # marcador de fim de sequência
                continue
            rows.append(
                {
                    "position": i,
                    "aa": aa,
                    "cc": float(cc),
                    "raw_cc": float(raw_cc),
                    "prob_a": float(prob_a),
                    "prob_d": float(prob_d),
                }
            )
    return pd.DataFrame(rows)


def find_segments(df, cc_threshold, min_len):
    """Identifica segmentos contíguos com cc > cc_threshold e comprimento
    mínimo min_len. Retorna lista de dicts com estatísticas do segmento.

    Implementação vetorizada: constrói um array booleano (cc > threshold)
    e varre os índices onde ele é True, agrupando os que são consecutivos.
    Evita guardar "última posição vista" em uma variável de estado separada
    do laço, que é frágil e propensa a erro caso a lógica do laço mude.
    """
    above = (df["cc"] > cc_threshold).to_numpy()
    positions = df["position"].to_numpy()

    segments = []
    i = 0
    n = len(above)
    while i < n:
        if not above[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and above[j + 1]:
            j += 1
        segments.append((positions[i], positions[j]))
        i = j + 1

    results = []
    for start, end in segments:
        length = end - start + 1
        if length < min_len:
            continue
        sub = df[(df["position"] >= start) & (df["position"] <= end)]
        results.append(
            {
                "start": start,
                "end": end,
                "length": length,
                "mean_cc": round(sub["cc"].mean(), 3),
                "max_cc": round(sub["cc"].max(), 3),
                "mean_raw_cc": round(sub["raw_cc"].mean(), 3),
                "mean_prob_a": round(sub["prob_a"].mean(), 3),
                "mean_prob_d": round(sub["prob_d"].mean(), 3),
            }
        )
    return results


def build_dataframe(input_dir, cc_threshold, min_len):
    rows = []

    for fname in sorted(os.listdir(input_dir)):
        if not fname.endswith(".out"):
            continue
        seq_id = fname[:-4]  # remove ".out" -> nome do arquivo = SEQUENCE_ID
        path = os.path.join(input_dir, fname)

        df = parse_out_file(path)
        if df.empty:
            continue

        segments = find_segments(df, cc_threshold, min_len)
        if not segments:
            continue  # sem coiled-coil predito -> não entra no dataframe final

        row = {
            "SEQUENCE_ID": seq_id,
            "Length": len(df),
            "has_coiled_coil": 1,
            "N_domains": len(segments),
        }
        for i, seg in enumerate(segments, start=1):
            for key, val in seg.items():
                row[f"domain_{i}_{key}"] = val

        rows.append(row)

    result_df = pd.DataFrame(rows)
    return result_df

predictions_deepcoil = build_dataframe(DEEPCOIL_RESULTS_DIR, CC_THRESHOLD, MIN_SEGMENT_LEN)

# InterProScan
INTERPROSCAN_COLUMNS = [
    'SEQUENCE_ID',
    'sequence_md5',
    'sequence_length',
    'analysis',
    'signature_accession',
    'signature_description',
    'start_location',
    'stop_location',
    'score',
    'status',
    'date',
    'interpro_accession',
    'interpro_description',
    'go_annotations',
    'pathway_annotations'
]
predictions_interproscan = pd.read_table('../data/rgas/InterProScan/r570_interpro.tsv',
                                         names=INTERPROSCAN_COLUMNS)

# SignalP6
"""
Antes do código, vale entender as classes que o SignalP 6.0 distingue (relevante pro seu objetivo de features de domínio):

SP (Sec/SPI) — peptídeo sinal clássico, clivado pela peptidase sinal I (via geral de secreção)
LIPO (Sec/SPII) — peptídeo sinal de lipoproteína, clivado pela peptidase sinal II (mais comum em bactérias, mas o modelo cobre todos os tipos)
TAT (Tat/SPI) — peptídeo sinal da via Tat (transporta proteínas já dobradas), clivado pela peptidase I
TATLIPO (Tat/SPII) — combinação Tat + lipoproteína
PILIN (Sec/SPIII) — peptídeo sinal de pilina/pseudopilina
OTHER — nenhum peptídeo sinal predito

Como você está trabalhando com plantas (cana-de-açúcar), na prática quase tudo que não for OTHER deve cair em SP, já que TAT/LIPO/PILIN são majoritariamente de procariotos — mas vale manter todas as colunas de probabilidade para não perder informação.

A posição do CS (cleavage site) é o dado mais importante para você: indica onde termina o peptídeo sinal e começa a proteína madura — é o limite real do domínio.

Ran 2 commands, created a file · 1 note

Funcionou: das 4 sequências, só SoffiXsponR570.7os1g018900.1.p tem peptídeo sinal predito (SP, prob 0.999), com sítio de clivagem entre as posições 30-31 (probabilidade 0.9323). As outras 3 (OTHER) ficaram de fora, como esperado.

Presented file

Pontos que valem sua atenção:

SEQUENCE_ID bate com o TMHMM, mas não necessariamente com os arquivos do DeepCoil2. Aqui o ID vem do cabeçalho FASTA original (com pontos, ex. SoffiXsponR570.7os1g018900.1.p), igual ao arquivo do TMHMM. Isso facilita o merge entre esses dois, mas reforça o problema que mencionei antes: você vai precisar de uma tabela de correspondência para juntar com os resultados do DeepCoil2 (cujo nome de arquivo tinha os pontos removidos).
CS_start/CS_end marcam o fim do peptídeo sinal, não o início de um domínio funcional. Para features de classificação, o mais útil costuma ser: (a) o dummy has_signal_peptide, (b) o comprimento do peptídeo sinal (CS_end), e (c) a probabilidade de confiança (CS_prob e prob_SP) — combine isso com a informação do TMHMM: uma proteína com peptídeo sinal E sem hélices TM é candidata a proteína secretada; com peptídeo sinal E hélices TM pode ser proteína de membrana com sinal N-terminal clivado.
Cuidado com sobreposição TMHMM x SignalP. É comum o TMHMM confundir peptídeo sinal com uma hélice transmembrana (você já viu isso no primeiro arquivo que me mandou — o registro com signal 1-30 do TMHMM bate exatamente com o CS pos: 30-31 do SignalP aqui, mesma proteína). Ou seja, esses dois parsers vão gerar informação redundante/confirmatória para o N-terminal — pode valer a pena usar o SignalP como fonte de verdade para peptídeo sinal (é mais especializado nisso) e o TMHMM só para as hélices internas (TMhelix).
Só o campo Prediction já resolve o filtro que você pediu (!= OTHER), mas se quiser ser mais rigoroso, pode também exigir prob_SP > 0.5 (ou o threshold que preferir) como segunda camada de confiança, do mesmo jeito que discutimos threshold no DeepCoil2.
"""
SIGNALP_INPUT_PATH = '../data/rgas/SignalP6/prediction_results.txt'
CS_PATTERN = re.compile(r"CS pos:\s*(\d+)-(\d+)\.\s*Pr:\s*([\d.]+)")

def parse_signalp(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 8:
                continue

            fasta_header = fields[0]
            seq_id = fasta_header.split(" ")[0]  # primeiro token = SEQUENCE_ID

            prediction = fields[1]
            prob_other = float(fields[2])
            prob_sp = float(fields[3])
            prob_lipo = float(fields[4])
            prob_tat = float(fields[5])
            prob_tatlipo = float(fields[6])
            prob_pilin = float(fields[7])
            cs_field = fields[8] if len(fields) > 8 else ""

            cs_start = cs_end = cs_prob = None
            m = CS_PATTERN.search(cs_field)
            if m:
                cs_start, cs_end, cs_prob = int(m.group(1)), int(m.group(2)), float(m.group(3))

            rows.append(
                {
                    "SEQUENCE_ID": seq_id,
                    "Prediction": prediction,
                    "prob_OTHER": prob_other,
                    "prob_SP": prob_sp,
                    "prob_LIPO": prob_lipo,
                    "prob_TAT": prob_tat,
                    "prob_TATLIPO": prob_tatlipo,
                    "prob_PILIN": prob_pilin,
                    "CS_start": cs_start,
                    "CS_end": cs_end,
                    "CS_prob": cs_prob,
                }
            )

    return pd.DataFrame(rows)

predictions_signalp = parse_signalp(SIGNALP_INPUT_PATH)

# DeepLoc2
DEEPLOC_COLUMNS = ['SEQUENCE_ID', 'Localizations', 'Signals', 'Membrane types', 'Cytoplasm',
       'Nucleus', 'Extracellular', 'Cell membrane', 'Mitochondrion', 'Plastid',
       'Endoplasmic reticulum', 'Lysosome/Vacuole', 'Golgi apparatus',
       'Peroxisome', 'Peripheral', 'Transmembrane', 'Lipid anchor', 'Soluble']
predictions_deeploc = pd.read_csv('../data/rgas/DeepLoc2/results_20260812-223612.csv',
                                  names=DEEPLOC_COLUMNS,
                                  header=0)

############
# RGA Rule
RGA_PATTERNS_INTERPROSCAN = {
    "LRR":            (r"leucine[- ]rich repeat|\blrr\b", None),
    "NBS":            (r"nb-arc|nucleotide[- ]binding", {"IPR002182"}),
    "TIR":            (r"toll/interleukin-1 receptor|\btir\b", {"IPR000157"}),
    "Kinase":         (r"protein kinase|receptor-like kinase|ser/thr kinase|"
                        r"serine/threonine-protein kinase", None),
    "CC":             (r"coiled[- ]coil|\bcoil\b", None),
    "LysM":           (r"\blysm\b", None),
    "B3_DNA_binding": (r"b3 dna binding|b3 domain", {"IPR003340"}),
    "WRKY":           (r"\bwrky\b", {"IPR003657"}),
    "RPW8":           (r"rpw8", None),  # marks helper NLRs (RNL-type, e.g. ADR1/NRG1)
}


def _match_any_column(df, columns, pattern):
    """OR a regex match across whichever of `columns` actually exist."""
    hit = pd.Series(False, index=df.index)
    for col in columns:
        if col in df.columns:
            hit = hit | df[col].str.contains(
                pattern, regex=True, na=False, flags=re.IGNORECASE
            )
    return hit


def _match_accession(df, accessions, acc_col):
    if accessions and acc_col in df.columns:
        return df[acc_col].isin(accessions)
    return pd.Series(False, index=df.index)



DESC_COLS_DEFAULT = ("interpro_description", "signature_description")
FEATURE_COLS = list(RGA_PATTERNS_INTERPROSCAN.keys())

description_cols=DESC_COLS_DEFAULT
accession_col="interpro_accession"
feature_series = {"SEQUENCE_ID": predictions_interproscan["SEQUENCE_ID"]}

for feature, (pattern, accessions) in RGA_PATTERNS_INTERPROSCAN.items():
    text_hit = _match_any_column(predictions_interproscan, description_cols, pattern)
    acc_hit = _match_accession(predictions_interproscan, accessions, accession_col)
    feature_series[feature] = (text_hit | acc_hit).astype(int)

features = pd.concat(feature_series, axis=1)

features_unified = (
    features.groupby("SEQUENCE_ID", as_index=False)[FEATURE_COLS].max()
)
