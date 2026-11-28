
from __future__ import annotations
from collections import Counter, defaultdict
from dataclasses import dataclass
import argparse
import json
import random
import statistics
import sys
import time
import os
import threading
from pathlib import Path

DNA = 'ACGT'
COMPLEMENT = str.maketrans('ACGT', 'TGCA')
# Gjenomi real fillestar: ADN mitokondriale e njeriut (NCBI NC_012920.1, 16.569 bp).
# Skedari FASTA duhet të jetë në të njëjtën dosje me këtë program.
DEFAULT_GENOME = Path(__file__).resolve().parent / 'mtADN_njeriut.fasta'
# Mbulimi i rekomanduar: sa herë lexohet mesatarisht çdo bazë e gjenomit.
RECOMMENDED_COVERAGE = 25
# Cilësia (Phred) që u jepet bazave kur reads vijnë pa cilësi (skedar FASTA).
DEFAULT_QUALITY = 30


def reverse_complement(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def write_fasta(path: Path, records):
    with path.open('w') as f:
        for name, seq in records:
            f.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + '\n')


def write_fastq(path: Path, records):
    """Shkruan reads në FASTQ: @emri, sekuenca, '+', cilësia e bazave (Phred + 33)."""
    with path.open('w') as f:
        for name, seq, quality in records:
            f.write(f'@{name}\n{seq}\n+\n{quality}\n')


def read_fasta(path: Path):
    records = []
    title, parts = None, []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if title is not None:
                records.append((title, ''.join(parts).upper()))
            title, parts = line[1:], []
        elif title is not None:
            parts.append(line)
        else:
            raise ValueError('FASTA must start with a >header')
    if title is not None:
        records.append((title, ''.join(parts).upper()))
    if not records or any(not s or set(s) - set(DNA) for _, s in records):
        raise ValueError('Need non-empty DNA reads containing only A,C,G,T')
    return records


def read_reads(path: Path):
    """Lexon reads nga FASTA ose FASTQ dhe kthen listën e (sekuenca, cilësia).
    FASTQ ka cilësi (Phred) për çdo bazë; te FASTA cilësia mungon (None)."""
    text = path.read_text()
    if not text.lstrip().startswith('@'):
        return [(seq, None) for _, seq in read_fasta(path)]
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) % 4 != 0:
        raise ValueError('FASTQ must have exactly 4 lines per read')
    reads = []
    for i in range(0, len(lines), 4):
        seq, quality = lines[i + 1].upper(), lines[i + 3]
        if (not lines[i].startswith('@') or not lines[i + 2].startswith('+')
                or not seq or len(seq) != len(quality) or set(seq) - set(DNA)):
            raise ValueError(f'Invalid FASTQ record at line {i + 1}')
        reads.append((seq, quality))
    return reads


def read_genome(path: Path):
    """Lexon gjenomin real nga skedari FASTA (merret vetëm rekordi i parë).
    Mbahen vetëm bazat A, C, G, T; shkronjat e tjera (p.sh. 'N') hiqen dhe numërohen.
    Te NC_012920.1 ka vetëm një 'N' (pozicioni 3107): ajo nuk është bazë e vërtetë,
    por vend bosh që ruan numërimin historik të referencës mitokondriale."""
    name = None
    bases = []
    removed = 0
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name is not None:
                break  # mjafton rekordi i parë
            name = line[1:].strip()
        elif name is None:
            raise ValueError('FASTA must start with a >header')
        else:
            for base in line.upper():
                if base in DNA:
                    bases.append(base)
                else:
                    removed += 1
    genome = ''.join(bases)
    if not genome:
        raise ValueError(f'No A/C/G/T bases found in {path}')
    return name, genome, removed


def add_read_errors(fragment, rng, error_rate, indel_rate):
    """Shton gabime leximi në një copë ADN dhe i jep çdo baze një cilësi Phred.
    Si te sekuenatorët e vërtetë, bazat e gabuara zakonisht kanë cilësi të ulët (3-20),
    ndërsa bazat e sakta cilësi të lartë (15-40). Kthen (sekuenca, cilësia si tekst)."""
    bases = []
    quality = []
    for base in fragment:
        r = rng.random()
        if r < indel_rate / 2:  # fshirje e shkaktuar nga gabimi i leximit
            continue
        if r < indel_rate:      # shtim i shkaktuar nga gabimi i leximit
            bases.append(rng.choice(DNA))
            quality.append(chr(rng.randint(3, 20) + 33))
        if rng.random() < error_rate:
            bases.append(rng.choice([b for b in DNA if b != base]))
            quality.append(chr(rng.randint(3, 20) + 33))
        else:
            bases.append(base)
            quality.append(chr(rng.randint(15, 40) + 33))
    return ''.join(bases), ''.join(quality)


def simulate(reference, n_reads=3600, read_length=115, error_rate=.01,
             indel_rate=.0015, seed=7, length_spread=15, duplicate_rate=.05,
             coverage_variation=.2):
    """Merr gjenomin REAL si referencë, i vendos disa mutacione (kjo është mostra)
    dhe pastaj krijon reads nga mostra, si te një sekuencim i vërtetë:
      - gabime leximi (zëvendësime, shtime, fshirje) dhe cilësi Phred për çdo bazë,
      - gjatësi të ndryshme të reads (read_length ± length_spread),
      - mbulim të ndryshëm sipas zonave (coverage_variation),
      - reads të dyfishuara (duplicate_rate) dhe reads në drejtimin e kundërt.
    Pozicionet e vërteta dhe ADN-ja e mostrës nuk i jepen algoritmit të rindërtimit."""
    rng = random.Random(seed)
    genome_length = len(reference)
    # Mutacionet biologjike (8 SNP, 1 shtim dhe 1 fshirje) dallohen nga
    # gabimet e leximit. Pozicionet mbahen larg për teste të qarta.
    sites = []
    possible = list(range(60, genome_length - 60))
    rng.shuffle(possible)
    for pos in possible:
        if all(abs(pos-other) > 25 for other in sites):
            sites.append(pos)
        if len(sites) == min(10, max(3, genome_length // 140)):
            break
    if len(sites) < 3:
        raise ValueError('Genome must have at least 200 bases')
    ins_pos, del_pos = sites[:2]
    snp_positions = set(sites[2:])
    true_variants = []
    sample_chars = []
    inserted = rng.choice([b for b in DNA if b != reference[ins_pos] and b != reference[ins_pos + 1]])
    for i, old in enumerate(reference):
        if i == del_pos:
            true_variants.append({'type':'DEL', 'position':i+1, 'ref':old, 'alt':'-'})
            continue
        base = old
        if i in snp_positions:
            base = rng.choice([b for b in DNA if b != old])
            true_variants.append({'type':'SNP', 'position':i+1, 'ref':old, 'alt':base})
        sample_chars.append(base)
        if i == ins_pos:
            sample_chars.append(inserted)
            true_variants.append({'type':'INS', 'position':i+2, 'ref':'-', 'alt':inserted})
    sample = ''.join(sample_chars)
    true_variants.sort(key=lambda v:v['position'])
    if len(sample) <= read_length + length_spread:
        raise ValueError('Read length must be shorter than genome length')
    # Mbulim i ndryshëm: ADN-ja ndahet në blloqe 1000 bazash dhe çdo blloku i jepet
    # një peshë midis (1 - coverage_variation) dhe (1 + coverage_variation).
    block = 1000
    blocks = list(range(0, len(sample), block))
    weights = [1 + rng.uniform(-coverage_variation, coverage_variation) for _ in blocks]
    # Sigurojmë mbulim në të dy skajet (deri në 5 reads për skaj), që ADN-ja të
    # rindërtohet plotësisht dhe edhe bazat e skajeve të kenë shumicë të qartë votash.
    edge_reads = max(1, min(5, n_reads // 4))
    records = []
    for ix in range(n_reads):
        name = f'read_{ix+1:06d}'
        if ix >= 2 * edge_reads and rng.random() < duplicate_rate:
            # Read i dyfishuar (p.sh. nga PCR): kopje e saktë e një read-i të mëparshëm.
            copy_name, seq, quality = rng.choice(records)
            records.append((name, seq, quality))
            continue
        length = rng.randint(read_length - length_spread, read_length + length_spread)
        max_start = len(sample) - length
        if ix < edge_reads:
            start = 0
        elif ix < 2 * edge_reads:
            start = max_start
        else:
            while True:  # zgjedhim një bllok sipas peshës dhe një pozicion brenda tij
                block_start = rng.choices(blocks, weights)[0]
                start = rng.randint(block_start, block_start + block - 1)
                if start <= max_start:
                    break
        seq, quality = add_read_errors(sample[start:start + length], rng, error_rate, indel_rate)
        if rng.random() < .5:
            seq, quality = reverse_complement(seq), quality[::-1]
        records.append((name, seq, quality))
    rng.shuffle(records)
    return sample, records, true_variants


@dataclass(frozen=True, slots=True)  # slots=True: pa fjalor për çdo objekt, më pak memorie
class Edge:
    src: int  # ID e orientimit: 2 * indeksi_i_read + (0 përpara / 1 mbrapsht)
    dst: int
    shift: int  # pozicioni ku fillon read-i i dytë në raport me të parin
    target_consumed: int  # pjesa fillestare e read-it B që mbulohet nga mbivendosja
    overlap: int
    identity: float
    votes: int

    @property
    def score(self):
        # Eksponenti i madh: çdo gabim në mbivendosje "kushton" rreth 20 baza mbivendosje.
        # Në zonat e përsëritura, mbivendosjet e rreme kanë zakonisht pak më shumë gabime
        # se ato të vërteta, prandaj DFS-ja i provon më vonë.
        return self.overlap * (self.identity ** 20) + min(self.votes, 10) * .3


def oriented_reads(reads):
    return [s for seq in reads for s in (seq, reverse_complement(seq))]


def remove_duplicates(reads):
    """Heq reads identike (p.sh. dyfishime nga PCR) dhe kthen (reads unike, sa u hoqën).
    Një read dhe komplementi i tij i kundërt konsiderohen i njëjti read."""
    seen = set()
    unique = []
    for seq, quality in reads:
        rc = reverse_complement(seq)
        key = seq if seq < rc else rc
        if key in seen:
            continue
        seen.add(key)
        unique.append((seq, quality))
    return unique, len(reads) - len(unique)


def approximate_overlap(a: str, b: str, estimated_shift: int,
                        min_overlap: int, max_error: float, padding=7, band=10):
    """Gjen mbivendosjen e përafërt me programim dinamik (edit distance).
    Fillimi i A dhe fundi i B mund të mbeten jashtë përputhjes pa penalizim.
    Votat e k-merëve tregojnë afërsisht ku fillon mbivendosja (diagonalja e pritur),
    prandaj llogarisim vetëm qelizat në një brez ±band rreth saj (banded DP).
    Qelizat larg diagonales nuk japin mbivendosje të pranueshme, kështu që rezultati
    praktikisht nuk ndryshon, por llogaritja bëhet disa herë më e shpejtë.
    Kthen (zhvendosja, baza të përdorura nga B, mbivendosja, ngjashmëria) ose None."""
    astart = max(0, estimated_shift - padding)
    x = a[astart:]
    # Read-i B duhet të vazhdojë përtej A; kontrollojmë vetëm pjesën që mund të mbivendoset.
    approx = len(a) - estimated_shift
    y = b[:min(len(b), approx + padding * 2 + 4)]
    if len(x) < min_overlap or len(y) < min_overlap:
        return None
    big = len(x) + len(y) + 1           # "pafundësi": më e madhe se çdo kosto e mundshme
    diagonal = estimated_shift - astart  # B pritet të fillojë te rreshti `diagonal` i x
    # Te dp[i][j], start[i][j] ruan fillimin e pjesës fundore të A që përputhet.
    prev = list(range(len(y) + 1))
    prev_start = [0] * (len(y) + 1)
    for i in range(1, len(x) + 1):
        curr = [big] * (len(y) + 1)
        curr_start = [0] * (len(y) + 1)
        curr[0] = 0                     # fillimi i A nuk penalizohet
        curr_start[0] = i
        ca = x[i - 1]
        low = max(1, i - diagonal - band)
        high = min(len(y), i - diagonal + band)
        for j in range(low, high + 1):
            # Zgjedhim veprimin më të lirë: përputhje/zëvendësim, shtim ose fshirje.
            cost = prev[j - 1] + (ca != y[j - 1])
            start = prev_start[j - 1]
            if prev[j] + 1 < cost:
                cost = prev[j] + 1
                start = prev_start[j]
            if curr[j - 1] + 1 < cost:
                cost = curr[j - 1] + 1
                start = curr_start[j - 1]
            curr[j] = cost
            curr_start[j] = start
        prev, prev_start = curr, curr_start
    best = None
    for j in range(min_overlap, len(y) + 1):
        if prev[j] >= big:  # qelizë jashtë brezit
            continue
        shift = astart + prev_start[j]
        left_ov = len(a) - shift
        if left_ov < min_overlap or shift < 3:
            continue
        if abs(shift - estimated_shift) > padding + 2:
            continue
        if len(b) - j < 3:  # read-i i dytë duhet të zgjatet në të djathtë
            continue
        overlap = min(left_ov, j)
        identity = 1 - prev[j] / max(left_ov, j)
        if identity < 1 - max_error:
            continue
        rank = (overlap * identity ** 6, identity, overlap)
        if best is None or rank > best[0]:
            best = (rank, shift, j, overlap, identity)
    return best[1:] if best else None

def build_overlap_graph(sequences, k=11, min_overlap=42, max_error=.09,
                        max_candidates=8, max_edges=4, stride=4):
    """Indekson k-merët e reads në të dy orientimet (k-mer hashing, indeks i kundërt).
    Në indeks ruhet vetëm çdo i `stride`-ti k-mer i çdo read-i; kërkimi bëhet me të
    gjithë k-merët e read-it tjetër. Çdo mbivendosje prej >= min_overlap bazash përmban
    gjithsesi disa k-merë të ruajtur, ndërsa memoria dhe votimi bien rreth `stride` herë.
    Kontrollon në detaje vetëm një pjesë të vogël të çifteve kandidate.
    Për çdo read mjaftojnë `max_edges` lidhje të mira (alternativat për DFS-në);
    kandidatët e tjerë nuk krahasohen, që të shmangim punë të panevojshme.
    Grafi është i drejtuar; lidhjet vlerësohen nga gjatësia dhe cilësia e përputhjes."""
    index = defaultdict(list)
    for node, s in enumerate(sequences):
        for p in range(0, len(s) - k + 1, stride):
            index[s[p:p+k]].append((node, p))
        # Heqim k-merët shumë të shpeshtë: një k-mer që shfaqet më shumë se 10 herë mbi
    # mesataren vjen nga një zonë e përsëritur (p.sh. elementet REP te E. coli, me qindra
    # kopje) dhe do të krijonte shumë çifte kandidate të rreme.
    total = 0
    for positions in index.values():
        total += len(positions)
    frequency_cap = max(60, 10 * total // max(1, len(index)))
    edges = defaultdict(list)
    comparisons = 0
    candidates_seen = 0
    for node, s in enumerate(sequences):
        votes = defaultdict(Counter)
        for p in range(len(s) - k + 1):
            occurrences = index.get(s[p:p + k], [])
            if len(occurrences) > frequency_cap:
                continue
            for other, q in occurrences:
                if other // 2 == node // 2:
                    continue
                shift = p - q
                if 3 <= shift <= len(s) - min_overlap + 4 and shift + len(sequences[other]) > len(s) + 2:
                    votes[other][shift] += 1
        candidates_seen += len(votes)
        # Për secilin pasues zgjedhim zhvendosjen me më shumë vota dhe i rendisim kandidatët.
        ranked = []
        for other, counts in votes.items():
            # Shtimet/fshirjet mund të zhvendosin përputhjen 1-2 vende; bashkojmë votat afër.
            best_shift = max(counts, key=lambda off: sum(n for st, n in counts.items() if abs(st-off) <= 2))
            support = sum(n for st, n in counts.items() if abs(st-best_shift) <= 2)
            ranked.append((support, other, best_shift))
        ranked.sort(reverse=True)
        for support, other, shift in ranked[:max_candidates]:
            if len(edges[node]) == max_edges:
                break
            comparisons += 1
            result = approximate_overlap(s, sequences[other], shift, min_overlap, max_error)
            if result:
                real_shift, b_used, ov, identity = result
                edges[node].append(Edge(node, other, real_shift, b_used, ov, identity, support))
        edges[node].sort(key=lambda e: e.score, reverse=True)
    return edges, {'candidate_pairs': candidates_seen, 'alignments': comparisons,
                   'overlaps': sum(len(v) for v in edges.values()),
                   'naive_oriented_comparisons': len(sequences) * (len(sequences) - 2)}

def dfs_backtracking(sequences, edges, tolerance=10, max_steps=None):
    """Gjen rrugën (backbone) në grafin e mbivendosjeve me DFS dhe BACKTRACKING.

    Problemi: duam një rrugë në graf (read -> read që e vazhdon) që jep sekuencën
    më të gjatë, pa përdorur të njëjtin read dy herë (si te rruga Hamiltoniane).
      1. Zgjedhim: nga read-i aktual provojmë me radhë çdo read që e vazhdon
         (lidhjet më të mira provohen të parat).
      2. Shkojmë thellë (DFS) derisa rruga të mos vazhdojë më dhe e krahasojmë
         me rrugën më të mirë të gjetur deri tani.
      3. Kthehemi mbrapa (backtracking): heqim read-in e fundit dhe provojmë
         alternativën tjetër. Ndryshe nga greedy, nëse një zgjedhje të çon në
         rrugë pa dalje, algoritmi nuk ngec, por provon rrugët e tjera.

    Krasitja (pruning): brenda të njëjtit kërkim (nga e njëjta pikë fillimi) çdo read
    vizitohet vetëm një herë, si te DFS-ja klasike, sepse vazhdimet prej tij janë
    provuar tashmë. Kjo e mban kohën lineare dhe pengon që rruga të kalojë dy herë
    nëpër një zonë të përsëritur (gjë që do ta zgjaste ADN-në artificialisht).
    Një read i vizituar nga një kërkim i mëparshëm rivizitohet vetëm nëse tani
    arrijmë me rrugë qartësisht më të gjatë (më shumë se `tolerance` baza; dallimet
    e vogla vijnë nga gabimet e leximit), p.sh. kur fillimi i parë s'ishte fillimi i ADN-së.
    `max_steps` siguron që kërkimi të përfundojë gjithmonë."""
    if max_steps is None:
        max_steps = 100 * len(sequences)  # kufi sigurie, shumë më i madh se nevoja
    # Numërojmë lidhjet hyrëse: reads pa lidhje hyrëse janë fillime të mundshme të ADN-së.
    incoming = [0] * len(sequences)
    for node in edges:
        for e in edges[node]:
            incoming[e.dst] += 1
    starts = []
    for node in range(len(sequences)):
        starts.append((incoming[node], node))
    starts.sort()  # provojmë së pari reads me më pak lidhje hyrëse

    best = {'span': 0, 'path': [], 'links': []}  # rruga më e mirë e gjetur deri tani
    best_at = {}  # nyja -> (gjatësia me të cilën arritëm tek ajo, kërkimi që e vizitoi)
    path = []     # nyjet e rrugës aktuale
    links = []    # lidhjet (Edge) e rrugës aktuale
    used = set()  # reads që janë tashmë në rrugë (në cilindo drejtim)
    stats = {'dfs_steps': 0, 'dfs_backtracks': 0, 'dfs_pruned': 0,
             'dfs_dead_ends': 0, 'dfs_start_nodes': 0, 'dfs_complete': True}

    def dfs(node, span, search):
        stats['dfs_steps'] += 1
        if stats['dfs_steps'] > max_steps:
            stats['dfs_complete'] = False
            return
        extended = False      # a e vazhduam rrugën te ndonjë read tjetër?
        can_continue = False  # a ka të paktën një read të lirë që e vazhdon?
        for e in edges.get(node, []):
            read = e.dst // 2
            if read in used:
                continue  # read-i është tashmë në rrugë
            can_continue = True
            new_span = span + len(sequences[e.dst]) - e.target_consumed
            if e.dst in best_at:
                old_span, old_search = best_at[e.dst]
                if old_search == search or new_span <= old_span + tolerance:
                    stats['dfs_pruned'] += 1
                    continue  # krasitje: këtë vazhdim e kemi provuar më parë
            best_at[e.dst] = (new_span, search)
            extended = True
            # 1) Zgjedhim lidhjen.
            used.add(read)
            path.append(e.dst)
            links.append(e)
            # 2) Vazhdojmë më thellë (DFS).
            dfs(e.dst, new_span, search)
            # 3) Backtracking: e anulojmë zgjedhjen dhe provojmë lidhjen tjetër.
            used.remove(read)
            path.pop()
            links.pop()
            stats['dfs_backtracks'] += 1
        if not can_continue:
            stats['dfs_dead_ends'] += 1  # asnjë read i lirë nuk e vazhdon këtë rrugë
        # Rruga mbaron këtu: e krahasojmë me më të mirën. E kopjojmë vetëm kur është
        # më e gjatë (jo në çdo hap), përndryshe koha do të rritej në katror.
        if not extended and span > best['span']:
            best['span'] = span
            best['path'] = path[:]
            best['links'] = links[:]

    # Rekursioni shkon aq thellë sa reads ka rruga; kufiri standard i Python-it (1000)
    # mund të mos mjaftojë, prandaj e rrisim sipas numrit të reads.
    if sys.getrecursionlimit() < len(sequences) + 200:
        sys.setrecursionlimit(len(sequences) + 200)

    for count, start in starts:
        if start in best_at:
            continue  # ky read është arritur tashmë nga një fillim më i mirë
        if stats['dfs_steps'] > max_steps:
            break
        stats['dfs_start_nodes'] += 1
        best_at[start] = (len(sequences[start]), start)
        used.add(start // 2)
        path.append(start)
        dfs(start, len(sequences[start]), start)  # numri i kërkimit = nyja e fillimit
        used.remove(start // 2)
        path.pop()

    if not best['path']:
        raise ValueError('No reads supplied')
    parts = [sequences[best['path'][0]]]
    for e in best['links']:
        parts.append(sequences[e.dst][e.target_consumed:])
    return ''.join(parts), best['path'], best['links'], stats

def local_alignment(read: str, reference_window: str, diagonal=12, band=10):
    """Përputh një read me një zonë të sekuencës duke përdorur edit distance.
    Skajet e zonës mund të mbeten të papërdorura, pa penalizim.
    Read-i pritet të fillojë rreth kolonës `diagonal` të zonës, prandaj llogarisim
    vetëm qelizat në brezin ±band rreth asaj diagonaleje (banded DP).
    Kthen koston, pozicionet dhe veprimet (veprimi, pozicioni në zonë, pozicioni në read)
    për bazat (B), shtimet (I) dhe fshirjet (D)."""
    n, m = len(read), len(reference_window)
    big = n + m + 1
    prev = [0] * (m + 1)
    back = [bytearray(m + 1) for _ in range(n + 1)]
    for i, base in enumerate(read, 1):
        curr = [big] * (m + 1)
        curr[0] = i
        back[i][0] = 1
        low = max(1, i + diagonal - band)
        high = min(m, i + diagonal + band)
        for j in range(low, high + 1):
            diag = prev[j-1] + (base != reference_window[j-1])
            up = prev[j] + 1  # shtim në read
            left = curr[j-1] + 1  # fshirje në read
            if diag <= up and diag <= left:
                curr[j] = diag
                back[i][j] = 0
            elif up <= left:
                curr[j] = up
                back[i][j] = 1
            else:
                curr[j] = left
                back[i][j] = 2
        prev = curr
    end = prev.index(min(prev))
    cost = prev[end]
    ops = []
    i, j = n, end
    while i > 0:
        operation = back[i][j]
        if operation == 0:
            ops.append(('B', j - 1, i - 1))
            i -= 1
            j -= 1
        elif operation == 1:
            ops.append(('I', j, i - 1))
            i -= 1
        else:
            ops.append(('D', j - 1, None))
            j -= 1
    ops.reverse()
    return cost, j, end, ops

def polish(backbone, reads, k=12, max_error=.15):
    """Rivendos TË GJITHA reads mbi backbone dhe ndërton konsensusin.
    Për çdo read zgjedhim drejtimin (përpara ose mbrapsht) me më shumë k-merë të
    përbashkët dhe e përputhim vetëm një herë.
    Vota e çdo baze peshohet me cilësinë e saj (Phred): bazat me cilësi të ulët,
    që shpesh janë gabime leximi, ndikojnë më pak se bazat e sigurta.
    Kështu korrigjohen gabimet e leximit dhe ruhen mutacionet e vërteta."""
    n = len(backbone)
    idx = defaultdict(list)
    for j in range(n - k + 1):
        idx[backbone[j:j + k]].append(j)
    # Për çdo pozicion: shuma e cilësive për secilën bazë (4 lista të thjeshta, jo një
    # Counter për pozicion, që të kursejmë memorie te gjenomet e mëdha).
    weight = {'A': [0] * n, 'C': [0] * n, 'G': [0] * n, 'T': [0] * n}
    cover = [0] * n          # sa reads kanë një bazë në këtë pozicion
    del_votes = [0] * n
    insert_votes = {}        # hapësira -> Counter me sekuencat e shtuara
    gap_cover = [0] * (n + 1)
    mapped = 0
    edits = []
    for read, quality in reads:
        if quality is None:  # FASTA: të gjitha bazat kanë të njëjtën peshë
            quality = chr(DEFAULT_QUALITY + 33) * len(read)
        # Fillojmë pranë pozicionit ku k-merët identikë japin provën më të fortë.
        best_read, best_quality, best_offset, best_support = None, None, 0, 0
        for orientation, orientation_quality in ((read, quality),
                                                 (reverse_complement(read), quality[::-1])):
            positions = Counter()
            for p in range(len(orientation) - k + 1):
                for j in idx.get(orientation[p:p + k], [])[:15]:
                    positions[j-p] += 1
            if not positions:
                continue
            offset, support = positions.most_common(1)[0]
            if support > best_support:
                best_read, best_quality = orientation, orientation_quality
                best_offset, best_support = offset, support
        if best_read is None:
            continue
        left = max(0, best_offset - 12)
        right = min(n, best_offset + len(best_read) + 12)
        if right <= left:
            continue
        cost, rel_start, rel_end, ops = local_alignment(best_read, backbone[left:right],
                                                        best_offset - left)
        if cost > max(3, round(len(best_read) * max_error)):
            continue
        start, end = left + rel_start, left + rel_end
        mapped += 1
        edits.append(cost)
        ins_per_gap = {}
        for action, local_j, read_i in ops:
            pos = left + local_j
            if action == 'B':
                weight[best_read[read_i]][pos] += ord(best_quality[read_i]) - 33
                cover[pos] += 1
            elif action == 'D':
                del_votes[pos] += 1
            else:
                ins_per_gap[pos] = ins_per_gap.get(pos, '') + best_read[read_i]
        for gap in range(max(0, start + 1), min(n, end) + 1):
            gap_cover[gap] += 1
        for gap in ins_per_gap:
            if gap not in insert_votes:
                insert_votes[gap] = Counter()
            insert_votes[gap][ins_per_gap[gap]] += 1
    output = []
    corrections = 0
    for i in range(n + 1):
        # Reads që kalojnë mbi hapësirën i mund të mbështesin një sekuencë të shtuar.
        if i in insert_votes:
            inserted, count = insert_votes[i].most_common(1)[0]
            if count > gap_cover[i] / 2 and (i < n or count > 1):
                output.append(inserted)
                corrections += len(inserted)
        if i == n:
            break
        if del_votes[i] > cover[i] and del_votes[i] > 1:
            corrections += 1
            continue
        old = backbone[i]
        consensus = old
        best_weight = 0
        for base in DNA:  # baza me peshën (cilësinë) më të madhe fiton
            if weight[base][i] > best_weight:
                consensus = base
                best_weight = weight[base][i]
        output.append(consensus)
        if consensus != old:
            corrections += 1
    return ''.join(output), {'reads_mapped': mapped, 'consensus_changes': corrections,
                             'mean_alignment_edits': round(statistics.mean(edits), 3) if edits else None}
