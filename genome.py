# PROJEKTI: Rindërtimi i ADN-së nga fragmente me gabime
#
# Rrjedha kryesore e programit:
#   1. Simulon ADN-në dhe krijon reads me gabime (ose lexon një skedar FASTA).
#   2. Përdor k-mer hashing për të gjetur çiftet që mund të mbivendosen.
#   3. Kontrollon mbivendosjet me edit distance dhe ndërton grafin.
#   4. Me metodën greedy krijon sekuencën fillestare (backbone).
#   5. Me konsensusin e reads korrigjon gabimet.
#   6. Krahason rezultatin me ADN-në e njohur vetëm për verifikim.
#


from __future__ import annotations
from collections import Counter, defaultdict
from dataclasses import dataclass
from difflib import SequenceMatcher
import argparse
import json
import random
import statistics
import time
import os
import threading
from pathlib import Path

DNA = 'ACGT'
COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def reverse_complement(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def write_fasta(path: Path, records):
    with path.open('w') as f:
        for name, seq in records:
            f.write(f'>{name}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + '\n')


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


def simulate(genome_length=1800, n_reads=350, read_length=115, error_rate=.01,
             indel_rate=.0015, seed=7):
    """Krijon një ADN me mutacione dhe pastaj reads me gabime të rastësishme.
    Pozicionet e vërteta dhe ADN-ja origjinale nuk i jepen algoritmit të rindërtimit."""
    rng = random.Random(seed)
    reference = ''.join(rng.choices(DNA, k=genome_length))
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
        raise ValueError('Choose genome length >= 200')
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
    # Sigurojmë mbulim në të dy skajet, që ADN-ja të rindërtohet plotësisht.
    max_start = len(sample) - read_length
    if max_start < 0:
        raise ValueError('Read length must be shorter than genome length')
    starts = [0, max_start] + [rng.randint(0, max_start) for _ in range(max(0, n_reads - 2))]
    records = []
    for ix, start in enumerate(starts):
        seq = sample[start:start + read_length]
        noisy = []
        for base in seq:
            r = rng.random()
            if r < indel_rate / 2:  # fshirje e shkaktuar nga gabimi i leximit
                continue
            if r < indel_rate:      # shtim i shkaktuar nga gabimi i leximit
                noisy.append(rng.choice(DNA))
            if rng.random() < error_rate:
                base = rng.choice([b for b in DNA if b != base])
            noisy.append(base)
        seq = ''.join(noisy)
        if rng.random() < .5:
            seq = reverse_complement(seq)
        records.append((f'read_{ix+1:04d}', seq))
    rng.shuffle(records)
    return reference, sample, records, true_variants
@dataclass(frozen=True)
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
        return self.overlap * (self.identity ** 5) + min(self.votes, 10) * .3


def oriented_reads(reads):
    return [s for seq in reads for s in (seq, reverse_complement(seq))]


def approximate_overlap(a: str, b: str, estimated_shift: int,
                        min_overlap: int, max_error: float, padding=7):
    """Gjen mbivendosjen e përafërt me programim dinamik (edit distance).
    Fillimi i A dhe fundi i B mund të mbeten jashtë përputhjes pa penalizim.
    Votat e k-merëve e kufizojnë zonën ku kërkohet mbivendosja.
    Kthen (zhvendosja, baza të përdorura nga B, mbivendosja, ngjashmëria) ose None."""
    astart = max(0, estimated_shift - padding)
    x = a[astart:]
    # Read-i B duhet të vazhdojë përtej A; kontrollojmë vetëm pjesën që mund të mbivendoset.
    approx = len(a) - estimated_shift
    y = b[:min(len(b), approx + padding * 2 + 4)]
    if len(x) < min_overlap or len(y) < min_overlap:
        return None
    # Te dp[i][j], start[i][j] ruan fillimin e pjesës fundore të A që përputhet.
    prev = list(range(len(y) + 1))
    prev_start = [0] * (len(y) + 1)
    for i, ca in enumerate(x, 1):
        curr = [0] * (len(y) + 1)   # fillimi i A nuk penalizohet
        curr_start = [i] + [0] * len(y)
        for j, cb in enumerate(y, 1):
            choices = ((prev[j - 1] + (ca != cb), prev_start[j - 1]),
                       (prev[j] + 1, prev_start[j]),
                       (curr[j - 1] + 1, curr_start[j - 1]))
            cost, start = min(choices, key=lambda z: z[0])
            curr[j], curr_start[j] = cost, start
        prev, prev_start = curr, curr_start
    best = None
    for j in range(min_overlap, len(y) + 1):
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
                        max_candidates=8):
    """Indekson k-merët e reads në të dy orientimet.
    Kontrollon në detaje vetëm një pjesë të vogël të çifteve kandidate.
    Grafi është i drejtuar; lidhjet vlerësohen nga gjatësia dhe cilësia e përputhjes."""
    index = defaultdict(list)
    for node, s in enumerate(sequences):
        for p in range(len(s) - k + 1):
            index[s[p:p+k]].append((node, p))
    # Heqim k-merët shumë të shpeshtë (zona të përsëritura) për të shmangur shumë krahasime.
    frequency_cap = max(60, len(sequences) // 5)
    edges = defaultdict(list)
    comparisons = 0
    candidates_seen = 0
    for node, s in enumerate(sequences):
        votes = defaultdict(Counter)
        for p in range(len(s) - k + 1):
            occurrences = index[s[p:p + k]]
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
            comparisons += 1
            result = approximate_overlap(s, sequences[other], shift, min_overlap, max_error)
            if result:
                real_shift, b_used, ov, identity = result
                edges[node].append(Edge(node, other, real_shift, b_used, ov, identity, support))
        edges[node].sort(key=lambda e: e.score, reverse=True)
    return edges, {'candidate_pairs': candidates_seen, 'alignments': comparisons,
                   'overlaps': sum(len(v) for v in edges.values()),
                   'naive_oriented_comparisons': len(sequences) * (len(sequences) - 2)}
