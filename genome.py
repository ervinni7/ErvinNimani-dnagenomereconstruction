
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