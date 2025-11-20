#!/usr/bin/env python3
from __future__ import annotations
import os, sys, time, random, argparse
from collections import defaultdict, Counter
from typing import List
from bisect import bisect_left
import itertools

try:
    from delta_debugging.DD import DD
except ImportError:
    print("Could not import delta_debugging.DD. Check PYTHONPATH.")
    sys.exit(1)

# ==================================================
# CONFIG
# ==================================================
MAX_DEPTH = 2
BEAM_W = 2000
TOP_PATHS_FOR_SURVIVAL = 1000
USE_LCS_UB = True
PREFIX_LEN = 600
PAIR_SAMPLE = 1000
random.seed(0)

# ==================================================
# HELPERS
# ==================================================
def is_valid_subsequence(subseq: str, seq: str) -> bool:
    it = iter(seq)
    return all(ch in it for ch in subseq)

def validate_mlcs(candidate: str, sequences: list[str]) -> bool:
    return all(is_valid_subsequence(candidate, s) for s in sequences)

# ==================================================
# APPROX-LCS PRIMARY RANKING
# ==================================================
def build_suffix_lcs_table(a: str, b: str):
    n, m = len(a), len(b)
    L = [[0]*(m+1) for _ in range(n+1)]
    for i in range(n-1, -1, -1):
        for j in range(m-1, -1, -1):
            L[i][j] = L[i+1][j+1]+1 if a[i]==b[j] else max(L[i+1][j], L[i][j+1])
    return L

def approximate_lcs_primary(sequences: List[str], prefix_len=PREFIX_LEN, pair_sample=PAIR_SAMPLE):
    n = len(sequences)
    scores = [0.0]*n
    seqs = [s[:prefix_len] for s in sequences]

    pairs = list(itertools.combinations(range(n), 2))
    if len(pairs) > pair_sample:
        pairs = random.sample(pairs, pair_sample)

    for (i, j) in pairs:
        s1, s2 = seqs[i], seqs[j]
        L = build_suffix_lcs_table(s1, s2)
        lcs_len = L[0][0]
        sim = lcs_len / min(len(s1), len(s2))
        scores[i] += sim
        scores[j] += sim

    avg_scores = [s / (n - 1) for s in scores]
    ranked = sorted(range(n), key=lambda i: -avg_scores[i])
    return ranked, avg_scores

# ==================================================
# BEAM-GUIDED DELTA REORDERING
# ==================================================
def _next_ge(lst, lo):
    i = bisect_left(lst, lo)
    return lst[i] if i < len(lst) else -1

def _build_pos_lists(strings: List[str]):
    pos_all = []
    for s in strings:
        d = defaultdict(list)
        for j, ch in enumerate(s):
            d[ch].append(j)
        pos_all.append(d)
    return pos_all

def build_all_suffix_lcs_tables(primary: str, others: List[str]):
    return [build_suffix_lcs_table(primary, s) for s in others]

def beam_rank_deltas(primary: str, strings: List[str]) -> List[int]:
    pidx = strings.index(primary)
    others = [s for i, s in enumerate(strings) if i != pidx]
    if not others:
        return list(range(len(primary)))

    tables = build_all_suffix_lcs_tables(primary, others) if USE_LCS_UB else None
    pos_all = _build_pos_lists(others)
    lens_b = [len(s) for s in others]
    n = len(primary)

    def ub_remaining(i: int, curs) -> int:
        if not tables: return 0
        vals = []
        for k, T in enumerate(tables):
            j = curs[k] + 1
            if j <= lens_b[k]:
                vals.append(T[i][j])
        return int(sum(vals)/len(vals)*1.2) if vals else 0

    start_cursors = tuple([-1]*len(others))
    beam = [(0, ub_remaining(0, start_cursors), 0, start_cursors, [])]
    best_ln, best_len_at = 0, {}
    surv = Counter()

    while beam:
        nxt = []
        for (ln, f, i, curs, picked) in beam:
            best_ln = max(best_ln, ln)
            if i >= n:
                continue

            ch = primary[i]
            ok, newc = True, []
            for sj, dpos in enumerate(pos_all):
                np = _next_ge(dpos.get(ch, []), curs[sj]+1)
                if np == -1:
                    ok = False
                    break
                newc.append(np)

            if ok:
                i2, ln2 = i+1, ln+1
                curs2 = tuple(newc)
                f2 = ln2 + ub_remaining(i2, curs2)
                if f2 >= best_ln - 10 and best_len_at.get((i2, curs2), -1) < ln2:
                    best_len_at[(i2, curs2)] = ln2
                    nxt.append((ln2, f2, i2, curs2, picked+[i]))

            i1 = i+1
            f1 = ln + ub_remaining(i1, curs)
            if f1 >= best_ln - 10 and best_len_at.get((i1, curs), -1) < ln:
                best_len_at[(i1, curs)] = ln
                nxt.append((ln, f1, i1, curs, picked))

        if not nxt:
            break

        nxt.sort(key=lambda s:(-s[1], -s[0], s[2], tuple(s[3]), random.random()))
        beam = nxt[:BEAM_W]

        if all(s[2] >= n for s in beam):
            break

    for (_, _, _, _, picked) in beam[:min(TOP_PATHS_FOR_SURVIVAL, len(beam))]:
        for i in picked:
            surv[i] += 1

    for k in surv:
        surv[k] = int(surv[k] * 5)

    additions = list(range(len(primary)))
    additions.sort(key=lambda i: (-surv.get(i, 0), i))
    return additions

class MLCS_Addition(DD):
    def __init__(self, strings: list[str], primary_string=None):
        super().__init__()
        self.strings = strings
        self.primary_string = primary_string or min(strings, key=len)

    @staticmethod
    def select(string: str, indices: list[int]) -> list[str]:
        return [string[i] for i in sorted(indices)]

    def _test(self, additions: list[int]):
        subseq = ''.join(self.select(self.primary_string, additions))
        for s in self.strings:
            it = iter(s)
            if not all(ch in it for ch in subseq):
                return self.PASS
        return self.FAIL

    def __listminus(self, c1, c2):
        return [x for x in c1 if x not in set(c2)]

    def ddmin_add(self, c, r=list(), *, _depth=0):
        set_c, list_c = set(c), list(c)
        solns = []
        while self.test(list(set(c)-set_c)+r) == self.FAIL:
            soln = set(self.vanilla_add_dd(list_c, n=2, r=list(set(c)-set_c)+r))
            soln |= (set(c)-set_c)
            set_c &= soln
            list_c = list(set_c)
            solns.append(soln)
        if _depth < MAX_DEPTH:
            all_soln = set(c).intersection(*solns) if solns else set()
            set_c = set(c) - all_soln
            for soln in solns:
                assert self.test(list(soln & set_c)+r) == self.FAIL
                irreplaceable = self.vanilla_add_dd(list(set(c)-soln), n=2, r=list(soln & set_c)+r)
                if not irreplaceable: continue
                irreplaceable = list((soln & set_c)|set(irreplaceable))
                result = self.ddmin_add(all_soln, irreplaceable+r, _depth=_depth+1)
                if len(result) > len(soln)-len(irreplaceable):
                    soln = set(result + irreplaceable)
            solns.append(soln)
        if not solns:
            return []
        return list(max(solns, key=len))

    def vanilla_add_dd(self, c, n=2, r=list()):
        cur_soln = []
        cbar_offset = 0
        while True:
            if n > len(c)-len(cur_soln): return cur_soln
            cs = self.split(self.__listminus(c, cur_soln), n)
            c_failed = cbar_failed = False
            next_n = n
            for i in range(n):
                subset = cur_soln + cs[i]
                if self.test(subset+r) == self.FAIL:
                    c_failed=True; cur_soln=subset; next_n=2; break
            if not c_failed:
                for j in range(n):
                    i=int((j+cbar_offset)%n)
                    complement=self.__listminus(c, cs[i])
                    if self.test(complement+r) == self.FAIL:
                        cbar_failed=True; cur_soln=complement; next_n-=1; cbar_offset=i; break
            if not c_failed and not cbar_failed:
                if n>=len(c)-len(cur_soln): return cur_soln
                next_n=min(len(c)-len(cur_soln), n*2)
            n=next_n


def load_sequences_from_file(path):
    return [line.strip() for line in open(path) if line.strip()]

def dd_MLCS_addition(sequences: List[str]):
    ranked, scores = approximate_lcs_primary(sequences)

    print(f"\nLoaded {len(sequences)} sequences\n")

    # ---- SHOW TOP 5 ----
    print(" Top-5 predicted primaries:")
    top5 = ranked[:5]
    for r, i in enumerate(top5, 1):
        print(f"  {r}. Seq #{i+1} | score={scores[i]:.4f}")
    print()

    # ---- RUN ARP ON TOP-5 ----
    print("Running ARP on Top-5 primaries:\n")

    results = []
    for rank, i in enumerate(top5, 1):
        s = sequences[i]
        print(f"[Run {rank}] Primary #{i+1} | len={len(s)} | score={scores[i]:.4f}")

        t0 = time.time()
        deltas = beam_rank_deltas(s, sequences)
        dd = MLCS_Addition(sequences, primary_string=s)
        result = dd.ddmin_add(deltas)

        mlcs = ''.join(MLCS_Addition.select(s, result))
        valid = validate_mlcs(mlcs, sequences)
        elapsed = round(time.time() - t0, 2)

        print(f"MLCS Len = {len(mlcs)}, Valid={valid}, Runtime={elapsed}s\n")

        results.append((len(mlcs), i, mlcs))

    # ---- PICK BEST ----
    best = max(results, key=lambda x: x[0])

    print("==============================")
    print(" BEST MLCS FOUND:")
    print(f" Length: {best[0]}")
    print(f" Primary used: #{best[1]+1}")
    print(f" MLCS: {best[2][:100]}...")
    print("==============================")

    return best[2]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to sequence file")
    args = parser.parse_args()

    sequences = load_sequences_from_file(args.input)
    start = time.time()
    final_mlcs = dd_MLCS_addition(sequences)
    end = time.time()

    print(f"\nTotal Runtime: {end-start:.2f}s")