#!/usr/bin/env python3
from __future__ import annotations
import os, sys, time, random, resource
from collections import defaultdict, Counter
from typing import List
from bisect import bisect_left
import argparse

try:
    from delta_debugging.DD import DD
except ImportError:
    print("Could not import delta_debugging.DD. Check PYTHONPATH.")
    sys.exit(1)

MAX_DEPTH = 2
BEAM_W = 2000          
TOP_PATHS_FOR_SURVIVAL = 1000
USE_LCS_UB = True
random.seed(0)


def is_valid_subsequence(subseq: str, seq: str) -> bool:
    it = iter(seq)
    return all(ch in it for ch in subseq)

def validate_mlcs(candidate: str, sequences: list[str]) -> bool:
    return all(is_valid_subsequence(candidate, s) for s in sequences)

def build_subseq(primary: str, idxs: List[int]) -> str:
    return "".join(primary[i] for i in sorted(idxs))

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


def build_suffix_lcs_table(a: str, b: str):
    n, m = len(a), len(b)
    L = [[0]*(m+1) for _ in range(n+1)]
    for i in range(n-1, -1, -1):
        for j in range(m-1, -1, -1):
            if a[i] == b[j]:
                L[i][j] = L[i+1][j+1] + 1
            else:
                L[i][j] = max(L[i+1][j], L[i][j+1])
    return L

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
        if not tables:
            return 0
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
    beam_seed_counter = Counter()

    while beam:
        nxt = []
        for (ln, f, i, curs, picked) in beam:
            best_ln = max(best_ln, ln)
            if i >= n:
                continue

            ch = primary[i]

            # Try including character
            ok = True
            newc = []
            for sj, dpos in enumerate(pos_all):
                np = _next_ge(dpos.get(ch, []), curs[sj] + 1)
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

            # Try skipping character
            i1 = i + 1
            f1 = ln + ub_remaining(i1, curs)
            if f1 >= best_ln - 10 and best_len_at.get((i1, curs), -1) < ln:
                best_len_at[(i1, curs)] = ln
                nxt.append((ln, f1, i1, curs, picked))

        if not nxt:
            break

        nxt.sort(key=lambda s: (-s[1], -s[0], s[2], tuple(s[3]), random.random()))
        beam = nxt[:BEAM_W]

        if all(s[2] >= n for s in beam):
            break

    # Accumulate survival scores
    for (_, _, _, _, picked) in beam[:min(TOP_PATHS_FOR_SURVIVAL, len(beam))]:
        for i in picked:
            surv[i] += 1

    for k in list(surv.keys()):
        surv[k] = int(surv[k] * 5)

    additions = list(range(len(primary)))
    additions.sort(key=lambda i: (-surv.get(i, 0), i))
    return additions


class MLCS_Addition(DD):
    def __init__(self, strings: list[str], primary_string=None):
        super().__init__()
        self.strings = strings

        if primary_string is None:
            self.primary_string = min(strings, key=len)
            charset = set(self.primary_string)
            for s in strings:
                charset &= set(s)
            self.primary_string = "".join([c for c in self.primary_string if c in charset])
        else:
            self.primary_string = primary_string

    @staticmethod
    def select(string: str, indices: list[int]) -> list[str]:
        return [string[i] for i in sorted(indices)]

    def _test(self, additions: list[int]):
        subseq = "".join(self.select(self.primary_string, additions))
        for s in self.strings:
            it = iter(s)
            if not all(ch in it for ch in subseq):
                return self.PASS
        return self.FAIL

    def __listminus(self, c1, c2):
        return [x for x in c1 if x not in set(c2)]

   
    def ddmin_add(self, c, r=list(), *, _depth=0):
        set_c = set(c)
        list_c = list(c)
        solns = []

        
        while self.test(list(set(c) - set_c) + r) == self.FAIL:
            soln = set(self.vanilla_add_dd(list_c, n=2, r=list(set(c) - set_c) + r))
            soln |= (set(c) - set_c)
            set_c &= soln
            list_c = list(set_c)
            solns.append(soln)

       
        if _depth < MAX_DEPTH:
            all_soln = set(c).intersection(*solns) if solns else set()
            set_c = set(c) - all_soln

            for soln in solns:
                assert self.test(list(soln & set_c)+r) == self.FAIL

                irreplaceable = self.vanilla_add_dd(
                    list(set(c) - soln), n=2, r=list(soln & set_c) + r
                )

                if not irreplaceable:
                    continue

                irreplaceable = list((soln & set_c) | set(irreplaceable))
                result = self.ddmin_add(all_soln, irreplaceable+r, _depth=_depth+1)

                if len(result) > len(soln) - len(irreplaceable):
                    soln = set(result + irreplaceable)

            solns.append(soln)

        if not solns:
            return []
        return list(max(solns, key=len))

    
    def vanilla_add_dd(self, c, n=2, r=list()):
        cur_soln = []
        cbar_offset = 0

        while True:
            if n > len(c) - len(cur_soln):
                return cur_soln

            cs = self.split(self.__listminus(c, cur_soln), n)
            c_failed = False
            cbar_failed = False
            next_n = n

            # try subsets
            for i in range(n):
                subset = cur_soln + cs[i]
                if self.test(subset+r) == self.FAIL:
                    c_failed = True
                    cur_soln = subset
                    next_n = 2
                    break

            if not c_failed:
                for j in range(n):
                    i = (j + cbar_offset) % n
                    complement = self.__listminus(c, cs[i])
                    if self.test(complement+r) == self.FAIL:
                        cbar_failed = True
                        cur_soln = complement
                        next_n -= 1
                        cbar_offset = i
                        break

            if not c_failed and not cbar_failed:
                if n >= len(c) - len(cur_soln):
                    return cur_soln
                next_n = min(len(c)-len(cur_soln), n*2)

            n = next_n


def load_sequences_from_file(path: str) -> List[str]:
    #with open(path) as f:
    with open(path, "r", encoding="latin-1") as f:
        return [line.strip() for line in f if line.strip()]


def dd_MLCS_addition(sequences: List[str]):
    results = []
    for i, s in enumerate(sequences):
        print(f"\n[Run] Primary {i+1}/{len(sequences)} | len={len(s)}", flush=True)

        deltas = beam_rank_deltas(s, sequences)
        dd = MLCS_Addition(sequences, primary_string=s)
        result = dd.ddmin_add(deltas)

        mlcs = "".join(MLCS_Addition.select(s, result))
        ok = validate_mlcs(mlcs, sequences)

        print(f"[Run] Primary {i+1} done | MLCS len={len(mlcs)} | valid={ok}\n")
        results.append((len(mlcs), mlcs))

    best = max(results, key=lambda x: x[0])

    print("\n===== Best MLCS Found =====")
    print(f"Length: {best[0]}")
    print(f"MLCS: {best[1][:100]}{'...' if len(best[1])>100 else ''}")
    return best[1]

"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to input sequence file")
    parser.add_argument("--beam", type=int, default=2000,
                        help="Beam width for beam-guided ordering (default=2000)")

    args = parser.parse_args()

    # overwrite config
    BEAM_W = args.beam

    sequences = load_sequences_from_file(args.input)
    print(f" Loaded {len(sequences)} sequences.")
    print(f" Lengths: {[len(s) for s in sequences]}")

    start = time.time()
    final_mlcs = dd_MLCS_addition(sequences)
    end = time.time()

    print(f"\n Final Total Runtime: {end - start:.4f} seconds")
"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to input sequence file")
    parser.add_argument("--beam", type=int, default=2000,
                        help="Beam width for beam-guided ordering (default=2000)")
    args = parser.parse_args()

    
    BEAM_W = args.beam

    sequences = load_sequences_from_file(args.input)
    print(f" Loaded {len(sequences)} sequences.")
    print(f" Lengths: {[len(s) for s in sequences]}")

    start = time.time()

    # ============================================================
    # *** ONLY RUN FOR PRIMARY SEQUENCE 1 ***
    # ============================================================
    primary = sequences[0]                      # primary string #1
    print(f"\n[Run] Primary 1/1 | len={len(primary)}", flush=True)

    deltas = beam_rank_deltas(primary, sequences)
    dd = MLCS_Addition(sequences, primary_string=primary)
    result = dd.ddmin_add(deltas)

    mlcs = "".join(MLCS_Addition.select(primary, result))
    ok = validate_mlcs(mlcs, sequences)

    print(f"\n[Run] Primary 1 done | MLCS len={len(mlcs)} | valid={ok}")
    print(f"MLCS: {mlcs[:120]}{'...' if len(mlcs)>120 else ''}")

    end = time.time()
    print(f"\n Final Total Runtime: {end - start:.4f} seconds")
