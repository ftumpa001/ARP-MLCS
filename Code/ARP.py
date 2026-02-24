#!/usr/bin/env python3

from __future__ import annotations
import time
import argparse
import resource
from delta_debugging.DD import DD

MAX_DEPTH = 2   # default value


def is_valid_subsequence(subseq: str, seq: str) -> bool:
    it = iter(seq)
    return all(ch in it for ch in subseq)


def validate_mlcs(candidate: str, sequences: list[str]) -> bool:
    return all(is_valid_subsequence(candidate, s) for s in sequences)


class MLCS_Addition(DD):

    def __init__(self, strings: list[str], primary_string=None):
        super().__init__()
        self.strings = strings

        if primary_string is None:
            base = min(strings, key=len)
            charset = set(base)
            for s in strings:
                charset &= set(s)
            self.primary_string = "".join([c for c in base if c in charset])
        else:
            self.primary_string = primary_string

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
        global MAX_DEPTH

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
                assert self.test(list(soln & set_c) + r) == self.FAIL

                irreplaceable = self.vanilla_add_dd(
                    list(set(c) - soln), n=2,
                    r=list(soln & set_c) + r
                )
                if not irreplaceable:
                    continue

                irreplaceable = list((soln & set_c) | set(irreplaceable))
                result = self.ddmin_add(all_soln, irreplaceable + r, _depth=_depth + 1)

                if len(result) > len(soln) - len(irreplaceable):
                    soln = set(result + irreplaceable)

            solns.append(soln)

        if not solns:
            return []

        return list(max(solns, key=len))

    def vanilla_add_dd(self, c, n=2, r=list()):
        cur_soln = []
        cbar_offset = 0
        c.sort()

        while True:
            if n > len(c) - len(cur_soln):
                return cur_soln

            cs = self.split(self.__listminus(c, cur_soln), n)
            c_failed = False
            cbar_failed = False
            next_n = n

            # subsets
            for i in range(n):
                subset = cur_soln + cs[i]
                if self.test(subset + r) == self.FAIL:
                    c_failed = True
                    cur_soln = subset
                    next_n = 2
                    break

            # complements
            if not c_failed:
                for j in range(n):
                    i = int((j + cbar_offset) % n)
                    complement = self.__listminus(c, cs[i])
                    if self.test(complement + r) == self.FAIL:
                        cbar_failed = True
                        cur_soln = complement
                        next_n -= 1
                        break

            # increase n
            if not c_failed and not cbar_failed:
                if n >= len(c) - len(cur_soln):
                    return cur_soln
                next_n = min(len(c) - len(cur_soln), n * 2)

            n = next_n

def dd_MLCS_addition(strings: list[str]):
    best_len = 0
    best_mlcs = ""
    best_primary = -1

    print("\n=== Running ARP ===\n")

    for idx, primary in enumerate(strings):
        print(f"Primary {idx+1}/{len(strings)} (len={len(primary)})")

        dd = MLCS_Addition(strings, primary_string=primary)
        deltas = list(range(len(primary)))

        start = time.time()
        result = dd.ddmin_add(deltas)
        end = time.time()

        mlcs = ''.join(MLCS_Addition.select(primary, result))
        length = len(mlcs)
        print(f"  MLCS = {length}, time={end-start:.2f}s")

        if length > best_len:
            best_len = length
            best_mlcs = mlcs
            best_primary = idx

    print("\n=== Best MLCS ===")
    print(best_mlcs)
    print(f"Length = {best_len}")
    print(f"Primary = {best_primary + 1}")

    return best_mlcs

# ============================================================
# Load sequences
# ============================================================

def load_sequences_from_file(filepath: str) -> list[str]:
    #with open(filepath, "r") as f:
    with open(filepath, "r", encoding="latin-1") as f:
        return [line.strip() for line in f if line.strip()]

# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to sequence file")
    parser.add_argument("--maxdepth", type=int, default=2,
                        help="Maximum recursion depth for ddmin_add")
    args = parser.parse_args()

   
    MAX_DEPTH = args.maxdepth

    seqs = load_sequences_from_file(args.input)
    print(f"Loaded {len(seqs)} sequences.")

    start = time.time()
    final_mlcs = dd_MLCS_addition(seqs)
    end = time.time()

    print(f"\nFinal Runtime = {end-start:.2f}s")
