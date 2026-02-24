#!/usr/bin/env python3
from __future__ import annotations
import time
import sys
from delta_debugging.DD import DD

# ==================================================
# CONFIG
# ==================================================
def is_valid_subsequence(subseq: str, seq: str) -> bool:
    it = iter(seq)
    return all(ch in it for ch in subseq)

def validate_mlcs(candidate: str, sequences: list[str]) -> bool:
    return all(is_valid_subsequence(candidate, s) for s in sequences)

# -------------------------------------------------------------
class MLCS_Addition(DD):
    """
    Pure Vanilla Addition DD (NO ddmin_add).
    """

    def __init__(self, strings: list[str], primary_string=None):
        super().__init__()
        self.strings = strings

        if primary_string is None:
            # default = shortest intersection string
            self.primary_string = min(self.strings, key=len)
            charset = set(self.primary_string)
            for s in self.strings:
                charset &= set(s)
            self.primary_string = "".join([c for c in self.primary_string if c in charset])
        else:
            self.primary_string = primary_string

    @staticmethod
    def select(string: str, indices: list[int]) -> list[str]:
        return [string[i] for i in sorted(indices)]

    def _test(self, additions: list[int]):
        """FAIL = subsequence is common, PASS = not common."""
        subseq = ''.join(self.select(self.primary_string, additions))
        for s in self.strings:
            it = iter(s)
            if not all(ch in it for ch in subseq):
                return self.PASS
        return self.FAIL

    def __listminus(self, c1, c2):
        return [x for x in c1 if x not in set(c2)]

    # -------------------------------------------------------------
    # ★ ONLY THIS FUNCTION IS USED (vanilla-add algorithm)
    # -------------------------------------------------------------
    def vanilla_add_dd(self, c, n=2, r=list()):
        """
        Pure addition DD:
        Start with empty solution and try adding subsets until FAIL holds.
        """
        cur_soln = []
        cbar_offset = 0
        c = sorted(c)

        while True:
            if n > len(c) - len(cur_soln):
                return cur_soln

            cs = self.split(self.__listminus(c, cur_soln), n)
            c_failed = False
            next_n = n

            # Try subsets
            for i in range(n):
                subset = cur_soln + cs[i]
                if self.test(subset + r) == self.FAIL:
                    c_failed = True
                    cur_soln = subset
                    next_n = 2
                    break

            # Try complements
            if not c_failed:
                for j in range(n):
                    i = int((j + cbar_offset) % n)
                    complement = self.__listminus(c, cs[i])
                    if self.test(complement + r) == self.FAIL:
                        cur_soln = complement
                        next_n -= 1
                        cbar_offset = i
                        break

            # If no progress, double granularity
            if n >= len(c) - len(cur_soln):
                return cur_soln
            n = min(len(c) - len(cur_soln), next_n * 2)

# -------------------------------------------------------------
def dd_MLCS_addition(strings: list[str]):
    print("\nRunning ADDITION-BASED MLCS...\n")

    best_len = 0
    best_mlcs = ""
    best_primary = -1

    for i, primary in enumerate(strings):
        print(f"[Primary {i+1}/{len(strings)}] len={len(primary)}")

        dd = MLCS_Addition(strings, primary_string=primary)
        deltas = list(range(len(primary)))

        start = time.time()
        result = dd.vanilla_add_dd(deltas)
        end = time.time()

        mlcs = ''.join(MLCS_Addition.select(primary, result))
        mlcs_len = len(mlcs)
        ok = validate_mlcs(mlcs, strings)

        print(f"  MLCS Len = {mlcs_len},  Valid = {ok},  Runtime = {end-start:.2f}s")

        if ok and mlcs_len > best_len:
            best_len = mlcs_len
            best_mlcs = mlcs
            best_primary = i

    print("\n==============================")
    print(" BEST MLCS FOUND:")
    print(f" Length: {best_len}")
    print(f" Primary used: #{best_primary+1}")
    print(f" MLCS: {best_mlcs[:100]}{'...' if len(best_mlcs)>100 else ''}")
    print("==============================")

    return best_mlcs

# -------------------------------------------------------------
def load_sequences_from_file(path: str) -> list[str]:
    #with open(path) as f:
    with open(path, "r", encoding="latin-1") as f:
        return [line.strip() for line in f if line.strip()]

# -------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Vanilla Addition-based MLCS")
    parser.add_argument("--input", required=True, help="Path to sequence file")
    args = parser.parse_args()

    sequences = load_sequences_from_file(args.input)
    print(f"Loaded {len(sequences)} sequences.")

    start = time.time()
    dd_MLCS_addition(sequences)
    end = time.time()

    print(f"\nTotal Runtime: {end-start:.2f}s")
