#!/usr/bin/env python3

from __future__ import annotations
import time
import argparse
import resource
from delta_debugging.DD import DD


MAX_DEPTH = 2   # recursion depth controlling refinement


# ============================================================
# Utility Functions
# ============================================================

def is_valid_subsequence(subseq: str, seq: str) -> bool:
    """Return True if subseq is a subsequence of seq."""
    it = iter(seq)
    return all(ch in it for ch in subseq)


def validate_mlcs(candidate: str, sequences: list[str]) -> bool:
    """Return True if candidate is a common subsequence for all strings."""
    return all(is_valid_subsequence(candidate, s) for s in sequences)


# ============================================================
# MLCS Addition
# ============================================================

class MLCS_Addition(DD):
    """
    Addition-based Delta Debugging for MLCS.

    Starts from an empty subsequence and *adds characters* from the primary
    string, using recursive refinement to maintain validity across all strings.
    """

    def __init__(self, strings: list[str], primary_string=None):
        super().__init__()
        self.strings = strings

        # If primary not given → choose the shortest, filtered by common charset
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
        """Return selected characters in sorted index order."""
        return [string[i] for i in sorted(indices)]

    # --------------------------------------------------------

    def _test(self, additions: list[int]):
        """
        Delta Debugging _test():
        FAIL → subsequence is still common → expand further
        PASS → subsequence invalid → stop this branch
        """
        subseq = ''.join(self.select(self.primary_string, additions))

        for s in self.strings:
            it = iter(s)
            if not all(ch in it for ch in subseq):
                return self.PASS
        return self.FAIL

    

    def __listminus(self, c1, c2):
        """Return elements of c1 not in c2."""
        return [x for x in c1 if x not in set(c2)]

    # --------------------------------------------------------
    # Recursive ADDITION: ddmin_add()
    # --------------------------------------------------------

    def ddmin_add(self, c, r=list(), *, _depth=0):
        set_c = set(c)
        list_c = list(c)
        solns = []

        # Outer loop: expand while additions remain valid
        while self.test(list(set(c) - set_c) + r) == self.FAIL:
            soln = set(self.vanilla_add_dd(list_c, n=2, r=list(set(c) - set_c) + r))
            soln |= (set(c) - set_c)
            set_c &= soln
            list_c = list(set_c)
            solns.append(soln)

        # Recursive refinement
        if _depth < MAX_DEPTH:
            all_soln = set(c).intersection(*solns) if solns else set()
            set_c = set(c) - all_soln
            list_c = list(set_c)

            for i, soln in enumerate(solns):
                assert self.test(list(soln & set_c) + r) == self.FAIL

                irreplaceable = self.vanilla_add_dd(
                    list(set(c) - soln), n=2,
                    r=list(soln & set_c) + r
                )

                if len(irreplaceable) == 0:
                    continue

                irreplaceable = list((soln & set_c) | set(irreplaceable))
                result = self.ddmin_add(all_soln, irreplaceable + r,
                                        _depth=_depth + 1)

                if len(result) > len(soln) - len(irreplaceable):
                    solns[i] = result + irreplaceable

        if not solns:
            return []

        return list(max(solns, key=len))

    # --------------------------------------------------------
    # Base ADDITION primitive
    # --------------------------------------------------------

    def vanilla_add_dd(self, c, n=2, r=list()):
        """
        Core ADDITION primitive (DD).
        Adds characters while keeping subsequence valid.
        """

        cur_soln = []
        run = 1
        cbar_offset = 0
        c.sort()

        while True:
            tc = self.test(cur_soln + r)
            assert tc in [self.FAIL, self.UNRESOLVED]

            if n > len(c) - len(cur_soln):
                return cur_soln

            cs = self.split(self.__listminus(c, cur_soln), n)
            c_failed = False
            cbar_failed = False
            next_n = n

            # Try adding subsets
            for i in range(n):
                subset = cur_soln + cs[i]
                t = self.test(subset + r)
                if t == self.FAIL:
                    c_failed = True
                    cur_soln = subset
                    next_n = 2
                    cbar_offset = 0
                    break

            # Try complements (remove sets)
            if not c_failed:
                for j in range(n):
                    i = int((j + cbar_offset) % n)
                    complement = self.__listminus(c, cs[i])
                    t = self.test(complement + r)
                    if t == self.FAIL:
                        cbar_failed = True
                        cur_soln = complement
                        next_n -= 1
                        cbar_offset = i
                        break

            # If nothing failed, increase granularity
            if not c_failed and not cbar_failed:
                if n >= len(c) - len(cur_soln):
                    return cur_soln

                next_n = min(len(c) - len(cur_soln), n * 2)
                cbar_offset = int((cbar_offset * next_n) / n)

            n = next_n
            run += 1


# ============================================================
# MLCS Runner
# ============================================================

def dd_MLCS_addition(strings: list[str]):
    best_len = 0
    best_mlcs = ""
    best_primary = -1
    best_runtime = 0.0

    print("\n=== Running ARP (Addition-based Recursive Pruning) ===\n")
    total_start = time.time()

    for idx, primary in enumerate(strings):
        print(f"> Primary {idx+1}/{len(strings)} (len={len(primary)})")

        dd_obj = MLCS_Addition(strings, primary_string=primary)
        deltas = list(range(len(primary)))

        start = time.time()
        result = dd_obj.ddmin_add(deltas)
        end = time.time()

        mlcs = ''.join(MLCS_Addition.select(primary, result))
        mlcs_len = len(mlcs)
        runtime = end - start

        valid = validate_mlcs(mlcs, strings)
        status = "VALID" if valid else "INVALID"

        print(f"  MLCS length = {mlcs_len}, time = {runtime:.2f}s, {status}")

        if valid and mlcs_len > best_len:
            best_len = mlcs_len
            best_mlcs = mlcs
            best_primary = idx
            best_runtime = runtime

    total = time.time() - total_start

    print("\n========================================")
    print(" Best MLCS Across All Primaries")
    print("----------------------------------------")
    print(f"  MLCS:    {best_mlcs}")
    print(f"  Length:  {best_len}")
    print(f"  Primary: Sequence #{best_primary+1}")
    print(f"  Runtime: {best_runtime:.2f}s")
    print(f"  Total Runtime (all primaries): {total:.2f}s")
    print("========================================\n")

    return best_mlcs


# ============================================================
# Load sequences
# ============================================================

def load_sequences_from_file(filepath: str) -> list[str]:
    with open(filepath, "r") as f:
        return [line.strip() for line in f if line.strip()]


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ARP: Addition-based Recursive Pruning for MLCS")
    parser.add_argument("--input", type=str, required=True,
                        help="Path to sequence file (one sequence per line).")

    args = parser.parse_args()

    seqs = load_sequences_from_file(args.input)
    print(f"Loaded {len(seqs)} sequences.")
    print(f"Lengths: {[len(s) for s in seqs]}")

    start = time.time()
    mlcs = dd_MLCS_addition(seqs)
    end = time.time()

    print(f"Final Total Runtime: {end - start:.4f}s")
