from __future__ import annotations
import os, sys, time, random, argparse
from collections import defaultdict, Counter
from typing import List
from bisect import bisect_left
import itertools
import numpy as np

try:
    from delta_debugging.DD import DD
except ImportError:
    print("Could not import delta_debugging.DD. Check PYTHONPATH.")
    sys.exit(1)


MAX_DEPTH = 2
BEAM_W = 2000
TOP_PATHS_FOR_SURVIVAL = 1000
USE_LCS_UB = True
PREFIX_LEN = 600
PAIR_SAMPLE = 1000
random.seed(0)

_HAS_NUMBA = False
try:
    import numba
    @numba.njit
    def _build_lcs_numba(a_arr, b_arr):
        n = len(a_arr)
        m = len(b_arr)
        L = np.zeros((n + 1, m + 1), dtype=np.int32)
        for i in range(n - 1, -1, -1):
            ai = a_arr[i]
            for j in range(m - 1, -1, -1):
                if ai == b_arr[j]:
                    L[i, j] = L[i+1, j+1] + 1
                else:
                    v1 = L[i+1, j]
                    v2 = L[i, j+1]
                    L[i, j] = v1 if v1 > v2 else v2
        return L
    _dummy = _build_lcs_numba(np.array([65,66], dtype=np.int8),
                               np.array([65,67], dtype=np.int8))
    _HAS_NUMBA = True
    print("[INFO] Numba JIT enabled", flush=True)
except Exception:
    print("[WARN] numba not available — suffix LCS will be slow")
    print("       pip install numba --break-system-packages --user")

def build_suffix_lcs_table(a, b):
    if _HAS_NUMBA:
        a_arr = np.frombuffer(a.encode('latin-1'), dtype=np.int8)
        b_arr = np.frombuffer(b.encode('latin-1'), dtype=np.int8)
        return _build_lcs_numba(a_arr, b_arr)
    else:
        n, m = len(a), len(b)
        L = [[0]*(m+1) for _ in range(n+1)]
        for i in range(n-1, -1, -1):
            for j in range(m-1, -1, -1):
                L[i][j] = L[i+1][j+1]+1 if a[i]==b[j] else max(L[i+1][j], L[i][j+1])
        return L

def build_all_suffix_lcs_tables(primary, others):
    return [build_suffix_lcs_table(primary, s) for s in others]

def is_valid_subsequence(subseq, seq):
    it = iter(seq)
    return all(ch in it for ch in subseq)

def validate_mlcs(candidate, sequences):
    return all(is_valid_subsequence(candidate, s) for s in sequences)

def _build_pos_lists(strings):
    pos_all = []
    for s in strings:
        d = defaultdict(list)
        for j, ch in enumerate(s):
            d[ch].append(j)
        pos_all.append(d)
    return pos_all

def _next_ge(lst, lo):
    i = bisect_left(lst, lo)
    return lst[i] if i < len(lst) else -1


def approximate_lcs_primary(sequences, prefix_len=PREFIX_LEN, pair_sample=PAIR_SAMPLE):
    n = len(sequences)
    scores = [0.0] * n
    seqs = [s[:prefix_len] for s in sequences]

    pairs = list(itertools.combinations(range(n), 2))
    if len(pairs) > pair_sample:
        pairs = random.sample(pairs, pair_sample)

    for (i, j) in pairs:
        s1, s2 = seqs[i], seqs[j]
        L = build_suffix_lcs_table(s1, s2)
        if _HAS_NUMBA:
            lcs_len = int(L[0, 0])
        else:
            lcs_len = L[0][0]
        sim = lcs_len / min(len(s1), len(s2))
        scores[i] += sim
        scores[j] += sim

    avg_scores = [s / (n - 1) for s in scores]
    ranked = sorted(range(n), key=lambda i: -avg_scores[i])
    return ranked, avg_scores

def beam_rank_deltas(primary, strings, beam_w=BEAM_W, top_paths=TOP_PATHS_FOR_SURVIVAL, primary_idx=0):
    others = [s for i, s in enumerate(strings) if i != primary_idx]
    if not others:
        return list(range(len(primary)))

    K = len(others)
    n = len(primary)
    lens_b = [len(s) for s in others]

    t0 = time.time()
    tables = build_all_suffix_lcs_tables(primary, others)
    t1 = time.time()
    

    
    tables = [np.array(t, dtype=np.int32) if not isinstance(t, np.ndarray) else t for t in tables]

   
    max_rows = max(t.shape[0] for t in tables)
    max_cols = max(t.shape[1] for t in tables)
    padded = []
    for t in tables:
        if t.shape[0] < max_rows or t.shape[1] < max_cols:
            p = np.zeros((max_rows, max_cols), dtype=np.int32)
            p[:t.shape[0], :t.shape[1]] = t
            padded.append(p)
        else:
            padded.append(t)
    all_tables = np.stack(padded)
    max_m = max_cols - 1  # max valid column index in padded tables
    lens_b_np = np.array(lens_b, dtype=np.int32)
    k_range = np.arange(K)

    pos_all = _build_pos_lists(others)
    charset = set(primary)
    pos_np = {}
    for ch in charset:
        pos_np[ch] = [
            np.array(pos_all[k].get(ch, []), dtype=np.int32)
            for k in range(K)
        ]

    def batch_ub(i_next, curs_2d):
        B_local = curs_2d.shape[0]
        if i_next > n or B_local == 0:
            return np.zeros(B_local, dtype=np.int32)
        j_vals = curs_2d + 1
        valid = j_vals <= lens_b_np[None, :]
        j_clip = np.clip(j_vals, 0, max_m)
        raw = all_tables[k_range[None, :], i_next, j_clip]
        raw = np.where(valid, raw, 0)
        counts = valid.sum(axis=1)
        sums = raw.sum(axis=1, dtype=np.float64)
        return np.where(counts > 0, (sums / counts * 1.2).astype(np.int32), 0)

    def batch_advance(ch, curs_2d):
        B_local = curs_2d.shape[0]
        if ch not in pos_np:
            return np.zeros(B_local, dtype=bool), None
        ok = np.ones(B_local, dtype=bool)
        new_curs = np.empty((B_local, K), dtype=np.int32)
        for k in range(K):
            parr = pos_np[ch][k]
            if len(parr) == 0:
                ok[:] = False
                return ok, None
            targets = curs_2d[:, k] + 1
            idxs = np.searchsorted(parr, targets)
            found = idxs < len(parr)
            ok &= found
            safe = np.minimum(idxs, len(parr) - 1)
            new_curs[:, k] = parr[safe]
        return ok, new_curs

    # Parent pool 
    pool = [(-1, -1)]
    pool_counter = 1

    # Initial beam
    B = 1
    beam_ln = np.zeros(1, dtype=np.int32)
    beam_curs = np.full((1, K), -1, dtype=np.int32)
    beam_pool = np.array([0], dtype=np.int32)

    best_ln = 0
    best_len_at = {}

    for i in range(n):
        if B == 0:
            break

        ch = primary[i]
        i_next = i + 1

        # BATCH phase: pre-compute all UBs and cursor advances
        inc_ok, inc_curs_all = batch_advance(ch, beam_curs)

        inc_which = np.where(inc_ok)[0]
        inc_ub_vals = np.zeros(B, dtype=np.int32)
        if len(inc_which) > 0:
            inc_ub_vals[inc_which] = batch_ub(i_next, inc_curs_all[inc_which])

        skip_ub_vals = batch_ub(i_next, beam_curs)

      
        candidates = []
        for bi in range(B):
            best_ln = max(best_ln, int(beam_ln[bi]))

           
            if inc_ok[bi]:
                ln_val = int(beam_ln[bi]) + 1
                f_val = ln_val + int(inc_ub_vals[bi])
                curs_tup = tuple(inc_curs_all[bi].tolist())
                key = (i_next, curs_tup)
                if f_val >= best_ln - 10 and best_len_at.get(key, -1) < ln_val:
                    best_len_at[key] = ln_val
                    pool.append((int(beam_pool[bi]), i))
                    candidates.append((f_val, ln_val, curs_tup, pool_counter))
                    pool_counter += 1

            
            ln_val = int(beam_ln[bi])
            f_val = ln_val + int(skip_ub_vals[bi])
            curs_tup = tuple(beam_curs[bi].tolist())
            key = (i_next, curs_tup)
            if f_val >= best_ln - 10 and best_len_at.get(key, -1) < ln_val:
                best_len_at[key] = ln_val
                pool.append((int(beam_pool[bi]), -1))
                candidates.append((f_val, ln_val, curs_tup, pool_counter))
                pool_counter += 1

        if not candidates:
            break

        candidates.sort(key=lambda s: (-s[0], -s[1], s[2], random.random()))
        candidates = candidates[:beam_w]

        B = len(candidates)
        beam_ln = np.array([c[1] for c in candidates], dtype=np.int32)
        beam_curs = np.array([list(c[2]) for c in candidates], dtype=np.int32)
        beam_pool = np.array([c[3] for c in candidates], dtype=np.int32)

    surv = Counter()
    n_top = min(top_paths, B)
    for idx in range(n_top):
        pid = int(beam_pool[idx])
        while pid >= 0:
            parent_pid, added = pool[pid]
            if added >= 0:
                surv[added] += 1
            pid = parent_pid

    for k in list(surv.keys()):
        surv[k] = int(surv[k] * 5)

    additions = list(range(n))
    additions.sort(key=lambda idx: (-surv.get(idx, 0), idx))

    t2 = time.time()
   
    return additions

class MLCS_Addition(DD):
    def __init__(self, strings, primary_string=None):
        super().__init__()
        self.strings = strings
        self.primary_string = primary_string or min(strings, key=len)

    @staticmethod
    def select(string, indices):
        return [string[i] for i in sorted(indices)]

    def _test(self, additions):
        subseq = "".join(self.select(self.primary_string, additions))
        for s in self.strings:
            it = iter(s)
            if not all(ch in it for ch in subseq):
                return self.PASS
        return self.FAIL

    def __listminus(self, c1, c2):
        s2 = set(c2)
        return [x for x in c1 if x not in s2]

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
                assert self.test(list(soln & set_c) + r) == self.FAIL
                irreplaceable = self.vanilla_add_dd(
                    list(set(c) - soln), n=2, r=list(soln & set_c) + r)
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
        while True:
            if n > len(c) - len(cur_soln):
                return cur_soln
            cs = self.split(self.__listminus(c, cur_soln), n)
            c_failed = False
            cbar_failed = False
            next_n = n
            for i in range(n):
                subset = cur_soln + cs[i]
                if self.test(subset + r) == self.FAIL:
                    c_failed = True
                    cur_soln = subset
                    next_n = 2
                    break
            if not c_failed:
                for j in range(n):
                    i = (j + cbar_offset) % n
                    complement = self.__listminus(c, cs[i])
                    if self.test(complement + r) == self.FAIL:
                        cbar_failed = True
                        cur_soln = complement
                        next_n -= 1
                        cbar_offset = i
                        break
            if not c_failed and not cbar_failed:
                if n >= len(c) - len(cur_soln):
                    return cur_soln
                next_n = min(len(c) - len(cur_soln), n * 2)
            n = next_n


def load_sequences_from_file(path):
    with open(path, "r", encoding="latin-1", errors="ignore") as f:
        return [line.strip() for line in f if line.strip()]

def dd_MLCS_addition(sequences, beam_w=BEAM_W):
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
    top1_runtime = None
    for rank, idx in enumerate(top5, 1):
        s = sequences[idx]

        
        charset = set(s)
        for seq in sequences:
            charset &= set(seq)
        filtered = "".join(c for c in s if c in charset)

        print(f"[Run {rank}] Primary #{idx+1} | len={len(s)} |  score={scores[idx]:.4f}")

        t0 = time.time()
        deltas = beam_rank_deltas(filtered, sequences, beam_w=beam_w, primary_idx=idx)

        t_dd = time.time()
        dd = MLCS_Addition(sequences, primary_string=filtered)
        result = dd.ddmin_add(deltas)

        mlcs = "".join(MLCS_Addition.select(filtered, result))
        valid = validate_mlcs(mlcs, sequences)
        elapsed = round(time.time() - t0, 2)

        #print(f"    DD addition: {time.time() - t_dd:.2f}s")
        print(f"    MLCS Len = {len(mlcs)}, Valid={valid}, Runtime={elapsed}s\n")

        results.append((len(mlcs), idx, mlcs))

        # Track top-1 ranked primary runtime
        if rank == 1:
            top1_runtime = elapsed

    # ---- PICK BEST ----
    best = max(results, key=lambda x: x[0])

    print("==============================")
    print(" BEST MLCS FOUND:")
    print(f" Length: {best[0]}")
    print(f" Primary used: #{best[1]+1}")
    print(f" MLCS: {best[2][:100]}...")
    print(f" Top-1 Runtime: {top1_runtime}s")
    print("==============================")

    return best[2], top1_runtime

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to sequence file")
    parser.add_argument("--beam", type=int, default=2000,
                        help="Beam width for beam-guided delta ordering")
    args = parser.parse_args()

    BEAM_W = args.beam

    sequences = load_sequences_from_file(args.input)
    start = time.time()
    final_mlcs, top1_runtime = dd_MLCS_addition(sequences, beam_w=args.beam)
    end = time.time()

    print(f"\nTop-1 Primary Runtime: {top1_runtime}s")
    print(f"Total Runtime (all 5): {end - start:.2f}s")