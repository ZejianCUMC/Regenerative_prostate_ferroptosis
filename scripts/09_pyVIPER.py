#!/usr/bin/env python
# coding: utf-8
# =============================================================================
# Wouter (Science) luminal -- pyVIPER pipeline (summarised archive)
#
# Sections:
#   1) Build meta cells
#   2) Load, merge and prune ARACNe-AP networks
#   3) Use aREA approach to get NES scores
#   4) Find Intact -> castration -> regeneration most changed MRs
#      ranked by absolute NES change:
#         score = |Intact - Cast28| + |Reg28 - Cast28|
#
# =============================================================================

import os
import glob
import warnings
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import pyviper

warnings.filterwarnings(
    "ignore", message=".*The 'nopython' keyword.*"
)  # jit decorator issue with sc.pp.neighbors

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
PROJECT_DIR = "/path/to/project/"
ARACNE_DIR  = os.path.join(PROJECT_DIR, "pyVIPER", "ARACNe_out")
PYVIPER_DIR = os.path.join(PROJECT_DIR, "pyVIPER")

RAW_COUNTS_FILE = "Wouter_science_luminal_raw_clean.csv"
META_FILE       = "Wouter_science_luminal_metainfo.csv"

# Timepoints to drop (T00 epi/non-epi sub-splits we don't keep)
EXCLUDED_TIMEPOINTS = {"T00_Epi", "T00_NonEpi"}

# Metacell parameters
METACELL_KEY            = "metacells_approach_2"   # n_cells_per_metacell variant
METACELL_SIZE           = 200
N_CELLS_PER_METACELL    = 20
CLUSTERS_SLOT           = "timepoint"

# ARACNe network combination
NETWORK_COLS = ["regulator", "target", "mor", "likelihood"]

# pyVIPER / aREA settings
N_JOBS              = 24
PRUNE_MAX_TARGETS   = 200
GROUPING_COLUMN     = "timepoint"

AREA_RAW_OUTPUT_FILE     = os.path.join(
    PYVIPER_DIR, "Wounter_Sci_aREA_estimated_protein_activity.csv"
)
AREA_GROUPED_OUTPUT_FILE = os.path.join(
    PYVIPER_DIR, "Wounter_Sci_area_NES_by_timepoint.csv"
)


# =============================================================================
# 1) Build meta cells
# =============================================================================
def build_metacells():
    """Load raw counts + metadata, run preprocessing, build metacells per
    timepoint, and dump the per-timepoint log2(CPM+1) matrices needed by
    ARACNe-AP."""
    os.chdir(PROJECT_DIR)

    # ---- load raw data + metadata --------------------------------------------
    raw_counts = pd.read_csv(RAW_COUNTS_FILE, index_col=0)
    metadata   = pd.read_csv(META_FILE,       index_col=0)

    adata = sc.AnnData(X=raw_counts.T)
    adata.obs = pd.merge(
        adata.obs, metadata, how="left", left_index=True, right_index=True
    )

    # drop unwanted timepoints
    adata = adata[~adata.obs["timepoint"].isin(EXCLUDED_TIMEPOINTS)].copy()

    # ---- basic QC + normalization --------------------------------------------
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_counts=10)
    adata.raw = adata

    sc.pp.normalize_total(adata, inplace=True, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, flavor="seurat", n_top_genes=3000, inplace=True
    )
    sc.pp.scale(adata)

    # ---- PCA + correlation distance for metacell construction -----------------
    sc.tl.pca(adata, svd_solver="arpack", random_state=0)
    pyviper.pp.corr_distance(adata)

    # ---- metacells: fix N cells per metacell ---------------------------------
    pyviper.pp.repr_metacells(
        adata,
        counts=None,
        pca_slot="X_pca",
        dist_slot="corr_dist",
        size=METACELL_SIZE,
        n_cells_per_metacell=N_CELLS_PER_METACELL,
        min_median_depth=None,
        clusters_slot=CLUSTERS_SLOT,
        key_added=METACELL_KEY,
    )

    # ---- concat all per-timepoint metacell matrices --------------------------
    metacell_dfs = []
    for tp in adata.obs["timepoint"].unique():
        df = adata.uns[f"{METACELL_KEY}_{tp}"]
        df.index = f"{tp}_" + df.index
        metacell_dfs.append(df)
    metacells_df = pd.concat(metacell_dfs, ignore_index=False)
    metacells_df.to_csv(
        "Wounter_Sci_Lum_viper_metacell_raw.tsv",
        sep="\t",
        index_label="CellName",
    )

    # ---- export per-timepoint log2(CPM+1) matrices for ARACNe-AP -------------
    for tp in adata.obs["timepoint"].unique():
        key = f"{METACELL_KEY}_{tp}"
        df = adata.uns[key].T                        # genes x metacells
        cpm       = df.div(df.sum(axis=0), axis=1) * 1_000_000
        log2_cpm  = np.log2(cpm + 1).round(3)
        log2_cpm.to_csv(f"{tp}_metaCells.tsv", sep="\t", index_label="gene")

    return adata, metadata


# =============================================================================
# 2) Load, merge and prune ARACNe-AP networks
# =============================================================================
def merge_networks_for_timepoint(file_list):
    """Concatenate per-bootstrap ARACNe networks for one timepoint,
    take the median mor/likelihood per regulator-target pair, and
    convert ARACNe's likelihood (a p-value-like quantity) to 1 - likelihood
    (a confidence score) which pyVIPER expects."""
    dfs = []
    for f in file_list:
        try:
            d = pd.read_csv(f, sep="\t", header=0)
            d.columns = NETWORK_COLS
            dfs.append(d)
        except Exception as e:
            print(f"  Could not read {f}: {e}")
    if not dfs:
        return None

    combined = pd.concat(dfs, ignore_index=True)
    combined = combined.groupby(["regulator", "target"], as_index=False).median()
    combined["likelihood"] = 1 - combined["likelihood"]
    return combined


def load_and_prune_networks(adata):
    """Find every network.txt under ARACNE_DIR, group by timepoint
    prefix of the parent folder, build merged TSVs, then construct
    pyVIPER Interactomes filtered to genes in adata and pruned to
    PRUNE_MAX_TARGETS targets per regulator."""
    os.chdir(ARACNE_DIR)
    unique_timepoints = adata.obs["timepoint"].unique()

    all_files = glob.glob(
        os.path.join(ARACNE_DIR, "**/network.txt"), recursive=True
    )
    print(f"Found {len(all_files)} network.txt files")

    # group by timepoint prefix on the parent folder name
    timepoint_files = {tp: [] for tp in unique_timepoints}
    for fp in all_files:
        parent = os.path.basename(os.path.dirname(fp))
        for tp in unique_timepoints:
            if parent.startswith(tp):
                timepoint_files[tp].append(fp)
                break

    # merge bootstrap networks per timepoint
    for tp, flist in timepoint_files.items():
        if not flist:
            print(f"  No networks for {tp}")
            continue
        print(f"Merging {len(flist)} networks for {tp}")
        merged = merge_networks_for_timepoint(flist)
        if merged is not None:
            merged.to_csv(f"combined_network_{tp}.tsv", index=False, sep="\t")

    # build + filter + prune Interactomes
    interactomes = {}
    for tp in unique_timepoints:
        net_file = f"combined_network_{tp}.tsv"
        if not os.path.exists(net_file):
            continue
        name = tp.replace("_", " ").title().replace(" ", "")
        net = pyviper.Interactome(name, net_table=net_file)
        net.filter_targets(adata.var_names)
        net.prune(max_targets=PRUNE_MAX_TARGETS, eliminate=True)
        interactomes[tp] = net
        print(f"  Loaded + pruned interactome for {tp}")

    return interactomes


# =============================================================================
# 3) Use aREA approach to get NES scores
# =============================================================================
def run_area(adata, interactomes, metadata):
    """Run pyVIPER (metaVIPER, aREA enrichment) across every cell,
    write per-cell NES, and the per-timepoint-mean NES matrix."""
    pruned_list = list(interactomes.values())

    adata_pa_area = pyviper.viper(
        gex_data=adata,
        interactome=pruned_list,
        eset_filter=False,
        njobs=N_JOBS,
        verbose=False,
    )

    # per-cell NES
    adata_pa_area_res = adata_pa_area.to_df()
    adata_pa_area_res.to_csv(AREA_RAW_OUTPUT_FILE, index=True)

    # cells x (MRs + timepoint)
    combined_df_area = pd.concat(
        [adata_pa_area_res, metadata[GROUPING_COLUMN]], axis=1
    )

    # per-timepoint mean NES (regulators x timepoints)
    average_nes_area_T = (
        combined_df_area.groupby(GROUPING_COLUMN).mean().transpose()
    )
    average_nes_area_T.to_csv(AREA_GROUPED_OUTPUT_FILE, index=True)

    return adata_pa_area, combined_df_area


# =============================================================================
# 4) Rank MRs by absolute NES change across castration & regeneration
#    score = |Intact - Cast28| + |Reg28 - Cast28|
#
#  Refactored from `find_v_shape_mrs` / `find_v_shape_mrs_v2`, but stripped
#  to only what we need: absolute amplitude of the V across the three
#  anchor timepoints. Direction-agnostic -- captures both V and inverted-V
#  responses.
# =============================================================================
def find_absolute_change_mrs(
    merged_count_mtx,
    intact_tp="T00_Intact",
    cast28_tp="T04_Cast_Day28",
    reg28_tp="T10_Regen_Day28",
    timepoint_col="timepoint",
    aggregator="median",
):
    """
    Identify MRs with the largest combined absolute NES change between:
        (i)  Intact   -> Cast Day 28      (castration response)
        (ii) Cast Day 28 -> Regen Day 28  (regeneration response)

    score = |Intact - Cast28| + |Reg28 - Cast28|

    Args:
        merged_count_mtx (pd.DataFrame): rows = cells (or metacells),
            columns = MRs + a `timepoint` column (default name).
        intact_tp, cast28_tp, reg28_tp (str): timepoint labels.
        timepoint_col (str): name of the grouping column.
        aggregator (str): "median" (default, robust) or "mean".

    Returns:
        pd.DataFrame sorted by `abs_change_score` desc, with columns:
            MR, intact_NES, cast28_NES, reg28_NES,
            cast_delta, regen_delta, abs_change_score
        where:
            cast_delta  = Cast28  - Intact   (signed)
            regen_delta = Reg28   - Cast28   (signed)
            abs_change_score = |cast_delta| + |regen_delta|
    """
    if timepoint_col not in merged_count_mtx.columns:
        raise ValueError(f"'{timepoint_col}' column not in input DataFrame")

    # aggregate per timepoint
    if aggregator == "median":
        agg = merged_count_mtx.groupby(timepoint_col).median()
    elif aggregator == "mean":
        agg = merged_count_mtx.groupby(timepoint_col).mean()
    else:
        raise ValueError("aggregator must be 'median' or 'mean'")

    for tp in (intact_tp, cast28_tp, reg28_tp):
        if tp not in agg.index:
            raise ValueError(
                f"timepoint '{tp}' missing from data; "
                f"available: {sorted(agg.index.tolist())}"
            )

    intact = agg.loc[intact_tp]
    cast28 = agg.loc[cast28_tp]
    reg28  = agg.loc[reg28_tp]

    cast_delta  = cast28 - intact
    regen_delta = reg28  - cast28
    score       = cast_delta.abs() + regen_delta.abs()

    out = pd.DataFrame({
        "MR":               score.index,
        "intact_NES":       intact.values,
        "cast28_NES":       cast28.values,
        "reg28_NES":        reg28.values,
        "cast_delta":       cast_delta.values,
        "regen_delta":      regen_delta.values,
        "abs_change_score": score.values,
    })

    # drop rows that became NaN (e.g. constant zero columns)
    out = out.dropna(subset=["abs_change_score"])

    out = out.sort_values("abs_change_score", ascending=False).reset_index(
        drop=True
    )
    return out


# =============================================================================
# Main
# =============================================================================
def main():
    # 1) metacells
    adata, metadata = build_metacells()

    # 2) ARACNe networks -> filtered + pruned interactomes
    interactomes = load_and_prune_networks(adata)

    # 3) aREA -> per-cell NES + per-timepoint mean NES
    adata_pa_area, combined_df_area = run_area(adata, interactomes, metadata)

    # 4) rank MRs by |Intact - Cast28| + |Reg28 - Cast28|
    abs_change_mrs = find_absolute_change_mrs(combined_df_area)
    abs_change_mrs.to_csv(
        os.path.join(
            PROJECT_DIR,
            "Wounter_sci_lum1_area_abs_change_mrs.csv",
        ),
        index=False,
        header=True,
    )
    print(abs_change_mrs.head(20))


if __name__ == "__main__":
    main()