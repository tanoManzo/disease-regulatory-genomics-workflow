#!/usr/bin/env python3
"""Create the model-performance table, experimental audit, and THE LIST."""

from pathlib import Path
import argparse
import pandas as pd


TRAIN = {f"chr{x}" for x in [1, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, "X", "Y"]}
VAL = {"chr7"}
TEST = {"chr8", "chr9"}


def bed_counts(path):
    chroms = pd.read_csv(path, sep="\t", header=None, usecols=[0])[0]
    return int(chroms.isin(TRAIN).sum()), int(chroms.isin(VAL).sum()), int(chroms.isin(TEST).sum())


def read_metric(path):
    return float(Path(path).read_text().split(":", 1)[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-dir", required=True)
    ap.add_argument("--positive-bed", required=True)
    ap.add_argument("--control-bed", required=True)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--population", required=True)
    ap.add_argument("--targets", required=True)
    ap.add_argument("--tfbs", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    model = Path(args.model_dir)
    thresholds = {}
    for line in (model / "fpr_threshold_scores.txt").read_text().splitlines():
        fpr, score = line.split()
        thresholds[int(fpr)] = float(score)
    pt, pv, pe = bed_counts(args.positive_bed)
    nt, nv, ne = bed_counts(args.control_bed)
    performance = pd.DataFrame([{
        "model": "MotorNeuron_Enhancer_DHS_x2", "regulatory_class": "enhancer",
        "biosample": "ALS donor iPSC-derived motor neuron", "donor": "ENCDO689YDN",
        "accessibility_experiment": "ENCSR131HOY", "histone_experiment": "ENCSR489LNU",
        "positive_train": pt, "positive_validation": pv, "positive_test": pe,
        "control_train": nt, "control_validation": nv, "control_test": ne,
        "test_ROC_AUC": read_metric(model / "auc.txt"),
        "test_PR_AUC": read_metric(model / "prc.txt"),
        "threshold_FPR10": thresholds[10], "threshold_FPR5": thresholds[5],
        "threshold_FPR3": thresholds[3], "threshold_FPR1": thresholds[1],
    }])
    performance.to_csv(out / "ALS-S4-R001_model_performance.tsv", sep="\t", index=False)

    variants = pd.read_csv(args.variants, sep="\t")
    causal = variants[variants.causal_regulatory.astype(str).str.lower().eq("true")].copy()
    pop = pd.read_csv(args.population, sep="\t")[["name", "EUR_AF", "EAS_AF", "frequency_class",
        "population_specific", "ancestral_allele", "derived_allele", "EUR_DAF", "EAS_DAF", "DAF_gt_0_5"]]
    targets = pd.read_csv(args.targets, sep="\t").rename(columns={"variant": "name"})
    keep_targets = ["name", "nearest_gene", "nearest_tss_distance_bp", "hic_target_genes",
                    "hic_supporting_experiments", "mapping_methods"]
    causal = causal.merge(pop, on="name", how="left").merge(targets[keep_targets], on="name", how="left")
    tfbs = pd.read_csv(args.tfbs, sep="\t")
    tfbs = tfbs[tfbs.variant_in_TFBS | tfbs.disrupts_TFBS | tfbs.creates_TFBS]
    tf_summary = tfbs.groupby("variant").apply(lambda g: pd.Series({
        "key_TFBS": ";".join(sorted(set(g.loc[g.variant_in_TFBS, "TF"]))),
        "disrupted_TFBS": ";".join(sorted(set(g.loc[g.disrupts_TFBS, "TF"]))),
        "created_TFBS": ";".join(sorted(set(g.loc[g.creates_TFBS, "TF"])))
    }), include_groups=False).reset_index().rename(columns={"variant": "name"})
    causal = causal.merge(tf_summary, on="name", how="left")
    causal["experimentally_profiled"] = False
    causal["experimental_evidence"] = "No exact-variant functional study identified (search date 2026-09-23)"
    causal["priority_rank"] = range(1, len(causal) + 1)
    causal.to_csv(out / "ALS-S4-R009_THE_LIST.tsv", sep="\t", index=False)

    unc = variants[variants.name.isin(["rs12608932", "rs12973192"])].copy()
    unc["experimentally_profiled"] = True
    unc["experimental_evidence"] = "TDP-43-dependent UNC13A cryptic-exon inclusion (PMID:35197626; PMID:35197628)"
    audit = pd.concat([causal, unc], ignore_index=True, sort=False)
    audit.to_csv(out / "ALS-S4-R008_experimental_audit.tsv", sep="\t", index=False)
    print(performance.to_string(index=False))
    print(f"THE LIST: {len(causal)} variants; experimental causal hits: 0/{len(causal)}")


if __name__ == "__main__":
    main()
