import type { Decision, PredictedVariant } from "./api/client";

export type Direction = "pathogenic" | "benign";

export type ReasonKey =
  | "conserved"
  | "notConserved"
  | "absent"
  | "common"
  | "esmUnusual"
  | "esmTolerated"
  | "chemDifferent"
  | "chemSimilar"
  | "inRegion";

export type Reason = { key: ReasonKey; direction: Direction };

export type Change = {
  chrom: string;
  position: string;
  ref: string;
  alt: string;
  residue: string;
  aaRef: string;
  aaAlt: string;
};

const VARIANT = /^(\w+):(\d+) ([ACGT])>([ACGT])$/;
const PROTEIN = /^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$/;

export const THRESHOLDS = {
  conserved: 4,
  notConserved: 1,
  common: 0.001,
  esmUnusual: -7,
  esmTolerated: -3,
  chemDifferent: 100,
  chemSimilar: 50,
};

export function parseChange(variant: string, hgvs: string): Change | null {
  const genomic = VARIANT.exec(variant);
  const protein = PROTEIN.exec(hgvs);
  if (!genomic || !protein) {
    return null;
  }
  const [, chrom = "", position = "", ref = "", alt = ""] = genomic;
  const [, aaRef = "", residue = "", aaAlt = ""] = protein;
  return { chrom, position, ref, alt, residue, aaRef, aaAlt };
}

export function reasons(evidence: PredictedVariant["evidence"]): Reason[] {
  const found: Reason[] = [];
  const { phylop, gnomad_af, esm_llr, grantham, protein_regions } = evidence;
  if (phylop !== null && phylop >= THRESHOLDS.conserved) {
    found.push({ key: "conserved", direction: "pathogenic" });
  } else if (phylop !== null && phylop <= THRESHOLDS.notConserved) {
    found.push({ key: "notConserved", direction: "benign" });
  }
  if (gnomad_af === 0) {
    found.push({ key: "absent", direction: "pathogenic" });
  } else if (gnomad_af !== null && gnomad_af >= THRESHOLDS.common) {
    found.push({ key: "common", direction: "benign" });
  }
  if (esm_llr !== null && esm_llr <= THRESHOLDS.esmUnusual) {
    found.push({ key: "esmUnusual", direction: "pathogenic" });
  } else if (esm_llr !== null && esm_llr >= THRESHOLDS.esmTolerated) {
    found.push({ key: "esmTolerated", direction: "benign" });
  }
  if (grantham !== null && grantham >= THRESHOLDS.chemDifferent) {
    found.push({ key: "chemDifferent", direction: "pathogenic" });
  } else if (grantham !== null && grantham <= THRESHOLDS.chemSimilar) {
    found.push({ key: "chemSimilar", direction: "benign" });
  }
  if (protein_regions.length > 0) {
    found.push({ key: "inRegion", direction: "pathogenic" });
  }
  return found;
}

export function shownDirections(decision: Decision): Direction[] {
  if (decision === "likely_pathogenic") {
    return ["pathogenic"];
  }
  return decision === "likely_benign" ? ["benign"] : ["pathogenic", "benign"];
}
