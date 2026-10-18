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
  | "inRegion"
  | "nmdTriggered"
  | "nmdEscapes"
  | "mostlyLost"
  | "mostlyKept";

export type Reason = { key: ReasonKey; direction: Direction };

export type Genomic = { chrom: string; position: string; ref: string; alt: string };

export type Substitution = { residue: string; aaRef: string; aaAlt: string };

const VARIANT = /^(\w+):(\d+) ([ACGT]+)>([ACGT]+)$/;
const SUBSTITUTION = /^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$/;
const RESIDUE = /(\d+)/;

export const THRESHOLDS = {
  conserved: 4,
  notConserved: 1,
  common: 0.001,
  esmUnusual: -7,
  esmTolerated: -3,
  chemDifferent: 100,
  chemSimilar: 50,
  mostlyLost: 0.5,
  mostlyKept: 0.9,
};

export function parseGenomic(variant: string): Genomic | null {
  const match = VARIANT.exec(variant);
  if (!match) {
    return null;
  }
  const [, chrom = "", position = "", ref = "", alt = ""] = match;
  return { chrom, position, ref, alt };
}

export function parseSubstitution(hgvs: string | null): Substitution | null {
  const match = hgvs ? SUBSTITUTION.exec(hgvs) : null;
  if (!match) {
    return null;
  }
  const [, aaRef = "", residue = "", aaAlt = ""] = match;
  return { residue, aaRef, aaAlt };
}

export function residueOf(hgvs: string | null): string {
  return (hgvs && RESIDUE.exec(hgvs)?.[1]) ?? "";
}

export function reasons(evidence: PredictedVariant["evidence"]): Reason[] {
  const found: Reason[] = [];
  const { phylop, gnomad_af, esm_llr, grantham, protein_regions } = evidence;
  const nmd = evidence.nmd_escape ?? null;
  const retained = evidence.protein_retained ?? null;
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
  if (nmd === 0) {
    found.push({ key: "nmdTriggered", direction: "pathogenic" });
  } else if (nmd === 1) {
    found.push({ key: "nmdEscapes", direction: "benign" });
  }
  if (retained !== null && retained < THRESHOLDS.mostlyLost) {
    found.push({ key: "mostlyLost", direction: "pathogenic" });
  } else if (retained !== null && retained >= THRESHOLDS.mostlyKept) {
    found.push({ key: "mostlyKept", direction: "benign" });
  }
  return found;
}

export function shownDirections(decision: Decision): Direction[] {
  if (decision === "likely_pathogenic") {
    return ["pathogenic"];
  }
  return decision === "likely_benign" ? ["benign"] : ["pathogenic", "benign"];
}
