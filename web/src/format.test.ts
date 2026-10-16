import { describe, expect, it } from "vitest";
import {
  barShare,
  collapseChromosomes,
  formatFrequency,
  formatPercent,
  strongestSignals,
  type Contribution,
} from "./format";

describe("formatPercent", () => {
  it("formats probabilities with the locale's conventions", () => {
    expect(formatPercent(0.954, "en")).toBe("95.4%");
    expect(formatPercent(0.954, "tr")).toBe("%95,4");
  });

  it("never rounds a tiny or near-certain probability to an exact 0% or 100%", () => {
    expect(formatPercent(0.00004, "en")).toBe("<0.1%");
    expect(formatPercent(0.99995, "en")).toBe(">99.9%");
    expect(formatPercent(0, "en")).toBe("0%");
    expect(formatPercent(1, "en")).toBe("100%");
  });
});

describe("formatFrequency", () => {
  it("separates not measured from not observed", () => {
    expect(formatFrequency(null, "en")).toEqual({ kind: "not_measured" });
    expect(formatFrequency(0, "en")).toEqual({ kind: "not_observed" });
  });

  it("uses percentages for common variants and '1 in N' for rare ones", () => {
    expect(formatFrequency(0.34, "en")).toEqual({ kind: "percent", value: "34%" });
    expect(formatFrequency(8.20928e-6, "en")).toEqual({ kind: "one_in", value: "121,813" });
    expect(formatFrequency(8.20928e-6, "tr")).toEqual({ kind: "one_in", value: "121.813" });
  });
});

const CONTRIBUTIONS: Contribution[] = [
  { group: "esm_scores", log_odds: 3.5 },
  { group: "gene_constraint", log_odds: -1.75 },
  { group: "protein_regions", log_odds: 0.05 },
];

describe("contribution helpers", () => {
  it("scales bars against the largest absolute contribution", () => {
    expect(barShare(3.5, CONTRIBUTIONS)).toBe(1);
    expect(barShare(-1.75, CONTRIBUTIONS)).toBe(0.5);
  });

  it("keeps tiny effects on a fixed scale when every contribution is small", () => {
    expect(barShare(0.05, [{ group: "conservation", log_odds: 0.05 }])).toBe(0.05);
  });

  it("names the strongest signals and ignores negligible ones", () => {
    expect(strongestSignals(CONTRIBUTIONS).map((c) => c.group)).toEqual([
      "esm_scores",
      "gene_constraint",
    ]);
    expect(strongestSignals([{ group: "protein_regions", log_odds: 0.01 }])).toEqual([]);
  });
});

describe("collapseChromosomes", () => {
  it("merges per-chromosome sources that share a version", () => {
    expect(
      collapseChromosomes({
        mane: "v1.5",
        gnomad_exomes_chr1: "v4.1",
        gnomad_exomes_chrX: "v4.1",
        phylop_chr2: "a",
        phylop_chr3: "b",
      }),
    ).toEqual([
      ["mane", "v1.5"],
      ["gnomad_exomes_chr*", "v4.1"],
      ["phylop_chr*", "a"],
      ["phylop_chr3", "b"],
    ]);
  });
});
