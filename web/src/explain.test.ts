import { describe, expect, it } from "vitest";
import { parseChange, reasons, shownDirections } from "./explain";

const NEUTRAL = { gnomad_af: null, phylop: 2, esm_llr: -5, grantham: 70, protein_regions: [] };

describe("parseChange", () => {
  it("reads the genomic and protein change", () => {
    expect(parseChange("17:43063931 G>A", "p.Arg1699Trp")).toEqual({
      chrom: "17",
      position: "43063931",
      ref: "G",
      alt: "A",
      residue: "1699",
      aaRef: "Arg",
      aaAlt: "Trp",
    });
  });

  it("returns null for shapes it cannot explain", () => {
    expect(parseChange("17:43124027 ACT>A", "p.Glu23fs")).toBeNull();
  });
});

describe("reasons", () => {
  it("lists harmful evidence for a conserved, absent, unusual change in a domain", () => {
    const found = reasons({
      gnomad_af: 0,
      phylop: 7.6,
      esm_llr: -13.9,
      grantham: 101,
      protein_regions: ["domain"],
    });
    expect(found.map((r) => r.key)).toEqual([
      "conserved",
      "absent",
      "esmUnusual",
      "chemDifferent",
      "inRegion",
    ]);
    expect(found.every((r) => r.direction === "pathogenic")).toBe(true);
  });

  it("lists harmless evidence for a common, tolerated, similar change", () => {
    const found = reasons({
      gnomad_af: 0.02,
      phylop: 0.3,
      esm_llr: -1,
      grantham: 30,
      protein_regions: [],
    });
    expect(found.map((r) => r.key)).toEqual([
      "notConserved",
      "common",
      "esmTolerated",
      "chemSimilar",
    ]);
  });

  it("stays silent on unmeasured or middling values", () => {
    expect(reasons(NEUTRAL)).toEqual([]);
  });
});

describe("shownDirections", () => {
  it("shows both directions only for uncertain results", () => {
    expect(shownDirections("likely_pathogenic")).toEqual(["pathogenic"]);
    expect(shownDirections("likely_benign")).toEqual(["benign"]);
    expect(shownDirections("uncertain")).toEqual(["pathogenic", "benign"]);
  });
});
