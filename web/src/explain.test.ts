import { describe, expect, it } from "vitest";
import { parseGenomic, parseSubstitution, reasons, residueOf, shownDirections } from "./explain";

const NEUTRAL = { gnomad_af: null, phylop: 2, esm_llr: -5, grantham: 70, protein_regions: [] };

describe("parsing", () => {
  it("reads SNV and indel genomic changes", () => {
    expect(parseGenomic("17:43063931 G>A")).toEqual({
      chrom: "17",
      position: "43063931",
      ref: "G",
      alt: "A",
    });
    expect(parseGenomic("17:43124027 ACT>A")?.ref).toBe("ACT");
    expect(parseGenomic("BRCA1 p.Arg1699Trp")).toBeNull();
  });

  it("reads substitutions and residue numbers", () => {
    expect(parseSubstitution("p.Arg1699Trp")).toEqual({ residue: "1699", aaRef: "Arg", aaAlt: "Trp" });
    expect(parseSubstitution("p.Glu23fs")).toBeNull();
    expect(residueOf("p.Glu23ValfsTer17")).toBe("23");
    expect(residueOf(null)).toBe("");
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

  it("reads NMD and retained protein for truncating variants", () => {
    const early = reasons({ ...NEUTRAL, nmd_escape: 0, protein_retained: 0.1 });
    expect(early.map((r) => r.key)).toEqual(["nmdTriggered", "mostlyLost"]);
    const late = reasons({ ...NEUTRAL, nmd_escape: 1, protein_retained: 0.98 });
    expect(late.map((r) => r.key)).toEqual(["nmdEscapes", "mostlyKept"]);
    expect(late.every((r) => r.direction === "benign")).toBe(true);
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
