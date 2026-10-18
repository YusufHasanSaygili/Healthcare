import { describe, expect, it } from "vitest";
import type { BatchItem } from "./api/client";
import { batchCsv, CSV_COLUMNS, decisionCounts, MAX_UPLOAD_BYTES, uploadProblem } from "./batch";

const PREDICTED: BatchItem = {
  line: 3,
  result: {
    status: "predicted",
    variant: "17:43063931 G>A",
    gene: "BRCA1",
    transcript_id: "ENST00000357654.9",
    hgvs_p: "p.Arg1699Trp",
    consequence: "missense",
    probability: 0.954,
    decision: "likely_pathogenic",
    lightgbm: 0.894,
    dual_tower: 0.911,
    contributions: [],
    evidence: { gnomad_af: null, phylop: null, esm_llr: null, grantham: null, protein_regions: [] },
  },
};

const OUT_OF_SCOPE: BatchItem = {
  line: 7,
  result: {
    status: "out_of_scope",
    variant: "1:10500 A>G",
    gene: null,
    consequence: "intergenic",
    reason: 'no MANE Select transcript overlaps it, "intergenic"',
  },
};

describe("uploadProblem", () => {
  it("rejects empty, oversized and compressed files before upload", () => {
    expect(uploadProblem("a.vcf", 0)).toBe("empty");
    expect(uploadProblem("a.vcf", MAX_UPLOAD_BYTES + 1)).toBe("too_large");
    expect(uploadProblem("a.vcf.gz", 100)).toBe("compressed");
    expect(uploadProblem("a.vcf", MAX_UPLOAD_BYTES)).toBeNull();
  });
});

describe("batchCsv", () => {
  it("writes one row per result under a fixed header", () => {
    const lines = batchCsv([PREDICTED, OUT_OF_SCOPE]).trimEnd().split("\r\n");
    expect(lines[0]).toBe(CSV_COLUMNS.join(","));
    expect(lines).toHaveLength(3);
    expect(lines[1]).toBe(
      "3,17:43063931 G>A,predicted,BRCA1,ENST00000357654.9,p.Arg1699Trp,0.954000,likely_pathogenic,0.894000,0.911000,missense,",
    );
  });

  it("quotes cells with commas or quotes", () => {
    const last = batchCsv([OUT_OF_SCOPE]).trimEnd().split("\r\n")[1];
    expect(last).toBe(
      '7,1:10500 A>G,out_of_scope,,,,,,,,intergenic,"no MANE Select transcript overlaps it, ""intergenic"""',
    );
  });

  it("neutralises spreadsheet formulas", () => {
    const item: BatchItem = {
      ...OUT_OF_SCOPE,
      result: {
        ...OUT_OF_SCOPE.result,
        status: "out_of_scope",
        reason: "=HYPERLINK(1)",
      } as BatchItem["result"],
    };
    expect(batchCsv([item])).toContain(",'=HYPERLINK(1)");
  });
});

describe("decisionCounts", () => {
  it("counts each decision and out-of-scope results", () => {
    expect(decisionCounts([PREDICTED, OUT_OF_SCOPE, PREDICTED])).toEqual({
      likely_pathogenic: 2,
      uncertain: 0,
      likely_benign: 0,
      out_of_scope: 1,
    });
  });
});
