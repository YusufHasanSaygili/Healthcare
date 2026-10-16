import { describe, expect, it } from "vitest";
import schema from "../openapi.json";
import { MESSAGES, fill } from "./i18n";

const BACKEND_CONSEQUENCES = [
  "intergenic",
  "multiple_genes",
  "multi_nucleotide",
  "synonymous",
  "stop_gained",
  "stop_lost",
  "start_lost",
  "canonical_splice",
  "intronic",
  "utr",
];

const UNIPROT_REGIONS = [
  "domain",
  "active_site",
  "binding_site",
  "motif",
  "region",
  "coiled_coil",
  "transmembrane",
  "signal_or_propeptide",
  "disulfide",
  "ptm_site",
  "helix",
  "strand",
  "turn",
];

function keys(value: unknown, prefix = ""): string[] {
  if (typeof value !== "object" || value === null) {
    return [prefix];
  }
  return Object.entries(value).flatMap(([key, child]) =>
    keys(child, prefix ? `${prefix}.${key}` : key),
  );
}

describe("translations", () => {
  it("English and Turkish define exactly the same messages", () => {
    expect(keys(MESSAGES.tr).sort()).toEqual(keys(MESSAGES.en).sort());
  });

  it("cover every out-of-scope reason and protein region the backend can return", () => {
    for (const messages of Object.values(MESSAGES)) {
      expect(Object.keys(messages.consequences).sort()).toEqual([...BACKEND_CONSEQUENCES].sort());
      expect(Object.keys(messages.regions).sort()).toEqual([...UNIPROT_REGIONS].sort());
    }
  });

  it("cover every evidence group in the API schema", () => {
    const groups: string[] = schema.components.schemas.Contribution.properties.group.enum;
    for (const messages of Object.values(MESSAGES)) {
      expect(Object.keys(messages.groups).sort()).toEqual([...groups].sort());
    }
  });

  it("keep every placeholder in both languages", () => {
    const placeholders = (text: string) => [...text.matchAll(/\{(\w+)\}/g)].map((m) => m[1]).sort();
    const english = keys(MESSAGES.en);
    for (const key of english) {
      const pick = (messages: object) =>
        key.split(".").reduce<unknown>((node, part) => (node as Record<string, unknown>)[part], messages);
      const en = pick(MESSAGES.en);
      const tr = pick(MESSAGES.tr);
      if (typeof en === "string" && typeof tr === "string") {
        expect(placeholders(tr), key).toEqual(placeholders(en));
      }
    }
  });

  it("fills placeholders and leaves unknown ones visible", () => {
    expect(fill("1 in {value} alleles", { value: "42" })).toBe("1 in 42 alleles");
    expect(fill("{missing}", {})).toBe("{missing}");
  });
});
