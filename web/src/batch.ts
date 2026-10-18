import type { BatchItem, Decision } from "./api/client";

export const MAX_UPLOAD_BYTES = 2_000_000;

export type UploadProblem = "empty" | "too_large" | "compressed" | "over_limit" | "not_text";

export const API_PROBLEMS: Partial<Record<number, UploadProblem>> = {
  413: "over_limit",
  415: "compressed",
  422: "not_text",
};

export function uploadProblem(name: string, size: number): UploadProblem | null {
  if (size === 0) {
    return "empty";
  }
  if (size > MAX_UPLOAD_BYTES) {
    return "too_large";
  }
  if (/\.(gz|bgz|zip)$/i.test(name)) {
    return "compressed";
  }
  return null;
}

export const CSV_COLUMNS = [
  "line",
  "variant",
  "status",
  "gene",
  "transcript",
  "hgvs_p",
  "probability",
  "decision",
  "lightgbm",
  "dual_tower",
  "consequence",
  "reason",
] as const;

function cell(value: string | number | null): string {
  if (value === null) {
    return "";
  }
  if (typeof value === "number") {
    return value.toFixed(6);
  }
  const safe = /^[=+\-@\t\r]/.test(value) ? `'${value}` : value;
  return /[",\r\n]/.test(safe) ? `"${safe.replaceAll('"', '""')}"` : safe;
}

function row(item: BatchItem): (string | number | null)[] {
  const result = item.result;
  if (result.status === "out_of_scope") {
    return [
      String(item.line),
      result.variant,
      "out_of_scope",
      result.gene,
      null,
      null,
      null,
      null,
      null,
      null,
      result.consequence,
      result.reason,
    ];
  }
  return [
    String(item.line),
    result.variant,
    "predicted",
    result.gene,
    result.transcript_id,
    result.hgvs_p,
    result.probability,
    result.decision,
    result.lightgbm,
    result.dual_tower,
    result.consequence,
    null,
  ];
}

export function batchCsv(items: BatchItem[]): string {
  const lines = [CSV_COLUMNS.join(","), ...items.map((item) => row(item).map(cell).join(","))];
  return `${lines.join("\r\n")}\r\n`;
}

export function decisionCounts(items: BatchItem[]): Record<Decision | "out_of_scope", number> {
  const counts = { likely_pathogenic: 0, uncertain: 0, likely_benign: 0, out_of_scope: 0 };
  for (const { result } of items) {
    counts[result.status === "out_of_scope" ? "out_of_scope" : result.decision] += 1;
  }
  return counts;
}
