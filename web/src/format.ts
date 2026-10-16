import type { EvidenceGroup } from "./api/client";

export const LIKELY_BENIGN_MAX = 0.1;
export const LIKELY_PATHOGENIC_MIN = 0.9;

export type Locale = "en" | "tr";

export type Contribution = { group: EvidenceGroup; log_odds: number };

export type Frequency =
  | { kind: "not_measured" }
  | { kind: "not_observed" }
  | { kind: "percent"; value: string }
  | { kind: "one_in"; value: string };

const LOCALE_TAGS: Record<Locale, string> = { en: "en-US", tr: "tr-TR" };

export function formatPercent(probability: number, locale: Locale): string {
  if (probability > 0 && probability < 0.001) {
    return `<${formatPercent(0.001, locale)}`;
  }
  if (probability < 1 && probability > 0.999) {
    return `>${formatPercent(0.999, locale)}`;
  }
  return new Intl.NumberFormat(LOCALE_TAGS[locale], {
    style: "percent",
    maximumFractionDigits: 1,
  }).format(probability);
}

export function formatDecimal(value: number | null, locale: Locale, digits = 2): string | null {
  if (value === null) {
    return null;
  }
  return new Intl.NumberFormat(LOCALE_TAGS[locale], { maximumFractionDigits: digits }).format(
    value,
  );
}

export function formatFrequency(frequency: number | null, locale: Locale): Frequency {
  if (frequency === null) {
    return { kind: "not_measured" };
  }
  if (frequency === 0) {
    return { kind: "not_observed" };
  }
  if (frequency >= 0.01) {
    return { kind: "percent", value: formatPercent(frequency, locale) };
  }
  const oneIn = Math.round(1 / frequency);
  return { kind: "one_in", value: new Intl.NumberFormat(LOCALE_TAGS[locale]).format(oneIn) };
}

export function barShare(value: number, contributions: Contribution[]): number {
  const largest = Math.max(...contributions.map((c) => Math.abs(c.log_odds)), 1);
  return Math.min(Math.abs(value) / largest, 1);
}

export function strongestSignals(contributions: Contribution[], count = 2): Contribution[] {
  return [...contributions]
    .filter((c) => Math.abs(c.log_odds) >= 0.1)
    .sort((a, b) => Math.abs(b.log_odds) - Math.abs(a.log_odds))
    .slice(0, count);
}
