import { useId, type ReactNode } from "react";
import type { Decision, OutOfScopeVariant, PredictedVariant, VariantResult } from "../api/client";
import { formatDecimal, formatFrequency, formatPercent, type Frequency } from "../format";
import { fill } from "../i18n";
import { useSettings } from "../settings";
import { ContributionChart } from "./ContributionChart";
import { BlockedIcon, ClearIcon, HazardIcon, UnknownIcon } from "./Icons";
import { ProbabilityGauge } from "./ProbabilityGauge";

const DECISION_ICONS: Record<Decision, typeof HazardIcon> = {
  likely_pathogenic: HazardIcon,
  uncertain: UnknownIcon,
  likely_benign: ClearIcon,
};

function Measurement({ label, hint, children }: { label: string; hint?: string; children: ReactNode }) {
  return (
    <div className="measurement">
      <dt>{label}</dt>
      <dd>
        <span className="measurement-value">{children}</span>
        {hint && <span className="measurement-hint">{hint}</span>}
      </dd>
    </div>
  );
}

function describeFrequency(
  frequency: Frequency,
  messages: ReturnType<typeof useSettings>["messages"],
): string {
  switch (frequency.kind) {
    case "not_measured":
      return messages.frequencyNotMeasured;
    case "not_observed":
      return messages.frequencyNotObserved;
    case "percent":
      return frequency.value;
    case "one_in":
      return fill(messages.frequencyOneIn, { value: frequency.value });
  }
}

function Measurements({ result }: { result: PredictedVariant }) {
  const { locale, messages } = useSettings();
  const { evidence } = result;
  const frequencyText = describeFrequency(formatFrequency(evidence.gnomad_af, locale), messages);
  const number = (value: number | null, digits?: number) =>
    formatDecimal(value, locale, digits) ?? messages.notAvailable;
  return (
    <dl className="measurements">
      <Measurement label={messages.frequencyLabel}>{frequencyText}</Measurement>
      <Measurement label={messages.conservationLabel} hint={messages.conservationHint}>
        {number(evidence.phylop)}
      </Measurement>
      <Measurement label={messages.esmLabel} hint={messages.esmHint}>
        {number(evidence.esm_llr, 1)}
      </Measurement>
      <Measurement label={messages.granthamLabel} hint={messages.granthamHint}>
        {number(evidence.grantham, 0)}
      </Measurement>
      <Measurement label={messages.regionsLabel}>
        {evidence.protein_regions.length > 0
          ? evidence.protein_regions.map((region) => messages.regions[region] ?? region).join(", ")
          : messages.regionsNone}
      </Measurement>
    </dl>
  );
}

function PredictedCard({ result }: { result: PredictedVariant }) {
  const { locale, messages } = useSettings();
  const DecisionIcon = DECISION_ICONS[result.decision];
  const whyId = useId();
  const measurementsId = useId();
  return (
    <article className="card hud" data-decision={result.decision}>
      <header className="card-header">
        <p className="specimen">
          {messages.specimen} · <span translate="no">{result.variant} · {result.transcript_id}</span>
        </p>
        <h2 className="card-title" translate="no">
          {result.gene} <span className="hgvs">{result.hgvs_p}</span>
        </h2>
      </header>
      <div className="verdict">
        <p className="decision" data-decision={result.decision}>
          <DecisionIcon width={28} height={28} />
          {messages.decisions[result.decision]}
        </p>
        <p className="probability">
          <span className="probability-value">{formatPercent(result.probability, locale)}</span>
          <span className="probability-label">{messages.probabilityLabel}</span>
        </p>
      </div>
      <ProbabilityGauge probability={result.probability} decision={result.decision} />
      <section className="card-section" aria-labelledby={whyId}>
        <h3 id={whyId}>{messages.whyTitle}</h3>
        <ContributionChart contributions={result.contributions} />
      </section>
      <section className="card-section" aria-labelledby={measurementsId}>
        <h3 id={measurementsId}>{messages.measurementsTitle}</h3>
        <Measurements result={result} />
      </section>
      <p className="note">
        {fill(messages.members, {
          lightgbm: formatPercent(result.lightgbm, locale),
          dualTower: formatPercent(result.dual_tower, locale),
        })}
      </p>
    </article>
  );
}

function OutOfScopeCard({ result }: { result: OutOfScopeVariant }) {
  const { messages } = useSettings();
  const showGene = result.gene && !result.variant.toUpperCase().startsWith(result.gene.toUpperCase());
  return (
    <article className="card hud" data-decision="out_of_scope">
      <header className="card-header">
        <p className="specimen">
          {messages.specimen} ·{" "}
          <span translate="no">
            {result.variant}
            {showGene ? ` · ${result.gene}` : ""}
          </span>
        </p>
        <h2 className="card-title decision" data-decision="out_of_scope">
          <BlockedIcon width={28} height={28} />
          {messages.outOfScopeTitle}
        </h2>
      </header>
      <p>{messages.consequences[result.consequence] ?? result.reason}</p>
    </article>
  );
}

export function ResultCard({ result }: { result: VariantResult }) {
  return result.status === "predicted" ? (
    <PredictedCard result={result} />
  ) : (
    <OutOfScopeCard result={result} />
  );
}
