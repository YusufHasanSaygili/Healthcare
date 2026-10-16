import type { CSSProperties } from "react";
import type { Decision } from "../api/client";
import { LIKELY_BENIGN_MAX, LIKELY_PATHOGENIC_MIN, formatPercent } from "../format";
import { fill } from "../i18n";
import { useSettings } from "../settings";

const TICKS = [0, LIKELY_BENIGN_MAX, LIKELY_PATHOGENIC_MIN, 1];

export function ProbabilityGauge({
  probability,
  decision,
}: {
  probability: number;
  decision: Decision;
}) {
  const { locale, messages } = useSettings();
  const label = fill(messages.gaugeLabel, {
    value: formatPercent(probability, locale),
    band: messages.decisions[decision].toLocaleLowerCase(locale),
  });
  return (
    <figure className="gauge" role="img" aria-label={label}>
      <div
        className="gauge-track"
        style={
          {
            "--benign-end": `${LIKELY_BENIGN_MAX * 100}%`,
            "--pathogenic-start": `${LIKELY_PATHOGENIC_MIN * 100}%`,
          } as CSSProperties
        }
      >
        <span className="gauge-marker" data-decision={decision} style={{ left: `${probability * 100}%` }} />
      </div>
      <div className="gauge-ticks" aria-hidden="true">
        {TICKS.map((tick) => (
          <span key={tick} style={{ left: `${tick * 100}%` }}>
            {formatPercent(tick, locale)}
          </span>
        ))}
      </div>
      <figcaption className="gauge-legend" aria-hidden="true">
        <span data-decision="likely_benign">{messages.gaugeBenign}</span>
        <span data-decision="uncertain">{messages.gaugeUncertain}</span>
        <span data-decision="likely_pathogenic">{messages.gaugePathogenic}</span>
      </figcaption>
    </figure>
  );
}
