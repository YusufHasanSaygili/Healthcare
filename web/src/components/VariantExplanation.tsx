import { useId } from "react";
import type { PredictedVariant } from "../api/client";
import { type Direction, parseChange, type Reason, reasons, shownDirections } from "../explain";
import { formatFrequency } from "../format";
import { fill } from "../i18n";
import { useSettings } from "../settings";

function ReasonList({ items, values }: { items: Reason[]; values: Record<string, string> }) {
  const { messages } = useSettings();
  return (
    <ul className="explain-reasons">
      {items.map((reason) => (
        <li key={reason.key} data-direction={reason.direction}>
          {fill(messages.explainReasons[reason.key], values)}
        </li>
      ))}
    </ul>
  );
}

export function VariantExplanation({ result }: { result: PredictedVariant }) {
  const { locale, messages } = useSettings();
  const titleId = useId();
  const change = parseChange(result.variant, result.hgvs_p);
  const found = reasons(result.evidence);
  const directions = shownDirections(result.decision);
  const frequency = formatFrequency(result.evidence.gnomad_af, locale);
  const values = {
    frequency:
      frequency.kind === "one_in"
        ? fill(messages.frequencyOneIn, { value: frequency.value })
        : frequency.kind === "percent"
          ? frequency.value
          : "",
    regions: result.evidence.protein_regions
      .map((region) => (messages.regions[region] ?? region).toLocaleLowerCase(locale))
      .join(", "),
  };
  const groups = directions.map((direction): [Direction, Reason[]] => [
    direction,
    found.filter((reason) => reason.direction === direction),
  ]);
  const shown = groups.some(([, items]) => items.length > 0);
  const name = (code: string) => messages.aminoAcids[code] ?? code;
  return (
    <section className="card-section explain" aria-labelledby={titleId}>
      <h3 id={titleId}>{messages.explainTitle}</h3>
      {change && (
        <p>
          {fill(messages.explainVariant, {
            chrom: change.chrom,
            position: Number(change.position).toLocaleString(locale),
            ref: change.ref,
            alt: change.alt,
          })}{" "}
          {fill(messages.explainProtein, {
            gene: result.gene,
            residue: change.residue,
            aaRef: name(change.aaRef),
            aaAlt: name(change.aaAlt),
          })}
        </p>
      )}
      <h4>{messages.explainWhy[result.decision]}</h4>
      {!shown && <p>{messages.explainNoReasons}</p>}
      {shown &&
        groups.map(([direction, items]) =>
          items.length === 0 ? null : (
            <div key={direction}>
              {directions.length > 1 && (
                <p className="explain-toward" data-direction={direction}>
                  {messages.explainToward[direction]}
                </p>
              )}
              <ReasonList items={items} values={values} />
            </div>
          ),
        )}
      {result.gene_function && (
        <div className="gene-function">
          <h4>{fill(messages.geneFunctionTitle, { gene: result.gene })}</h4>
          <p lang="en">{result.gene_function}</p>
          <p className="field-hint">{messages.geneFunctionSource}</p>
        </div>
      )}
    </section>
  );
}
