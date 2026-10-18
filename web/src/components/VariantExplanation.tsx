import { useId } from "react";
import type { PredictedVariant } from "../api/client";
import {
  type Direction,
  type Genomic,
  parseGenomic,
  parseSubstitution,
  type Reason,
  reasons,
  residueOf,
  shownDirections,
} from "../explain";
import { formatFrequency } from "../format";
import { fill, type Messages } from "../i18n";
import { useSettings } from "../settings";

function dnaSentence(change: Genomic, messages: Messages, locale: string): string {
  const position = Number(change.position).toLocaleString(locale);
  if (change.ref.length === 1 && change.alt.length === 1) {
    return fill(messages.explainVariant, {
      chrom: change.chrom,
      position,
      ref: change.ref,
      alt: change.alt,
    });
  }
  const shared = change.ref[0] === change.alt[0] ? 1 : 0;
  const deleted = change.ref.length - shared;
  const inserted = change.alt.length - shared;
  const text =
    inserted === 0
      ? fill(messages.indelDeleted, { count: String(deleted) })
      : deleted === 0
        ? fill(messages.indelInserted, { count: String(inserted) })
        : fill(messages.indelReplaced, { deleted: String(deleted), inserted: String(inserted) });
  return fill(messages.explainIndel, { chrom: change.chrom, position, change: text });
}

function proteinSentence(result: PredictedVariant, messages: Messages): string | null {
  if (result.consequence === "missense") {
    const substitution = parseSubstitution(result.hgvs_p);
    if (!substitution) {
      return null;
    }
    return fill(messages.explainProtein, {
      gene: result.gene,
      residue: substitution.residue,
      aaRef: messages.aminoAcids[substitution.aaRef] ?? substitution.aaRef,
      aaAlt: messages.aminoAcids[substitution.aaAlt] ?? substitution.aaAlt,
    });
  }
  return fill(messages.explainType[result.consequence], {
    gene: result.gene,
    residue: residueOf(result.hgvs_p),
    hgvs: result.hgvs_p ?? "",
  });
}

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
  const change = parseGenomic(result.variant);
  const protein = proteinSentence(result, messages);
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
  return (
    <section className="card-section explain" aria-labelledby={titleId}>
      <h3 id={titleId}>{messages.explainTitle}</h3>
      {change && (
        <p>
          {dnaSentence(change, messages, locale)} {protein}
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
