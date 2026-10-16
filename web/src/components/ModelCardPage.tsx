import { useEffect, useState } from "react";
import { ApiError, fetchModelCard, type ModelCard } from "../api/client";
import { collapseChromosomes, formatDecimal, formatPercent } from "../format";
import { fill } from "../i18n";
import { useSettings } from "../settings";
import { HazardIcon } from "./Icons";

type CardState = { card: ModelCard } | { error: string } | null;

function Evaluation({ card }: { card: ModelCard }) {
  const { locale, messages } = useSettings();
  const decimal = (value: number) => formatDecimal(value, locale, 4);
  return (
    <section className="card-section" aria-labelledby="evaluation-title">
      <h3 id="evaluation-title">
        {messages.evaluationTitle} · <span translate="no">{card.variant_types.join(", ")}</span>
      </h3>
      <p className="field-hint">{messages.evaluationNote}</p>
      <div className="table-scroll">
        <table className="data-table">
          <thead>
            <tr>
              <th scope="col">{messages.metrics.set}</th>
              <th scope="col">{messages.metrics.variants}</th>
              <th scope="col">{messages.metrics.rocAuc}</th>
              <th scope="col">{messages.metrics.prAuc}</th>
              <th scope="col">{messages.metrics.ece}</th>
              <th scope="col">{messages.metrics.uncertain}</th>
            </tr>
          </thead>
          <tbody>
            {card.evaluation.map((row) => (
              <tr key={row.name}>
                <th scope="row">{row.name}</th>
                <td className="numeric">{row.variants.toLocaleString(locale)}</td>
                <td className="numeric">
                  {decimal(row.roc_auc)}{" "}
                  <span className="muted">
                    [{decimal(row.roc_auc_ci[0])}–{decimal(row.roc_auc_ci[1])}]
                  </span>
                </td>
                <td className="numeric">{decimal(row.pr_auc)}</td>
                <td className="numeric">{decimal(row.ece)}</td>
                <td className="numeric">{formatPercent(row.uncertain_share, locale)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </section>
  );
}

function Ablations({ card }: { card: ModelCard }) {
  const { locale, messages } = useSettings();
  return (
    <section className="card-section" aria-labelledby="ablations-title">
      <h3 id="ablations-title">{messages.ablationsTitle}</h3>
      <p className="field-hint">{messages.ablationsNote}</p>
      <div className="table-scroll">
        <table className="data-table">
          <thead>
            <tr>
              <th scope="col">{messages.ablation}</th>
              <th scope="col">{messages.metrics.rocAuc}</th>
            </tr>
          </thead>
          <tbody>
            {card.ablations.map((row) => (
              <tr key={row.name}>
                <th scope="row">{row.name}</th>
                <td className="numeric">{formatDecimal(row.roc_auc, locale, 4)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </section>
  );
}

function CardBody({ card }: { card: ModelCard }) {
  const { locale, messages } = useSettings();
  return (
    <article className="card hud">
      <header className="card-header">
        <p className="specimen">
          {messages.modelFile} · <span translate="no">{card.model}</span>
        </p>
        <h2 className="card-title">{messages.modelTitle}</h2>
        {messages.sourceLanguageNote && <p className="field-hint">{messages.sourceLanguageNote}</p>}
      </header>
      <section className="card-section" aria-labelledby="thresholds-title">
        <h3 id="thresholds-title">{messages.thresholdsTitle}</h3>
        <p>
          {fill(messages.thresholdsText, {
            high: formatPercent(card.thresholds.likely_pathogenic_min, locale),
            low: formatPercent(card.thresholds.likely_benign_max, locale),
          })}
        </p>
      </section>
      <Evaluation card={card} />
      <Ablations card={card} />
      <section className="card-section" aria-labelledby="data-title">
        <h3 id="data-title">{messages.dataTitle}</h3>
        <dl className="data-versions">
          {collapseChromosomes(card.data_versions).map(([source, version]) => (
            <div key={source}>
              <dt translate="no">{source}</dt>
              <dd translate="no">{version}</dd>
            </div>
          ))}
        </dl>
      </section>
      <section className="card-section" aria-labelledby="limitations-title">
        <h3 id="limitations-title">{messages.limitationsTitle}</h3>
        <ul className="limitations">
          {card.limitations.map((limitation) => (
            <li key={limitation}>{limitation}</li>
          ))}
        </ul>
      </section>
      <p className="disclaimer">
        <HazardIcon width={18} height={18} />
        {messages.disclaimer}
      </p>
    </article>
  );
}

export function ModelCardPage() {
  const { messages } = useSettings();
  const [state, setState] = useState<CardState>(null);

  useEffect(() => {
    let active = true;
    fetchModelCard()
      .then((card) => active && setState({ card }))
      .catch((error: unknown) => {
        if (active) {
          setState({ error: error instanceof ApiError && error.status !== 0 ? error.message : "" });
        }
      });
    return () => {
      active = false;
    };
  }, []);

  return (
    <main className="layout" aria-busy={state === null}>
      {state === null && <p className="awaiting">{messages.modelLoading}</p>}
      {state && "error" in state && (
        <div className="alert" role="alert">
          {state.error || messages.networkError}
        </div>
      )}
      {state && "card" in state && <CardBody card={state.card} />}
    </main>
  );
}
