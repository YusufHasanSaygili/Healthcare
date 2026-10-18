import { useActionState, useId } from "react";
import { ApiError, type BatchItem, type BatchResponse, predictBatch } from "../api/client";
import {
  API_PROBLEMS,
  batchCsv,
  decisionCounts,
  type UploadProblem,
  uploadProblem,
} from "../batch";
import { formatPercent } from "../format";
import { fill } from "../i18n";
import { useSettings } from "../settings";
import { HazardIcon } from "./Icons";

type BatchError =
  | { kind: "upload"; problem: UploadProblem }
  | { kind: "api"; detail: string }
  | { kind: "network" };

type BatchState = { response: BatchResponse | null; error: BatchError | null };

const INITIAL_STATE: BatchState = { response: null, error: null };

async function runBatch(_: BatchState, form: FormData): Promise<BatchState> {
  const file = form.get("vcf");
  if (!(file instanceof File)) {
    return { response: null, error: { kind: "upload", problem: "empty" } };
  }
  const problem = uploadProblem(file.name, file.size);
  if (problem) {
    return { response: null, error: { kind: "upload", problem } };
  }
  try {
    return { response: await predictBatch(file), error: null };
  } catch (error) {
    if (error instanceof ApiError && error.status !== 0) {
      const known = API_PROBLEMS[error.status];
      return {
        response: null,
        error: known ? { kind: "upload", problem: known } : { kind: "api", detail: error.message },
      };
    }
    return { response: null, error: { kind: "network" } };
  }
}

function downloadCsv(items: BatchItem[]) {
  const url = URL.createObjectURL(new Blob([batchCsv(items)], { type: "text/csv;charset=utf-8" }));
  const link = document.createElement("a");
  link.href = url;
  link.download = "variant-predictions.csv";
  link.click();
  URL.revokeObjectURL(url);
}

function UploadError({ error }: { error: BatchError }) {
  const { messages } = useSettings();
  const text =
    error.kind === "upload"
      ? messages.uploadProblems[error.problem]
      : error.kind === "api"
        ? error.detail
        : messages.networkError;
  return (
    <div className="alert" role="alert">
      <strong>{messages.batchRejectedTitle}</strong>
      <span>{text}</span>
    </div>
  );
}

function ResultRow({ item }: { item: BatchItem }) {
  const { locale, messages } = useSettings();
  const { result } = item;
  return (
    <tr>
      <td className="numeric">{item.line}</td>
      <td translate="no">{result.variant}</td>
      <td translate="no">{result.gene ?? "—"}</td>
      <td translate="no">{result.status === "predicted"
          ? (result.hgvs_p ?? messages.typeNames[result.consequence])
          : "—"}</td>
      <td className="numeric">
        {result.status === "predicted" ? formatPercent(result.probability, locale) : "—"}
      </td>
      <td>
        {result.status === "predicted" ? (
          <span className="tag" data-decision={result.decision}>
            {messages.decisions[result.decision]}
          </span>
        ) : (
          <span className="tag" data-decision="out_of_scope" title={result.reason}>
            {messages.consequences[result.consequence] ?? messages.outOfScope}
          </span>
        )}
      </td>
    </tr>
  );
}

function BatchResults({ response }: { response: BatchResponse }) {
  const { locale, messages } = useSettings();
  const counts = decisionCounts(response.results);
  const hidden = response.error_count - response.errors.length;
  return (
    <section className="batch-results hud" aria-label={messages.batchTitle}>
      <div className="batch-summary">
        <p>
          {fill(messages.batchSummary, {
            variants: String(response.variants),
            scored: String(response.results.length - counts.out_of_scope),
            skipped: String(counts.out_of_scope),
            errors: String(response.error_count),
            seconds: (response.milliseconds / 1000).toLocaleString(locale, {
              maximumFractionDigits: 1,
            }),
          })}
        </p>
        <ul className="batch-counts">
          {(["likely_pathogenic", "uncertain", "likely_benign"] as const).map((decision) => (
            <li key={decision} className="tag" data-decision={decision}>
              {messages.decisions[decision]}: {counts[decision]}
            </li>
          ))}
        </ul>
        <button
          type="button"
          className="secondary-button"
          disabled={response.results.length === 0}
          onClick={() => downloadCsv(response.results)}
        >
          {messages.downloadCsv}
        </button>
      </div>
      {response.results.length === 0 ? (
        <p className="batch-empty">{messages.noResults}</p>
      ) : (
        <div className="table-scroll">
          <table className="data-table">
            <thead>
              <tr>
                <th scope="col">{messages.columns.line}</th>
                <th scope="col">{messages.columns.variant}</th>
                <th scope="col">{messages.columns.gene}</th>
                <th scope="col">{messages.columns.protein}</th>
                <th scope="col">{messages.columns.probability}</th>
                <th scope="col">{messages.columns.decision}</th>
              </tr>
            </thead>
            <tbody>
              {response.results.map((item) => (
                <ResultRow key={`${item.line}-${item.result.variant}`} item={item} />
              ))}
            </tbody>
          </table>
        </div>
      )}
      {response.errors.length > 0 && (
        <details className="line-errors">
          <summary>
            {messages.lineErrorsTitle} ({response.error_count})
          </summary>
          <ul>
            {response.errors.map((error) => (
              <li key={`${error.line}-${error.message}`}>
                <span className="numeric">{error.line}</span> {error.message}
              </li>
            ))}
          </ul>
          {hidden > 0 && <p>{fill(messages.lineErrorsHidden, { count: String(hidden) })}</p>}
        </details>
      )}
    </section>
  );
}

export function BatchPage() {
  const { messages } = useSettings();
  const [state, submit, pending] = useActionState(runBatch, INITIAL_STATE);
  const hintId = useId();
  return (
    <main className="layout">
      <section className="panel hud" aria-labelledby="batch-title">
        <form action={submit} className="query-form">
          <h2 id="batch-title" className="field-label">
            {messages.batchTitle}
          </h2>
          <div className="query-row">
            <input
              type="file"
              name="vcf"
              accept=".vcf,.txt,text/plain"
              className="file-input"
              aria-label={messages.batchFileLabel}
              aria-describedby={hintId}
              required
            />
            <button type="submit" className="primary-button" disabled={pending}>
              {pending ? messages.batchSubmitting : messages.batchSubmit}
            </button>
          </div>
          <p id={hintId} className="field-hint">
            {messages.batchHint} {messages.batchPrivacy}{" "}
            <a href="/sample.vcf" download className="text-link">
              {messages.batchSample}
            </a>
          </p>
        </form>
        <p className="disclaimer">
          <HazardIcon width={18} height={18} />
          {messages.disclaimer}
        </p>
      </section>
      <div className="results" aria-live="polite" aria-busy={pending}>
        {pending && <div className="scan" aria-hidden="true" />}
        {state.error && <UploadError error={state.error} />}
        {state.response && <BatchResults response={state.response} />}
      </div>
    </main>
  );
}
