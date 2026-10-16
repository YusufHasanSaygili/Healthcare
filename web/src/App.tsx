import { useActionState } from "react";
import {
  ApiError,
  predictVariant,
  type Decision,
  type EvidenceGroup,
  type OutOfScopeVariant,
  type PredictResponse,
  type PredictedVariant,
  type VariantResult,
} from "./api/client";

const DISCLAIMER =
  "For research and educational use only; not for clinical diagnosis or treatment decisions.";

const DECISION_LABELS: Record<Decision, string> = {
  likely_pathogenic: "Likely pathogenic",
  uncertain: "Uncertain",
  likely_benign: "Likely benign",
};

const GROUP_LABELS: Record<EvidenceGroup, string> = {
  allele_frequency: "Allele frequency",
  conservation: "Conservation",
  esm_scores: "Protein language model (ESM-2)",
  gene_constraint: "Gene constraint",
  protein_regions: "Protein regions",
  amino_acid_change: "Amino acid change",
  protein_position: "Position in protein",
};

type QueryState = {
  query: string;
  response: PredictResponse | null;
  error: string | null;
};

const INITIAL_STATE: QueryState = { query: "", response: null, error: null };

async function runQuery(_previous: QueryState, form: FormData): Promise<QueryState> {
  const query = String(form.get("variant") ?? "").trim();
  if (!query) {
    return { query, response: null, error: "Enter a variant." };
  }
  try {
    return { query, response: await predictVariant(query), error: null };
  } catch (error) {
    const message = error instanceof ApiError ? error.message : "Unexpected error.";
    return { query, response: null, error: message };
  }
}

function formatNumber(value: number | null, digits = 3): string {
  return value === null ? "n/a" : value.toPrecision(digits);
}

function PredictedResult({ result }: { result: PredictedVariant }) {
  const { evidence } = result;
  return (
    <article>
      <h2>
        {result.gene} {result.hgvs_p}
      </h2>
      <p>
        {result.variant} ({result.transcript_id})
      </p>
      <p>
        <strong>{DECISION_LABELS[result.decision]}</strong>, probability{" "}
        {result.probability.toFixed(3)}
      </p>
      <h3>Evidence that moved the prediction most</h3>
      <ol>
        {result.contributions.slice(0, 4).map((contribution) => (
          <li key={contribution.group}>
            {GROUP_LABELS[contribution.group]}: {contribution.log_odds > 0 ? "+" : ""}
            {contribution.log_odds.toFixed(2)} log-odds (
            {contribution.log_odds > 0 ? "toward pathogenic" : "toward benign"})
          </li>
        ))}
      </ol>
      <h3>Measurements</h3>
      <dl>
        <dt>gnomAD allele frequency</dt>
        <dd>{formatNumber(evidence.gnomad_af)}</dd>
        <dt>phyloP conservation</dt>
        <dd>{formatNumber(evidence.phylop)}</dd>
        <dt>ESM-2 log-likelihood ratio</dt>
        <dd>{formatNumber(evidence.esm_llr)}</dd>
        <dt>Grantham distance</dt>
        <dd>{formatNumber(evidence.grantham)}</dd>
        <dt>Protein regions</dt>
        <dd>{evidence.protein_regions.join(", ") || "none"}</dd>
      </dl>
    </article>
  );
}

function OutOfScopeResult({ result }: { result: OutOfScopeVariant }) {
  return (
    <article>
      <h2>{result.variant}</h2>
      <p>
        Out of scope ({result.consequence.replaceAll("_", " ")}): {result.reason}.
      </p>
    </article>
  );
}

function Result({ result }: { result: VariantResult }) {
  return result.status === "predicted" ? (
    <PredictedResult result={result} />
  ) : (
    <OutOfScopeResult result={result} />
  );
}

export function App() {
  const [state, submit, pending] = useActionState(runQuery, INITIAL_STATE);
  return (
    <main>
      <h1>Variant pathogenicity</h1>
      <p>{DISCLAIMER}</p>
      <form action={submit}>
        <label htmlFor="variant">Variant</label>{" "}
        <input
          id="variant"
          name="variant"
          defaultValue={state.query}
          placeholder="17:43063931 G>A or BRCA1 p.Arg1699Trp"
          maxLength={100}
          autoComplete="off"
          spellCheck={false}
        />{" "}
        <button type="submit" disabled={pending}>
          {pending ? "Predicting…" : "Predict"}
        </button>
      </form>
      <section aria-live="polite">
        {state.error && <p role="alert">{state.error}</p>}
        {state.response?.results.map((result) => (
          <Result key={result.variant} result={result} />
        ))}
      </section>
    </main>
  );
}
