import { startTransition, useActionState, useState } from "react";
import { ApiError, predictVariant, type PredictResponse } from "./api/client";
import { BioBackground, type Lab } from "./components/BioBackground";
import { BatchPage } from "./components/BatchPage";
import { Header, PageNav } from "./components/Header";
import { ModelCardPage } from "./components/ModelCardPage";
import { QueryPanel } from "./components/QueryPanel";
import { ResultCard } from "./components/ResultCard";
import { useRoute } from "./route";
import { useSettings } from "./settings";

type QueryError = { kind: "invalid"; detail: string } | { kind: "network" } | { kind: "unexpected" };

type QueryState = {
  response: PredictResponse | null;
  error: QueryError | null;
};

const INITIAL_STATE: QueryState = { response: null, error: null };

async function runQuery(previous: QueryState, form: FormData): Promise<QueryState> {
  const query = String(form.get("variant") ?? "").trim();
  if (!query) {
    return previous;
  }
  try {
    return { response: await predictVariant(query), error: null };
  } catch (error) {
    if (error instanceof ApiError) {
      return {
        response: null,
        error: error.status === 0 ? { kind: "network" } : { kind: "invalid", detail: error.message },
      };
    }
    return { response: null, error: { kind: "unexpected" } };
  }
}

function labState(pending: boolean, state: QueryState): Lab {
  if (pending) {
    return { status: "analyzing", label: null, probability: null };
  }
  if (state.error) {
    return { status: "error", label: null, probability: null };
  }
  const results = state.response?.results ?? [];
  const predicted = results.find((result) => result.status === "predicted");
  if (predicted?.status === "predicted") {
    return {
      status: predicted.decision,
      label: `${predicted.gene} ${predicted.hgvs_p}`,
      probability: predicted.probability,
    };
  }
  const first = results[0];
  return first
    ? { status: "out_of_scope", label: first.variant, probability: null }
    : { status: "idle", label: null, probability: null };
}

function ErrorBanner({ error }: { error: QueryError }) {
  const { messages } = useSettings();
  if (error.kind === "invalid") {
    return (
      <div className="alert" role="alert">
        <strong>{messages.invalidInputTitle}</strong>
        <span>{error.detail}</span>
      </div>
    );
  }
  return (
    <div className="alert" role="alert">
      {error.kind === "network" ? messages.networkError : messages.unexpectedError}
    </div>
  );
}

function AwaitingPanel() {
  const { messages } = useSettings();
  return (
    <div className="awaiting">
      <p className="awaiting-title">{messages.awaitingTitle}</p>
      <p>{messages.awaitingHint}</p>
    </div>
  );
}

export function App() {
  const { messages } = useSettings();
  const [value, setValue] = useState("");
  const [state, submit, pending] = useActionState(runQuery, INITIAL_STATE);

  function runExample(example: string) {
    setValue(example);
    const form = new FormData();
    form.set("variant", example);
    startTransition(() => submit(form));
  }

  const route = useRoute();
  const results = state.response?.results ?? [];
  const lab: Lab =
    route === "single"
      ? labState(pending, state)
      : { status: "idle", label: null, probability: null };
  return (
    <>
      <BioBackground status={lab.status} label={lab.label} probability={lab.probability} />
      <div className="page">
        <Header />
        <PageNav route={route} />
        {route === "batch" && <BatchPage />}
        {route === "model" && <ModelCardPage />}
        {route === "single" && (
          <main className="layout">
            <QueryPanel
              action={submit}
              value={value}
              pending={pending}
              onValueChange={setValue}
              onExample={runExample}
            />
            <section className="results" aria-live="polite" aria-busy={pending}>
              {pending && <div className="scan" aria-hidden="true" />}
              {state.error && <ErrorBanner error={state.error} />}
              {!state.error && results.length === 0 && !pending && <AwaitingPanel />}
              {results.map((result) => (
                <ResultCard key={result.variant} result={result} />
              ))}
            </section>
          </main>
        )}
        <footer className="footer">
          <p>{messages.disclaimer}</p>
          <p className="footer-meta">{messages.footer}</p>
        </footer>
      </div>
    </>
  );
}
