import { useSettings } from "../settings";
import { HazardIcon } from "./Icons";

const EXAMPLES = [
  "BRCA1 p.Arg1699Trp",
  "TP53 p.Arg175His",
  "BRCA1 p.Ser1613Gly",
  "13:32316497 G>T",
  "17:43124027 ACT>A",
  "7:117559592 CTTT>C",
  "17:7669612 G>A",
];

export function QueryPanel({
  action,
  value,
  pending,
  onValueChange,
  onExample,
}: {
  action: (form: FormData) => void;
  value: string;
  pending: boolean;
  onValueChange: (value: string) => void;
  onExample: (example: string) => void;
}) {
  const { messages } = useSettings();
  return (
    <section className="panel hud" aria-labelledby="query-label">
      <form action={action} className="query-form">
        <label id="query-label" htmlFor="variant" className="field-label">
          {messages.inputLabel}
        </label>
        <div className="query-row">
          <input
            id="variant"
            name="variant"
            className="query-input"
            value={value}
            onChange={(event) => onValueChange(event.target.value)}
            placeholder="BRCA1 p.Arg1699Trp…"
            aria-describedby="variant-hint"
            maxLength={100}
            autoComplete="off"
            autoCapitalize="off"
            spellCheck={false}
            required
          />
          <button type="submit" className="primary-button" disabled={pending}>
            {pending ? messages.submitting : messages.submit}
          </button>
        </div>
        <p id="variant-hint" className="field-hint">
          {messages.inputHint}
        </p>
      </form>
      <div className="examples">
        <span className="examples-label">{messages.examplesLabel}</span>
        {EXAMPLES.map((example) => (
          <button
            key={example}
            type="button"
            className="chip"
            translate="no"
            disabled={pending}
            onClick={() => onExample(example)}
          >
            {example}
          </button>
        ))}
      </div>
      <p className="disclaimer">
        <HazardIcon width={18} height={18} />
        {messages.disclaimer}
      </p>
    </section>
  );
}
