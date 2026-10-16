import {
  barShare,
  formatDecimal,
  strongestSignals,
  type Contribution,
  type Locale,
} from "../format";
import { fill } from "../i18n";
import { useSettings } from "../settings";

function signed(value: number, text: string | null): string {
  return value > 0 ? `+${text}` : `${text}`;
}

function lowerFirst(text: string, locale: Locale): string {
  return text.charAt(0).toLocaleLowerCase(locale) + text.slice(1);
}

function summarize(
  signals: Contribution[],
  messages: ReturnType<typeof useSettings>["messages"],
  locale: Locale,
): string | null {
  const [first, second] = signals;
  if (!first) {
    return null;
  }
  const name = (c: Contribution) => lowerFirst(messages.groups[c.group], locale);
  const direction = (c: Contribution) =>
    c.log_odds > 0 ? messages.towardPathogenic : messages.towardBenign;
  if (second && first.log_odds > 0 !== second.log_odds > 0) {
    const sentence = fill(messages.strongestMixed, {
      first: name(first),
      firstDirection: direction(first),
      second: name(second),
      secondDirection: direction(second),
    });
    return sentence.charAt(0).toLocaleUpperCase(locale) + sentence.slice(1);
  }
  const list = signals.map(name).join(` ${messages.and} `);
  return fill(messages.strongestOneWay, { direction: direction(first), list });
}

export function ContributionChart({ contributions }: { contributions: Contribution[] }) {
  const { locale, messages } = useSettings();
  const direction = (value: number) =>
    value > 0 ? messages.towardPathogenic : messages.towardBenign;
  const summary = summarize(strongestSignals(contributions), messages, locale);
  return (
    <div className="contributions">
      {summary && <p className="contributions-summary">{summary}</p>}
      <div className="contributions-axis" aria-hidden="true">
        <span>← {messages.towardBenign}</span>
        <span>{messages.towardPathogenic} →</span>
      </div>
      <ul className="contributions-list">
        {contributions.map((contribution) => {
          const value = signed(
            contribution.log_odds,
            formatDecimal(contribution.log_odds, locale),
          );
          const side = contribution.log_odds > 0 ? "pathogenic" : "benign";
          return (
            <li
              key={contribution.group}
              className="contribution"
              aria-label={`${messages.groups[contribution.group]}: ${value}, ${direction(contribution.log_odds)}`}
            >
              <span className="contribution-label" aria-hidden="true">
                {messages.groups[contribution.group]}
              </span>
              <span className="contribution-track" aria-hidden="true">
                <span
                  className="contribution-bar"
                  data-side={side}
                  style={{ width: `${barShare(contribution.log_odds, contributions) * 50}%` }}
                />
              </span>
              <span className="contribution-value" data-side={side} aria-hidden="true">
                {value}
              </span>
            </li>
          );
        })}
      </ul>
      <p className="note">{messages.whyNote}</p>
    </div>
  );
}
