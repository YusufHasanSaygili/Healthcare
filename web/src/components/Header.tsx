import type { Locale } from "../format";
import { ROUTE_HASHES, type Route } from "../route";
import { useSettings } from "../settings";
import { HelixMark, MoonIcon, SunIcon } from "./Icons";

const LOCALES: Locale[] = ["en", "tr"];

export function Header() {
  const { theme, toggleTheme, locale, setLocale, messages } = useSettings();
  return (
    <header className="topbar">
      <div className="brand">
        <HelixMark className="brand-mark" width={34} height={34} />
        <div>
          <h1 className="brand-title">{messages.title}</h1>
          <p className="brand-tagline">{messages.tagline}</p>
        </div>
      </div>
      <div className="controls">
        <div className="segmented" role="group" aria-label={messages.languageLabel}>
          {LOCALES.map((option) => (
            <button
              key={option}
              type="button"
              lang={option}
              aria-pressed={locale === option}
              onClick={() => setLocale(option)}
            >
              {option.toUpperCase()}
            </button>
          ))}
        </div>
        <button
          type="button"
          className="icon-button"
          onClick={toggleTheme}
          aria-label={theme === "dark" ? messages.themeToLight : messages.themeToDark}
        >
          {theme === "dark" ? <SunIcon /> : <MoonIcon />}
        </button>
      </div>
    </header>
  );
}

const ROUTES: Route[] = ["single", "batch", "model"];

export function PageNav({ route }: { route: Route }) {
  const { messages } = useSettings();
  return (
    <nav className="page-nav" aria-label={messages.navLabel}>
      {ROUTES.map((option) => (
        <a
          key={option}
          href={ROUTE_HASHES[option]}
          aria-current={route === option ? "page" : undefined}
        >
          {messages.nav[option]}
        </a>
      ))}
    </nav>
  );
}
