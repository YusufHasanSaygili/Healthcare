import { createContext, use, useEffect, useState, type ReactNode } from "react";
import type { Locale } from "./format";
import { MESSAGES } from "./i18n";

export type Theme = "dark" | "light";

type Settings = {
  locale: Locale;
  theme: Theme;
  messages: (typeof MESSAGES)[Locale];
  setLocale: (locale: Locale) => void;
  toggleTheme: () => void;
};

const SettingsContext = createContext<Settings | null>(null);

function readStored<T extends string>(key: string, allowed: readonly T[]): T | null {
  try {
    const value = localStorage.getItem(key);
    return allowed.find((item) => item === value) ?? null;
  } catch {
    return null;
  }
}

function store(key: string, value: string): void {
  try {
    localStorage.setItem(key, value);
  } catch {
    return;
  }
}

function initialTheme(): Theme {
  const stored = readStored<Theme>("theme", ["dark", "light"]);
  if (stored) {
    return stored;
  }
  return window.matchMedia("(prefers-color-scheme: light)").matches ? "light" : "dark";
}

function initialLocale(): Locale {
  return readStored<Locale>("locale", ["en", "tr"]) ?? "en";
}

export function SettingsProvider({ children }: { children: ReactNode }) {
  const [theme, setTheme] = useState<Theme>(initialTheme);
  const [locale, setLocale] = useState<Locale>(initialLocale);

  useEffect(() => {
    document.documentElement.dataset.theme = theme;
    store("theme", theme);
  }, [theme]);

  useEffect(() => {
    document.documentElement.lang = locale;
    store("locale", locale);
  }, [locale]);

  const settings: Settings = {
    locale,
    theme,
    messages: MESSAGES[locale],
    setLocale,
    toggleTheme: () => setTheme((current) => (current === "dark" ? "light" : "dark")),
  };
  return <SettingsContext value={settings}>{children}</SettingsContext>;
}

export function useSettings(): Settings {
  const settings = use(SettingsContext);
  if (!settings) {
    throw new Error("useSettings must be used inside SettingsProvider");
  }
  return settings;
}
