import { useSyncExternalStore } from "react";

export type Route = "single" | "batch" | "model";

export const ROUTE_HASHES: Record<Route, string> = {
  single: "#/",
  batch: "#/batch",
  model: "#/model",
};

export function routeFromHash(hash: string): Route {
  if (hash === ROUTE_HASHES.batch) {
    return "batch";
  }
  return hash === ROUTE_HASHES.model ? "model" : "single";
}

function subscribe(onChange: () => void): () => void {
  window.addEventListener("hashchange", onChange);
  return () => window.removeEventListener("hashchange", onChange);
}

export function useRoute(): Route {
  return useSyncExternalStore(subscribe, () => routeFromHash(window.location.hash));
}
