import type { components } from "./schema";

type Schemas = components["schemas"];
export type PredictResponse = Schemas["PredictResponse"];
export type VariantResult = PredictResponse["results"][number];
export type PredictedVariant = Schemas["PredictedVariant"];
export type OutOfScopeVariant = Schemas["OutOfScopeVariant"];
export type Decision = PredictedVariant["decision"];
export type EvidenceGroup = Schemas["Contribution"]["group"];
export type BatchResponse = Schemas["BatchResponse"];
export type BatchItem = BatchResponse["results"][number];
export type ModelCard = Schemas["ModelCard"];

const API_URL = import.meta.env.VITE_API_URL ?? "http://127.0.0.1:8000";

export class ApiError extends Error {
  constructor(
    message: string,
    readonly status: number,
  ) {
    super(message);
  }
}

async function errorMessage(response: Response): Promise<string> {
  try {
    const body: { detail?: unknown } = await response.json();
    if (typeof body.detail === "string") {
      return body.detail;
    }
    if (Array.isArray(body.detail)) {
      return body.detail.map((item: { msg?: string }) => item.msg ?? "invalid input").join("; ");
    }
  } catch {
    return response.statusText || `request failed (${response.status})`;
  }
  return `request failed (${response.status})`;
}

async function send<T>(path: string, init?: RequestInit): Promise<T> {
  let response: Response;
  try {
    response = await fetch(`${API_URL}${path}`, init);
  } catch {
    throw new ApiError(`cannot reach the API at ${API_URL}; is \`healthcare serve\` running?`, 0);
  }
  if (!response.ok) {
    throw new ApiError(await errorMessage(response), response.status);
  }
  return (await response.json()) as T;
}

export function predictVariant(variant: string): Promise<PredictResponse> {
  const request: Schemas["PredictRequest"] = { variant };
  return send("/api/predict", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(request),
  });
}

export function predictBatch(vcf: Blob): Promise<BatchResponse> {
  return send("/api/predict/batch", {
    method: "POST",
    headers: { "Content-Type": "text/plain" },
    body: vcf,
  });
}

export function fetchModelCard(): Promise<ModelCard> {
  return send("/api/model-card");
}
