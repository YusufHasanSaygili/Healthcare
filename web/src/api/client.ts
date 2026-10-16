import type { components } from "./schema";

type Schemas = components["schemas"];
export type PredictResponse = Schemas["PredictResponse"];
export type VariantResult = PredictResponse["results"][number];
export type PredictedVariant = Schemas["PredictedVariant"];
export type OutOfScopeVariant = Schemas["OutOfScopeVariant"];
export type Decision = PredictedVariant["decision"];
export type EvidenceGroup = Schemas["Contribution"]["group"];

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

export async function predictVariant(variant: string): Promise<PredictResponse> {
  const request: Schemas["PredictRequest"] = { variant };
  let response: Response;
  try {
    response = await fetch(`${API_URL}/api/predict`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(request),
    });
  } catch {
    throw new ApiError(`cannot reach the API at ${API_URL}; is \`healthcare serve\` running?`, 0);
  }
  if (!response.ok) {
    throw new ApiError(await errorMessage(response), response.status);
  }
  return (await response.json()) as PredictResponse;
}
