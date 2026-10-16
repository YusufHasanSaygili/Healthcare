import type { Decision } from "../api/client";
import { formatPercent } from "../format";
import { useSettings } from "../settings";

export type LabStatus = Decision | "idle" | "analyzing" | "out_of_scope" | "error";

const BASELINE = 92;
const LEFT = 16;
const RIGHT = 234;

function beat(x: number, height: number, width: number): string {
  const w = width / 60;
  return [
    `L${x + 8 * w},${BASELINE}`,
    `Q${x + 12 * w},${BASELINE - height * 0.12} ${x + 16 * w},${BASELINE}`,
    `L${x + 21 * w},${BASELINE}`,
    `L${x + 24 * w},${BASELINE + height * 0.12}`,
    `L${x + 28 * w},${BASELINE - height}`,
    `L${x + 32 * w},${BASELINE + height * 0.35}`,
    `L${x + 35 * w},${BASELINE}`,
    `L${x + 42 * w},${BASELINE}`,
    `Q${x + 48 * w},${BASELINE - height * 0.22} ${x + 54 * w},${BASELINE}`,
    `L${x + width},${BASELINE}`,
  ].join(" ");
}

function rhythm(beats: { width: number; height: number }[]): string {
  let x = LEFT;
  const parts = [`M${LEFT},${BASELINE}`];
  for (const item of beats) {
    parts.push(beat(x, item.height, item.width));
    x += item.width;
  }
  parts.push(`L${RIGHT},${BASELINE}`);
  return parts.join(" ");
}

function noise(amplitude: number, seed: number): string {
  const points: string[] = [];
  for (let x = LEFT; x <= RIGHT; x += 3) {
    const y = BASELINE + amplitude * Math.sin(x * 0.9 + seed) * Math.sin(x * 0.23 + seed * 2);
    points.push(`${x},${y.toFixed(1)}`);
  }
  return `M${points.join("L")}`;
}

function sweep(): string {
  const points: string[] = [];
  for (let x = LEFT; x <= RIGHT; x += 3) {
    points.push(`${x},${(BASELINE + 16 * Math.sin((x - LEFT) / 9)).toFixed(1)}`);
  }
  return `M${points.join("L")}`;
}

const WAVES: Record<LabStatus, string> = {
  idle: noise(0.8, 1),
  analyzing: sweep(),
  likely_benign: rhythm([
    { width: 72, height: 34 },
    { width: 72, height: 34 },
    { width: 72, height: 34 },
  ]),
  uncertain: rhythm([
    { width: 54, height: 30 },
    { width: 82, height: 22 },
    { width: 40, height: 38 },
    { width: 42, height: 18 },
  ]),
  likely_pathogenic: rhythm([
    { width: 26, height: 52 },
    { width: 30, height: 34 },
    { width: 22, height: 58 },
    { width: 34, height: 26 },
    { width: 24, height: 56 },
    { width: 28, height: 40 },
    { width: 24, height: 60 },
    { width: 30, height: 30 },
  ]),
  out_of_scope: noise(1.2, 4),
  error: noise(3, 7),
};

const SCANLINES = Array.from({ length: 38 }, (_, index) => 18 + index * 4);

function truncate(text: string, length: number): string {
  return text.length > length ? `${text.slice(0, length - 1)}…` : text;
}

export function VitalsMonitor({
  status,
  label,
  probability,
}: {
  status: LabStatus;
  label: string | null;
  probability: number | null;
}) {
  const { locale, messages } = useSettings();
  const readout = probability === null ? "--.-%" : formatPercent(probability, locale);
  return (
    <svg viewBox="0 0 250 190" className="monitor-svg" data-status={status}>
      <defs>
        <linearGradient id="monitor-case" x1="0" y1="0" x2="0" y2="1">
          <stop offset="0" stopColor="var(--metal-light)" />
          <stop offset="0.12" stopColor="var(--metal-dark)" />
          <stop offset="1" stopColor="var(--metal-dark)" />
        </linearGradient>
        <radialGradient id="monitor-phosphor" cx="0.5" cy="0.45" r="0.75">
          <stop offset="0" stopColor="var(--crt-tint)" stopOpacity="0.35" />
          <stop offset="0.7" stopColor="var(--crt-screen)" />
          <stop offset="1" stopColor="#000000" />
        </radialGradient>
        <filter id="monitor-glow" x="-10%" y="-40%" width="120%" height="180%">
          <feGaussianBlur stdDeviation="2.4" result="blur" />
          <feMerge>
            <feMergeNode in="blur" />
            <feMergeNode in="SourceGraphic" />
          </feMerge>
        </filter>
        <clipPath id="monitor-screen-clip">
          <rect x="12" y="12" width="226" height="150" rx="16" />
        </clipPath>
      </defs>
      <rect x="1" y="1" width="248" height="188" rx="14" fill="url(#monitor-case)" className="monitor-case" />
      <rect x="12" y="12" width="226" height="150" rx="16" fill="url(#monitor-phosphor)" />
      <g clipPath="url(#monitor-screen-clip)">
        <path d="M16 52H234M16 92H234M16 132H234M56 12V162M96 12V162M136 12V162M176 12V162M216 12V162" className="monitor-grid" />
        <path d={WAVES[status]} className="monitor-wave" pathLength={1} filter="url(#monitor-glow)" />
        {status === "analyzing" && <rect x="12" y="12" width="18" height="150" className="monitor-scanbar" />}
        <g className="monitor-scanlines">
          {SCANLINES.map((y) => (
            <line key={y} x1="12" x2="238" y1={y} y2={y} />
          ))}
        </g>
      </g>
      <text x="22" y="32" className="monitor-text">
        {truncate(label ?? "SPECIMEN --", 24)}
      </text>
      <circle cx="226" cy="28" r="4.5" className="monitor-led" />
      <text x="22" y="152" className="monitor-text monitor-status">
        {messages.monitor[status]}
      </text>
      <text x="228" y="152" textAnchor="end" className="monitor-text monitor-readout">
        {readout}
      </text>
      <circle cx="26" cy="176" r="3" className="monitor-screw" />
      <circle cx="224" cy="176" r="3" className="monitor-screw" />
      <rect x="96" y="171" width="58" height="9" rx="2" className="monitor-plate" />
    </svg>
  );
}
