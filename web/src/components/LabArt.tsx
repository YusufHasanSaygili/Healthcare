type Point = [number, number];

function seeded(index: number): number {
  const value = Math.sin(index * 12.9898) * 43758.5453;
  return value - Math.floor(value);
}

function round(value: number): number {
  return Math.round(value * 10) / 10;
}

function cubicPoint(points: Point[], t: number): Point {
  const [p0, p1, p2, p3] = points as [Point, Point, Point, Point];
  const u = 1 - t;
  return [
    u ** 3 * p0[0] + 3 * u * u * t * p1[0] + 3 * u * t * t * p2[0] + t ** 3 * p3[0],
    u ** 3 * p0[1] + 3 * u * u * t * p1[1] + 3 * u * t * t * p2[1] + t ** 3 * p3[1],
  ];
}

function cubicAngle(points: Point[], t: number): number {
  const [a, b] = [cubicPoint(points, Math.max(0, t - 0.01)), cubicPoint(points, Math.min(1, t + 0.01))];
  return (Math.atan2(b[1] - a[1], b[0] - a[0]) * 180) / Math.PI;
}

function curve(points: Point[]): string {
  const [p0, p1, p2, p3] = points as [Point, Point, Point, Point];
  return `M${p0.join(",")} C${p1.join(",")} ${p2.join(",")} ${p3.join(",")}`;
}

type Cable = {
  key: string;
  d: string;
  width: number;
  ties: { x: number; y: number; angle: number }[];
  plug: Point | null;
  duration: string;
  delay: string;
};

function hanging(start: number, span: number, sag: number): Point[] {
  return [
    [start, -12],
    [start + span * 0.15, sag],
    [start + span * 0.85, sag],
    [start + span, -12],
  ];
}

const FAR_CABLES: Cable[] = Array.from({ length: 9 }, (_, i) => {
  const start = round(20 + i * 165 + seeded(i + 3) * 70);
  const points = hanging(start, round(110 + seeded(i + 13) * 170), round(70 + seeded(i + 23) * 120));
  return {
    key: `far-${i}`,
    d: curve(points),
    width: round(2 + seeded(i + 33) * 2),
    ties: [],
    plug: null,
    duration: `${16 + seeded(i + 43) * 8}s`,
    delay: `${-seeded(i + 53) * 10}s`,
  };
});

const NEAR_CABLES: Cable[] = [
  { start: 90, span: 300, sag: 210, width: 11 },
  { start: 520, span: 240, sag: 150, width: 8 },
  { start: 980, span: 330, sag: 250, width: 12 },
  { start: 1230, span: 200, sag: 130, width: 7 },
].map((item, i) => {
  const points = hanging(item.start, item.span, item.sag);
  return {
    key: `near-${i}`,
    d: curve(points),
    width: item.width,
    ties: [],
    plug: null,
    duration: `${18 + i * 3}s`,
    delay: `${-i * 4}s`,
  };
});

const BUNDLES: Cable[] = [
  { start: 300, span: 420, sag: 280 },
  { start: 780, span: 360, sag: 210 },
].flatMap((item, b) => {
  const lanes = [-7, 0, 7];
  return lanes.map((offset, lane) => {
    const points = hanging(item.start + offset, item.span, item.sag + offset);
    return {
      key: `bundle-${b}-${lane}`,
      d: curve(points),
      width: 6,
      ties:
        lane === 1
          ? [0.25, 0.5, 0.75].map((t) => {
              const [x, y] = cubicPoint(points, t);
              return { x: round(x), y: round(y), angle: round(cubicAngle(points, t)) };
            })
          : [],
      plug: null,
      duration: `${20 + b * 4}s`,
      delay: `${-b * 6}s`,
    };
  });
});

const PENDANTS: Cable[] = [
  { x: 690, length: 230 },
  { x: 1120, length: 180 },
].map((item, i) => ({
  key: `pendant-${i}`,
  d: `M${item.x},-12 C${item.x + 8},${item.length * 0.4} ${item.x - 10},${item.length * 0.7} ${item.x},${item.length}`,
  width: 6,
  ties: [],
  plug: [item.x, item.length],
  duration: `${9 + i * 3}s`,
  delay: `${-i * 2}s`,
}));

function CableShape({ cable, near }: { cable: Cable; near: boolean }) {
  return (
    <g
      className={near ? "cable cable-near" : "cable cable-far"}
      style={{ animationDuration: cable.duration, animationDelay: cable.delay }}
    >
      <path d={cable.d} className="cable-body" strokeWidth={cable.width} />
      {near && (
        <path
          d={cable.d}
          className="cable-rim"
          strokeWidth={Math.max(1, cable.width * 0.28)}
          transform={`translate(0 ${-cable.width * 0.28})`}
        />
      )}
      {cable.ties.map((tie) => (
        <rect
          key={`${tie.x}-${tie.y}`}
          x={tie.x - 3}
          y={tie.y - 14}
          width="6"
          height="28"
          rx="1.5"
          className="cable-tie"
          transform={`rotate(${tie.angle} ${tie.x} ${tie.y})`}
        />
      ))}
      {cable.plug && (
        <g className="cable-plug">
          <rect x={cable.plug[0] - 8} y={cable.plug[1] - 2} width="16" height="24" rx="3" />
          <path d={`M${cable.plug[0] - 4} ${cable.plug[1] + 22}v8M${cable.plug[0] + 4} ${cable.plug[1] + 22}v8`} />
        </g>
      )}
    </g>
  );
}

export function Cables() {
  return (
    <svg className="bio-cables" viewBox="0 0 1440 340" preserveAspectRatio="xMidYMin slice">
      <g className="cables-far">
        {FAR_CABLES.map((cable) => (
          <CableShape key={cable.key} cable={cable} near={false} />
        ))}
      </g>
      {[...NEAR_CABLES, ...BUNDLES, ...PENDANTS].map((cable) => (
        <CableShape key={cable.key} cable={cable} near />
      ))}
    </svg>
  );
}


function normalAt(spine: Point[], i: number): Point {
  const last = spine.length - 1;
  const [ax, ay] = spine[Math.max(0, i - 1)] as Point;
  const [bx, by] = spine[Math.min(last, i + 1)] as Point;
  const length = Math.hypot(bx - ax, by - ay) || 1;
  return [-(by - ay) / length, (bx - ax) / length];
}

type Profile = (t: number) => number;

function ribbon(spine: Point[], width: Profile): string {
  const last = spine.length - 1;
  const left: string[] = [];
  const right: string[] = [];
  spine.forEach(([x, y], i) => {
    const [nx, ny] = normalAt(spine, i);
    const half = width(i / last) / 2;
    left.push(`${round(x + nx * half)},${round(y + ny * half)}`);
    right.push(`${round(x - nx * half)},${round(y - ny * half)}`);
  });
  return `M${left.join("L")}L${right.reverse().join("L")}Z`;
}

function shine(spine: Point[], width: Profile, from: number, to: number): string {
  const last = spine.length - 1;
  const indices = spine.map((_, i) => i).filter((i) => i >= from * last && i <= to * last);
  const [mx, my] = normalAt(spine, indices[Math.floor(indices.length / 2)] ?? 0);
  const side = mx * -0.6 + my * -0.8 >= 0 ? 1 : -1;
  const points = indices.map((i): Point => {
    const [x, y] = spine[i] as Point;
    const [nx, ny] = normalAt(spine, i);
    const offset = side * width(i / last) * 0.26;
    return [x + nx * offset, y + ny * offset];
  });
  const widths = indices.map((i) => width(i / last));
  return ribbon(points, (t) => (widths[Math.round(t * (widths.length - 1))] ?? 0) * 0.24 * Math.sin(Math.PI * t));
}

function crescent(cx: number, cy: number, rx: number, ry: number, start: number, span: number, thickness: number) {
  const arc = Array.from({ length: 17 }, (_, k): Point => {
    const angle = ((start + (span * k) / 16) * Math.PI) / 180;
    return [cx + rx * Math.cos(angle), cy + ry * Math.sin(angle)];
  });
  return ribbon(arc, (t) => thickness * Math.sin(Math.PI * t));
}

function blob(cx: number, cy: number, rx: number, ry: number, seed: number): string {
  const points = Array.from({ length: 48 }, (_, k) => {
    const angle = (k / 48) * 2 * Math.PI;
    const r =
      1 + 0.1 * Math.sin(3 * angle + seed) + 0.06 * Math.sin(5 * angle + 2 * seed) + 0.035 * Math.sin(9 * angle + 3 * seed);
    return `${round(cx + rx * r * Math.cos(angle))},${round(cy + ry * r * Math.sin(angle))}`;
  });
  return `M${points.join("L")}Z`;
}

function drip([x, y]: Point, length: number): string {
  const bottom = round(y + length);
  return [
    `M${round(x - 0.8)},${round(y)}`,
    `C${round(x - 0.7)},${round(y + length * 0.5)} ${round(x - 2.6)},${bottom} ${round(x)},${bottom}`,
    `C${round(x + 2.6)},${bottom} ${round(x + 0.7)},${round(y + length * 0.5)} ${round(x + 0.8)},${round(y)}Z`,
  ].join(" ");
}

type Limb = {
  x: number;
  y: number;
  angle: number;
  length: number;
  width: number;
  curl: number;
  wobble: number;
  waves: number;
  phase: number;
};

function grow(limb: Limb, steps: number): Point[] {
  const points: Point[] = [[limb.x, limb.y]];
  let [x, y] = [limb.x, limb.y];
  for (let i = 1; i <= steps; i += 1) {
    const t = i / steps;
    const heading = limb.angle + limb.curl * t * t + limb.wobble * Math.sin(t * limb.waves + limb.phase);
    x += (Math.cos((heading * Math.PI) / 180) * limb.length) / steps;
    y += (Math.sin((heading * Math.PI) / 180) * limb.length) / steps;
    points.push([x, y]);
  }
  return points;
}

function taper(limb: Limb): Profile {
  return (t) => limb.width * (1 - t) ** 1.25 * (1 + 0.2 * (1 - t) * Math.sin(t * 11 + limb.phase)) + 0.5;
}

function limbShapes(limb: Limb, spine: Point[]) {
  const width = taper(limb);
  return { body: ribbon(spine, width), shine: shine(spine, width, 0.06, 0.6) };
}

function branch(
  parent: Limb,
  spine: Point[],
  at: number,
  turn: number,
  share: number,
  thickness: number,
  curl: number,
): Limb {
  const i = Math.round(at * (spine.length - 1));
  const [x, y] = spine[i] as Point;
  const [ax, ay] = spine[i - 1] as Point;
  const heading = (Math.atan2(y - ay, x - ax) * 180) / Math.PI;
  return {
    x,
    y,
    angle: heading + turn,
    length: parent.length * share,
    width: taper(parent)(at) * thickness,
    curl,
    wobble: 14,
    waves: 12,
    phase: parent.phase + turn,
  };
}

function sample(points: Point[], steps: number): Point[] {
  return Array.from({ length: steps + 1 }, (_, i) => cubicPoint(points, i / steps));
}

type GooGroup = {
  key: string;
  origin: Point;
  motion: "goo-writhe" | "goo-hair";
  duration: number;
  delay: number;
  bodies: string[];
  shines: string[];
  rim: boolean;
};

const GOO_CORE: Point = [175, 450];

const HEROES: Limb[] = [
  { x: 110, y: 425, angle: -118, length: 250, width: 30, curl: 75, wobble: 12, waves: 8, phase: 0.5 },
  { x: 172, y: 410, angle: -90, length: 330, width: 34, curl: -85, wobble: 10, waves: 7, phase: 2.1 },
  { x: 238, y: 415, angle: -62, length: 280, width: 28, curl: -105, wobble: 14, waves: 9, phase: 4 },
  { x: 288, y: 445, angle: -28, length: 210, width: 22, curl: -125, wobble: 10, waves: 8, phase: 1.3 },
  { x: 62, y: 458, angle: -150, length: 160, width: 20, curl: 95, wobble: 12, waves: 8, phase: 3.3 },
  { x: 205, y: 432, angle: -74, length: 190, width: 15, curl: 150, wobble: 16, waves: 10, phase: 5.2 },
];

const HERO_SPINES = HEROES.map((limb) => grow(limb, 44));

const HERO_GROUPS: GooGroup[] = HEROES.map((limb, i) => {
  const spine = HERO_SPINES[i] as Point[];
  const main = limbShapes(limb, spine);
  const side = i % 2 === 0 ? 1 : -1;
  const fork = branch(limb, spine, 0.5, side * (34 + seeded(i * 10) * 20), 0.32 + seeded(i * 10 + 1) * 0.14, 0.6, side * 120);
  const branches = [limbShapes(fork, grow(fork, 28))];
  return {
    key: `hero-${i}`,
    origin: [limb.x, limb.y],
    motion: "goo-writhe",
    duration: 7 + i * 1.3,
    delay: -i * 1.7,
    bodies: [main.body, ...branches.map((branch) => branch.body)],
    shines: [main.shine, ...branches.map((branch) => branch.shine)],
    rim: true,
  };
});

const MEDIUM: Limb[] = Array.from({ length: 8 }, (_, i) => {
  const s = i / 7;
  return {
    x: round(70 + s * 230 + (seeded(i + 600) - 0.5) * 20),
    y: round(440 - Math.sin(s * Math.PI) * 30),
    angle: -150 + s * 120 + (seeded(i + 610) - 0.5) * 30,
    length: 80 + seeded(i + 620) * 70,
    width: 9 + seeded(i + 630) * 6,
    curl: (seeded(i + 640) - 0.5) * 200,
    wobble: 16,
    waves: 11,
    phase: seeded(i + 650) * 6,
  };
});

const HAIRS: Limb[] = Array.from({ length: 44 }, (_, i) => ({
  x: round(110 + seeded(i + 700) * 190),
  y: round(425 + seeded(i + 710) * 20),
  angle: -78 + (seeded(i + 720) - 0.5) * 36,
  length: 140 + seeded(i + 730) * 180,
  width: 1 + seeded(i + 740) * 1.6,
  curl: (seeded(i + 750) - 0.5) * 60,
  wobble: 14 + seeded(i + 760) * 10,
  waves: 10 + seeded(i + 770) * 6,
  phase: seeded(i + 780) * 6,
}));

function bundle(limbs: Limb[], count: number, steps: number) {
  return Array.from({ length: count }, (_, g) =>
    limbs.filter((_, i) => i % count === g).map((limb) => limbShapes(limb, grow(limb, steps))),
  );
}

const HAIR_GROUPS: GooGroup[] = bundle(HAIRS, 4, 48).map((shapes, g) => ({
  key: `hair-${g}`,
  origin: [190, 440],
  motion: "goo-hair",
  duration: 6 + g * 1.4,
  delay: -g * 2.2,
  bodies: shapes.map((shape) => shape.body),
  shines: [],
  rim: false,
}));

const MEDIUM_GROUPS: GooGroup[] = bundle(MEDIUM, 3, 30).map((shapes, g) => ({
  key: `medium-${g}`,
  origin: GOO_CORE,
  motion: "goo-writhe",
  duration: 9 + g * 2,
  delay: -g * 3,
  bodies: shapes.map((shape) => shape.body),
  shines: shapes.map((shape) => shape.shine),
  rim: false,
}));

const GOO_GROUPS = [...HAIR_GROUPS, ...MEDIUM_GROUPS, ...HERO_GROUPS];

const MEMBRANES = [
  [4, 0],
  [0, 1],
  [1, 5],
  [5, 2],
  [2, 3],
].map(([a, b]) => {
  const first = HERO_SPINES[a as number] as Point[];
  const second = HERO_SPINES[b as number] as Point[];
  const [ax, ay] = first[0] as Point;
  const [bx, by] = second[0] as Point;
  const [pax, pay] = first[8] as Point;
  const [pbx, pby] = second[8] as Point;
  const control = `${round((pax + pbx) * 0.3 + (ax + bx) * 0.2)},${round((pay + pby) * 0.3 + (ay + by) * 0.2)}`;
  return `M${ax},${ay}L${round(pax)},${round(pay)}Q${control} ${round(pbx)},${round(pby)}L${bx},${by}Z`;
});

const STRINGS = [
  [4, 0],
  [0, 1],
  [1, 5],
  [5, 2],
  [2, 3],
  [0, 1],
  [1, 5],
  [2, 3],
].map(([a, b], i) => {
  const at = (seed: number) => Math.round((0.22 + seeded(seed) * 0.28) * 44);
  const [ax, ay] = (HERO_SPINES[a as number] as Point[])[at(i + 950)] as Point;
  const [bx, by] = (HERO_SPINES[b as number] as Point[])[at(i + 960)] as Point;
  const sag = 16 + seeded(i + 970) * 18;
  const points: Point[] = [
    [ax, ay],
    [ax + (bx - ax) / 3, ay + (by - ay) / 3 + sag],
    [ax + ((bx - ax) * 2) / 3, ay + ((by - ay) * 2) / 3 + sag],
    [bx, by],
  ];
  return [
    ribbon(sample(points, 20), (t) => 0.7 + 4.4 * (2 * t - 1) ** 4),
    drip(cubicPoint(points, 0.5), 10 + seeded(i + 975) * 8),
  ];
}).flat();

function coreTop(x: number): number {
  return 470 - 68 * Math.sqrt(Math.max(0, 1 - ((x - 175) / 130) ** 2));
}

const GOO_CORE_SHAPES = [
  blob(175, 470, 130, 68, 1.1),
  blob(92, 448, 52, 38, 2.3),
  blob(258, 452, 64, 44, 3.7),
  blob(178, 430, 52, 32, 4.4),
  ...Array.from({ length: 7 }, (_, i) => {
    const x = 70 + i * 36 + seeded(i + 980) * 12;
    return blob(x, coreTop(x) + 4, 14 + seeded(i + 990) * 12, 11 + seeded(i + 995) * 8, i);
  }),
];

const CORE_SPECULARS = [
  crescent(175, 470, 112, 56, 200, 50, 4),
  crescent(175, 470, 100, 50, 262, 30, 2.5),
  crescent(92, 448, 40, 28, 200, 60, 3),
  crescent(258, 452, 50, 34, 235, 55, 3),
  crescent(178, 430, 40, 24, 215, 70, 2.4),
];

const DROPLETS = HERO_SPINES.map((spine, i) => {
  const [x, y] = spine[spine.length - 1] as Point;
  const [ax, ay] = spine[spine.length - 2] as Point;
  const distance = Math.hypot(x - ax, y - ay) || 1;
  const gap = 9 + seeded(i + 900) * 10;
  return {
    cx: round(x + ((x - ax) / distance) * gap),
    cy: round(y + ((y - ay) / distance) * gap),
    r: round(1.4 + seeded(i + 910) * 1.8),
  };
});

function GooLayer({ rimOnly }: { rimOnly: boolean }) {
  return GOO_GROUPS.filter((group) => !rimOnly || group.rim).map((group) => (
    <g
      key={group.key}
      className={group.motion}
      style={{
        transformOrigin: `${group.origin[0]}px ${group.origin[1]}px`,
        animationDuration: `${group.duration}s`,
        animationDelay: `${group.delay}s`,
      }}
    >
      {group.bodies.map((d) => (
        <path key={d} d={d} />
      ))}
      {!rimOnly &&
        group.shines.map((d) => <path key={d} d={d} className="goo-shine" />)}
    </g>
  ));
}

export function Symbiote() {
  return (
    <svg viewBox="0 0 420 520" className="symbiote-svg">
      <defs>
        <radialGradient id="goo-backlight">
          <stop offset="0" stopColor="var(--goo-backlight)" />
          <stop offset="1" stopColor="var(--goo-backlight)" stopOpacity="0" />
        </radialGradient>
      </defs>
      <ellipse cx="190" cy="330" rx="230" ry="210" fill="url(#goo-backlight)" />
      <g className="goo-rim" transform="translate(-1.3 -1.7)">
        {STRINGS.map((d) => (
          <path key={d} d={d} />
        ))}
        <GooLayer rimOnly />
        {[...MEMBRANES, ...GOO_CORE_SHAPES].map((d) => (
          <path key={d} d={d} />
        ))}
      </g>
      <g className="goo-body">
        {STRINGS.map((d) => (
          <path key={d} d={d} />
        ))}
        <GooLayer rimOnly={false} />
        {[...MEMBRANES, ...GOO_CORE_SHAPES].map((d) => (
          <path key={d} d={d} />
        ))}
        {CORE_SPECULARS.map((d) => (
          <path key={d} d={d} className="goo-shine" />
        ))}
        {DROPLETS.map((drop) => (
          <circle key={`${drop.cx}-${drop.cy}`} cx={drop.cx} cy={drop.cy} r={drop.r} />
        ))}
      </g>
    </svg>
  );
}

const KNOT: Point = [35, 214];
const ANCHORED: Profile = (t) => 1.6 + 5 * (1 - t) ** 2.2 + 3.5 * t ** 6;
const SPANNING: Profile = (t) => 1.2 + 3.6 * (2 * t - 1) ** 4;
const STEM: Profile = (t) => 5 + 4 * (1 - t) ** 2 + 7 * t ** 5;

const WEB_STRANDS: { points: Point[]; profile: Profile }[] = [
  { points: [KNOT, [26, 196], [12, 168], [3, 150]], profile: ANCHORED },
  { points: [KNOT, [46, 190], [58, 150], [67, 128]], profile: ANCHORED },
  { points: [KNOT, [48, 222], [58, 240], [67, 252]], profile: ANCHORED },
  { points: [KNOT, [28, 236], [12, 262], [3, 284]], profile: ANCHORED },
  { points: [KNOT, [44, 250], [60, 300], [67, 334]], profile: ANCHORED },
  { points: [KNOT, [32, 170], [24, 110], [22, 72]], profile: ANCHORED },
  { points: [KNOT, [40, 160], [50, 110], [52, 72]], profile: ANCHORED },
  { points: [KNOT, [54, 262], [14, 320], [35, 372]], profile: STEM },
  { points: [[3, 104], [20, 136], [46, 178], [67, 190]], profile: SPANNING },
  { points: [[3, 206], [22, 208], [46, 196], [67, 168]], profile: SPANNING },
  { points: [[3, 318], [24, 324], [46, 312], [67, 292]], profile: SPANNING },
  { points: [[4, 244], [16, 290], [40, 340], [64, 360]], profile: SPANNING },
];

const WEB_MEMBRANES = [
  [0, 5],
  [1, 6],
  [2, 4],
  [3, 7],
];

const WEB_NODES = [
  { strand: 8, at: 0.5, size: 4 },
  { strand: 11, at: 0.35, size: 3.5 },
  { strand: 10, at: 0.55, size: 3 },
  { strand: 1, at: 0.6, size: 2.6 },
  { strand: 3, at: 0.55, size: 2.8 },
];

const WEB_DRIPS = [
  { strand: 8, at: 0.45, length: 9 },
  { strand: 9, at: 0.6, length: 7 },
  { strand: 4, at: 0.55, length: 10 },
  { strand: 10, at: 0.4, length: 8 },
];

function cup([, y]: Point): { rx: number; ry: number } {
  return y <= 75 ? { rx: 4.5, ry: 1.4 } : { rx: 1.4, ry: 4.5 };
}

function GooWeb({ x }: { x: number }) {
  const strands = WEB_STRANDS.map(({ points, profile }) => {
    const shifted = points.map(([px, py]): Point => [px + x, py]);
    return { points: shifted, spine: sample(shifted, 24), profile };
  });
  const [kx, ky] = [KNOT[0] + x, KNOT[1]];
  const ends = strands.flatMap((strand) =>
    strand.profile === SPANNING ? [strand.points[0], strand.points[3]] : [strand.points[3]],
  ) as Point[];
  return (
    <g className="goo-web" style={{ transformOrigin: `${kx}px ${ky}px` }}>
      <ellipse cx={kx} cy="232" rx="36" ry="150" fill="url(#web-backlight)" />
      {WEB_MEMBRANES.map(([a, b]) => {
        const [pax, pay] = cubicPoint(strands[a as number]?.points as Point[], 0.5);
        const [pbx, pby] = cubicPoint(strands[b as number]?.points as Point[], 0.5);
        const control = `${round((pax + pbx) * 0.3 + kx * 0.4)},${round((pay + pby) * 0.3 + ky * 0.4)}`;
        return (
          <path
            key={`${a}-${b}`}
            d={`M${kx},${ky}L${round(pax)},${round(pay)}Q${control} ${round(pbx)},${round(pby)}Z`}
            className="web-goo web-membrane"
          />
        );
      })}
      {strands.map((strand) => (
        <path key={strand.points.join()} d={ribbon(strand.spine, strand.profile)} className="web-goo" />
      ))}
      {ends.map((end) => (
        <ellipse key={end.join()} cx={end[0]} cy={end[1]} {...cup(end)} className="web-goo" />
      ))}
      {WEB_NODES.map((node) => {
        const [nx, ny] = cubicPoint(strands[node.strand]?.points as Point[], node.at);
        return (
          <path
            key={`${node.strand}-${node.at}`}
            d={blob(nx, ny, node.size * 1.3, node.size, node.strand)}
            className="web-goo"
          />
        );
      })}
      <g className="web-knot" style={{ transformOrigin: `${kx}px ${ky}px` }}>
        <path d={blob(kx, ky, 10, 8, 2)} className="web-goo" />
      </g>
      <path d={blob(x + 35, 376, 18, 6, 5)} className="web-goo" />
      {WEB_DRIPS.map((item, i) => (
        <path
          key={`${item.strand}-${item.at}`}
          d={drip(cubicPoint(strands[item.strand]?.points as Point[], item.at), item.length)}
          className="web-goo web-drip"
          style={{ animationDuration: `${3.2 + i * 0.7}s`, animationDelay: `${-i * 1.3}s` }}
        />
      ))}
      {strands.map((strand) => (
        <path key={strand.points.join()} d={shine(strand.spine, strand.profile, 0.1, 0.55)} className="web-shine" />
      ))}
    </g>
  );
}

const FLESH_LOBES = [
  { cx: 35, cy: 216, rx: 21, ry: 28 },
  { cx: 22, cy: 196, rx: 12, ry: 13 },
  { cx: 49, cy: 232, rx: 13, ry: 15 },
  { cx: 27, cy: 244, rx: 12, ry: 10 },
  { cx: 47, cy: 192, rx: 10, ry: 11 },
];

const FLESH_VEINS: Limb[] = Array.from({ length: 12 }, (_, i) => ({
  x: 35 + (seeded(i + 1300) - 0.5) * 36,
  y: 216 + (seeded(i + 1310) - 0.5) * 56,
  angle: (i / 12) * 360 + seeded(i + 1320) * 30,
  length: 22 + seeded(i + 1330) * 18,
  width: 1 + seeded(i + 1340) * 0.6,
  curl: (seeded(i + 1350) - 0.5) * 140,
  wobble: 28,
  waves: 9,
  phase: seeded(i + 1360) * 6,
}));

const FLESH_TUBES: Point[][] = [
  [
    [24, 72],
    [18, 118],
    [16, 160],
    [23, 186],
  ],
  [
    [46, 72],
    [53, 126],
    [55, 160],
    [47, 184],
  ],
];

const FLESH_TENDRILS: Limb[] = [
  { x: 24, y: 248, angle: 100, length: 92, width: 5.5, curl: 40, wobble: 16, waves: 9, phase: 0.4 },
  { x: 34, y: 252, angle: 88, length: 110, width: 6, curl: -30, wobble: 14, waves: 8, phase: 2.2 },
  { x: 45, y: 246, angle: 78, length: 84, width: 5, curl: -50, wobble: 18, waves: 10, phase: 4.1 },
  { x: 54, y: 234, angle: 60, length: 64, width: 4, curl: 70, wobble: 16, waves: 9, phase: 1.5 },
];

const FLESH_PARTICLES = Array.from({ length: 22 }, (_, i) => ({
  key: i,
  cx: round(6 + seeded(i + 1400) * 58),
  cy: round(96 + seeded(i + 1410) * 270),
  r: round(0.5 + seeded(i + 1420) * 1.1),
}));

function FleshCore({ x }: { x: number }) {
  const shift = (limb: Limb): Limb => ({ ...limb, x: limb.x + x });
  const veins = FLESH_VEINS.map((limb) => {
    const vein = shift(limb);
    return ribbon(grow(vein, 18), taper(vein));
  });
  const capillaries = FLESH_VEINS.map((limb, i) => {
    const vein = shift(limb);
    const twig = branch(vein, grow(vein, 18), 0.45, (i % 2 === 0 ? 1 : -1) * 40, 0.6, 0.5, 60);
    return ribbon(grow(twig, 12), taper(twig));
  });
  return (
    <g className="flesh">
      <ellipse cx={x + 35} cy="216" rx="34" ry="80" fill="url(#flesh-backlight)" />
      <ellipse cx={x + 30} cy="300" rx="22" ry="14" className="flesh-cloud" />
      <ellipse cx={x + 44} cy="150" rx="16" ry="10" className="flesh-cloud" />
      {FLESH_PARTICLES.map((particle) => (
        <circle key={particle.key} cx={x + particle.cx} cy={particle.cy} r={particle.r} className="flesh-particle" />
      ))}
      {FLESH_TENDRILS.map((limb) => {
        const tendril = shift(limb);
        const spine = grow(tendril, 30);
        const width = taper(tendril);
        return (
          <g
            key={limb.phase}
            className="goo-writhe"
            style={{
              transformOrigin: `${tendril.x}px ${tendril.y}px`,
              animationDuration: `${6 + limb.phase}s`,
              animationDelay: `${-limb.phase}s`,
            }}
          >
            <path d={ribbon(spine, width)} className="flesh-tendril" />
            <path d={shine(spine, width, 0.05, 0.6)} className="flesh-shine" />
          </g>
        );
      })}
      {FLESH_TUBES.map((points) => {
        const spine = sample(
          points.map(([px, py]): Point => [px + x, py]),
          24,
        );
        return (
          <g key={points.join()}>
            <path d={ribbon(spine, () => 4)} className="flesh-tube" />
            <path d={shine(spine, () => 4, 0.02, 0.98)} className="flesh-tube-shine" />
            {[4, 8, 12, 16, 20].map((i) => {
              const [rx, ry] = spine[i] as Point;
              return <ellipse key={i} cx={round(rx)} cy={round(ry)} rx="2.8" ry="1" className="flesh-tube-rib" />;
            })}
          </g>
        );
      })}
      <g className="flesh-pulse" style={{ transformOrigin: `${x + 35}px 216px` }}>
        {FLESH_LOBES.map((lobe, i) => (
          <path key={lobe.cy} d={blob(x + lobe.cx, lobe.cy, lobe.rx, lobe.ry, i + 1)} fill="url(#flesh-body)" />
        ))}
        <clipPath id="flesh-clip">
          {FLESH_LOBES.map((lobe, i) => (
            <path key={lobe.cy} d={blob(x + lobe.cx, lobe.cy, lobe.rx, lobe.ry, i + 1)} />
          ))}
        </clipPath>
        <g clipPath="url(#flesh-clip)">
          {veins.map((d) => (
            <path key={d} d={d} className="flesh-vein" />
          ))}
          {capillaries.map((d) => (
            <path key={d} d={d} className="flesh-capillary" />
          ))}
        </g>
        {FLESH_LOBES.map((lobe) => (
          <path
            key={lobe.cy}
            d={crescent(x + lobe.cx, lobe.cy, lobe.rx * 0.72, lobe.ry * 0.72, 200, 70, Math.max(1.2, lobe.rx * 0.12))}
            className="flesh-gloss"
          />
        ))}
      </g>
    </g>
  );
}

const SERUM_TOP = 118;

const VIRIONS = [
  { cx: 34, cy: 214, r: 11, far: false },
  { cx: 18, cy: 150, r: 5.5, far: false },
  { cx: 52, cy: 170, r: 7, far: false },
  { cx: 22, cy: 270, r: 6.5, far: false },
  { cx: 50, cy: 298, r: 8, far: false },
  { cx: 16, cy: 320, r: 3.5, far: false },
  { cx: 46, cy: 246, r: 4, far: true },
  { cx: 27, cy: 190, r: 3, far: true },
  { cx: 56, cy: 134, r: 3, far: true },
  { cx: 38, cy: 326, r: 3.5, far: true },
];

const SERUM_PARTICLES = Array.from({ length: 26 }, (_, i) => ({
  key: i,
  cx: round(6 + seeded(i + 1000) * 58),
  cy: round(SERUM_TOP + 10 + seeded(i + 1010) * 250),
  r: round(0.5 + seeded(i + 1020) * 0.9),
}));

function Virion({ cx, cy, r }: { cx: number; cy: number; r: number }) {
  const spikes = Math.max(10, Math.round(r * 1.7));
  return (
    <>
      <g className="virion-spin" style={{ animationDuration: `${10 + r * 1.5}s` }}>
        {Array.from({ length: spikes }, (_, k) => {
          const angle = (k / spikes) * 2 * Math.PI;
          const [dx, dy] = [Math.cos(angle), Math.sin(angle)];
          return (
            <g key={k}>
              <line
                x1={round(cx + dx * r * 0.9)}
                y1={round(cy + dy * r * 0.9)}
                x2={round(cx + dx * r * 1.42)}
                y2={round(cy + dy * r * 1.42)}
                className="virion-spike"
                strokeWidth={round(Math.max(0.5, r * 0.11))}
              />
              <circle
                cx={round(cx + dx * r * 1.5)}
                cy={round(cy + dy * r * 1.5)}
                r={round(Math.max(0.7, r * 0.14))}
                className="virion-knob"
              />
            </g>
          );
        })}
      </g>
      <circle cx={cx} cy={cy} r={r} fill="url(#virion-body)" />
      {r > 5 &&
        [0, 1, 2, 3, 4].map((k) => (
          <circle
            key={k}
            cx={round(cx + Math.cos(k * 1.3 + r) * r * 0.5)}
            cy={round(cy + Math.sin(k * 1.3 + r) * r * 0.5)}
            r={round(r * 0.12)}
            className="virion-pit"
          />
        ))}
    </>
  );
}

function Serum({ x }: { x: number }) {
  return (
    <g className="serum">
      <rect x={x} y={SERUM_TOP} width="70" height="280" fill="url(#serum-fluid)" />
      <ellipse cx={x + 35} cy={SERUM_TOP} rx="33" ry="3.5" className="serum-meniscus" />
      <g className="virion-far">
        {VIRIONS.filter((virion) => virion.far).map((virion) => (
          <Virion key={virion.cy} cx={x + virion.cx} cy={virion.cy} r={virion.r} />
        ))}
      </g>
      {SERUM_PARTICLES.map((particle) => (
        <circle
          key={particle.key}
          cx={x + particle.cx}
          cy={particle.cy}
          r={particle.r}
          className="serum-particle"
          style={{
            animationDuration: `${1.6 + seeded(particle.key + 1030) * 2.4}s`,
            animationDelay: `${-seeded(particle.key + 1040) * 4}s`,
          }}
        />
      ))}
      {VIRIONS.filter((virion) => !virion.far).map((virion, i) => (
        <g
          key={virion.cy}
          className="virion-drift"
          style={{ animationDuration: `${5 + seeded(i + 1100) * 4}s`, animationDelay: `${-seeded(i + 1110) * 6}s` }}
        >
          <Virion cx={x + virion.cx} cy={virion.cy} r={virion.r} />
        </g>
      ))}
      <ellipse cx={x + 35} cy="381" rx="30" ry="4" className="serum-sediment" />
    </g>
  );
}

type Ring = { y: number; height: number; inset: number; finish: "metal" | "dark" | "rubber" | "light" };

const TOP_RINGS: Ring[] = [
  { y: 18, height: 6, inset: 20, finish: "dark" },
  { y: 24, height: 10, inset: 3, finish: "metal" },
  { y: 34, height: 24, inset: -8, finish: "metal" },
  { y: 58, height: 8, inset: -4, finish: "dark" },
  { y: 66, height: 6, inset: -1, finish: "rubber" },
];

const BASE_RINGS: Ring[] = [
  { y: 386, height: 6, inset: -1, finish: "rubber" },
  { y: 392, height: 6, inset: -4, finish: "light" },
  { y: 398, height: 24, inset: -8, finish: "metal" },
  { y: 422, height: 14, inset: -13, finish: "dark" },
];

const SAG = 2.5;

function arc(x1: number, x2: number, y: number): string {
  return `Q${round((x1 + x2) / 2)},${y + SAG} ${x2},${y}`;
}

function ringPath(x: number, ring: Ring): string {
  const [x1, x2] = [x + ring.inset, x + 70 - ring.inset];
  const bottom = ring.y + ring.height;
  return `M${x1},${ring.y} ${arc(x1, x2, ring.y)} L${x2},${bottom} Q${round((x1 + x2) / 2)},${bottom + SAG} ${x1},${bottom}Z`;
}

const RING_FILL = {
  metal: "url(#pod-metal)",
  dark: "url(#pod-metal-dark)",
  rubber: "var(--pod-rubber)",
  light: "var(--pod)",
};

function Housing({ x, y, height, led }: { x: number; y: number; height: number; led: boolean }) {
  const center = x + 35;
  return (
    <g>
      <path d={`M${x - 8},${y + 5} ${arc(x - 8, x + 78, y + 5)}`} className="pod-groove" />
      <path d={`M${x - 8},${y + height - 5} ${arc(x - 8, x + 78, y + height - 5)}`} className="pod-groove" />
      <rect x={center - 6} y={y + 2} width="12" height={height - 2} rx="1.5" className="pod-clamp" />
      <rect x={x - 6} y={y + 2} width="4" height={height - 2} rx="1" className="pod-clamp" />
      <rect x={x + 72} y={y + 2} width="4" height={height - 2} rx="1" className="pod-clamp" />
      {led && <rect x={x + 50} y={y + height / 2 - 1.5} width="9" height="3" rx="1" className="pod-led" />}
      {led &&
        [0, 1, 2].map((k) => (
          <rect key={k} x={x + 12 + k * 4} y={y + height / 2 - 3} width="1.6" height="6" className="pod-vent" />
        ))}
    </g>
  );
}

function Bolts({ x, y }: { x: number; y: number }) {
  return [6, 20, 35, 50, 64].map((dx) => (
    <circle key={dx} cx={x + dx} cy={y + (dx === 35 ? SAG : dx === 20 || dx === 50 ? SAG * 0.7 : SAG * 0.3)} r="1.3" className="pod-bolt" />
  ));
}

function Hose({ x }: { x: number }) {
  const center = x + 35;
  return (
    <g>
      <rect x={center - 7} y="-4" width="14" height="24" fill="url(#pod-metal-dark)" />
      {Array.from({ length: 7 }, (_, k) => (
        <g key={k}>
          <line x1={center - 7} x2={center + 7} y1={-2 + k * 3} y2={-2 + k * 3} className="pod-hose-ridge" />
          <line x1={center - 7} x2={center + 7} y1={-1 + k * 3} y2={-1 + k * 3} className="pod-hose-light" />
        </g>
      ))}
    </g>
  );
}

function PodFrame({ x }: { x: number }) {
  return (
    <g className="pod-frame">
      <Hose x={x} />
      {[...TOP_RINGS, ...BASE_RINGS].map((ring) => (
        <path
          key={ring.y}
          d={ringPath(x, ring)}
          fill={RING_FILL[ring.finish]}
          className={ring.finish === "light" ? "pod-ring pod-ring-light" : "pod-ring"}
        />
      ))}
      <Housing x={x} y={34} height={24} led />
      <Housing x={x} y={398} height={24} led={false} />
      <Bolts x={x} y={61} />
      <Bolts x={x} y={428} />
    </g>
  );
}

const PODS = [
  { x: 40, tone: "red", Virus: FleshCore },
  { x: 130, tone: "purple", Virus: GooWeb },
  { x: 220, tone: "green", Virus: Serum },
];

const POD_BUBBLES = PODS.flatMap(({ x }, pod) =>
  [0, 1, 2].map((i) => ({
    key: `${pod}-${i}`,
    cx: x + 14 + seeded(pod * 3 + i + 60) * 42,
    r: 1.5 + seeded(pod * 3 + i + 70) * 2,
    duration: `${6 + seeded(pod * 3 + i + 80) * 5}s`,
    delay: `${-seeded(pod * 3 + i + 90) * 8}s`,
  })),
);

export function SpecimenPods() {
  return (
    <svg viewBox="0 -6 300 460" className="pods-svg">
      <defs>
        {PODS.map(({ tone }) => (
          <linearGradient key={tone} id={`pod-fluid-${tone}`} x1="0" y1="0" x2="0" y2="1">
            <stop offset="0" stopColor={`var(--pod-${tone})`} stopOpacity="0.1" />
            <stop offset="1" stopColor={`var(--pod-${tone})`} stopOpacity="0.34" />
          </linearGradient>
        ))}
        {PODS.map(({ tone }) => (
          <radialGradient key={tone} id={`pod-underlight-${tone}`} cy="1" r="1">
            <stop offset="0" stopColor={`var(--pod-${tone})`} stopOpacity="0.55" />
            <stop offset="1" stopColor={`var(--pod-${tone})`} stopOpacity="0" />
          </radialGradient>
        ))}
        <linearGradient id="pod-metal" x1="0" x2="1">
          <stop offset="0" stopColor="var(--metal-dark)" />
          <stop offset="0.18" stopColor="var(--metal-light)" />
          <stop offset="0.26" stopColor="var(--metal-shine)" />
          <stop offset="0.36" stopColor="var(--metal-mid)" />
          <stop offset="0.84" stopColor="var(--metal-dark)" />
          <stop offset="0.94" stopColor="var(--metal-light)" />
          <stop offset="1" stopColor="var(--metal-dark)" />
        </linearGradient>
        <linearGradient id="pod-metal-dark" x1="0" x2="1">
          <stop offset="0" stopColor="var(--pod-rubber)" />
          <stop offset="0.24" stopColor="var(--metal-mid)" />
          <stop offset="0.5" stopColor="var(--metal-dark)" />
          <stop offset="0.92" stopColor="var(--metal-mid)" />
          <stop offset="1" stopColor="var(--pod-rubber)" />
        </linearGradient>
        <linearGradient id="pod-volume" x1="0" x2="1">
          <stop offset="0" stopColor="#000000" stopOpacity="0.5" />
          <stop offset="0.28" stopColor="#000000" stopOpacity="0" />
          <stop offset="0.72" stopColor="#000000" stopOpacity="0" />
          <stop offset="1" stopColor="#000000" stopOpacity="0.55" />
        </linearGradient>
        <radialGradient id="web-backlight">
          <stop offset="0" stopColor="var(--pod-purple)" stopOpacity="0.6" />
          <stop offset="1" stopColor="var(--pod-purple)" stopOpacity="0" />
        </radialGradient>
        <radialGradient id="flesh-backlight">
          <stop offset="0" stopColor="#ff4a2a" stopOpacity="0.3" />
          <stop offset="1" stopColor="#ff4a2a" stopOpacity="0" />
        </radialGradient>
        <radialGradient id="flesh-body" cx="0.4" cy="0.35" r="0.8">
          <stop offset="0" stopColor="#ff9c80" />
          <stop offset="0.3" stopColor="#c8321f" />
          <stop offset="0.7" stopColor="#6e0b0c" />
          <stop offset="1" stopColor="#2a0304" />
        </radialGradient>
        <linearGradient id="serum-fluid" x1="0" y1="0" x2="0" y2="1">
          <stop offset="0" stopColor="var(--pod-green)" stopOpacity="0.62" />
          <stop offset="0.45" stopColor="var(--pod-green)" stopOpacity="0.3" />
          <stop offset="1" stopColor="var(--pod-green)" stopOpacity="0.55" />
        </linearGradient>
        <radialGradient id="virion-body" cx="0.36" cy="0.32" r="0.75">
          <stop offset="0" stopColor="#f1fff4" />
          <stop offset="0.35" stopColor="var(--pod-green)" />
          <stop offset="1" stopColor="var(--virion-shadow)" />
        </radialGradient>
        <filter id="pod-glow" x="-60%" y="-30%" width="220%" height="160%">
          <feGaussianBlur stdDeviation="14" />
        </filter>
        <filter id="virus-glow" x="-30%" y="-30%" width="160%" height="160%">
          <feGaussianBlur stdDeviation="1.6" result="blur" />
          <feMerge>
            <feMergeNode in="blur" />
            <feMergeNode in="SourceGraphic" />
          </feMerge>
        </filter>
        <clipPath id="pod-inside">
          {PODS.map(({ x }) => (
            <rect key={x} x={x + 2} y="70" width="66" height="318" rx="2" />
          ))}
        </clipPath>
      </defs>
      {PODS.map(({ x, tone }) => (
        <g key={x} className={`pod pod-${tone}`}>
          <rect x={x} y="70" width="70" height="320" rx="4" className="pod-halo" filter="url(#pod-glow)" />
          <rect x={x} y="68" width="70" height="322" rx="3" fill={`url(#pod-fluid-${tone})`} className="pod-glass" />
        </g>
      ))}
      <g clipPath="url(#pod-inside)">
        {PODS.map(({ x, tone }) => (
          <ellipse key={x} cx={x + 35} cy="388" rx="40" ry="90" fill={`url(#pod-underlight-${tone})`} />
        ))}
        {PODS.map(({ x, tone, Virus }) => (
          <g key={x} className={`pod-${tone}`} filter="url(#virus-glow)">
            <Virus x={x} />
          </g>
        ))}
        {POD_BUBBLES.map((bubble) => (
          <circle
            key={bubble.key}
            cx={bubble.cx}
            cy="380"
            r={bubble.r}
            className="pod-bubble"
            style={{ animationDuration: bubble.duration, animationDelay: bubble.delay }}
          />
        ))}
      </g>
      {PODS.map(({ x, tone }) => (
        <g key={x} className={`pod-${tone}`}>
          <rect x={x} y="68" width="70" height="322" fill="url(#pod-volume)" className="pod-volume" />
          <path d={`M${x + 1.5},70V388M${x + 68.5},70V388`} className="pod-edge" />
          <rect x={x + 12} y="76" width="5" height="306" rx="2.5" className="pod-shine" />
          <rect x={x + 56} y="76" width="2" height="306" rx="1" className="pod-shine pod-shine-soft" />
          <PodFrame x={x} />
        </g>
      ))}
    </svg>
  );
}
