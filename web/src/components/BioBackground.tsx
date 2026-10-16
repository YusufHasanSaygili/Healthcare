import { memo } from "react";
import { Cables, SpecimenPods, Symbiote } from "./LabArt";
import { type LabStatus, VitalsMonitor } from "./VitalsMonitor";

export type Lab = { status: LabStatus; label: string | null; probability: number | null };

function Background({ status, label, probability }: Lab) {
  return (
    <div className="bio-background" aria-hidden="true" data-status={status}>
      <div className="bio-furnace" />
      <Cables />
      <div className="bio-symbiote">
        <Symbiote />
      </div>
      <div className="bio-pods">
        <SpecimenPods />
      </div>
      <div className="bio-monitor">
        <VitalsMonitor status={status} label={label} probability={probability} />
      </div>
      <div className="bio-vignette" />
      <div className="bio-grain" />
      <div className="bio-tape bio-tape-top" />
      <div className="bio-tape bio-tape-bottom" />
    </div>
  );
}

export const BioBackground = memo(Background);
