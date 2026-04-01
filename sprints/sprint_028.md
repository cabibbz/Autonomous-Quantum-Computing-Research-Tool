# Sprint 028 — Repeated Syndrome Measurement with Majority Vote

**Date:** 2026-03-31
**Goal:** Test whether repeated syndrome extraction (3+ rounds) with majority vote provides enough syndrome reliability for active correction to beat passive encoding on [[5,1,3]].

**Literature check:**
- Searched: "repeated syndrome measurement majority vote quantum error correction [[5,1,3]] fault tolerant threshold 2025 2026"
- Key findings: Standard FT protocol requires d rounds of syndrome measurement for distance-d code. Google's below-threshold result (Nature 2024) uses repeated measurement on surface codes. Majority voting is the simplest decoder for repeated syndromes. No specific study on whether this is sufficient at [[5,1,3]] scale.
- Gap: Whether repeated measurement can overcome the single-round limitation proven in Sprints 026-027 for the [[5,1,3]] code specifically.

**Hypothesis:** 3 rounds of syndrome measurement with majority vote should dramatically improve syndrome reliability. The key question is whether the additional gate overhead (3× more 2Q gates) cancels the reliability gain.

## Experiments

### 28a: Repeated syndrome circuit — build and verify
**Status:** pending

### 28b: Repeated syndrome vs single-round vs passive under noise
**Status:** pending

### 28c: Number of rounds sweep — how many rounds needed?
**Status:** pending
