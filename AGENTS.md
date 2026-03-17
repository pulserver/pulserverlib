# AGENTS

This repository is the active `pulserverlib` codebase.

## Context continuity
- Refer to `AGENT.md` for domain architecture and implementation guidance.
- Refer to `.github/prompts/` for planning and implementation prompt context.

## Note
- `migration/` contains historical split scaffolds and planning artifacts.
- Keep new changes aligned with the current top-level repository layout.

## Current workflow
- For multipass prep/cooldown safety semantics, first implement and validate behavior in MATLAB truth tooling (`tests/generators/TruthBuilder.m` and `tests/generators/+testutils/`) before touching C library safety checks.
- For off-isocenter frequency-modulation updates, first extend MATLAB truth artifacts and plotting (`tests/generators/TruthBuilder.m`, `tests/generators/+testutils/`), then update C tests to compile against the new truth schema, and only afterwards tighten runtime assertions and modify C library behavior.
