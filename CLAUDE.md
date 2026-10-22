# CLAUDE.md

## BEHAVIOR

- Be concise. Sacrifice grammar for concision.
- Do not be agreeable by default. If something seems wrong or unclear, stop and say so.
- If requirements are missing or ambiguous, ask before proceeding — not after.
- Flag bad ideas. A short "this seems off because X — confirm?" is better than silent compliance.
- Never assume. Verify.
- If an approach fails twice, stop and reassess — don't retry the same thing.

## WORKFLOW

- Read `docs/README.md` for a summary of available docs, then read only what's relevant to the task
- If anything is unclear, ask before proceeding
- Before a multi-file change, surface a plan — list the files with a one-line summary each — and ask: "Do you want this plan in chat output, saved as a file (e.g. `plans/<slug>.md`), or skip and proceed directly?" Then act accordingly.
- **Batch independent reads/searches in parallel** — calls with no dependency between them should go out together in one response. Only serialize when one call depends on another's output.
- **But edit files one at a time, sequentially — never batched.** Issue the next edit only after the previous one resolves. If I reject an early edit in a batch, every edit after it is built on rejected state and has to be thrown out too — wasted work.

## VALIDATION

- After making changes, run the relevant validation (lint, type checks, build, tests, code generation). Read failures, fix them, and iterate until green, then report results.
- Pause and ask first before running anything slow, expensive, destructive, or outward-facing (deploys, migrations, mass deletes).

## PLANS

- See `docs/plans.md` for format and session behavior.

## MEMORY

- Do not rely on internal/implicit memory between sessions. Treat each conversation as a fresh start.
- All persistent context lives in `docs/` — check `docs/README.md` for what exists.
- Never assume you "remember" something — if it's not in the repo, you don't know it.

## RULES

- Doc filenames must be self-describing (e.g., `auth.md`, `api-patterns.md`, `conventions.md`)
- `docs/` = all project docs — both short context notes and full feature docs, readable by Claude and engineers
- `docs/` holds prescriptive conventions and patterns (rules to follow), not descriptive state (claims about current code that go stale). The *why*, not a snapshot of the *what*. Never put source/business logic in `docs/`
- DO NOT re-read CLAUDE.md mid-session — read once, then follow the flow
- Do NOT read `.claude/prompts/` unless explicitly asked to (e.g., "fix this prompt" or "add a prompt file")