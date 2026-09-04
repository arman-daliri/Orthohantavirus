# Orthohantavirus pipeline — working notes for Claude sessions

## RESUME HERE

**READ `docs/PROJECT-STATE.md` FIRST.** It is the single source of truth for resuming:
current branch, decisions of record, measured performance numbers, what was completed,
and the exact next action. This file only carries the guidance that outlives any one
session — repo policy, environment quirks, and a recurring defect class to watch for.

`docs/PROJECT-STATE.md` itself may lag the actual `git log` on a fast-moving branch;
when in doubt, trust `git log --oneline -20` and the current branch's own commits over
a stale "current state" paragraph, and update `docs/PROJECT-STATE.md` once you've
re-established the real state.

## Repository policy

**This repo stays PRIVATE.** It is the single working repo: code, data, docs, notes and
session backups all live here together. Do not propose splitting them out.

**Never rewrite history. Never force-push.** `git filter-repo`, `filter-branch`,
`rebase` onto pushed commits, and `push --force` are all off the table. The remote is
not the primary user's account, so a rewrite would break other people's clones. If
something sensitive lands in history, the answer is "the repo is private," not a
rewrite.

**Public release is a SEPARATE, CLEAN repository**, built by copying an explicit file
list into a fresh repo. Never by flipping this repo's visibility. The list lives in
`PUBLIC_RELEASE.md` — read it before any release work, and keep it current as the
project grows.

**Saving session backups and working notes under `docs/` is fine** while the repo is
private. They are excluded from public release by `PUBLIC_RELEASE.md`, not by deletion.

**If unsure whether a file belongs in the public release, ask.** Do not guess.

## Environment

Full detail in `tools/VERSIONS.txt` and `requirements-pip.txt`; `environment.yml`
describes the Linux/conda machine. Check both before assuming which interpreter and
package set a session is running under — this has diverged between machines before
(see `docs/PROJECT-STATE.md` for the currently-reconciled state).

## A recurring defect class in this codebase

Multiple reviews have caught the same shape: **an invariant asserted on an intermediate
value while the value that actually ships quietly violates it.** Examples so far: a
status field reported "scored" for a constant placeholder; a sensitivity sweep reported
rank correlations that were arithmetic identities; a total rounded to a literal `0.0` in
the shipped row while a test asserted on the raw dict instead.

**Assert on the row that ships, not on the function that feeds it.** When adding a test
for a guarantee, check it against the CSV/JSON artifact, not the intermediate return.
Also add a guard-the-guard assertion where a fixture could drift into making the test
trivially true.

## Debugging habits worth stealing

From an early defect-fixing session (commit `b8f5f0c`, 2026-07-30) — see the
"Debugging habits" subsection of the README's Known Limitations section for the full
writeup:

- **Distributions catch what code review misses.** If a function branches on data, log
  how often each branch fires, not just that the branching logic reads correctly.
- **Measure a fix before committing to its direction.** The "obvious" fix has twice
  turned out wrong in this codebase. Let a measurement, not which side looks more like a
  bug, decide.
