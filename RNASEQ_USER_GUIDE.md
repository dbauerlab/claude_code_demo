---
marp: true
theme: default
paginate: true
title: "Hands-on Guide — Bulk RNA-seq with Claude Code"
style: |
  section { font-size: 21px; padding: 34px 46px; justify-content: flex-start; }
  h1 { font-size: 32px; }
  h2 { font-size: 26px; }
  h3 { font-size: 22px; }
  pre, code { font-size: 0.78em; }
  pre { line-height: 1.28; }
  table { font-size: 17px; }
  blockquote { font-size: 0.95em; }
  section.dense { font-size: 17px; padding: 22px 40px; }
  section.dense pre { font-size: 12.5px; line-height: 1.22; }
---

<!--
Marp deck. View in VS Code with the "Marp for VS Code" extension (marp-team.marp-vscode):
open this file and click the Marp preview icon (top-right). Slides are separated by `---`.
Export via Command Palette → "Marp: Export slide deck…".
-->

# Hands-on Guide — Bulk RNA-seq Differential Expression with Claude Code

Work through this at your own pace. By the end you'll have used Claude Code to produce a runnable,
reproducible **Quarto (`.qmd`)** document that performs a standard DESeq2 differential expression
workflow in **R** on the **airway** bulk RNA-seq dataset.

> **Important:** Claude does *not* run the analysis for you. It **writes a `.qmd` file** that *you*
> render locally afterwards to produce all the plots and tables. The file is yours to keep, re-run,
> and adapt.

---

## 0. Before you start

- [ ] **Claude desktop app** installed and signed in (Mac/Windows — see the setup deck)
- [ ] A terminal available for `git` commands (cloning + branching)
- [ ] Data present: `data/rsem.merged.gene_counts.tsv` and `data/experiment_table.csv`
- [ ] To render the `.qmd` at the end: a local R install with `DESeq2`, `tidyverse`, `pheatmap`,
      `EnhancedVolcano` (or equivalents) and `quarto`. *(Not needed just to generate the file.)*

---

### The dataset

The **airway** dataset: primary human airway smooth muscle cells from **4 donors (cell lines)**, each
sequenced **untreated** and **treated with dexamethasone** (a glucocorticoid). It's a **paired
design** — the same donor appears in both conditions.

| File | What it is |
|------|-----------|
| `data/rsem.merged.gene_counts.tsv` | Gene-level counts from **RSEM**. Columns: `gene_id`, `transcript_id(s)`, then 8 sample columns (`SRR1039508` … `SRR1039521`). ~58,700 genes. |
| `data/experiment_table.csv` | Sample sheet: `ID`, `cellLine` (N61311, N052611, N080611, N061011), `treatment` (Untreated / Dexamethasone) |

**One thing to watch:** RSEM produces **estimated, non-integer counts** (e.g. `284.50`). DESeq2
requires integers, so the analysis must round the matrix (or use `tximport`). Keep an eye on whether
this is handled — it's a common trip-up.

---

## Set up Git first: clone, branch, commit

We'll treat this like real project work: get your own copy of the repo, do everything on a branch, and
commit as you go so nothing is ever lost. Where possible, just **ask Claude to run the Git for you**.

1. **Clone the repo** to your machine (get the URL from the demonstrator):
   ```bash
   git clone https://github.com/dbauerlab/claude_code_demo
   cd claude_code_demo
   ```
2. **Create a branch to work on** — never work directly on `main`. Name it after yourself:
   ```bash
   git checkout -b demo/your-name
   ```
   Or, once you've opened the project in the **Claude desktop app**, just ask:
   ```
   Create a new git branch called demo/your-name and switch to it.
   ```

---

3. **Commit as you go.** Each time you finish a step (added `CLAUDE.md`, added the skill, generated the
   `.qmd`), save a checkpoint. You can ask Claude to do it:
   ```
   Commit the current changes with a sensible message.
   ```
   Claude will stage the changes, write a message, and commit. Doing this regularly means you always
   have a clear history and can roll back if something goes wrong.

> You do all your work on your own branch; `main` stays clean. Keeping your history on the branch is
> plenty for this exercise (opening a pull request later is optional).

---

## What you'll do

You've already set up Git (clone + branch, above) — commit as you go from here on. Then the five
steps of the analysis itself:

1. **Start the project** — open the repo in the Claude desktop app and let it orient itself.
2. **Add context** — a `CLAUDE.md` file that teaches Claude about this project.
3. **Add a skill** — reusable, packaged expertise (ClawBio or a custom RNA-seq skill).
4. **Plan first** — use **plan mode** so Claude proposes an approach *before* writing code.
5. **Execute** — approve the plan and let Claude write the `.qmd`.

---

## 1. Start the project

Open the **Claude desktop app**, then open the cloned `claude_code_demo` folder as the working
directory.

As your first message, get Claude to orient itself and generate a project memory file:

```
Have a look around this repository and tell me what's here. Then run /init to
create a CLAUDE.md that captures what this project is.
```

Or run the built-in command directly:

```
/init
```

Claude will list the two data files, inspect the counts matrix and sample sheet, identify this as
bulk RNA-seq, and draft a `CLAUDE.md`. Notice it reads the *actual data* — it sees the non-integer
counts and the sample-sheet columns rather than guessing from the folder name.

---

## 2. Add context with `CLAUDE.md`

`CLAUDE.md` is a file Claude automatically reads every time it starts in this folder — your project's
standing instructions (conventions, gotchas, the analysis design). `/init` gives you a starting
point; extend it with the project-specific detail below.

Ask Claude to update it:

```
Update CLAUDE.md so it captures the experimental design and analysis conventions
for this project. Include: the paired design (~ cellLine + treatment), that
Untreated is the reference level, that RSEM counts are non-integer and must be
rounded for DESeq2, that we work in R, and that deliverables are Quarto .qmd files.
```

You're aiming for something like this (you can also paste it in yourself):

---

<!-- _class: dense -->

```markdown
# Project: Airway Bulk RNA-seq — Dexamethasone Differential Expression

## What this is
Bulk RNA-seq of primary human airway smooth muscle cells. 4 donors (cell lines),
each with an untreated and a dexamethasone-treated sample (8 samples total).
Goal: identify genes differentially expressed on dexamethasone treatment.

## Data
- `data/rsem.merged.gene_counts.tsv` — RSEM gene-level counts. Columns: `gene_id`,
  `transcript_id(s)`, then 8 sample columns. Gene IDs are Ensembl (ENSG...).
- `data/experiment_table.csv` — sample sheet: `ID`, `cellLine`, `treatment`.

## Experimental design (IMPORTANT)
- Paired design: model as `~ cellLine + treatment` so donor-to-donor variation is
  controlled for.
- `treatment` reference level MUST be `Untreated` (so log2FC is treated-vs-untreated).
- Sample order in the counts columns must be matched to the sample sheet before building
  the DESeqDataSet.

## Data gotchas
- RSEM counts are ESTIMATES and NON-INTEGER. DESeq2 needs integer counts — round the
  matrix with `round()` before `DESeqDataSetFromMatrix()`.
- Drop the `transcript_id(s)` column; set `gene_id` as rownames.

## Conventions
- Language: R.
- Deliverables are Quarto documents (`.qmd`) that render end-to-end with no manual steps.
- Use DESeq2 for DE; save results tables as CSV and figures as PNG/PDF under `results/`.
- Keep everything reproducible: set a seed where relevant, print `sessionInfo()`.
```

---

Everything in here is context Claude would otherwise have to guess at or ask about every session.
Written once, it applies to every future session — and it's version-controlled alongside your code.

---

## 3. Add a skill

Where `CLAUDE.md` is per-project, **skills** are reusable, packaged expertise you can carry across
projects. A skill is a folder with a `SKILL.md` that Claude loads *only when it's relevant*.

### Option A — install the ClawBio skill library

[ClawBio](https://github.com/ClawBio/ClawBio) is "the first bioinformatics-native AI agent skill
library" — genomics, GWAS, scRNA-seq, and bulk RNA-seq skills. Install it from inside Claude Code:

```
/plugin marketplace add ClawBio/ClawBio
/plugin install clawbio
```

> ClawBio's RNA-seq DE skill is built around **pydeseq2** (Python). It's excellent, but this exercise
> uses **R / DESeq2**, so we'll add a small *custom* skill tailored to exactly that. Both approaches
> are valid — the point is that skills are composable and you can write your own.

---

### Option B — a custom R/DESeq2 skill (used in this exercise)

Create the folder:

```bash
mkdir -p .claude/skills/bulk-rnaseq-deseq2
```

Then create `.claude/skills/bulk-rnaseq-deseq2/SKILL.md`. A skill is a Markdown file with **YAML
frontmatter** (`name` + `description` — the `description` is how Claude decides when to use it) and a
body of instructions:

---

<!-- _class: dense -->

```markdown
---
name: bulk-rnaseq-deseq2
description: >
  Use when performing bulk RNA-seq differential expression in R with DESeq2 from a
  gene-level counts matrix and a sample sheet. Covers reading counts (including
  non-integer RSEM/tximport counts), building a DESeqDataSet, running the standard
  DESeq2 workflow, and producing QC, results tables, and figures. Deliverables are
  reproducible Quarto (.qmd) documents.
---

# Bulk RNA-seq Differential Expression (R / DESeq2)

## When to use
A gene-level counts matrix + a sample sheet, and the user wants differential expression.

## Standard workflow
1. **Load inputs.** Read the counts matrix; drop non-count columns (e.g. `transcript_id(s)`);
   set gene IDs as rownames. Read the sample sheet.
2. **Handle non-integer counts.** RSEM/tximport counts are estimates — `round()` the matrix
   to integers before building the DESeqDataSet (or import via `tximport`).
3. **Align samples.** Reorder counts columns to match the sample-sheet row order; assert they
   are identical before proceeding.
4. **Factors & reference levels.** Convert design variables to factors. Explicitly set the
   control/reference level with `relevel()` (e.g. `Untreated`).
5. **Build & run.** `DESeqDataSetFromMatrix(countData, colData, design)` → `DESeq()`.
   Respect a paired/blocked design if donors/batches exist (e.g. `~ subject + condition`).
6. **QC.** Pre-filter low-count genes. `vst()`/`rlog` transform, then PCA and a
   sample-distance heatmap to check clustering by condition.
7. **Results.** `results()` with a stated alpha (e.g. 0.05); `lfcShrink()` for ranking/plots.
   Order by adjusted p-value; report how many genes pass FDR.
8. **Figures.** MA plot, volcano plot, and a heatmap of the top variable/DE genes.
9. **Export.** Save the full results table as CSV and figures under `results/`.
10. **Reproducibility.** End with `sessionInfo()`.

## Output contract
Produce a single Quarto (`.qmd`) document that renders end-to-end with no manual steps,
with clear section headings and a short interpretation under each figure.

## Common pitfalls
- Forgetting to round non-integer counts → DESeq2 errors.
- Reference level defaults to alphabetical → wrong sign on log2FC. Always `relevel()`.
- Counts columns not matched to the sample sheet → silently wrong results.
```

---

The `description` line does the important work — it's how Claude knows to reach for this skill when
you say "run DESeq2", without you having to name it. Skills keep your best practices encoded so every
analysis follows them.

### Option C — run with no skill at all (a useful baseline)

You don't *need* a skill — Claude already knows DESeq2. Skip Options A and B (or disable the skill for
a moment) and go straight to planning in §4, relying only on your `CLAUDE.md` context plus Claude's
general knowledge.

Try all three and compare — it's the most instructive part of the exercise:

---

| Setup | What tends to happen |
|-------|----------------------|
| **No skill, no `CLAUDE.md`** | Reasonable generic DESeq2, but more likely to miss project specifics — e.g. forget to round the RSEM counts, default the reference level alphabetically, or ignore the paired design. |
| **No skill, with `CLAUDE.md`** | Usually gets the project-specific details right because they're written down; workflow structure is up to Claude. |
| **With the skill (+ `CLAUDE.md`)** | Most consistent — the skill enforces the full workflow (QC, `lfcShrink`, exports, `sessionInfo()`) *and* the context supplies the project facts. |

The takeaway: the base model is capable, but context (`CLAUDE.md`) and packaged expertise (skills) are
what make it *reliable* and *repeatable*. Ask yourself which output you'd trust in a paper.

---

## 4. Plan before writing code

For anything beyond a one-liner, you don't want the agent to just start typing — you want it to
*think first* and show you the plan. That's **plan mode**: Claude investigates and proposes an
approach but can't edit files until you approve.

**Enter plan mode:** switch the mode to **plan mode** before sending — in the desktop app use the mode
selector next to the message box (**Shift+Tab** also cycles the modes: normal → auto-accept → plan
mode). Then send your request:

```
Plan a bulk RNA-seq differential expression analysis for this project using R and
DESeq2. The end product should be a Quarto .qmd file I can render locally that runs
the whole workflow and produces all outputs.

Read data/rsem.merged.gene_counts.tsv and data/experiment_table.csv first so the
plan reflects the real data. Follow the design and gotchas in CLAUDE.md. Show me the
plan before writing anything.
```

---

Claude will read both files, then lay out a step-by-step plan — load/clean counts, round non-integer
values, match samples, `~ cellLine + treatment` with `Untreated` as reference, filtering, VST + PCA
QC, `DESeq()`, `lfcShrink`, results table, MA/volcano/heatmap, exports, `sessionInfo()` — and
describe the `.qmd` it will create. **It stops and asks for your approval.**

Read the plan like a lab protocol from a new student: Does it round the RSEM counts? Is the reference
level right (`Untreated`)? Is it a paired model over `cellLine`? If anything's off, refine in plain
English (e.g. *"make sure Untreated is the reference level and it's a paired design over cellLine"*)
and let it re-plan.

---

## 5. Execute the plan

Once the plan looks right, let Claude build. Approve the plan (accept the prompt / choose "Yes,
proceed"). If you'd exited plan mode, this nudge works too:

```
The plan looks good. Go ahead and create the .qmd file. Don't run it — I'll render
it locally. Make sure it renders end-to-end with no manual edits.
```

Claude writes something like `airway_deseq2_analysis.qmd` with a YAML header, a setup chunk
(libraries), data loading and rounding, `DESeqDataSet` construction with the paired design, QC, DE,
figures, and exports — each chunk annotated. It may also create a `results/` folder.

---

### Render it locally

If your R + Quarto environment is ready:

```bash
quarto render airway_deseq2_analysis.qmd
```

Open the resulting HTML and look at the PCA (samples should separate by treatment), the volcano plot,
and the top DE genes. Classic dexamethasone-responsive genes such as **DUSP1, PER1, KLF15, ZBTB16,
CRISPLD2** should surface — a good sign the biology is real.

If your environment isn't set up, that's fine — the deliverable is a **reproducible file** you can
render whenever you like.

---

## Quick reference — the whole flow in prompts

1. `git clone https://github.com/dbauerlab/claude_code_demo` then `cd claude_code_demo`
2. `git checkout -b demo/your-name` *(or ask Claude to do it once launched)*
3. Open the `claude_code_demo` folder in the **Claude desktop app**
4. `/init` → *"Commit the current changes with a sensible message."*
5. *"Update CLAUDE.md with the paired design, Untreated as reference, RSEM non-integer rounding, R, and .qmd deliverables."* → commit
6. `/plugin marketplace add ClawBio/ClawBio` then `/plugin install clawbio` *(optional)* + create `.claude/skills/bulk-rnaseq-deseq2/SKILL.md` (block in §3) → commit
7. **Plan mode** (mode selector / Shift+Tab), then the planning prompt in §4
8. Approve → *"Create the .qmd, don't run it, must render end-to-end."* → commit
9. `quarto render airway_deseq2_analysis.qmd`

---

## FAQ

- **Why not just let Claude run the analysis?** Reproducibility and trust. A `.qmd` is auditable,
  version-controllable, and re-runnable by anyone.
- **How does it know DESeq2 needs integers?** Because you told it once, in `CLAUDE.md` and the skill.
  That's the whole idea: context turns a general model into a competent lab member.
- **`CLAUDE.md` vs a skill — when do I use which?** `CLAUDE.md` = facts about *this* project. A skill
  = reusable *how-to* you want across many projects.
- **Does my data leave the machine?** The analysis runs locally in R. Skills like ClawBio are
  local-first by design.

---

*Sources: [ClawBio on GitHub](https://github.com/ClawBio/ClawBio) · [ClawBio site](https://clawbio.ai/)*
