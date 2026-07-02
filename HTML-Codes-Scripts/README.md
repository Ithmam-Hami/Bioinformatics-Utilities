# HTML Codes and Scripts

Small browser-based utilities for bioinformatics workflows.

## Mito SVG Review Dashboard

`mito-svg-review-dashboard.html` is a single-file browser dashboard for reviewing mitochondrial genome coverage map SVGs.

It lets users load local SVG files, review them one by one, assign formal coverage categories, add notes, and export the review table as a TSV file.

### Categories

- `ACCEPTABLE_COVERAGE`
- `REVIEW_REQUIRED`
- `EXCLUDE_FROM_ANALYSIS`

### Shortcuts

- `A`: Acceptable Coverage
- `R`: Review Required
- `E`: Exclude From Analysis
- Left/right arrows: move between samples

The dashboard runs locally in the browser and does not upload data anywhere.
