# Mitochondrial Coverage Map Review Dashboard

A single-file browser dashboard for reviewing mitochondrial genome coverage map SVGs.

The dashboard was built for rapid visual screening of IGV-generated mitochondrial coverage snapshots across many samples. It lets the user load a local folder of SVG files, review them one by one, assign formal calls, add notes, and export the results as a TSV file.

## Review Categories

- `ACCEPTABLE_COVERAGE`
- `REVIEW_REQUIRED`
- `EXCLUDE_FROM_ANALYSIS`

## Usage

1. Open `index.html` in a web browser.
2. Click **Load SVGs** and select one or more `.svg` files.
3. Review each mitochondrial coverage map.
4. Assign calls using the buttons or keyboard shortcuts:
   - `A`: Acceptable Coverage
   - `R`: Review Required
   - `E`: Exclude From Analysis
5. Add optional notes for each sample.
6. Click **Export TSV** to download the review table.

## Notes

- The dashboard runs locally in the browser.
- No data are uploaded anywhere.
- Calls are temporarily saved in browser `localStorage` for convenience.
- Exported TSV columns are `sample`, `file`, `call`, and `note`.
