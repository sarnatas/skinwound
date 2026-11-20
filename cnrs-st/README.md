# CNRS professional training: spatial transcriptomics (Visium)

**CNRS Formation Entreprise**, 24–25 November 2026, ENS de Lyon / Spatial-Cell-ID. Instructors: Sergio and Marco.
**Audience:** bioinformaticians (engineers, post-docs, PhD students).
**Data:** Visium (GSE241124), 16 sections; scRNA-seq reference from GSE241132 for the sc–st integration on Day 2.

Online: the course pages are part of the site of this repository (menu *CNRS – Spatial transcriptomics*).

## Contents

| Path | What it contains |
|---|---|
| `index.qmd` | Course page: access to IFB Biosphere, permanent storage, Python environment, schedule of the two days |
| `practical.qmd` | The practical work: what the participants do, step by step, and the questions to discuss |
| `CNRS_ST_participant_guide.pdf` | The two pages above in one PDF (sent by e-mail to the participants) |
| `environment.yml` | Conda environment `cnrs-st` (versions tested together; also on the course storage) |
| `scripts/` | Preparation of `DATA_CNRS` from the GEO downloads |
| `slides/` | Lecture slides |

Kept locally, not in the repository (see `.gitignore`): the instructor notebook (`notebooks/instructor/`), the PowerPoint sources, the planning spreadsheet, `resources/`.

## Computing resources

IFB Biosphere VMs (Jupyter appliance, group *Spatial-Cell-ID*). Data, environment file and one folder per participant are on the permanent storage
`/home/jovyan/ifbdata/spatial_cell_id/CNRS_Formation_Entreprise/` (`DATA_CNRS/`, `environment.yml`, `trainees/<first letter of first name + last name>/`).
