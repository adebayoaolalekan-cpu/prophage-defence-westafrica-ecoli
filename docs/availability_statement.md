# Declarations for BMC Genomics

Zenodo mints the DOI at publication, so the three placeholders marked
`[ZENODO DOI]` are the only thing to fill in. Everything else is final.

---

## Availability of data and materials

All data analysed in this study are publicly available and were obtained
entirely from open archives. No new sequence data were generated.

Genome assemblies and their sample metadata were retrieved from the European
Nucleotide Archive through its Portal API
(https://www.ebi.ac.uk/ena/portal/api/) on 25 November 2026, using result type
assembly under NCBI Taxonomy ID 562 filtered to sixteen West African countries.
Assembly statistics and the completeness reconciliation came from the NCBI
Datasets API (https://api.ncbi.nlm.nih.gov/datasets/v2alpha/) and NCBI Entrez
E-utilities (https://eutils.ncbi.nlm.nih.gov/entrez/eutils/). Every request was
an unauthenticated HTTP GET. No authentication, application or data use
agreement is required to repeat the retrieval, and the executed query strings
are reproduced in full in the study protocol together with the row counts each
returned.

Individual records are identified by INSDC assembly accession. The complete
accession list, with the ecological source, country, originating study and
quality band assigned to each, is deposited at Zenodo under
https://doi.org/[ZENODO DOI].

The study protocol, including the prespecified eligibility rule, the provenance
adjudication and the analysis decision rules, was registered before any genomic
analysis was carried out and is available at
https://doi.org/10.17605/OSF.IO/FXMEQ.

## Availability of computer code

Project name: Prophage and anti-phage defence genomic pilot workflow

Project home page: https://doi.org/[ZENODO DOI]

Archived version: https://doi.org/[ZENODO DOI]
(SHA-256 of the archive 4e4db7497479b3c3d246cbc7a3d3d685ece2150f9cf49ae39f28a9c8bbd13852,
also recorded in the registered protocol)

Operating system: Linux. Not tested on macOS or Windows.

Programming language: Python 3 and the Snakemake workflow language

Other requirements: Snakemake 7 or later, Conda or Mamba, and the tool
databases, which the workflow does not download. Database paths are given as
absolute paths in the configuration so that the release in use is pinned and
reportable rather than resolved at run time. The workflow calls geNomad, an
independent prophage predictor, CheckV, DefenseFinder, CRISPRCasFinder, Bakta,
CheckM2, seqkit and PopPUNK, each pinned through its own conda environment
specification.

License: MIT for the software. The selection tables and retrieval metadata are
dedicated to the public domain under CC0 1.0, since they derive from ENA and
NCBI records that are themselves in the public domain.

Any restrictions to use by non-academics: none.

## Competing interests

The author declares no competing interests.

## Funding

[State the funding or declare that no funding was received.]

## Authors' contributions

[Single-author study: state that the author designed the study, carried out the
analysis and wrote the manuscript.]

## Ethics approval and consent to participate

Not applicable. The study analyses genome assemblies and metadata already
deposited in public archives. No new human or animal material was collected, no
clinical intervention was performed, and no individually identifiable
information was accessed.

## Consent for publication

Not applicable.

---

# Shorter version, for a journal that wants one paragraph

All data are publicly available. Genome assemblies and metadata were retrieved
from the European Nucleotide Archive Portal API and reconciled against NCBI
Datasets and Entrez E-utilities on 25 November 2026. The full accession list,
the analysis workflow and the study protocol are deposited at Zenodo
(https://doi.org/[ZENODO DOI]) under the MIT License for software and CC0 for
data. The protocol was registered before any genomic analysis at
https://doi.org/10.17605/OSF.IO/FXMEQ.

---

# A note on two things reviewers may ask

The registration DOI and the Zenodo DOI point at the same study and should both
appear. The registration fixes the design and the decision rules. The Zenodo
deposit holds the artefacts those rules operate on, and its checksums are
recorded inside the registered protocol, so a reader can confirm that the
deposited genome set is the one the registration describes.

The deposit postdates the registration by hours rather than preceding it. That
ordering is deliberate and is explained in the protocol: the pilot's inference
block is a subset of the primary inference set, so running it before
registering would have exposed the primary contrast on roughly half the genomes
it is later tested on. Nothing in the deposit contains a result. At the time of
deposit no prophage had been predicted, no defence system called and no genome
sequence downloaded.
