These files are needed to run the RevBayes scripts located in `../code/`. The file `../code/run_inf.Rev` controls how the input files are used. A brief description of the data files:

- `viburnum.tre` contains a time-calibrated MCC tree for 119 (of ~165) Viburnum species.
- `viburnum.nex` contains the coded biome-region states for 173 Viburnum taxa. Note, the analysis itself uses only the 119 taxa present in the provided phylogeny (`viburnum.tre`).
- `atlas/` is a subdirectory that contains the following files/folders that describe the biome structure model
    - `atlas/epoch_times.txt` gives the time intervals for the modeled epochs
    - `atlas/epoch_names.txt` gives the indices, time intervals, and text names for the modeled epochs
    - `atlas/null/` contains the null/uniform adjacency matrices (A_G(m)) where all regions have strong (2) features
    - `atlas/land/` contains the geography adjacency matrices (A_G(m)) for marginal (0), weak (1), and strong (2) features
    - `atlas/tropical/` contains the tropical adjacency matrices (A_G(m)) for marginal (0), weak (1), and strong (2) features
    - `atlas/warm/` contains the warm temperate adjacency matrices (A_G(m)) for marginal (0), weak (1), and strong (2) features
    - `atlas/cold/` contains the cold temperate adjacency matrices (A_G(m)) for marginal (0), weak (1), and strong (2) features
