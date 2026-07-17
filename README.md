# MSc-Project
## Repository structure

```text
MSc-Project/
│
├── SimBu.R
│   └── Prepares the single-cell datasets and generates the primary
│       breast and lung pseudo-bulk samples.
│
├── InstaPrism.R
│   └── Builds InstaPrism references, performs deconvolution, and
│       evaluates the resulting estimates.
│
├── DWLS.R
│   └── Builds or loads DWLS signature matrices, performs deconvolution,
│       and evaluates the estimates.
│
├── CIBERSORTx.R
│   └── Prepares CIBERSORTx inputs, processes deconvolution outputs,
│       and evaluates the resulting estimates.
│
├── Compare_all_methods.R
│   └── Combines and compares the performance results from InstaPrism,
│       DWLS, and CIBERSORTx.
│
├── Prepare_Controlled_Simulations.R
│   └── Generates the controlled CAF-abundance scenarios and associated
│       pseudo-bulk simulations.
│
├── Controlled_Simulation_Functions.R
│   └── Contains reusable functions for scenario generation, simulation,
│       deconvolution, alignment, evaluation, plotting, and result storage.
│
├── Controlled_InstaPrism.R
│   └── Runs the complete controlled sensitivity workflow with InstaPrism.
│
├── Controlled_DWLS.R
│   └── Runs the complete controlled sensitivity workflow with DWLS.
│
├── cibersortx/
│   └── CIBERSORTx input files, container-related files, and outputs.
│
├── instaprism_results/
│   └── InstaPrism reference objects, estimates, metrics, and figures.
│
├── dwls_results/
│   └── DWLS signature matrices, estimates, metrics, and figures.
│
├── combined_results/
│   └── Combined benchmarking metrics and method-comparison results.
│
└── figures/
    └── controlled_analysis/
        └── Selected figures from the controlled sensitivity experiments.
