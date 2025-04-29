# Supplementary Material for PhD Thesis

This repository contains supplementary material for the PhD thesis:

**"Inference for meaningful estimands in factorial survival designs and competing risks settings"** by **Merle Munko**

## Contents

The repository includes the following folders and files:

1. **Tables.pdf** - A PDF file with detailed simulation results of the simulation studies.
2. **RMST_Simu** - R code for the simulation studies in Section 4.
3. **RMST_data_example** - R code for the data example in Section 4.
4. **RMTL_Simu** - R code for the simulation study in Section 5.
5. **RMTL_data_example** - R code for the data example in Section 5.

## Running the Simulations

- The entire simulation can be started by running the files `reproduceSimulationresults.R` in **2. (RMST_Simu)** and **4. (RMTL_Simu)**. 
- **Note:** Running the full simulations is not recommended as they are time-consuming. Instead, specific scenarios can be run using `Startasimulation.R`.
- Simulations were performed on the Linux Server **'fourier'** of Otto-von-Guericke University Magdeburg.
- Results are saved in the `results` folders, with one file per simulation scenario.
- Use `load_results.R` and `load_results_multiple.R` in the results folders to load and summarize results, as well as create plots.

## Running the Data Examples

- Run `data_example.R` in **3. (RMST_data_example)** and **5. (RMTL_data_example)**.
- Other files are dependencies and do not need to be executed separately.

## Data Availability Statement

- The original dataset **"gabriel.csv"**, loaded in line 73 of `RMST_data_example/data_example.R` and line 7 of `RMST_Simu/data_data` is **not uploaded** due to ethical restrictions and informed consent limitations.
- However, we uploaded a similar, artificially generated dataset for illustrative reasons.
- Data access for the original dataset can be requested from **Jon Genuneit** for collaborative research efforts.

## Licensing

This repository contains content under different licenses:

- The **code** is licensed under **GPL-3.0-or-later**.
- The file **Tables.pdf** is licensed under **CC-BY 4.0**.

You may copy, distribute, and modify the content as long as you comply with the terms of the respective licenses:

- The **code** must be used and redistributed under **GPL-3.0-or-later**.
- The **tables** must be used and redistributed under **CC-BY 4.0**.

### Full License Descriptions:
- **GPL-3.0-or-later**: [GNU GPL-3.0 License](https://www.gnu.org/licenses/gpl-3.0.html)
- **CC-BY 4.0**: [Creative Commons BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Credits

- The file **RMST_data_example/functions_factorial.R** was written by **Marc Ditzhaus**.
- Some implemented functions are based on those by **Marc Ditzhaus**.
- **Tables.pdf** includes tables from the supplement of **Munko et al. (2024, Statistics in Medicine)**.
- The R code contains parts of the source code from the CRAN packages **GFDrmst** and **GFDrmtl**.
