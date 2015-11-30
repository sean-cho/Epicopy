# Epicopy

Epicopy is an R package that obtains copy number variation (CNV) information from Illumina Human Methylation 450K microarrays. If the users have their own reference normal samples or decide to not use any normals for reference (see vignette after installation), Epicopy can be used as a standalone package.

If not, the users may choose to install its companion package **"EpicopyData"** which contains a series of normal tissues compiled from the Cancer Genome Atlas (TCGA) effort. It contains normal tissues from three sources; thyroid, breast, and lung. For more information on EpicopyData please visit the [EpicopyData page](https://github.com/sean-cho/EpicopyData).

---

## Installation

The preferred method for installation is through Hadley Wickham's `install_github` function through the package [devtools](https://github.com/hadley/devtools).

If users have their own normal samples, it is recommended that they install the standalone version as it is much smaller in size.

### Standalone

To install the standalone package use the following lines:
```
devtools::install_github('sean-cho/Epicopy')
```

### Complete

**Warning:** The EpicopyData file contains a large data file and requires the manual installation of EpicopyData.

The complete version includes all the normals compiled from TCGA.

#### Step 1: Install EpicopyData

Please download the binary release [here](https://github.com/sean-cho/EpicopyData/releases/download/v1.0.1/EpicopyData_1.0.1.tar.gz) and run the following line of code:

```
path <- "path.to.file/EpicopyData_1.0.1.tar.gz"
install.packages(path, repos = NULL, type = 'source')
```

#### Step 2a: Install Epicopy without vignette
```
devtools::install_github('sean-cho/Epicopy')
```

#### Step 2b: Install Epicopy with vignette
**Warning:** Please note that the vignette includes the complete processing of three example raw data files and may take a while to build.
```
devtools::install_github('sean-cho/Epicopy', build_vignette = TRUE)
```

---

## Usage

### Standalone

Users can view the vignette for examples [here](https://github.com/sean-cho/Epicopy/blob/master/vignettes/Epicopy.Rmd). Do not build it as EpicopyNormals is a dependency.

The key to using Epicopy as a standalone package would be to specify normals in the samplesheet or use the argument `Normals = NA` in the `epicopy` or `getLRR` functions. This is important, as the default argument uses all normals included in the EpicopyData package.

### With EpicopyNormals

Install and build the vignette. Users may specify any of the normals included in with the EpicopyData package.
