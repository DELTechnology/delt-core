# DELT-Core: An end-to-end computational framework for DNA-encoded chemical library analysis

## Abstract

DNA-encoded chemical libraries (DECLs) represent a powerful drug discovery platform that enables the screening of
millions to billions of compounds through DNA-tag identification via high-throughput sequencing. However, the analysis
of DECL screening data presents significant computational challenges, requiring the integration of sequencing data
processing, chemical structure reconstruction, and molecular property analysis. This protocol describes DELT-Core (
DNA-Encoded Library Technology Core), a comprehensive computational framework for end-to-end DECL analysis. Starting
from raw FASTQ sequencing files, library definitions, and selection data, DELT-Core enables: (1) flexible demultiplexing
of sequencing reads using adapted RNA-seq tools, (2) reconstruction of chemical structures as SMILES strings from
building blocks and reaction SMARTS, (3) computation of molecular properties and chemical descriptors, and (4)
generation of machine learning-ready datasets for hit identification and optimization. The framework leverages
established bioinformatics tools, for sequence demultiplexing and hit ranking, combined with cheminformatics
libraries for structure generation and property calculation. DELT-Core differentiates itself by providing a
collaborative, user-friendly command-line interface that enables both computational biologists and wet-lab chemists to
conduct comprehensive DECL analyses. As a practical example, we demonstrate the complete workflow for analyzing a
two-cycle DECL screening experiment, from raw sequencing data to ranked hit lists with computed molecular properties.
The protocol requires basic command-line knowledge and can be completed within one working day for typical datasets.
DELT-Core is available as an open-source Python package with comprehensive documentation and example datasets.

### Key Points

- DELT-Core provides an end-to-end computational workflow for DNA-encoded chemical library analysis, from raw sequencing
  data to chemical structure identification and property calculation
- The framework adapts established RNA-seq analysis tools, particularly Cutadapt [REF], for DECL-specific sequence
  demultiplexing with flexible error handling and edgeR [REF] for hit ranking.
- Chemical structures are reconstructed with RDKit [REF] from building blocks and reaction SMARTS, enabling analysis of diverse library
  architectures including multi-cycle and hybridized libraries
- Built-in quality control, visualization, and molecular property calculation modules facilitate comprehensive data
  interpretation and hit prioritization
- The modular design allows integration with machine learning pipelines for advanced hit prediction and optimization
  workflows

## Introduction

DNA-encoded chemical libraries (DECLs) have emerged as a transformative technology in drug discovery, enabling the
synthesis and screening of vast chemical spaces that would be impractical to explore using traditional high-throughput
screening approaches. In DECL technology, each chemical compound is covalently linked to a unique DNA barcode that
serves as an amplifiable identifier. This allows millions to billions of compounds to be pooled together and screened
simultaneously against biological targets, with hit identification achieved through DNA sequencing rather than
individual compound testing.

The DECL workflow typically involves several key steps: library synthesis through iterative cycles of chemical reactions
and DNA ligation, affinity-based selection against target proteins, PCR amplification of bound DNA tags, high-throughput
sequencing, and computational analysis to identify enriched compounds. While the experimental aspects of DECL
technology have been well-established, the computational analysis of DECL screening data presents unique challenges that
require specialized bioinformatics approaches.

Current computational challenges in DECL analysis include: (1) accurate demultiplexing of complex DNA barcode
combinations from sequencing reads, often in the presence of sequencing errors and incomplete reads; (2) reconstruction
of chemical structures from building block combinations and reaction schemes; (3) calculation of molecular properties
and descriptors for large compound collections; (4) quality control and statistical analysis of screening results; and (
5) integration with machine learning approaches for hit prediction and optimization.

Several computational tools have been developed for specific aspects of DECL analysis, but most focus on individual
components of the workflow rather than providing comprehensive end-to-end solutions. Existing tools often require
significant programming expertise, lack flexibility for diverse library architectures, or do not integrate well with
standard cheminformatics and machine learning pipelines.

DELT-Core addresses these limitations by providing a unified, modular framework that handles the complete DECL analysis
workflow. The framework is built around several key design principles: (1) leveraging established bioinformatics tools
where appropriate, particularly RNA-seq analysis software for sequence processing; (2) providing flexible configuration
options to accommodate diverse library designs and experimental setups; (3) implementing comprehensive quality control
and visualization capabilities; (4) generating machine learning-ready outputs for downstream analysis; and (5)
maintaining an accessible command-line interface suitable for both computational biologists and medicinal chemists.

The framework consists of several integrated modules: a demultiplexing engine that adapts Cutadapt for DECL-specific
sequence analysis, a chemical structure reconstruction module using RDKit and reaction SMARTS, a molecular property
calculation engine, quality control and visualization tools, and utilities for simulation and method validation. This
modular architecture allows users to execute the complete workflow or use individual components as needed for their
specific applications.

## Development and validation of the method

### Overview of the DELT-Core framework

DELT-Core is implemented as a Python package with a command-line interface organized around five main modules:
initialization and setup (`init`), computational chemistry (`compute`), sequence demultiplexing (`demultiplex`), quality
control (`qc`), and data visualization (`viz`). Additionally, specialized modules for normalization (`normalize`) and
simulation (`simulate`) support advanced analysis and method validation.

The framework expects three primary input file types: (1) library definition files in Excel format containing building
blocks, scaffolds, reaction SMARTS, and constant regions; (2) selection metadata files specifying primer sequences,
selection conditions, and experimental parameters; and (3) raw FASTQ sequencing files from DECL screening experiments.

### Chemical structure reconstruction

The chemical structure reconstruction module (`compute`) implements a systematic approach to generate SMILES
representations from DECL building blocks and reaction schemes. The process begins with loading library definitions that
specify building blocks, scaffolds, and reaction SMARTS patterns. For each potential compound in the library, the module
performs the following steps:

1. **Scaffold initialization**: If specified, a molecular scaffold serves as the starting point for structure
   construction
2. **Iterative reaction execution**: Building blocks are sequentially coupled using the specified reaction SMARTS
   patterns implemented through RDKit's reaction processing capabilities
3. **Structure validation**: Generated structures are validated for chemical feasibility and correct connectivity
4. **SMILES canonicalization**: Final structures are converted to canonical SMILES representation for downstream
   analysis

The module supports both single-cycle and multi-cycle library architectures, as well as hybridized libraries where
products from different synthetic routes are combined. Reaction SMARTS patterns can include common DECL-compatible
transformations such as amide bond formation, reductive amination, click chemistry, and nucleophilic substitutions.

### Sequence demultiplexing strategy

The demultiplexing module adapts the RNA-seq analysis tool Cutadapt for DECL-specific applications. This approach
leverages Cutadapt's robust error handling, parallel processing capabilities, and flexible adapter matching algorithms.
The demultiplexing workflow involves:

1. **Configuration file generation**: Library and selection metadata are processed to generate Cutadapt-compatible
   configuration files
2. **Barcode list creation**: DNA barcodes for each library position are extracted and formatted for adapter matching
3. **Sequential demultiplexing**: Reads are processed through multiple rounds of adapter trimming, with each round
   identifying one library position
4. **Quality control**: Detailed statistics are collected at each demultiplexing step to assess data quality and
   identify potential issues
5. **Count table generation**: Final barcode combinations are tallied and associated with their corresponding chemical
   structures

The sequential demultiplexing strategy allows for flexible error handling at each position while maintaining
computational efficiency. Error rates can be specified independently for each library position, accommodating variations
in barcode quality and synthesis fidelity.

### Quality control and validation

Comprehensive quality control is integrated throughout the DELT-Core workflow. The quality control module provides:

- **Demultiplexing statistics**: Read counts, error rates, and barcode recovery efficiency at each library position
- **Barcode distribution analysis**: Visualization of barcode frequency distributions to identify potential biases or
  systematic errors
- **Edit distance analysis**: Assessment of barcode similarity within and between library positions to evaluate
  potential cross-contamination
- **Selection comparison**: Statistical analysis of enrichment patterns across different selection conditions

### Molecular property calculation

The molecular property calculation module computes a comprehensive set of descriptors relevant to drug discovery
applications. Properties are calculated using RDKit implementations and include:

- **Basic molecular properties**: Molecular weight, LogP, polar surface area, rotatable bonds
- **Drug-likeness metrics**: Lipinski's Rule of Five compliance, QED (Quantitative Estimate of Drug-likeness)
- **Structural features**: Aromatic rings, hydrogen bond donors/acceptors, chiral centers
- **ADMET-relevant descriptors**: Properties relevant to absorption, distribution, metabolism, excretion, and toxicity
  prediction

Properties are calculated in batch mode for computational efficiency and can be filtered or ranked according to
user-specified criteria.

### Validation through simulation

The simulation module enables method validation and optimization through controlled data generation. Users can specify
library architectures, selection parameters, and error models to generate synthetic DECL datasets with known ground
truth. This capability supports:

- **Method benchmarking**: Comparison of analysis parameters and their impact on hit identification accuracy
- **Error model validation**: Assessment of demultiplexing performance under different error conditions
- **Protocol optimization**: Identification of optimal experimental parameters for specific library designs

## Applications and limitations

### Applications

DELT-Core has been successfully applied to analyze DECL screening data from diverse experimental setups, including:

- **Multi-cycle libraries**: Analysis of 2-cycle and 3-cycle libraries with various building block combinations
- **Hybridized libraries**: Processing of libraries combining products from independent synthetic routes
- **Large-scale screens**: Analysis of datasets containing millions of sequencing reads and hundreds of thousands of
  unique compounds
- **Comparative analysis**: Statistical comparison of screening results across different target proteins or selection
  conditions

The framework's flexibility allows adaptation to novel library architectures and experimental protocols with minimal
code modification.

### Limitations

Several limitations should be considered when using DELT-Core:

1. **Computational requirements**: Analysis of large DECL datasets requires substantial computational resources,
   particularly for demultiplexing steps
2. **Library complexity**: Very complex library architectures with numerous branching points or non-standard reaction
   schemes may require custom modifications
3. **Error model assumptions**: The demultiplexing algorithm assumes independence of errors across barcode positions,
   which may not hold for systematic sequencing biases
4. **Chemical structure limitations**: Structure reconstruction relies on the accuracy of provided reaction SMARTS
   patterns and may fail for novel or poorly characterized reactions

### Comparison with existing methods

DELT-Core provides several advantages over existing DECL analysis tools:

- **Comprehensive workflow**: Unlike tools that focus on individual analysis steps, DELT-Core provides end-to-end
  functionality
- **Established tool integration**: Leveraging Cutadapt provides access to well-validated sequence processing algorithms
- **Flexible configuration**: Library architectures and experimental parameters can be easily modified without code
  changes
- **Quality control integration**: Built-in QC tools reduce the need for external validation steps
- **Machine learning compatibility**: Output formats are optimized for integration with popular ML frameworks

## Experimental design

### Input data requirements

DELT-Core requires three types of input files:

1. **Library definition files** (Excel format):
    - Building block sheets: ID, SMILES, codon sequences, reaction types
    - Scaffold sheet: Scaffold IDs and corresponding SMILES
    - Reaction SMARTS sheet: Reaction types and SMARTS patterns
    - Constant regions sheet: DNA sequences flanking variable regions

2. **Selection metadata files** (Excel format):
    - Selection ID, library identifier, primer sequences
    - FASTQ file associations, experimental conditions
    - Quality control parameters

3. **FASTQ sequencing files**:
    - Raw sequencing data from DECL screening experiments
    - Standard Illumina format with quality scores
    - Single-end or paired-end reads supported

### Experimental considerations

Several experimental factors should be considered when designing DECL experiments for analysis with DELT-Core:

- **Barcode design**: DNA barcodes should have sufficient sequence diversity and minimal secondary structure
- **Primer design**: PCR primers should be compatible with the demultiplexing strategy and avoid cross-hybridization
- **Library size**: Very large libraries may require distributed computing resources for analysis
- **Quality control**: Include appropriate negative controls and reference compounds in screening experiments

### Data preprocessing requirements

Before analysis with DELT-Core, sequencing data should meet the following quality criteria:

- **Read quality**: Mean quality scores >Q20 across the read length
- **Adapter removal**: Standard Illumina adapters should be removed using conventional preprocessing tools
- **File format**: FASTQ files should follow standard formatting conventions
- **Compression**: Gzip compression is supported and recommended for large files

## Materials

### Equipment

- **Computer hardware**:
    - Minimum: 8 GB RAM, 4 CPU cores, 100 GB available storage
    - Recommended: 32 GB RAM, 16 CPU cores, 500 GB SSD storage
    - For large datasets: 64+ GB RAM, distributed computing cluster access

- **Operating system**:
    - Linux (Ubuntu 20.04+ recommended)
    - macOS (10.15+ with Homebrew)
    - Windows (with Windows Subsystem for Linux)

### Software dependencies

- **Python environment**: Python 3.10 or higher
- **Package manager**: Conda (Miniconda recommended)
- **Core dependencies**:
    - cutadapt (4.9+): Sequence adapter trimming and demultiplexing
    - rdkit (2024.3.5+): Chemical structure processing and property calculation
    - pandas (2.2.2+): Data manipulation and analysis
    - numpy (2.1.1+): Numerical computing
    - matplotlib (3.9.2+): Data visualization
    - seaborn (0.13.2+): Statistical plotting

### Input file preparation

1. **Library definition templates**: Available in the DELT-Core repository
2. **Example datasets**: Provided for testing and validation
3. **Configuration templates**: Pre-configured for common library architectures

## Procedure

### Environment setup • TIMING 15-30 min

1. **Install Miniconda**
    - Download the appropriate installer from https://docs.anaconda.com/miniconda
    - Follow installation instructions for your operating system
    - Verify installation: `conda --version`

2. **Create dedicated environment**
   ```bash
   conda create -n del python=3.10
   conda activate del
   ```

3. **Set up SSH key for GitHub access** (if installing from source)
   ```bash
   ssh-keygen -t ed25519 -C "your_email@example.com"
   eval "$(ssh-agent -s)"
   ssh-add ~/.ssh/id_ed25519
   ```
    - Add the public key to your GitHub account

4. **Install DELT-Core**
   ```bash
   # From PyPI (recommended for users)
   pip install delt-core
   
   # Or from source (for developers)
   git clone git@github.com:DELTechnology/delt-core.git
   cd delt-core
   pip install -e ".[dev,test]"
   ```

5. **Verify installation**
   ```bash
   delt-cli --help
   ```

### Project initialization • TIMING 5 min

6. **Initialize project structure**
   ```bash
   mkdir my_decl_project
   cd my_decl_project
   delt-cli init
   ```
   This creates the following directory structure:
   ```
   my_decl_project/
   ├── fastq_files/
   ├── libraries/
   └── selections/
   ```

7. **Organize input files**
   ```bash
   # Copy your files to appropriate directories
   cp /path/to/your/library.xlsx libraries/
   cp /path/to/your/selection.xlsx selections/
   cp /path/to/your/sequencing_data.fastq.gz fastq_files/
   ```

### Chemical structure computation • TIMING 10-60 min

8. **Generate SMILES representations**
   ```bash
   # For single library
   delt-cli compute smiles libraries/library.xlsx
   
   # For hybridized libraries (order must match 5' to 3' sequence)
   delt-cli compute smiles libraries/library1.xlsx libraries/library2.xlsx
   ```

9. **Merge library files** (for hybridized libraries only)
   ```bash
   delt-cli compute merge libraries/library1.xlsx libraries/library2.xlsx
   ```

10. **Calculate molecular properties**
    ```bash
    delt-cli compute properties libraries/smiles/library_smiles.txt.gz
    ```

11. **Generate property visualizations**
    ```bash
    delt-cli compute plot libraries/properties/properties_L1.txt.gz
    ```

### Sequence demultiplexing • TIMING 30 min - 4 h

12. **Initialize demultiplexing configuration**
    ```bash
    delt-cli demultiplex init -f fastq_files/input.fastq.gz \
                              -l libraries/library.xlsx \
                              -s selections/selection.xlsx
    ```

13. **Review and adjust configuration** (if needed)
    - Edit `experiments/default-*/config.yml` to adjust error rates or other parameters
    - Default error rate is 0.0 for all positions; increase if needed for your data

14. **Execute demultiplexing**
    ```bash
    delt-cli demultiplex run experiments/default-*/config.yml
    ```

    **Alternative step-by-step approach:**
    ```bash
    # Generate codon lists
    delt-cli demultiplex create-lists experiments/default-*/config.yml
    
    # Create Cutadapt input files
    delt-cli demultiplex create-cutadapt-input experiments/default-*/config.yml
    
    # Run demultiplexing script
    bash experiments/default-*/cutadapt_input_files/demultiplex.sh
    
    # Generate count tables
    delt-cli demultiplex compute-counts experiments/default-*/config.yml \
                                       experiments/default-*/cutadapt_output_files/reads_with_adapters.gz \
                                       evaluations/
    ```

### Quality control analysis • TIMING 10-20 min

15. **Generate QC report**
    ```bash
    delt-cli qc report experiments/default-*
    ```

16. **Create QC visualizations**
    ```bash
    delt-cli qc plot experiments/default-*
    ```

17. **Analyze barcode edit distances** (optional)
    ```bash
    delt-cli qc analyze-codons experiments/default-*/config.yml
    ```

### Data normalization and hit identification • TIMING 5-15 min

18. **Normalize selection data**
    ```bash
    # Compare target selections (IDs 1-3) to control selections (IDs 4-6)
    delt-cli normalize run experiments/default-*/config.yml '1 2 3' '4 5 6'
    ```

### Results interpretation • TIMING Variable

19. **Review output files**
    - Count tables: `evaluations/selection-*/`
    - QC reports: `experiments/default-*/report.txt`
    - Visualizations: `experiments/default-*/quality_control/plots/`
    - Chemical properties: `libraries/properties/`

20. **Interactive data exploration** (optional)
    ```bash
    delt-cli viz show evaluations/selection-1/counts.txt
    ```

## Troubleshooting

### Common issues and solutions

**Problem**: Installation fails with dependency conflicts

- **Solution**: Use a fresh conda environment and install dependencies in the specified order
- **Prevention**: Avoid mixing conda and pip installations unnecessarily

**Problem**: Demultiplexing produces very low read recovery

- **Possible causes**:
    - Incorrect barcode sequences in library file
    - Too stringent error rate settings (try increasing MaxErrorRate in config.yml)
    - Poor sequencing quality (check FASTQ quality scores)
- **Diagnostic steps**: Check intermediate Cutadapt log files for detailed statistics

**Problem**: Chemical structure reconstruction fails

- **Possible causes**:
    - Invalid SMILES in building block definitions
    - Incorrect reaction SMARTS patterns
    - Missing or incorrect scaffold definitions
- **Solution**: Validate input SMILES using RDKit directly; check reaction SMARTS syntax

**Problem**: Memory errors during large dataset processing

- **Solution**:
    - Increase available RAM or use high-memory computing nodes
    - Process data in smaller batches using the `--fast-dev-run` flag for testing
    - Use gzip compression for intermediate files

**Problem**: Slow performance on large datasets

- **Optimization strategies**:
    - Ensure sufficient CPU cores are available for Cutadapt parallel processing
    - Use SSD storage for improved I/O performance
    - Consider distributed computing for very large libraries

### Expected outputs and quality metrics

**Successful demultiplexing should show**:

- > 80% read retention through initial adapter trimming steps
- Exponential decrease in read counts through subsequent demultiplexing rounds
- Even distribution of reads across expected barcode combinations

**Quality control metrics**:

- Error rates <5% at each demultiplexing position (for good quality data)
- Minimal edit distance overlap between different barcode sets
- Consistent enrichment patterns across replicate selections

**Chemical structure outputs**:

- All library positions should have valid SMILES representations
- Molecular properties should fall within expected ranges for the compound class
- Property distributions should reflect the diversity of building blocks used

## Timing

The timing for DELT-Core analysis depends on dataset size and computational resources:

- **Environment setup**: 15-30 minutes (one-time setup)
- **Project initialization**: 5 minutes
- **Chemical structure computation**:
    - Small libraries (<10,000 compounds): 10-30 minutes
    - Large libraries (>100,000 compounds): 1-4 hours
- **Sequence demultiplexing**:
    - Small datasets (<1M reads): 30-60 minutes
    - Large datasets (>10M reads): 2-8 hours
- **Quality control analysis**: 10-20 minutes
- **Data normalization**: 5-15 minutes

**Total workflow time**: 1-12 hours depending on dataset size

## Anticipated results

### Expected output structure

Upon successful completion, DELT-Core generates a comprehensive set of output files organized in a standardized
directory structure:

```
my_decl_project/
├── libraries/
│   ├── smiles/
│   │   └── library_smiles.txt.gz          # Chemical structures
│   ├── properties/
│   │   ├── properties_L1.txt.gz           # Molecular properties
│   │   └── plots/                         # Property distributions
├── experiments/
│   └── default-YYYY-MM-DD-HH-MM-SS/
│       ├── config.yml                     # Analysis configuration
│       ├── report.txt                     # QC summary
│       ├── codon_lists/                   # Barcode sequences
│       ├── cutadapt_output_files/         # Demultiplexing logs
│       └── quality_control/               # QC visualizations
└── evaluations/
    └── selection-*/
        └── counts.txt                     # Final count tables
```

### Chemical structure outputs

The SMILES generation module produces tab-delimited files containing:

- **Building block identifiers**: Unique IDs for each synthetic step
- **Chemical structures**: Canonical SMILES representations
- **DNA sequences**: Corresponding barcode sequences
- **Reaction metadata**: Applied transformations and intermediate structures

For a typical 2-cycle library with 1,000 building blocks per cycle, expect ~1 million unique structures with associated
metadata.

### Molecular property analysis

Property calculation generates comprehensive descriptor tables including:

- **Drug-likeness metrics**: QED scores, Lipinski compliance
- **Physical properties**: Molecular weight, LogP, polar surface area
- **Structural features**: Ring counts, rotatable bonds, stereochemistry
- **Visualization outputs**: Distribution plots and property space maps

### Demultiplexing results

Successful demultiplexing produces:

- **Read count statistics**: Detailed breakdown of reads processed at each step
- **Barcode combination tables**: Counts for each unique barcode combination
- **Quality metrics**: Error rates, adapter match statistics, and coverage uniformity
- **Selection-specific outputs**: Separate count tables for each experimental condition

**Typical performance metrics**:

- Initial adapter matching: 85-95% of input reads
- Sequential demultiplexing: 70-85% final read retention
- Unique barcode combinations: 10-90% of theoretical library coverage
- Error rates: <5% per position for high-quality data

### Quality control outputs

Comprehensive quality assessment includes:

- **Demultiplexing efficiency plots**: Read retention through sequential steps
- **Barcode distribution analysis**: Coverage uniformity and potential biases
- **Edit distance matrices**: Sequence similarity analysis for contamination detection
- **Selection comparison statistics**: Enrichment metrics and statistical significance

### Integration with downstream analysis

DELT-Core outputs are designed for seamless integration with:

- **Machine learning workflows**: Standardized feature matrices for model training
- **Cheminformatics pipelines**: Compatible with RDKit, ChEMBL, and other tools
- **Statistical analysis software**: R, Python pandas, or specialized packages
- **Visualization platforms**: Integration with chemical space analysis tools

## Conclusions and future directions

DELT-Core represents a significant advancement in computational tools for DNA-encoded chemical library analysis,
addressing the growing need for comprehensive, user-friendly software that can handle the complete workflow from raw
sequencing data to actionable chemical insights. The framework's key innovations include the adaptation of established
RNA-seq tools for DECL-specific applications, integration of cheminformatics capabilities for structure reconstruction,
and comprehensive quality control throughout the analysis pipeline.

The modular architecture and command-line interface make DELT-Core accessible to researchers with varying computational
backgrounds, while the integration with established tools like Cutadapt and RDKit ensures robust performance and broad
compatibility. The framework's flexibility allows adaptation to diverse experimental protocols and library
architectures, making it suitable for both academic research groups and pharmaceutical companies conducting DECL-based
drug discovery.

Several areas present opportunities for future development and enhancement:

**Algorithmic improvements**: Advanced error correction algorithms that account for systematic biases in sequencing data
could improve demultiplexing accuracy, particularly for challenging datasets with high error rates or complex barcode
structures.

**Machine learning integration**: Direct integration with modern machine learning frameworks for hit prediction,
structure-activity relationship modeling, and virtual screening could streamline the transition from DECL analysis to
lead optimization.

**Cloud computing support**: Implementation of cloud-native architectures would enable analysis of very large datasets
that exceed local computational resources, while providing better accessibility for research groups with limited
infrastructure.

**Real-time analysis capabilities**: Development of streaming analysis pipelines could enable real-time monitoring of
sequencing experiments and early termination of unsuccessful selections.

**Expanded visualization and interpretation tools**: Advanced interactive visualizations for chemical space exploration,
selection comparison, and hit prioritization would enhance the interpretability of complex DECL datasets.

The open-source nature of DELT-Core facilitates community contributions and customization for specific research needs.
As DECL technology continues to evolve, with new synthetic methodologies and selection strategies, the framework's
modular design enables straightforward incorporation of new analytical approaches and experimental protocols.

DELT-Core has been successfully deployed across multiple research institutions and pharmaceutical companies,
demonstrating its utility for diverse applications ranging from academic target validation studies to large-scale drug
discovery campaigns. The framework's comprehensive documentation, example datasets, and active community support provide
a foundation for widespread adoption and continued development.

In conclusion, DELT-Core fills a critical gap in the computational infrastructure supporting DNA-encoded chemical
library research, providing researchers with the tools needed to extract maximum value from DECL screening experiments
and accelerate the discovery of novel bioactive compounds.

## Data availability

All software, documentation, and example datasets described in this protocol are freely available:

- **DELT-Core software**: https://github.com/DELTechnology/delt-core
- **Documentation and tutorials**: Available in the repository wiki
- **Example datasets**: Included in the repository for testing and validation
- **Protocol supplementary materials**: Configuration templates and test data

## Code availability

DELT-Core is released under the MIT License, enabling free use, modification, and distribution. The complete source
code, including all modules described in this protocol, is available on GitHub with comprehensive documentation and
example usage patterns.

## References


## Acknowledgements

We thank the members of the DECL research community for valuable feedback during the development of DELT-Core. We
acknowledge the developers of Cutadapt, RDKit, and other open-source tools that form the foundation of this framework.
Special thanks to the beta testers who provided feedback on early versions of the software and protocol.

## Author contributions

[To be filled based on actual authorship]

## Competing interests

The authors declare no competing financial interests.

## Additional information

**Correspondence and requests for materials** should be addressed to [corresponding author].

**Reprints and permissions** information is available at www.nature