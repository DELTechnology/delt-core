# Methods

## Input Data and Environment

All analyses were conducted on a high-performance computing (HPC) system running Linux. Parallel data compression and decompression were enabled via the pigz utility (v2.8) to accelerate file operations on large FASTQ datasets. The input consisted of raw sequencing reads in FASTQ format (*.fastq.gz), generated from DNA-encoded chemical library (DECL) selections.

## Adapter Trimming and Sequential Filtering

Adapter trimming was performed using Cutadapt (v4.9), executed in a sequential, multi-step fashion to match and remove a predefined set of barcode and linker sequences. Each adapter sequence was stored in an individual FASTQ file and passed to Cutadapt via the -g "^file:<path>" parameter to enable anchored matching from the start of the read. Matching was performed with an error tolerance of 0 (-e 0.0) and with indels disabled (--no-indels) to ensure exact sequence recognition. Reads lacking the expected adapter were discarded (--discard-untrimmed).

After each Cutadapt run:
	1.	The output FASTQ file was renamed as the input for the next step.
	2.	A new adapter file was supplied for matching.
	3.	The process was repeated until all expected adapters (S1, C1, B1, C2, B2, C3, B3, C4, S2) were sequentially matched.

Cutadapt logs (*.cutadapt.log) and structured JSON reports (*.cutadapt.json) were generated for each step to record trimming statistics, including read retention, adapter frequency, and quality metrics. Processing was parallelized using --cores=48.

## Extraction of Adapter-Matching Reads

Following the final adapter-matching step, all retained reads were extracted with zgrep by selecting entries with adapter annotations (@ in the FASTQ header). The resulting subset was compressed into a file (reads_with_adapters.gz) using pigz for downstream analysis.

## Demultiplexing and Read Counting

Demultiplexing and count computation were performed using the delt-cli demultiplex compute-counts command from the DELT-CLI toolkit. The process took as input:
	1.	A configuration file specifying barcode-to-construct mappings (config.yml).
	2.	The compressed set of adapter-containing reads (reads_with_adapters.gz).
	3.	An evaluation output directory for count tables and summary metrics.

This step generated read count matrices stratified by barcode identifiers, enabling quantitative analysis of the selection results.
