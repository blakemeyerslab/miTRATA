About Tailing Analysis:

miTRATA is the first web-based tool for microRNA Truncation and Tailing Analysis—the analysis of 3′ modifications of microRNAs including the loss or gain of nucleotides relative to the canonical sequence. miTRATA is implemented in Python and employs parallel processing modules to enhance its scalability when analyzing multiple small RNA (sRNA) sequencing datasets. It utilizes miRBase, currently version 21, as a source of known microRNAs for analysis. miTRATA notifies user(s) via email to download as well as visualize the results online. 

miTRATA’s strengths lie in (i) its biologist-focused web interface, (ii) improved scalability via parallel processing and (iii) its uniqueness as a webtool to perform microRNA truncation and tailing analysis.

The service consists of a website (under "Web") and a background job (under "Handler").
The service needs its own database and working directory.
The service reads from $ALLDATA.
The service reads from genome databases used by the "next_gen" service.
The genome databases are expected to reside on the DB server of the "next_gen" service.
The website is not associated with a particular genome, expression or methylation database.
The website is usually viewed as a child website of a "next_gen" website.

The "service_name" is "ta".  The "website_name" is "ta".

It requires the following inputs:
1. Small RNA sequence files. (tag count file)
2. Genome of interest. (aka the bowtie index)
3. List of mature miRNA sequence(s) from miRBase. (FASTA format)

Standalone Script parameters:

PATH to the INPUT FOLDER- Folder containing sequence file(s) in Tag_count format with “tag” preceding “count”; no blank lines separating sequences are permitted. Tag should be separated from count by only 1 TAB.

PATH to the GENOME INDEX- Bowtie index for the plant/animal of interest

miRNA file- Sequences must be in FASTA format with the ">name" proceeding each sequence; no blank lines separating sequences are permitted. Only A, U, C, G, T are allowed. Keep the name of the sequences short to retain optimal quality of the figures. For example, instead of ">osa-miR523f MIMAT0023516 Oryza sativa" as a name for a sequence in the FASTA file, use ">osa-miR818f"

PATH to the OUTPUT FOLDER- Output folder for results

Running the script
python3.4 ["PATH to the INPUT FOLDER"] ["PATH to the GENOME INDEX"] [miRNA file] ["PATH to the OUTPUT FOLDER"]

Example Usage python3.4 Tailing_Pipeline_v13.py "/Tailing /Input " "/GenomeIndex/MAIZE_AGPv2_genome" miRNA.fa "/Tailing /Output"
Output File
Results Folder
merged_full.pdf — a single pdf with all images
individual pdfs
text files stores tailing and truncation information (i.e. tail patterns and their lengths)
