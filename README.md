miTRATA: a tool for microRNA Truncation and Tailing Analysis
-----

Version: v1.3   
Updated: 08/27/2020    
Original author: Parth Patel (parth1415@gmail.com)    

Please refer to the [wiki](https://github.com/pupatel/miTRATA/wiki) page for the overview of this pipeline along with the detailed information about a standalone version, including scripts, installing dependencies, usage, and output files.


About the standalone version of this tool (this repository)
--

See information on [scripts and dependencies](https://github.com/pupatel/miTRATA/wiki/Scripts-and-Dependencies).

Script parameters:

1. PATH to the INPUT FOLDER- Folder containing sequence file(s) in Tag_count format with “tag” preceding “count”; no blank lines separating sequences are permitted. Tag should be separated from count by only 1 TAB.

2. PATH to the GENOME INDEX- Bowtie index for the plant/animal of interest

3. miRNA file- Sequences must be in FASTA format with the ">name" proceeding each sequence; no blank lines separating sequences are permitted. Only A, U, C, G, T are allowed. Keep the name of the sequences short to retain optimal quality of the figures. For example, instead of ">osa-miR523f MIMAT0023516 Oryza sativa" as a name for a sequence in the FASTA file, use ">osa-miR818f"

4. PATH to the OUTPUT FOLDER- Output folder for results

Running the script:   
python3.4 ["PATH to the INPUT FOLDER"] ["PATH to the GENOME INDEX"] [miRNA file] ["PATH to the OUTPUT FOLDER"]

Output File   
Results Folder   
merged_full.pdf — a single pdf with all images   
individual pdfs   
text files store tailing and truncation information (i.e. tail patterns and their lengths)   


About the web-based version of this tool 
--

URL:  https://wasabi.ddpsc.org/~apps/ta/   
published: [Bioinformatics](https://academic.oup.com/bioinformatics/article/32/3/450/1743711) 

miTRATA is the first web-based tool for microRNA Truncation and Tailing Analysis—the analysis of 3′ modifications of microRNAs including the loss or gain of nucleotides relative to the canonical sequence. miTRATA is implemented in Python and employs parallel processing modules to enhance its scalability when analyzing multiple small RNA (sRNA) sequencing datasets. It utilizes miRBase, currently version 21, as a source of known microRNAs for analysis. miTRATA notifies user(s) via email to download as well as visualize the results online. 

miTRATA’s strengths lie in 1) its biologist-focused web interface, 2) improved scalability via parallel processing, and 3) its uniqueness as a webtool to perform microRNA truncation and tailing analysis.

It requires the following inputs:
1. Small RNA sequence files. (tag count file)
2. Genome of interest (aka the bowtie index)
3. List of mature miRNA sequence(s) from miRBase (FASTA format)
