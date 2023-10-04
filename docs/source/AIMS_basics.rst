AIMS Basics
=====

In this section, we will focus on the absolute basics that you need to know about while running AIMS. This includes the :ref:`formatting` of the files used for the analysis, the :ref:`core` of the software, and the :ref:`bphysProp` that are the central pillars of the analysis. This information is key whether you're using the GUI, the CLI, or the Jupyter notebook. The hope is that if there are any fundamental questions with the meaning of a given output, or necessary troubleshooting with inputs, the users can find them here. While details are provided here, users can also check out the :doc:`Testing` section to learn how to download example data and use this to compare formatting.

.. _formatting:

Input Formatting
------------

At present, the AIMS software requires a specific input formatting to be read by the software. Future versions of the software will hopefully have more relaxed formatting requirements, assuming users request such functionality. Here we will go into specifics for each molecular species, but the general format maintained for all species involves a comma separated value (CSV) file with a header in the first row, and no associated metadata in the given file.

**Immunoglobulin (TCR/Ab) Formatting**

The immunoglobulin (Ig) input file formatting is the most flexible of the three, with options meant to satisfy the needs of any given AIMS input data. Specifically these data must contain only the complementarity determining region (CDR) loops, as AIMS excludes the framework regions of antibodies and TCRs in the analysis. An example of the proper input formatting for the first few lines of an input file can be seen below:

.. code-block:: python
    
    [Example file from AIMS/app/ab_testData/flu_poly.csv]
    l1,l2,l3,h1,h2,h3
    QSISSY,DAS,QHRSTWPPN,GGTFSSRA,IIPIFNTP,AREMATIFGRMDV
    ESLLHSDGKTY,EVS,MQTIQLPGT,GGIMRRNG,IIAIFGTP,VASSGYHLHRETWGY
    QDIKNY,HVS,HQCYNLPYT,GFIFGHFA,ISGGGLNT,ARFDSSGYNYVRGMVV
    .
    .
    .

The key features are that 
1. the general format must follow that of a comma separated value (csv) file. In this csv file each row represents a unique sequence, each column represents a given CDR loop, and each column is separated by a comma. 
2. Each column must have a header with no sequence information in it. The contents of this header are not critical, as the standard AIMS Ig file loader disregards this header and replaces it. Descriptive headers are preferred, if only for downstream transparency of the data. 
3. Single letter amino acid codes should be used, with only capitalized letters. The defined function "aimsLoad.convert_3let" can convert three letter codes to single letter codes [see API for more details... once I make this section]. 
4. The sequences do not have any extranneous characters or spaces. Currently AIMS is capable of identifying an "X" in a sequence and by default removes these sequences. Any other non-amino acid single letter characters will result in an error, and spaces will be encoded into the AIMS matrix, which could confound analysis.

At present, AIMS does NOT identify other issues with sequences. Missing CDR loops in a sequence, spaces included in a sequence, or other miscellaneous mistakes will not result in an error, and could lead to inaccuracies in the downstream analysis. Either visual inspection or other quality control steps should be taken to ensure proper analysis. If analyzing multiple files at once, all input files must have the same number of loops.

.. note::
    Currently, the AIMS GUI cannot analyze TCR/Ab inputs of 4 or 5 loops. For most repertoire analysis, the requirement to analyze such a dataset would be unexpected. Users should submit an issue on the GitHub if this analysis is needed for some reason.

**MHC/MHC-Like Formatting**

The current MHC/MHC-like analysis requires the most rigid and restrictive of all the input formatting, but this will hopefully change in future updates. Users can either input entire MHC sequences or just the platform domain sequences. Additionally, users can leverage these restrictive requirements to analyze other molecular species if needed. The initial input is a simple aligned FASTA, as seen below:

.. code-block:: python
    
    [Example file from AIMS/app/mhc_testData/hlaA_seqs.fasta]
    >3VJ6_A Chain A, H-2 Class I Histocompatibility Antigen, D-37 Alpha Chain [Mus musculus]
    ------------------------------------------------------------
    ------------------------------------------------------------
    -------------------------------------------------SPHSLRYFTTA
    VSRPGLGEPRFIIVGYVDDTQFVRFDSDAENPRMEPRARWIEQEGPEYWERETWKAR
    DMGRNFRVNLRTLLGYYNQSNDESHTLQWMYGCDVGPDGRLLRGYCQEAYDGQDYISLNE
    DLRSWTANDIASQISKHKSEAVDEAH-QQRAYLQGPCVEWLHRYLRLGNETLQRSDPPKA
    HVTHHPRSEDEVTLRCWALGFYPADITLTWQLNGEELTQDMELVETRPAGDGTFQKWAAV
    VVPLGKEQYYTCHVYHEGLPEPLTLRWEPP------------------------------
    -------------------------------------------------
    >5VCL_A Chain A, H2-t23 Protein [Mus musculus]
    ------------------------------------------------------------
    ------------------------------------------------------------
    ------------------------------------------------MSSHSLRYFHTA
    .
    .
    .

Importantly, each FASTA entry must be pre-aligned using BLAST or a similar alignment software. AIMS does not internally align the sequences, and requires that the inputs can be expected to be structurally very similar. For MHC and MHC-like molecules, this requirement is satisfied. If a subset of sequences align poorly for some reason, they can be included as a separate file. Each individual file will have its own user-specified region of the alignment that will ultimately be input into the analysis. The user specification can be done on-the-fly, or input as a separate file formatted as such:

.. code-block:: python
    
    [Example file from aims_immune/app_data/test_data/mhcs/ex_cd1_hla_uda_uaa.csv]
    Name,S1s,S1e/H1s,H1e/S2s,S2e/H2s,H2e
    cd1,124,167,209,262,303
    hla,170,210,260,306,348
    uda,2,49,93,152,193
    uaa,2,49,93,152,193

The above file is formatted again as a comma separated value (csv), with the first column giving the name of the dataset, and the remaining columns identifying the start and end point of four distinct structural features in the provided FASTA alignment. Specifically for the analysis of MHC and MHC-like molecules, these four structural features are the beta-strand of the alpha 1 domain, the alpha helix of the alpha 1 domain, the beta-strand of the alpha 2 domain, and the alpha helix of the alpha 2 domain. Each number represents either the start of one structural feature, the end of another structural feature, or both. In the example file, for the hla alignmemnt (corresponding to the FASTA above) the first beta strand starts at alignment position 170 and ends at position 210. Likewise, the first alpha helix starts at position 210 and ends at position 260. And so on.

Currently, the Phyre server (http://www.sbg.bio.ic.ac.uk/phyre2/html/page.cgi?id=index) is recommended to identify these structural features. Other software may be used to identify the key structural features for analysis, but the numbering provided in standard Phyre outputs makes translation to the above csv file easy. Generally only one sequence should be necessary to be used as input, as structural similarity is a requirement for comparable analysis using AIMS. Users can take advantage of this ambiguity in the software to analyze any four connected structual features in evolutionarily and structurally related molecules in the AIMS GUI. Users comfortable with the Jupyter Notebooks can instead follow the Multi-Sequence Alignment formatting instructions.

**Immunopeptidomics Formatting**

This is the first of two sections that are not yet implemented in the AIMS GUI, but can be analyzed using the AIMS notebook or CLI. Specifically, AIMS can be used to analyze immunopeptidomics data. Again the input is simply a comma separated value (csv) formatted file. However, since the input should only have one column, the precise format is a little less important. An example can be seen below:

.. code-block:: python
    
    [Example file from aims_immune/app_data/test_data/peptides/pancreas_hla_atlas.csv]
    sequence
    ALVSGNNTVPF
    TYRGVDLDQLL
    NYIDIVKYV
    SYIPIFPQ
    NYFPGGVALI
    .
    .
    .

Example data provided from the HLA Ligand Atlas (https://hla-ligand-atlas.org/welcome). In future releases, data related to mass spectrometry approaches used for the identification of these peptides will be included in the analysis. Metadata can be included in additional columns of a separate csv.

**Multi-Sequence Alignment Formatting**

Again, this multi-sequence alignment input is not yet available in the AIMS GUI, but is available in the notebook and CLI. As it turns out, the same file formatting that is used for loading MHC molecules works for the more general MSA input. The difference is that the old MHC module (and by extension, the GUI) required a subset of the MSA to be selected. This step is now optional, so if you'd like to pre-select certain regions of an MSA or input the entire MSA, you can do so! Careful though, very large sequences will likely process quite slowly in AIMS.

.. _core:

Core Functionalities
------------

Functionalities coming soon!

.. _bphysProp:

Biophysical Properties
------------

In generating the core biophysical property matrix of the AIMS analysis, the same 61 biophysical properties are used in all analyses, with an option to use fewer if the user decides to. The properties are listed in the table below:

.. list-table:: Table of AIMS Biophysical Properties
  :widths: 20 40 40
  :header-rows: 1

  * - Number
    - Property [Shorthand]
    - Decription
  * - 0
    - Hydrophobicity1 [Phob1]
    - Hydrophobicity Scale [-1,1]
  * - 1
    - Charge [Charge]
    - Charge [ec]
  * - 2
    - Hydrophobicity2 [Phob2]
    - Octanol-Interface Hydrophobicity Scale
  * - 3
    - Bulkiness [Bulk]
    - Side-Chain Bulkiness
  * - 4
    - Flexibility [Flex]
    - Side-Chain Flexibility
  * - 5 
    - Kidera 1 [KD1]
    - Helix/Bend Preference
  * - 6
    - Kidera 2 [KD2]
    - Side-Chain Size
  * - 7
    - Kidera 3 [KD3]
    - Extended Structure Preference
  * - 8
    - Kidera 4 [KD4]
    - Hydrophobicity
  * - 9 
    - Kidera 5 [KD5]
    - Double-bend Preference
  * - 10
    - Kidera 6 [KD6]
    - Flat Extended Preference
  * - 11
    - Kidera 7 [KD7]
    - Partial Specific Volume
  * - 12
    - Kidera 8 [KD8]
    - Occurrence in alpha-region
  * - 13
    - Kidera 9 [KD9]
    - pK-C
  * - 14
    - Kidera 10 [KD10]
    - Surrounding Hydrophobicity
  * - 15
    - Hotspot 1 [HS1]
    - Normalized Positional Residue Freq at Helix C-term
  * - 16
    - Hotspot 2 [HS2]
    - Normalized Positional Residue Freq at Helix C4-term
  * - 17
    - Hotspot 3 [HS3]
    - Spin-spin coupling constants
  * - 18
    - Hotspot 4 [HS4]
    - Random Parameter
  * - 19
    - Hotspot 5 [HS5]
    - pK-N
  * - 20
    - Hotspot 6 [HS6]
    - Alpha-Helix Indices for Beta-Proteins
  * - 21
    - Hotspot 7 [HS7]
    - Linker Propensity from 2-Linker Dataset
  * - 22
    - Hotspot 8 [HS8]
    - Linker Propensity from Long Dataset
  * - 23
    - Hotspot 9 [HS9]
    - Normalized Relative Freq of Helix End
  * - 24
    - Hotspot 10 [HS10]
    - Normalized Relative Freq of Double Bend
  * - 25
    - Hotspot 11 [HS11]
    - pK-COOH
  * - 26
    - Hotspot 12 [HS12]
    - Relative Mutability
  * - 27
    - Hotspot 13 [HS13]
    - Kerr-Constant Increments
  * - 28
    - Hotspot 14 [HS14]
    - Net Charge
  * - 29
    - Hotspot 15 [HS15]
    - Norm Freq Zeta-R
  * - 30
    - Hotspot 16 [HS16]
    - Hydropathy Scale
  * - 31
    - Hotspot 17 [HS17]
    - Ratio of Average Computed Composition
  * - 32
    - Hotspot 18 [HS18]
    - Intercept in Regression Analysis
  * - 33
    - Hotspot 19 [HS19]
    - Correlation coefficient in Reg Anal
  * - 34
    - Hotspot 20 [HS20]
    - Weights for Alpha-Helix at window pos
  * - 35
    - Hotspot 21 [HS21]
    - Weights for Beta-sheet at window pos -3
  * - 36
    - Hotspot 22 [HS22]
    - Weights for Beta-sheet at window pos 3
  * - 37
    - Hotspot 23 [HS23]
    - Weights for coil at win pos -5
  * - 38
    - Hotspot 24 [HS24]
    - Weights coil win pos -4
  * - 39
    - Hotspot 25 [HS25]
    - Weights coil win pos 6
  * - 40
    - Hotspot 26 [HS26]
    - Avg Rel Frac occur in AL
  * - 41
    - Hotspot 27 [HS27]
    - Avg Rel Frac occur in EL
  * - 42
    - Hotspot 28 [HS28]
    - Avg Rel Frac occur in A0
  * - 43
    - Hotspot 29 [HS29]
    - Rel Pref at N
  * - 44
    - Hotspot 30 [HS30]
    - Rel Pref at N1
  * - 45
    - Hotspot 31 [HS31]
    - Rel Pref at N2
  * - 46
    - Hotspot 32 [HS32]
    - Rel Pref at C1
  * - 47
    - Hotspot 33 [HS33]
    - Rel Pref at C
  * - 48
    - Hotspot 34 [HS34]
    - Information measure for extended without H-bond
  * - 49
    - Hotspot 35 [HS35]
    - Information measure for C-term turn
  * - 50
    - Hotspot 36 [HS36]
    - Loss of SC hydropathy by helix formation
  * - 51
    - Hotspot 37 [HS37]
    - Principal Component 4 (Sneath 1966)
  * - 52
    - Hotspot 38 [HS38]
    - Zimm-Bragg Parameter
  * - 53
    - Hotspot 39 [HS39]
    - Normalized Freq of ZetaR
  * - 54
    - Hotspot 40 [HS40]
    - Rel Pop Conformational State A
  * - 55
    - Hotspot 41 [HS41]
    - Rel Pop Conformational State C
  * - 56
    - Hotspot 42 [HS42]
    - Electron-Ion Interaction Potential
  * - 57
    - Hotspot 43 [HS43]
    - Free energy change of epsI to epsEx
  * - 58
    - Hotspot 44 [HS44]
    - Free energy change of alphaRI to alphaRH
  * - 59
    - Hotspot 45 [HS45]
    - Hydrophobicity coeff
  * - 60 
    - Hotspot 46 [HS46]
    - Principal Property Value z3 (Wold et. al. 1987)

The so-called Kidera factors are from the published work: 

Kidera et al. Statistical analysis of the physical properties of the 20 naturally occurring amino acids.
Journal of Protein Chemistry (1985)

While the hotspot variables mentioned above are from:

Liu et al. Hot spot prediction in protein-protein interactions by an ensemble system.
BMC Systems Biology (2018)