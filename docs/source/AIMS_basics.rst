AIMS Basics
=====

In this section, we will focus on the absolute basics that you need to know about while running AIMS. This includes the :ref:`formatting` of the files used for the analysis, the :ref:`core` of the software, and the :ref:`bphysProp` that are the central pillars of the analysis. This information is key whether users are interested in working with the AIMS GUI or with the AIMS Jupyter Notebooks. The hope is that if there are any fundamental questions with the meaning of a given output, or necessary troubleshooting with inputs, the users can find them here.

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

The key features are that 1. the general format must follow that of a comma separated value (csv) file. In this csv file each row represents a unique sequence, each column represents a given CDR loop, and each column is separated by a comma. 2. Each column must have a header with no sequence information in it. The contents of this header are not critical, as the standard AIMS Ig file loader disregards this header and replaces it. Descriptive headers are preferred, if only for downstream transparency of the data. 3. Single letter amino acid codes should be used, with only capitalized letters. The defined function "aimsLoad.convert_3let" can convert three letter codes to single letter codes [see API for more details... once I make this section]. 4. The sequences do not have any extranneous characters or spaces. Currently AIMS is capable of identifying an "X" in a sequence and by default removes these sequences. Any other non-amino acid single letter characters will result in an error, and spaces will be encoded into the AIMS matrix, which could confound analysis.

At present, AIMS does NOT identify other issues with sequences. Missing CDR loops in a sequence, spaces included in a sequence, or other miscellaneous mistakes will not result in an error, and could lead to inaccuracies in the downstream analysis. Either visual inspection or other quality control steps should be taken to ensure proper analysis. If analyzing multiple files at once, all input files must have the same number of loops.

.. note::
    Currently, AIMS cannot analyze TCR/Ab inputs of 4 or 5 loops. For most repertoire analysis, the requirement to analyze such a dataset would be unexpected. Users should submit an issue on the GitHub if this analysis is needed for some reason.

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
    
    [Example file from AIMS/app/mhc_testData/ex_cd1_hla_uda_uaa.csv]
    Name,S1s,S1e/H1s,H1e/S2s,S2e/H2s,H2e
    cd1,124,167,209,262,303
    hla,170,210,260,306,348
    uda,2,49,93,152,193
    uaa,2,49,93,152,193

The above file is formatted again as a comma separated value (csv), with the first column giving the name of the dataset, and the remaining columns identifying the start and end point of four distinct structural features in the provided FASTA alignment. Specifically for the analysis of MHC and MHC-like molecules, these four structural features are the beta-strand of the alpha 1 domain, the alpha helix of the alpha 1 domain, the beta-strand of the alpha 2 domain, and the alpha helix of the alpha 2 domain. Each number represents either the start of one structural feature, the end of another structural feature, or both. In the example file, for the hla alignmemnt (corresponding to the FASTA above) the first beta strand starts at alignment position 170 and ends at position 210. Likewise, the first alpha helix starts at position 210 and ends at position 260. And so on.

Currently, the Phyre server (http://www.sbg.bio.ic.ac.uk/phyre2/html/page.cgi?id=index) is recommended to identify these structural features. Other software may be used to identify the key structural features for analysis, but the numbering provided in standard Phyre outputs makes translation to the above csv file easy. Generally only one sequence should be necessary to be used as input, as structural similarity is a requirement for comparable analysis using AIMS. Users can take advantage of this ambiguity in the software to analyze any four connected structual features in evolutionarily and structurally related molecules in the AIMS GUI. Users comfortable with the Jupyter Notebooks can instead follow the Multi-Sequence Alignment formatting instructions.

**Immunopeptidomics Formatting**

This is the first of two sections that are not yet implemented in the AIMS GUI, but can be analyzed using the AIMS notebooks. Specifically, the AIMS_peptide.ipynb file can be used to analyze immunopeptidomics data. Again the input is simply a comma separated value (csv) formatted file. However, since the input should only have one column, the precise format is a little less important. An example can be seen below:

.. code-block:: python
    
    [Example file not yet on the GitHub]
    pep_input
    APATPAVVL
    PSPEAAVAV
    VALGGPHDP
    PAALPVPSL
    PTAPVTPSI
    CPGASQPIL
    DRGSCGVTV
    .
    .
    .

In future releases, data related to mass spectrometry approaches used for the identification of these peptides will be included in the analysis. Metadata will be included in additional columns of the csv.

**Multi-Sequence Alignment Formatting**

Again, this multi-sequence alignment input is not yet available in the AIMS GUI, is available in the notebook. Unlike the MHC input formatting described above, structural features should be identified before being loaded into AIMS. Each entry should have the exact same number of characters, and an associated sequence name, as seen below: 

.. code-block:: python
    
    [Example file not yet on the GitHub]
    Dpr,Sequence
    Dpr1,DKDVSWIRKRDLHILTAGGTTYTSD-----QINTEPKMSLSYTFNVVEL
    Dpr2,DKSVSWIRKRDLHILTAGILTYTSD-----QVNTEPKISMAFRLNVIVT
    Dpr3,DKSVSWIRKRDLHILTVGTATYTSD-----QVNTEPKMSMAFQLNIIEI
    .
    .
    .

The required sequence name is used in the downstream analysis in this case, to help identify potential clusters of biophysically similar sequences. Unlike in the analysis of TCR, MHC, or peptide sequences, the associated names or identifiers for each sequence are more likely to be clinically or evolutionarily important. 

.. _core:

Core Functionalities
------------

Functionalities coming soon!

.. _bphysProp:

Biophysical Properties
------------

Details on the BPHYS properties coming soon!