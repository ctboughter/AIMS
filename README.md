# AIMS - An Automated Immune Molecule Separator

# Quick Start
As of AIMS v0.9, everything should be nicely wrapped up as an installable pypi package. You can simply install the AIMS GUI, CLI, and notebook using pip:

```
pip install aims-immune
```

Whether you are a new or returning AIMS user, it is strongly recommended you check out the documentation (see below) to learn details about formatting and usage. For returning users especially, the way AIMS is called has changed completely.

# Description
The primary goal of AIMS is to identify discriminating factors between two distinct sets of immune molecules. As of versions 0.8 and later, the software is now capable of analyzing any set of sequences with general conservation and localized diversity. AIMS has specific analysis modes for Immunoglobulins (Ig - T Cell Receptors and Antibodies) and Peptides (Specifically those isolated from MHC), as well as a more general multi-sequence alignment analysis mode that has been used to characterize MHC molecules, MHC-like molecules, and the non-immunological Dpr-DIP proteins. 

AIMS is a python package distributed in a notebook, CLI, and GUI format. An example of an application of AIMS can be seen in
this peer-reviewed article: https://elifesciences.org/articles/61393

When publishing analysis from this software, please cite:

Boughter CT, Borowska MT, Guthmiller JJ, Bendelac A, Wilson PC, Roux B, Adams EJ. Biochemical Patterns of Antibody Polyreactivity Revealed Through a Bioinformatics-Based Analysis of CDR Loops. eLife. 2020. DOI: 10.7554/eLife.61393

&

Boughter CT, Meier-Schellersheim M. An Integrated Approach to the Characterization of Immune Repertoires Using AIMS: An Automated Immune Molecule Separator. BioRxiv. 2022. DOI: 10.1101/2022.12.07.519510

# Documentation
Rather than have all of the instructions on this GitHub page, all information on installation and usage (and more!) has been moved to a separate, more readable documentation page. Please follow this link:

https://aims-doc.readthedocs.io/en/latest/

For the comprehensive AIMS user guide.

# Reproduction of Published Results
As of versions 0.8 and later, the data necessary for reproducing data published thus far have been moved to a separate repository. This repository can be found here:

https://github.com/ctboughter/AIMS_manuscripts

The underlying code remains the same, and will continue to be updated. This has been done to keep the AIMS analysis software more streamlined and less cluttered with manuscript-specific analysis.

# Further Reading
Now that AIMS has been out and in the wild for around two years, there have been additional published peer-reviewed manuscripts or posted preprints that highlight the capabilities of AIMS! I'll try to keep this list relatively up to date, and if it ever gets lengthy will likely move it to the ReadTheDocs page. Manuscripts thus far include:

- An application of AIMS to non-immune molecules using multi-sequence alignment (MSA) encoding: https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.2c02173
- The AIMS bible, with a thorough explanation of the rationale behind the AIMS analysis: https://www.biorxiv.org/content/10.1101/2022.12.07.519510v1
- An investigation of the nature of the germline interactions between TCR CDR loops and MHC: https://www.biorxiv.org/content/10.1101/2022.12.07.519507v1