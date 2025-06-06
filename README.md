# Molecule-Database-Generator-and-Viewer

A Python-based GUI for parsing SDF files, calculating molecular properties, storing data in SQLite, visualizing structures, and filtering/exporting compounds. Developed as coursework for MSc Bioinformatics 2024/25, University of Glasgow.


## Overview

This script manages complex chemical compound datasets by automating data extraction, calculation of molecular descriptors, and interactive visualization. 
The software performs the following key functions:

### Data Extraction and Processing

* **Parsing SDF files:** Reads molecular information from Structure Data Format (SDF) files, a widely used format for storing molecular structures and metadata.
* **Extraction of molecular data:** Retrieves compound identifiers (e.g., compound names, IDs), atom and bond counts, and other structural details.
* **Calculation of molecular descriptors:** Computes important chemical properties such as:

  * Molecular Weight (MW)
  * LogP (lipophilicity)
  * Number of Hydrogen Bond Donors (HBD)
  * Number of Hydrogen Bond Acceptors (HBA)
  * Number of Rotatable Bonds
  * Topological Polar Surface Area (TPSA)

These descriptors help in assessing the drug-likeness and physicochemical properties of compounds.

### Data Storage

* **SQLite database integration:** All extracted and computed data is stored in a local SQLite database. This ensures:

  * Efficient data retrieval and querying
  * Persistent storage across sessions
  * Structured relational organization linking compounds, atoms, bonds, and calculated descriptors

### Graphical User Interface (GUI)

* **Interactive browsing:** Users can load datasets and view all compounds in a sortable, filterable table.
* **Dynamic filtering:** Apply filters on molecular properties (e.g., molecular weight ranges, LogP values, H-bond counts) to narrow down compounds of interest.
* **2D Molecular Visualization:** Upon selecting a compound, its 2D chemical structure is rendered in the GUI, providing intuitive insight into its molecular makeup.
* **Export Functionality:** Filtered or selected data can be exported as CSV files for downstream analysis or sharing.
* **External Links:** Direct integration with external databases such as PubChem allows users to quickly access additional compound information or references by redirecting to PubChem webpages based on compound IDs.

