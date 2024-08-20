# Integrating Spatial Transcriptomics for Intercellular Communication with NicheNet
A case study on the Kupffer cell niche of the liver

Master's thesis project at UGent 2023-2024.

Yunseol Park

Supervisors:
- Chananchida Sang-aram
- Ruth Seurinck
- Yvan Saeys

Files:
- `spatialCCCtools`: contains files for DeepCOLOR and NicheCompass analyses
  - `DeepCOLOR_Resolve.ipynb`: the main file for DeepCOLOR analysis using Resolve data (subcellular resolution)
  - `final_coexp_500epochs.csv`: the result from the above DeepCOLOR analysis
  - `DeepCOLOR_Visium.ipynb`: DeepCOLOR analysis usning Visium data (spot resolution)
  - `DeepCOLOR_colocMatrix.ipynb`: analysis on DeepCOLOR's colocalization matrix
  - `NicheCompass.ipynb`: 
- `spatialNicheNet`: contains files for spatially informed NicheNet analysis
  - `giotto_proximity_enrichment`: Giotto proximity enrichment results and visualizations
  - `spatialNicheNet`: the main spatial NicheNet analysis
  - `interaction_input.csv`: CellPhoneDB database for ligand-receptor pairs and their communication type
  - All R scripts, files with extension `.R`, contain functions necessary for the R markdown files
