# BrownDwarf-SpectralIndices
Method to determine if ONE brown dwarf is a candidate to be a photometric variable or non-variable, by means of the study of spectral indices.

We have two files:
- To analyze a single Brown Dwarf (modify the file name BD spectrum) [LTdwarfIndices_individualBD.py]
- To analyse a set sample indices of many Brown Dwarfs (modify the file name BD indices list) [LTdwarfIndices_Set-indices.py]

This method works in Python and the only modifications to be made at the beginning of the code are:
  - Path and name of the file containing the brown dwarf spectrum or BD indices list
      - The spectrum format must contain 3 columns with the names:
      'lambda', 'flux', 'eflux' and they must be separated with tab '\t'.
      - The indices list format must contain the number columns of the quantity of indices, depends on if L or T type, with the names like mention here:
        - L type: index_mostH, index_mostJ, index_less, index_Jcurve, index_H2OJ, index_CH4J
        - T type: index_H, index_J, index_HJ, index_J_H, index_Jslope, index_Jcurve
  - Define whether it is an L-type or T-type brown dwarf.
  - Whether or not to save a .pdf file with the index-index plot
  - Whether or not to save a .pdf file with the histogram plot (only for sample case)
Run code with: python LTdwarfIndices.py

Outputs:
- In how many graphs this object appears in the variable area
- whether it is a variable or non-variable candidate
- (optional) The index-index plot 
- (optional) The histogram plot
  
Attribution:
Please cite Oliveros-Gomez, N. et al. (2022) whenever results are used in a publication.

Future versions:
Decide to use the method for a single object or a list of objects in the same file.

Questions & feedback:
LTdwarfIndices has been developed and is maintained by Natalia Oliveros-Gomez (nl.oliverosgomez@ugto.mx). Feel free to send an email for questions, comments, or suggestions.

