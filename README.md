# BrownDwarf-SpectralIndices
Method to determine if ONE brown dwarf is a candidate to be a photometric variable or non-variable, by means of the study of spectral indices.

This first version works in Python and the only modifications to be made at the beginning of the code are:
  - Path and name of the file containing the brown dwarf spectrum.
      - The spectrum format must contain 3 columns with the names:
      'lambda', 'flux', 'eflux' and they must be separated with tab '\t'.
  - Define whether it is an L-type or T-type brown dwarf.
  - Whether or not to save a .pdf file with the index-index plot
Run code with: python LTdwarfIndices.py

Outputs:
- In how many graphs this object appears in the variable area
- whether it is a variable or non-variable candidate
- (optional) The index-index chart of the object

Attribution:
Please cite Oliveros-Gomez, N. et al. (2022) whenever results are used in a publication.

Future versions:
Decide to use the method for a single object or a list of objects.

Questions & feedback:
species has been developed and is maintained by Natalia Oliveros-Gomez (nl.oliverosgomez@ugto.mx). Feel free to send an email for questions, comments, or suggestions.

