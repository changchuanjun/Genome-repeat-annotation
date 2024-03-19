
For repeat annotation, This script implemented a comprehensive strategy to identify repetitive sequences within the genome. In the first round, the script utilized RepeatMasker (https://www.repeatmasker.org/RepeatMasker/) (Version 4.1.5) to mask the whole-genome sequences, using a combined database from RepBase[1] and Dfam[2]. In the second round, the script employed RepeatModeler(https://www.repeatmasker.org/RepeatModeler/) (Version 2.0.4) to generate a de novo repeat library. Then, the script re-masked the genome sequences that were masked in the first round with this new library using RepeatMasker. The final set of repetitive sequences was obtained by integrating the results from both rounds.

[1] Bao, W., Kojima, K.K. & Kohany, O. Repbase Update, a database of repetitive elements in eukaryotic genomes. Mobile DNA 6, 11 (2015). https://doi.org/10.1186/s13100-015-0041-9.

[2] Robert Hubley, Robert D. Finn, Jody Clements, Sean R. Eddy, Thomas A. Jones, Weidong Bao, Arian F.A. Smit, Travis J. Wheeler, The Dfam database of repetitive DNA families, Nucleic Acids Research, Volume 44, Issue D1, 4 January 2016, Pages D81–D89, https://doi.org/10.1093/nar/gkv1272.

Acknowledgments
This script was initially supported by genek(https://genek.cn/). I gratefully acknowledge X.D. Zhang for creating and integrating the pipeline. I have reorganized and modified it myself.
