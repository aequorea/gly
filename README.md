# gly

N-linked glycodesign software

<p align="center">
  <img src="Variety_of_glycans.svg.png" width="500">
  
<p align="center">  Varieties of glycosylation patterns --
  Image by Dna 621 - Own work, CC BY-SA 3.0, https://commons.wikimedia.org/w/index.php?curid=31825161
</p>
<p>
This software finds locations on a protein surface such that if we mutate the sites to have the N-linked glycosylation consensus sequence, it has a good chance of being efficiently expressed and glycosylated.
</p>
<p>
The algorithm implements a simple heuristic. If we want to mutate a protein residue, as long as the mutation doesn't disturb the molecule's hydrophobic core, and the mutation is to a hydrophilic residue, there is a good chance the mutation will be successful. The N-linked glycosylation "sequon" which is the attachment point for a glycan, involves pairs of hydrophilic residues in the mutant. So the algorithm looks for pairs of residues that have solvent exposed sidechains as potential mutation sites. Prolines are known to be "deal killers" for N-linked glycosylation, so prolines or other unknown residues in the neighborhood of a potential glycosylation site causes that mutation candidate to be excluded.
</p>
<p>
The program takes PDB files as inputs. There is a file with the complete structure, a file with solvent exposed atoms, a file with beta sheet residues and a file with alpha helix residues. I am currently using the pymol plugin findSurfaceResidues to get the solvent exposed atoms and pymol itself to get the secondary structure files.
</p>
<p>
This software has been tested against a number of data sets. For example 
[1cn4](aequorea/1cn4)
</p>
