# Amino Acid Network

A Python Script to generate and analyze an Amino Acid Network.

## Getting Started

To clone this repository, type

```
git clone https://github.com/lorpac/amino_acid_network.git
```

### Prerequisites

The required packages are listed in the file [requirements.txt](https://github.com/lorpac/amino_acid_network/blob/master/requirements.txt). To install the requirements with [pip](https://pypi.org/project/pip/), type:

```
pip install -r requirements.txt
```

You also need Rodigo Gilardi's class [biographs](https://github.com/rodogi/biographs). Clone it and do not forget to add the class folder (`biograps/biographs`) to you Python path. Alternatively, if you are using a pip virtual environment `env`, you can manually add it to your packages by saving a copy in the directory `env/Lib/site-packages`. To do so, move to the site-packages folder of your virtual environment or of your Python system install ad type
```
git clone https://github.com/rodogi/biographs .tmp;
mv .tmp/biographs .;
rm .tmp -f -r
```
If you are not familiar with virtual environments, [here](https://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/26/python-virtual-env/) you will find the detailed instructions on how to create and manage virtual environments with `pip`.

Finally, you need to have LaTex installed on your system in order to produce the plots with Matplotlib. On Ubuntu, you can install LaTex and the necessary extensions by running

```
sudo apt-get install dvipng texlive-latex-base texlive-latex-extra texlive-fonts-recommended
```
and if you still get an error, try installing `cm-super`:
```
sudo apt-get install cm-super
```

For other operating systems, or if you encounter problems, please follow the instructions in Matplotlib's [tutorial](https://matplotlib.org/3.1.0/tutorials/text/usetex.html).


## Set-up the analysis

To set up the parameters of you analysis, you need to create a configuration file following the template of [analysis_config_template.json](https://github.com/lorpac/amino_acid_network/blob/master/analysis_config_template.json) and save it as `analysis_config.json`. In the configuration file, the following parameters must be given (in any order):

- `name`: a name of your analysis
- `pdb_id`: either a `.pdb` or `.ent` file name found in `data/`, or a valid PDB identifier. Note: in the second casem a internet connection will be needed to download the required pdb structure from the [RCSB database](http://www.rcsb.org/).
- `cutoff`: cutoff distance to determine the connectivity between amino acids
- `dim`: type of links to consider. Valid options are: `all` (or an empty string, all links), `1D` (between first neighbors in the sequence), `2D` (between amino acids belonging to the same secondary structure element), `3D` (between amino acids belonging to the same chain but to the different secondary structure elements), `4D` (between amino acids belonging to different chains), `3-4D` (3D and 4D), `1-2D` (3D and 4D).
- `select positions`: boolean.
  - If `select positions` is true: `start` and `stop` are the sequence positions to start and stop the analysis
- `variant`: boolean. If `true`, the Perturbation Network is calculated
  - If `variant` is true, `threshold` is the difference threshold value to add links in the Peturbation Network and `reference folder` is the folder path where the WT Amino Acid Network is saved.
  - `single_mutation`: boolean.
    - If  `single_mutation` is true, `source` is the mutation position in the sequence.
- `subnet_highnw`: boolean. If true, the subnetwork of nodes with high Nw value are drawn.
- `subnet_highwij`: boolean. If true, the subnetwork of links with high wij value are drawn.
- `draw_neighborhoods`: boolean. If true, draw the ego networks of each node.
- `calculate_cliques`: boolean. If true, the cliques of the Amino Acid Network are calculated.
  - `draw_cliques`: boolean. If true, the cliques of the Amino Acid Network are drawn.
- `remove_hydrogen_atoms`: boolean. If true, pre-process the PDB file to remove the hydrogen atoms.

## Run the analysis

To run the analysis, type:

```
python analysis.py
```
You will find the results in the `results/` folder.
## Authors

Lorenza Pacini - [lorpac](https://github.com/lorpac)

## How to cite

If you use this code in you work, please cite:

- Dorantes-Gilardi, R., Bourgeat, L., Pacini, L., Vuillon, L., & Lesieur, C. (2018). In proteins, the structural responses of a position to mutation rely on the Goldilocks principle: not too many links, not too few. Physical Chemistry Chemical Physics, 20(39), 25399-25410.
  
## Licence

The Building network source code is available under the [CeCILL](http://cecill.info/) licence. Please see `LICENCE.txt` for details.

