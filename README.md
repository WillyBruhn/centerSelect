# centerSelect
Selects a box-shaped isosurface around the active center.

<img style="float: right;" src="https://github.com/WillyBruhn/centerSelect/blob/master/exampleImages/img.png">

The active center is approximated as a single point obtained as the geometric center of all the atoms from the active center. 

## 1. Finding the position of the active center
In the context of this implementation a motife is a vector consisting of multiple aminoacids e.g. *CYS GLY ALA*. A set of motives is specified by the user. Then the pqr-file is scanned for any potential motifes. A list of found motives is returned. When multiple motifes are found the first occurence is taken. Then the atom-coordinates from the found motive are extracted. The geometrical center is then calculated. This *geometric center* is used as the approximation for the active center.

## 2. Iso-surface in close proximity to the active center
Given a single point (*p.x*,*p.y*,*p.z*) in 3D representing the active center a subset of the iso-surface is calculated as follows:
Given a parameter *boxSize* and *depth* a square centered in *p* with height *boxSize*/2 and width *boxSize*/2 can be obtained. Then extruding this square along the z-axis from the front-side of the box to a position at *p.z* -*depth* one obtains a cuboid.
Then all points from the iso-surface that lie within this box are used for further investigation.

For our paper we used as parameters *boxSize* = 30 and *depth* = 15 angstrom.



Call:
```bash
./centerSelect.sh _parametersFile_
```
where _parametersFile_ specifies the file with the parameters. If you call 
```bash
./centerSelect.sh
```
without any parameter a file with the name "parameters.txt" in the same folder is read in.

## ParametersFile
The _parametersFile_ contains the following parameters:

| Parameter         | Value                                                      | Description                                                                                           |
|-------------------|------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| path2centerSelect | /home/willy/RedoxChallenges/centerSelect/                  | path to this repo on the local machine                                                                |
| path2files        | /home/willy/RedoxChallenges/centerSelect/Redox_old/Output/ | path to the Output-folder where all folders with the proteins are stored                              |
| boxSize           | 30                                                         | size of the square that is to be drawn around the active centre  parallel to the X- and Y-axis        |
| depth             | 10                                                         | the box is drawn along the Z-axis.                                                                    |
| eps               | 0.3                                                        | How much difference to the actual values -1.0 and 1.0 is tolerated.                                   |
| onlyDoMissing     | true                                                       | only select the active centres for the folders that do not already contain a file "activeCenter.csv". |
| SELECTBOX     | false                                                       | do you want to select a subset of points that are in the box? the pts-file will then contain these points. |

## Folder Structure

```
└─── Redox_old
│   │   file011.txt
│   │   file012.txt
│   │
│   └─── Input
│       │   └───  pdb
│       │         │   001.pdb
│       │         │   002.pdb
│       │         │   ...
│       │   └───  pqr
│       │         │   001.pqr
│       │         │   002.pqr
│       │         │   ...
│   └─── Output
│       │   └───  001
│       │         │   001.pqr
│       │         │   001_el.tga
│       │         │   001_es.tga
│       │         │   001_ss.tga
│       │         │   001_neg.pts
│       │         │   001_pos.pts
│       │         │   001_pot.dx
│       │         │   apbs.in
│       │         │   ...
│       │   └───  002
│       │         │   002.pqr
│       │         │   002_el.tga
│       │         │   ...
```
In each folder in Output **_centerSelecter.R_** creates a file **_active_center.csv_** that specifies the number of atoms that are within the active centre. The files **001_neg.pts** and **001_pos.pts** contain the points of the iso-surface around the active centre that can be used in the next step to compare the similarity of the iso-surfaces.

## AA_utility
### Active Center Motifes
The file **AA_utility/motifs3letterCode.txt** contains the motifes of the active centers that the script will look for in the pdb-files.

| CYS | ALA | HIS |
|-----|-----|-----|
| CYS | GLY | ALA |
| CYS | ASN | HIS |
| ... | ... | ... |

### AAs.csv
**AAs.csv** specifies the atom-chains that make up the aminoacids.

| ALA | ARG | ... |
|-----|-----|-----|
| N   | N   | ... |
| CA  | CA  | ... |
| C   | C   | ... |
| O   | O   | ... |
| CB  | CB  | ... |
| ... | ... | ... |

### AA_lengths.csv
**AA_lengths.csv** specifies the lenghts of the atom-chains of the amino-acids.

| ALA | ARG | ASN | ASP | ... |
|-----|-----|-----|-----|-----|
| 10  | 24  | 14  | 14  | ... |


## Missing motifes
The file **missing.txt** contains the names of the proteins for which the script did not find an acitve center.


(cool editor -> https://dillinger.io/)

(table -> https://www.tablesgenerator.com/markdown_tables)
