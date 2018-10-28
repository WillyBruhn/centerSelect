# centerSelect
Selects a box-shaped isosurface around the active center.
//![alt text](https://github.com/WillyBruhn/centerSelect/blob/master/exampleImages/img.png)

<img style="float: right;" src="https://github.com/WillyBruhn/centerSelect/blob/master/exampleImages/img.png">

Call:
```bash
./centerSelect.sh _parametersFile_
```
where _parametersFile_ specifies the file with the parameters. If you call 
```bash
./centerSelect.sh
```
without any parameter a file with the name "parameters.txt" in the same folder is read in.

# ParametersFile
The _parametersFile_ contains the following parameters:

| Parameter         | Value                                                      | Description                                                                                           |
|-------------------|------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| path2centerSelect | /home/willy/RedoxChallenges/centerSelect/                  | path to this repo on the local machine                                                                |
| path2files        | /home/willy/RedoxChallenges/centerSelect/Redox_old/Output/ | path to the Output-folder where all folders with the proteins are stored                              |
| boxSize           | 30                                                         | size of the square that is to be drawn around the active centre  parallel to the X- and Y-axis        |
| depth             | 10                                                         | the box is drawn along the Z-axis.                                                                    |
| eps               | 0.3                                                        | How much difference to the actual values -1.0 and 1.0 is tolerated.                                   |
| onlyDoMissing     | true                                                       | only select the active centres for the folders that do not already contain a file "activeCenter.csv". |

# Folder Structure

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


(cool editor -> https://dillinger.io/)

(table -> https://www.tablesgenerator.com/markdown_tables)
