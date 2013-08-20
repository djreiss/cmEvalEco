cmEvalEco
=========

Evaluate output of pyMonkey (in tsv format) on E. coli against RegulonDB.
This first comparison identifies how many RegulonDB motifs are recapitulated.

=========

To run, in R:

You will need the cMonkey package installed (and its dependencies), as well as having 

```
progs/meme
progs/mast
```

in the current directory (i.e. meme/mast in the "progs" directory; they can be symlinks to their actual locations).

Next, check out the "Defaults" in clusters.R; if they are not correct, then pre-set them in the R environment. Then:

```
source("clusters.R")
source("compare_regulondb_single.R")
```

the final number printed out will be the overlap.
A typical value for cMonkey (R) is about 50 (although for the example file, you will not get that value ;).

=========

I have added (gzipped) copies of a cMonkey tsv file and the full DISTILLER E. coli data set. You will want to gunzip them after cloning this repo.
