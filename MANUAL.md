# An Improved Probabilistic Simulation Framework for Gene Family Evolution

**Author:** *Soumya Kundu*

## HostTreeGen:

### Usage:

java -jar jprime-0.3.7.jar HostTreeGen [options] <time interval> <birth rate> <death rate> <out prefix>

### Required Arguments:

<time interval> <birth rate> <death rate> <output prefix>

### Options:

-h, --help

Display the help menu with all of the required and optional arguments.

-q, --quiet

Write pruned tree directly to standard output and suppress creation of any auxiliary files.

-s, --seed

Specify a seed for the pseudorandom number generator.

-min, --min-leaves

Enforce a minimum number of extant leaves on the pruned tree. Default: 2

-max, --max-leaves

Enforce a maximum number of extant leaves on the pruned tree. Default: 64

-nox, --no-auxiliary-tags

Exclude auxiliary PrIME tags in output trees.

-p, --leaf-sampling-probability

Set the probability of observing a leaf. Leaves that are unobserved will yield to their lineages being pruned away from the pruned tree. Default: 1.0

-bi, --start-with-bifurcation

Force the simulation process to start with a bifurcation in the tree.

-a, --max-attempts

Set the maximum number of attempts that can be made to generate a tree that meets the requirements. Default: 10000

-vp, --vertex-prefix

Set the prefix for the vertex label. Default: H

### Example:

java -jar jprime-0.3.7.jar HostTreeGen -min 100 -max 100 1.0 5.0 0.05 species

---

## GuestTreeGen:

### Usage:

java -jar jprime-0.3.7.jar GuestTreeGen [options] <host tree> <dup rate> <loss rate> <trans rate> <out prefix>

### Required arguments:

<host tree file or string> <duplication rate> <loss rate> <transfer rate> <output prefix>

### Options:

-h, --help

Display the help menu with all of the required and optional arguments.

-q, --quiet

Write pruned tree directly to standard output and suppress creation of any auxiliary files.

-s, --seed

Specify a seed for the pseudorandom number generator.

-min, --min-leaves

Enforce a minimum number of extant leaves on the pruned tree. Default: 2

-max, --max-leaves

Enforce a maximum number of extant leaves on the pruned tree. Default: 64

-minper, --min-leaves-per-host-leaf

Enforce a minimum number of extant guest leaves per host leaf. Default: 0

-maxper, --max-leaves-per-host-leaf

Enforce a maximum number of extant guest leaves per host leaf. Default: 10

-nox, --no-auxiliary-tags

Exclude auxiliary PrIME tags in output trees.

-p, --leaf-sampling-probability

Set the probability of observing a leaf. Leaves that are unobserved will yield to their lineages being pruned away from the pruned tree. Default: 1.0

-a, --max-attempts

Set the maximum number of attempts that can be made to generate a tree that meets the requirements. Default: 10000

-vp, --vertex-prefix

Set the prefix for the vertex label. Default: G

-rt, --replacing-transfers

Set the probability of a horizontal gene transfer event being a replacing transfer. Default: 0.5

-db, --distance-bias

Select the type of 
