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

Select the type of distance-bias to use when sampling transfer recipients. The three options are: none, simple, exponential. Default: none

None disables distance-bias and samples transfer recipients uniformly at random.

Simple implements a distance-bias that scales proportionally to the inverse of the phylogenetic distance.

Exponential implements a distance-bias that scales exponentially in relation to the phylogenetic distance.

-dbr, --distance-bias-rate

When using the exponential distance-bias, set the rate parameter of the exponential distribution to be used. Default: 1.0

A higher value more strongly biases the selection of the transfer recipient towards one that is phylogenetically closer to the origin of the transfer.

-gb, --gene-birth-sampling

Randomly sample the location of gene birth on the species tree.

-gbc, --gene-birth-coefficient

Set level of bias towards the root of the species tree for gene tree birth location. Default: 1.0

A higher value more strongly biases the location of gene birth towards the root of the species tree.

### Example:

java -jar jprime-0.3.7.jar GuestTreeGen -rt 0.6 -db exponential -dbr 1.5 -gb -gbc 2.1 species.pruned.tree 0.2 0.1 0.3 gene

---

## DomainTreeGen

### Usage:

java -jar jprime-0.3.7.jar DomainTreeGen [options] <host tree> <guest trees> <dup rate> <loss rate> <trans rate> <out prefix>

### Required arguments:

<host tree file or string> <guest tree folder> <duplication rate> <loss rate> <transfer rate> <output prefix>

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

Enforce a minimum number of extant domain leaves per host leaf. Default: 0

-maxper, --max-leaves-per-host-leaf

Enforce a maximum number of extant domain leaves per host leaf. Default: 10

-nox, --no-auxiliary-tags

Exclude auxiliary PrIME tags in output trees.

-p, --leaf-sampling-probability

Set the probability of observing a leaf. Leaves that are unobserved will yield to their lineages being pruned away from the pruned tree. Default: 1.0

-a, --max-attempts

Set the maximum number of attempts that can be made to generate a tree that meets the requirements. Default: 10000

-vp, --vertex-prefix

Set the prefix for the vertex label. Default: D

-rt, --replacing-transfers

Set the probability of a horizontal domain transfer event being a replacing transfer. Default: 0.5

-db, --distance-bias

Select the type of distance-bias to use when sampling transfer recipients. The three options are: none, simple, exponential. Default: none

None disables distance-bias and samples transfer recipients uniformly at random.

Simple implements a distance-bias that scales proportionally to the inverse of the phylogenetic distance.

Exponential implements a distance-bias that scales exponentially in relation to the phylogenetic distance.

-dbr, --distance-bias-rate

When using the exponential distance-bias, set the rate parameter of the exponential distribution to be used. Default: 1.0

A higher value more strongly biases the selection of the transfer recipient towards one that is phylogenetically closer to the origin of the transfer.

-dbs, --domain-birth-sampling

Randomly sample the location of domain birth on the first gene tree.

-dbc, --domain-birth-coefficient

Set level of bias towards the root of the gene tree for domain tree birth location. Default: 1.0

A higher value more strongly biases the location of domain birth towards the root of the gene tree.

-ig, --inter-gene-transfers

Set the probability of a horizontal domain transfer to occur across gene trees, if more than one gene tree is provided. Default: 0.5

-is, --inter-species-transfers

Set the probability of a horizontal domain transfer to occur across species. Default: 0.5

-all, --all-gene-trees

Enforce the evolution of the domain tree within all given gene trees.

### Example:

java -jar jprime-0.3.7.jar DomainTreeGen -rt 0.6 -db exponential -dbr 1.5 -dbs -dbc 2.1 -ig 0.9 -is 0.7 -all species.pruned.tree genes-directory 0.2 0.1 0.4 domain

---
