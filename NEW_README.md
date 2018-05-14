# Modified Simulation Framework for Accurate Modeling of Horizontal Gene Transfers

**Author:** *Soumya Kundu*

## Feature 1: Replacing Transfers

---

### Overview:

* The genphylodata simulation framework has been modified to support replacing transfers.

* The majority of the changes have been made to GuestTreeInHostTreeCreator.java, along with a few minor changes to GuestVertex.java.

---

### Method:

1. Every time a vertex is assigned a transfer event using a pseudorandom number generator, a new pseudorandom number generator is used to generate a random value that is used to choose between additive and replacing transfers.

2. The default probability of a replacing transfer is 0.5, but that number can be changed by the user.

3. If an additive transfer is chosen, then the simulation proceeds as usual.

4. If a replacing transfer is chosen, then the vertex with the replacing transfer event is added to a NavigableMap<Double, GuestVertex> where the key is the event time and the value is the vertex itself.

5. The NavigableMap is used to order the replacing transfer vertices in decreasing order of their event times.

6. After the rest of the guest tree is generated, for the first vertex in the NavigableMap, a contemporary guest tree lineage suitable to be the receipient of the transfer is sampled, such that the guest tree lineage cannot be a child of the replacing transfer lineage and that the guest tree lineage must have an event time lower than that for the replacing transfer lineage.

7. If such a guest tree lineage is found, then the event for that lineage is changed to a replacing loss, while its event time is changed to that of the replacing transfer.

8. After that, all of the other replacing transfers that have not been processed yet are checked to see if they are in the subtree of the replacing loss lineage; if any are found to be in its subtree, then those replacing transfers are removed as they are now in the subtree of a lost lineage.

9. Finally, the replacing loss lineage's children pointer is set to null and the subtree of the replacing transfer is generated using the birth-death process.

10. If no suitable recipient for the replacing transfer is found, then the replacing transfer's event is changed to additive transfer to better reflect the true nature of the event; this happens very rarely.

11. This process continues until there are no vertices left to process.

---

## Feature 2: Distance-biased Transfers

---

### Overview:

* The simulation framework has been modified to account for distance bias in horizontal gene transfers.

* Other than the original method for sampling eligible transfer recipients uniformly at random, two new methods that take into account the phylogenetic distance from the donor to the recipient have been incorporated.

* The majority of the changes have been made to RBTReeEpochDiscretiser.java and Epoch.java, along with a few minor changes to HostTreeGen.java, GuestTreeGenParameters.java, and GuestTreeInHostTreeCreator.java.

---

### Method:

1. Every time a transfer event is assigned, unless the distance bias model is turned off, in which case a transfer recipient is selected uniformly at random, one of the two distance bias models described below are used.

2. If the simple distance bias model is selected, then a NavigableMap of values corresponding to the inverse of the phylogenetic distances to all possible transfer recipients is populated, from which one receipient is randomly chosen, where the expectation of choosing any recipient is directly correlated to the inverse of its phylogentic distance from the donor.

3. If the exponential distance bias model is selected, then a NavigableMap of values corresponding to 10 raised to the power of the negative of the phylogenetic distances to all possible transfer recipients is populated, from which one receipient is randomly chosen, where the expectation of choosing any recipient is directly negatively correlated to the log of its phylogentic distance from the donor.

---

## Feature 3: Sampling Location of Gene Birth on Species Tree

---

### Overview:

* The gene tree simulation process has been modified to sample with biased probabilities a location on the species tree from where the gene tree will start to evolve.

* Previously, the gene tree was always evolved from the root of the species tree.

* We bias the selection gene birth locus towards the top of the species tree with a non-zero probability of not picking the species tree root.

---

### Method:

1. For k vertices in the species tree, construct a NavigableMap of values corresponding to Euler's number raised to the negative value of the index of the vertex divided by a constant such as 2 or 5.

2. The keys are the running total of the values and the values for the NavigableMap are the vertex numbers between 0 and k.

3. Choose a random value between 0 and the highest value in the NavigableMap and pick the vertex with a value immediately higher than the sampled value.

4. Return k - 1 - (sampled vertex) to return the sampled root, as the vertices of the tree are numbered in post-order with k being the value of the root and 0 being one of the leaves.

---
