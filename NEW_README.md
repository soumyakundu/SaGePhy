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

2. The default probability of a replacing transfer is 0.7, but that number can be changed by the user.

3. If an additive transfer is chosen, then the simulation proceeds as usual.

4. If a replacing transfer is chosen, then after performing all actions performed for the case of an additive transfer, the vertex with the replacing transfer event is added to a NavigableMap<Double, GuestVertex> where the key is the event time and the value is the vertex itself.

5. The NavigableMap is used to order the replacing transfer vertices by their event times.

6. After the guest tree is generated, the vertices in the NavigableMap are added to an ArrayList in decreasing order of event times in order to ensure that a replacing transfer that is processed later does not replace a previously processed replacing transfer lineage or its ancestor.

7. For each vertex in the ArrayList, a guest tree lineage in the host tree edge to which the transfer has occured is searched, such that the guest tree lineage cannot be a child of the replacing transfer lineage and that the guest tree lineage must have an event time lower than that for the replacing transfer lineage.

8. If such a guest tree lineage is found, then the event for that lineage is changed to a replacing loss, while its event time is changed to that of the replacing transfer.

9. After that, all of the other replacing transfers that have not been processed yet are checked to see if they are in the subtree of the replacing loss lineage; if any are found to be in its subtree, then those replacing transfers are added to an invalidated ArrayList that every processed replacing transfer is checked against to make sure that no replacing losses are invoked for a replacing transfer that could not have happened since it is now in the subtree of a lost lineage.

10. Finally, the replacing loss lineage's children pointer is set to null and the processed replacing transfer is removed from the ArrayList.

11. If a suitable lineage is not found for replacement, then the replacing transfer's event is changed to additive transfer to better reflect the true nature of the event.

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

3. If the exponential distance bias model is selected, then a NavigableMap of values corresponding to 10 raised to the power of the inverse of the phylogenetic distances to all possible transfer recipients is populated, from which one receipient is randomly chosen, where the expectation of choosing any recipient is directly correlated to the log of the inverse of its phylogentic distance from the donor.

---

