Modified Simulation Framework for Accurate Modeling of Horizontal Gene Transfers

Author: Soumya Kundu

Feature 1: Replacing Transfers

-------------------------------

Overview:

The genphylodata simulation framework has been modified to support replacing transfers.

The majority of the changes have been made to GuestTreeInHostTreeCreator.java, along with a few minor changes to GuestVertex.java.

-------------------------------------------------------------------------------

Method:

Every time a vertex is assigned a transfer event using a pseudorandom number generator, the generated value is used again to choose between additive and replacing transfers.

If an additive transfer is chosen, then the simulation proceeds as usual.

If a replacing transfer is chosen, then first, after the vertex and its sibling are added to the LinkedList "alive" of lineages being currently processed, a target host tree edge is randomly chosen from branches that are eligible to receive transfers from the current lineage.

Once this recipient branch is chosen, we search for any guest tree lineage that is being processed (currently in the LinkedList "alive") that is present inside the chosen host tree branch; the first one found that has not already been assigned a loss event is re-assigned to be a loss.

If none are found, then a guest tree lineage target is searched for again when the replacing transfer vertex is popped from the LinkedList.

If none are found again, then the integer corresponding to the target host tree branch is added to an ArrayList "replacements" that every new addition to the LinkedList "alive" is checked against to find any new guest tree vertex assigned to the target host tree branch.

If a match is found where the guest tree vertex event is not already a loss, then the guest tree vertex event is changed to a loss and the integer is removed from the ArrayList "replacements".

The replacing guest tree vertex is marked with a true value for the boolean "doNotReplace" to prevent it from replacing itself.

-------------------------------------------------------------------------------

