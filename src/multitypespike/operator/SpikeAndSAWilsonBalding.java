package multitypespike.operator;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import sa.evolution.operators.SAWilsonBalding;


public class SpikeAndSAWilsonBalding extends SAWilsonBalding {

    final public Input<RealParameter> spikesInput = new Input<>("spikes", "spikes associated with each branch on the tree",
            Input.Validate.REQUIRED);

    final public Input<RealParameter> spikeShapeInput = new Input<>("spikeShape", "shape parameter for the " +
            "gamma distribution of the spikes", Input.Validate.REQUIRED);



    // Needs testing
    /**
     * Implements spike-aware prune and regraft moves for sampled ancestor trees
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        Tree tree = treeInput.get();
        tree.startEditing(this);


        double oldMinAge, newMinAge, newRange, oldRange, newAge, fHastingsRatio, DimensionCoefficient;
        int newDimension, oldDimension;

        double pRev = 0;
        double pFwd = 0;

        //STEP 1: Select a subtree root 'i', avoiding root or sampled ancestor tips
        int nodeCount = tree.getNodeCount();
        Node i;
        do {
            i = tree.getNode(Randomizer.nextInt(nodeCount));
        } while (i.isRoot() || i.isDirectAncestor());

        if(i.isDirectAncestor()) System.out.println("i is SA");
        // Parent of subtree root
        Node iP = i.getParent();

        // Sibling of subtree root
        Node CiP = (iP.getLeft().getNr() == i.getNr()) ? iP.getRight() : iP.getLeft();


        // Make sure that there is at least one candidate edge to attach node iP to
        if (iP.getParent() == null && CiP.getHeight() <= i.getHeight()) {
            return Double.NEGATIVE_INFINITY;
        }

        //STEP 2: Select a random target node 'j' (an internal node or a leaf)
        Node j;
        Node jP;

        final int leafNodeCount = tree.getLeafNodeCount();

        if (leafNodeCount != tree.getExternalNodes().size()) {
            System.out.println("Node counts are incorrect: NodeCount = " + nodeCount + " leafNodeCount = " + leafNodeCount + " external node count = " + tree.getExternalNodes().size());
        }

        // Make sure that the target branch <jP, j> or target leaf j is above the subtree being moved

        int nodeNumber;
        double newParentHeight;
        boolean attachingToLeaf;
        boolean adjacentEdge;
        do {
            adjacentEdge = false;
            nodeNumber = Randomizer.nextInt(nodeCount + leafNodeCount);
            if (nodeNumber < nodeCount) {
                j = tree.getNode(nodeNumber);
                jP = j.getParent();
                if (jP != null)
                    newParentHeight = jP.getHeight();
                else newParentHeight = Double.POSITIVE_INFINITY;
                if (!CiP.isDirectAncestor())
                    adjacentEdge = (CiP.getNr() == j.getNr() || iP.getNr() == j.getNr());
                attachingToLeaf = false;
            } else {
                j = tree.getExternalNodes().get(nodeNumber - nodeCount);
                jP = j.getParent();
                newParentHeight = j.getHeight();
                attachingToLeaf = true;
            }
        } while (j.isDirectAncestor() || (newParentHeight <= i.getHeight()) || (i.getNr() == j.getNr()) || adjacentEdge /*|| adjacentLeaf */);

        // Reject if 'j' is parent of 'i'
        if (attachingToLeaf && iP.getNr() == j.getNr()) {
            System.out.println("Proposal failed because j = iP");
            return Double.NEGATIVE_INFINITY;
        }

        // Reject if move creates a cycle
        if (jP != null && jP.getNr() == i.getNr()) {
            System.out.println("Proposal failed because jP = i. Heights of i = " + i.getHeight() + " Height of jP = " + jP.getHeight());
            return Double.NEGATIVE_INFINITY;
        }

        oldDimension = nodeCount - tree.getDirectAncestorNodeCount() - 1;


        // Hastings numerator calculation + newAge of iP
        // Propose a new age for node iP based on where it's being attached

        if (attachingToLeaf) {
            // Attach to leaf 'j'
            newRange = 1;
            newAge = j.getHeight();

        } else {

            if (jP != null) {
                // Attach to edge between 'jP' and 'j'
                newMinAge = Math.max(i.getHeight(), j.getHeight());
                newRange = jP.getHeight() - newMinAge;
                newAge = newMinAge + (Randomizer.nextDouble() * newRange);

            } else {
                // Attach as new root
                double randomNumberFromExponential;
                randomNumberFromExponential = Randomizer.nextExponential(1);
                newRange = Math.exp(randomNumberFromExponential);
                newAge = j.getHeight() + randomNumberFromExponential;
            }
        }


        Node PiP = iP.getParent(); // grandparent of 'i'

        // Hastings denominator calculation
        // Compute old range for reverse move

        if (CiP.isDirectAncestor()) {
            oldRange = 1;
        }
        else {
            oldMinAge = Math.max(i.getHeight(), CiP.getHeight());
            if (PiP != null) {
                oldRange = PiP.getHeight() - oldMinAge;
            } else {
                oldRange = Math.exp(iP.getHeight() - oldMinAge);
            }
        }


        // Handle spikes

        // Detach subtree from sampled ancestor parent (fake) node
        if(CiP.isDirectAncestor()){
            double spikeShape = spikeShapeInput.get().getValue();
            GammaDistribution gamma = new GammaDistributionImpl(spikeShape, 1/spikeShape);
            RealParameter spikes = spikesInput.get();

            // Reattach subtree to a random non-sampled ancestor leaf
            if (attachingToLeaf){ // a.2 –> b.2 (see SA paper)
                // spike(j)>0, spike(CiP)=0 –> spike(j)=0, spike(CiP)>0

                double sjOld = spikes.getValue(j.getNr());
                double sCiPOld = spikes.getValue(CiP.getNr());

                // (sjNew = 0)
                double sCiPNew = Randomizer.nextGamma(spikeShape, 1/spikeShape);

                // Check
                if(sjOld == 0) {
                    System.out.println("Wilson-Balding reject: a.2 –> b.2 move: spike of j is zero");
                    return Double.NEGATIVE_INFINITY; // Reject
                }
                if (sCiPOld != 0) {
                    System.out.println("Wilson-Balding reject: a.2 –> b.2 move: spike of CiP is non-zero. sCiP = " + sCiPOld);
                    return Double.NEGATIVE_INFINITY; // Reject
                }

                pRev = gamma.logDensity(sjOld);
                pFwd = gamma.logDensity(sCiPNew);

                // Update spikes
                spikes.setValue(j.getNr(), 0.0);
                spikes.setValue(CiP.getNr(), sCiPNew);


            // Reattach subtree to origin branch, parent of 'i' becomes the new root
            } else if (j.isRoot()) { // a.2 –> b.1 (root)
                // spike(j)=0, spike(CiP)=0, spike(iP)>0 –> spike(j)>0, spike(CiP)>0, spike(iP)=0

                double sjOld = spikes.getValue(j.getNr());
                double sCiPOld = spikes.getValue(CiP.getNr());
                double siPOld = spikes.getValue(iP.getNr());

                double sjNew = Randomizer.nextGamma(spikeShape, 1/spikeShape);
                double sCiPNew = Randomizer.nextGamma(spikeShape, 1/spikeShape);
                // (siPNew = 0)

                // Check
                if(sjOld != 0) {
                    System.out.println("Wilson-Balding reject: a.2 –> b.1 (root) move: spike of j is non-zero. sj = " + sjOld);
                    return Double.NEGATIVE_INFINITY; // Reject
                }
                if (sCiPOld != 0) {
                    System.out.println("Wilson-Balding reject: a.2 –> b.1 (root) move: spike of CiP is non-zero. sCiP = " + sCiPOld);
                    return Double.NEGATIVE_INFINITY; // Reject
                }
                if(siPOld == 0) {
                    System.out.println("Wilson-Balding reject: a.2 –> b.1 (root) move: spike of iP is zero");
                    return Double.NEGATIVE_INFINITY; // Reject
                }

                pRev = gamma.logDensity(siPOld);
                pFwd = gamma.logDensity(sjNew) + gamma.logDensity(sCiPNew);

                // Update spikes
                spikes.setValue(j.getNr(), sjNew);
                spikes.setValue(CiP.getNr(), sCiPNew);
                spikes.setValue(iP.getNr(), 0.0);

                // Reattach subtree to the middle of a normal branch
            } else { // a.2 –> b.1
                // spike(CiP)=0 –> spike(CiP)>0

                double sCiPOld = spikes.getValue(CiP.getNr());
                double sCiPNew = Randomizer.nextGamma(spikeShape, 1/spikeShape);

                // Check
                if (sCiPOld != 0) {
                    System.out.println("Wilson-Balding reject: a.2 –> b.1 move: spike of CiP is non-zero. sCiP = " + sCiPOld);
                    return Double.NEGATIVE_INFINITY; // Reject
                }

                pRev = 0;
                pFwd = gamma.logDensity(sCiPNew);

                // Update spikes
                spikes.setValue(CiP.getNr(), sCiPNew);

            }
        }

        // Detach subtree from a non-fake node
        if(!CiP.isDirectAncestor()) {
            double spikeShape = spikeShapeInput.get().getValue();
            GammaDistribution gamma = new GammaDistributionImpl(spikeShape, 1 / spikeShape);
            RealParameter spikes = spikesInput.get();

            // Reattach subtree to a random non-sampled ancestor leaf
            if (attachingToLeaf) {

                // Detach subtree from the root/origin branch
                if (iP.isRoot()) { // a.1 (root) -> b.2
                    // spike(j)>0, spike(iP)=0 –> spike(j)=0, spike(iP)>0

                    double sjOld = spikes.getValue(j.getNr());
                    double siPOld = spikes.getValue(iP.getNr());

                    // (sjNew = 0)
                    double siPNew = Randomizer.nextGamma(spikeShape, 1/spikeShape);

                    // Check
                    if(sjOld == 0) {
                        System.out.println("Wilson-Balding reject: a.1 (root) -> b.2 move: spike of j is zero");
                        return Double.NEGATIVE_INFINITY; // Reject
                    }
                    if(siPOld != 0) {
                        System.out.println("Wilson-Balding reject: a.1 (root) -> b.2 move: spike of iP is non-zero. siP = " + siPOld);
                        return Double.NEGATIVE_INFINITY; // Reject
                    }

                    pRev = gamma.logDensity(sjOld);
                    pFwd = gamma.logDensity(siPNew);

                    // Update spikes
                    spikes.setValue(j.getNr(), 0.0);
                    spikes.setValue(iP.getNr(), siPNew);

                    // Detach subtree from a normal internal node
                } else { // a.1 –> b.2
                    // spike(j)>0 –> spike(j)=0

                    double sjOld = spikes.getValue(j.getNr());
                    // (sjNew = 0)

                    // Check
                    if(sjOld == 0) {
                        System.out.println("Wilson-Balding reject: a.1 (root) -> b.2 move: spike of j is zero");
                        return Double.NEGATIVE_INFINITY; // Reject
                    }

                    pRev = gamma.logDensity(sjOld);
                    pFwd = 0;

                    // Update spikes
                    spikes.setValue(j.getNr(), 0.0);

                }

            // Reattach subtree to origin branch, parent of 'i' becomes the new root
            } else if (j.isRoot()) { // a.1 –> b.1 (root)
//                 spike(j)=0, spike(iP)>0 –> spike(j)>0, spike(iP)=0

                double sjOld = spikes.getValue(j.getNr());
                double siPOld = spikes.getValue(iP.getNr());

                double sjNew = Randomizer.nextGamma(spikeShape, 1/spikeShape);
                // (siPNew = 0)

                // Check
                if(sjOld != 0) {
                    System.out.println("Wilson-Balding reject: a.1 –> b.1 (root) move: spike of j is non-zero. sj = " + sjOld);
                    return Double.NEGATIVE_INFINITY; // Reject
                }
                if(siPOld == 0) {
                    System.out.println("Wilson-Balding reject: a.1 –> b.1 (root) move: spike of iP is zero");
                    return Double.NEGATIVE_INFINITY; // Reject
                }

                pRev = gamma.logDensity(siPOld);
                pFwd = gamma.logDensity(sjNew);

                // Update spikes
                spikes.setValue(j.getNr(), sjNew);
                spikes.setValue(iP.getNr(), 0.0);

            }
            // else a.1 –> b.1: no spike adjustment required
        }


        // STEP 5: Update the topology
        if (iP.getNr() != j.getNr() && CiP.getNr() != j.getNr()) {

            // Detach iP from its old parent, reattach CiP
            iP.removeChild(CiP); //remove <iP, CiP>

            if (PiP != null) {
                PiP.removeChild(iP);   // remove <PiP,iP>
                PiP.addChild(CiP);   // add <PiP, CiP>
                PiP.makeDirty(Tree.IS_FILTHY);
                CiP.makeDirty(Tree.IS_FILTHY);
            } else {
                CiP.setParent(null); // completely remove <iP, CiP>
                tree.setRootOnly(CiP); // new root
            }

            // Attach iP to target location
            if (jP != null) {
                jP.removeChild(j);  // remove <jP, j>
                jP.addChild(iP);   // add <jP, iP>
                jP.makeDirty(Tree.IS_FILTHY);
            } else {
                iP.setParent(null); // completely remove <PiP, iP>
                tree.setRootOnly(iP); // new root
            }

            iP.addChild(j);  // graft subtree
            iP.makeDirty(Tree.IS_FILTHY);
            j.makeDirty(Tree.IS_FILTHY);
        }

        iP.setHeight(newAge); // update height for iP

//        System.out.println(count + "CiP is not a SA: " + !CiP.isDirectAncestor());


        //make sure that either there are no direct ancestors or r<1
        if ((rInput.get() != null) && (tree.getDirectAncestorNodeCount() > 0 && rInput.get().getValue() == 1))  {
            return Double.NEGATIVE_INFINITY;
        }

        // STEP 6: Compute Hastings ratio
        newDimension = nodeCount - tree.getDirectAncestorNodeCount() - 1;
        DimensionCoefficient = (double) oldDimension / newDimension;

        fHastingsRatio = Math.abs(DimensionCoefficient * newRange / oldRange);

        return Math.log(fHastingsRatio) + pRev - pFwd;
    }

}
