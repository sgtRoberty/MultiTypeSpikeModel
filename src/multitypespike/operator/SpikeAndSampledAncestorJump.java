package multitypespike.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import sa.evolution.operators.LeafToSampledAncestorJump;

@Description("Performs a joint spike and LeafToSampledAncestorJump topology adjustment")
public class SpikeAndSampledAncestorJump extends LeafToSampledAncestorJump {

    final public Input<RealParameter> spikesInput = new Input<>("spikes", "spikes associated with each branch on the tree",
            Input.Validate.REQUIRED);

    final public Input<RealParameter> spikeShapeInput = new Input<>("spikeShape", "shape parameter for the " +
            "gamma distribution of the spikes", Input.Validate.REQUIRED);


    /**
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        double newHeight, newRange, oldRange;

//        int categoryCount = 1;
//        if (categoriesInput.get() != null) {
//
//            categoryCount = categoriesInput.get().getUpper() - categoriesInput.get().getLower() +1;
//        }

        Tree tree = treeInput.get();

        int leafNodeCount = tree.getLeafNodeCount();

        // Select a leaf
        Node leaf = tree.getNode(Randomizer.nextInt(leafNodeCount));
        Node parent = leaf.getParent();

        double spikeShape = spikeShapeInput.get().getValue();
        GammaDistribution gamma = new GammaDistributionImpl(spikeShape, 1/spikeShape);
        RealParameter spikes = spikesInput.get();

        // Spike of the sampled ancestor
        double sOld = spikes.getValue(leaf.getNr());
        double sNew, pFwd, pRev;

        // Sampled ancestor to normal leaf move: (spike=0, length=0) –> (spike>0, length>0)
        if (leaf.isDirectAncestor()) {

            if (sOld != 0) {
                return Double.NEGATIVE_INFINITY; // Reject
            }

            oldRange = (double) 1;
            if (parent.isRoot()) {

                final double randomNumber = Randomizer.nextExponential(1);
                newHeight = parent.getHeight() + randomNumber;
                newRange = Math.exp(randomNumber);
            } else {
                newRange = parent.getParent().getHeight() - parent.getHeight();
                newHeight = parent.getHeight() + Randomizer.nextDouble() * newRange;
            }

            // Add a spike
            sNew = Randomizer.nextGamma(spikeShape, 1/spikeShape);
            pRev = 0;
            pFwd = gamma.logDensity(sNew);

//            if (categoriesInput.get() != null) {
//                int index = leaf.getNr();
//                int newValue = Randomizer.nextInt(categoryCount) + categoriesInput.get().getLower(); // from 0 to n-1, n must > 0,
//                categoriesInput.get().setValue(index, newValue);
//            }

        }
        // Normal leaf to sampled ancestor move:  (spike>0, length>0) –> (spike=0, length=0)
        else {
            newRange = (double) 1;
            //make sure that the branch where a new sampled node to appear is not above that sampled node
            if (getOtherChild(parent, leaf).getHeight() >= leaf.getHeight())  {
                return Double.NEGATIVE_INFINITY;
            }
            if (parent.isRoot()) {
                oldRange = Math.exp(parent.getHeight() - leaf.getHeight());
            } else {
                oldRange = parent.getParent().getHeight() - leaf.getHeight();
            }
            newHeight = leaf.getHeight();

            // Remove the spike
            sNew = 0;
            pFwd = 0;
            pRev = gamma.logDensity(sOld);

//            if  (categoriesInput.get() != null) {
//                int index = leaf.getNr();
//                categoriesInput.get().setValue(index, -1);
//            }
        }

        //make sure that either there are no direct ancestors or r<1
        if ((rInput.get() != null) && (tree.getDirectAncestorNodeCount() > 0 && rInput.get().getValue() == 1))  {
            return Double.NEGATIVE_INFINITY;
        }

        // Update state
        spikes.setValue(leaf.getNr(), sNew);
        parent.setHeight(newHeight);

        return Math.log(newRange/oldRange) + pRev - pFwd;
    }

}
