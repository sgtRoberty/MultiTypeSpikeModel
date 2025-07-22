package multitypespike.operator;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;



// DELETE - UNUSED

/**
 * Utility class for selecting valid branch indices in spike operators.
 */
public final class SpikeOperatorBase {
    /**
     * Returns a random valid branch index from the current tree state.
     * Filters out root and sampled ancestor branches.
     */
    static int getRandomValidBranch(Tree tree, RealParameter spikes, int maxTries) {

        for (int tries = 0; tries < maxTries; tries++) {
            int index = Randomizer.nextInt(tree.getNodeCount());
            Node node = tree.getNode(index);

            // Skip root
            if (node.isRoot()) continue;

            // If sampled ancestor: allow only if spike is non-zero
            if (node.isDirectAncestor()) {
                if (spikes == null || spikes.getValue(index) == 0.0) continue;
            }

            return index;
        }

        return -1;
    }
}
