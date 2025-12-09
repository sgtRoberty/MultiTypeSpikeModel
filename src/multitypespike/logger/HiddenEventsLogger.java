package multitypespike.logger;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.util.Randomizer;
import multitypespike.distribution.BranchSpikePrior;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math.special.Gamma;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;


public class HiddenEventsLogger extends CalculationNode implements Function, Loggable {
    final public Input<BranchSpikePrior> branchSpikePriorInput =
            new Input<>("branchSpikePrior", "Branch spike prior", Input.Validate.REQUIRED);
    final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "if false then no spikes are inferred", Input.Validate.OPTIONAL);


    private BranchSpikePrior bsp;
    private int nTypes, nodeCount;

    @Override
    public void initAndValidate() {
        bsp = branchSpikePriorInput.get();
        nTypes = bsp.nTypes;
        nodeCount = bsp.nodeCount;
        if (nTypes>1) throw new RuntimeException("Not implemented for multi-type yet.");
    }

    @Override
    public void log(long sample, PrintStream out) {
        for (int i = 0; i < this.getDimension(); i ++) {
            out.print(this.getArrayValue(i) + "\t");
        }
    }

    @Override
    public int getDimension() {
        if(nTypes == 1) return nTypes * nodeCount;
        else return nodeCount;
    }

    @Override
    public double getArrayValue(int dim) {
        return sampleHiddenEvent(dim);
    }

    @Override
    public void init(PrintStream out) {
        String id = this.getID();
        if (id == null || id.isEmpty()) id = "nHiddenEvents";

        for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
            out.print(id + ".node" + nodeNr + "\t");
        }
    }


    private int sampleHiddenEvent(int nodeNr) {
        Node node = bsp.treeInput.get().getNode(nodeNr);
        if (node.isRoot() || node.isDirectAncestor()) {
            return 0;
        }

        double[] cf = getCumulativeProbs(nodeNr);
        return Randomizer.randomChoice(cf);
    }


    final double MAX_CUM_SUM = 0.999;

    /**
     * Calculate cumulative probabilities of sampling a stub, conditional on the gamma distribution and tree prior (theta)
     * p(hiddenEvents | spike size, theta) = p(spike size | hiddenEvents, theta) x p (hiddenEvents | theta) / p (spike size | theta)
     */
    public double[] getCumulativeProbs(int nodeNr) {

        double expNrHiddenEvents = bsp.getExpectedHiddenEvents(nodeNr);
        Node node = bsp.treeInput.get().getNode(nodeNr);

        List<Double> probs = new ArrayList<>();

        // Check spikeShape is positive
        double spikeShape = bsp.spikeShapeInput.get().getValue();
        if (spikeShape <= 0) {
            throw new IllegalArgumentException("Cannot sample spikes because spikeShape is non-positive " + spikeShape);
        }

        if (node.isRoot() || node.isDirectAncestor()) {
            return new double[0];
        }

        // Spike size
        double branchSpike = bsp.spikesInput.get().getValue(nodeNr);

        int k = 0;
        double poissonCumSum = 0;

        double branchPSum = 0;
        while (poissonCumSum < MAX_CUM_SUM) {

            double branchP = 0;

            // P(k observations) under a Poisson(mu)
            // Probability of observing k hidden events P(k) under a Poisson(mu)
            double logpk = -expNrHiddenEvents + k*Math.log(expNrHiddenEvents) - Gamma.logGamma(k+1);
            double pk = Math.exp(logpk);

            if (indicatorInput.get() != null && indicatorInput.get().getValue()) {

            // Integrate across all possible values in poisson distribution
            int nSpikes = node.getParent().isFake() ? k : k + 1;

                if (nSpikes == 0) {
                    // Valid zero spike
                    if (branchSpike == 0) branchP += Math.exp(logpk);

                } else {
                    // Compute log-probability of observed spike under Gamma distribution
                    GammaDistribution gamma = new GammaDistributionImpl(
                            spikeShape * nSpikes, 1 / spikeShape);
                    double gammaLogP = gamma.logDensity(branchSpike);
                    if (branchSpike != 0 && Double.isFinite(gammaLogP)) {

                        branchP += Math.exp(logpk + gammaLogP);
                    }
                }
            } else {
                branchP += pk;
            }

            poissonCumSum += pk;
            branchPSum += branchP;
            probs.add(branchP);

            k++;
        }


        // Normalise to sum to 1
        double[] array = new double[probs.size()];
        for(int i = 0; i < probs.size(); i++) {
            array[i] = probs.get(i) / branchPSum;
        }


        // Convert into cumulative sum
        double cumsum = 0;
        for(int i = 0; i < probs.size(); i++) {
            double p = array[i];
            array[i] = p + cumsum;
            cumsum += p;
        }

        return array;

    }


    @Override
    public void close(PrintStream out) {
    }

}
