package multitypespike.logger;

import beast.base.core.Description;
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


@Description("Samples and logs the exact number of hidden speciation events per branch")
public class HiddenEventsLogger extends CalculationNode implements Function, Loggable {
    final public Input<BranchSpikePrior> branchSpikePriorInput =
            new Input<>("branchSpikePrior", "Branch spike prior", Input.Validate.REQUIRED);
    final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "if false then no spikes are inferred", Input.Validate.OPTIONAL);
    final public Input<Boolean> logPerTypeInput = new Input<>(
            "logPerType","If true, log hidden events of each type separately for multi-type models; " +
                    "if false, log totals per node (sum across types).",false); // default: sum across types


    private BranchSpikePrior bsp;
    private int nTypes, nodeCount;
    private boolean logPerType;

    @Override
    public void initAndValidate() {
        bsp = branchSpikePriorInput.get();
        nTypes = bsp.nTypes;
        nodeCount = bsp.nodeCount;
        logPerType = logPerTypeInput.get();

        if (nTypes == 1 && logPerType) throw new RuntimeException("logPerType cannot be true for single-type models.");
    }

    @Override
    public void log(long sample, PrintStream out) {
        for (int i = 0; i < this.getDimension(); i ++) {
            out.print(this.getArrayValue(i) + "\t");
        }
    }

    @Override
    public int getDimension() {
        if(logPerType) return nTypes * nodeCount;
        else return nodeCount;
    }

    @Override
    public double getArrayValue(int dim) {
        if(nTypes == 1) {
            return sampleHiddenEvent(dim);
        } else if(!logPerType){
            // Sum across all types for each node
            double sum = 0.0;
            for (int type = 0; type < nTypes; type++) {
                sum += sampleHiddenEvent(dim, type);
            }
            return sum;
        }
        else {
            int type = dim % nTypes;
            int nodeNr = dim / nTypes;
                return sampleHiddenEvent(nodeNr, type);
        }
    }

    @Override
    public void init(PrintStream out) {
        String id = this.getID();
        if (id == null || id.isEmpty()) id = "nHiddenEvents";

        if (logPerType) {
            for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
                for (int t = 0; t < nTypes; t++) {
                    out.print(id + ".node" + nodeNr + ".type" + t + "\t");
                }
            }
        } else {
            for (int nodeNr = 0; nodeNr < nodeCount; nodeNr++) {
                out.print(id + ".node" + nodeNr + "\t");
            }
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

    private int sampleHiddenEvent(int nodeNr, int type) {
        Node node = bsp.treeInput.get().getNode(nodeNr);
        if (node.isRoot() || node.isDirectAncestor()) {
            return 0;
        }
        double[] cf = getCumulativeProbs(nodeNr, type);

        return Randomizer.randomChoice(cf);
    }


    final double MAX_CUM_SUM = 0.999;

    /**
     * Calculate cumulative probabilities of sampling hidden events, conditional on the gamma distribution and tree prior (theta)
     * p(hiddenEvents | spike size, theta) = p(spike size | hiddenEvents, theta) x p (hiddenEvents | theta) / p (spike size | theta)
     */

    // Single-type version
    public double[] getCumulativeProbs(int nodeNr) {

        Node node = bsp.treeInput.get().getNode(nodeNr);

        if (node.isRoot() || node.isDirectAncestor()) {
            return new double[0];
        }

        List<Double> probs = new ArrayList<>();

        // Check spikeShape is positive
        double spikeShape = bsp.spikeShapeInput.get().getValue();
        if (spikeShape <= 0) {
            throw new IllegalArgumentException("Cannot sample spikes because spikeShape is non-positive " + spikeShape);
        }

        // Spike size
        double branchSpike = bsp.spikesInput.get().getValue(nodeNr);

        double expNrHiddenEvents = bsp.getExpectedHiddenEvents(nodeNr);

        int k = 0;
        double poissonCumSum = 0;

        double branchPSum = 0;
        while (poissonCumSum < MAX_CUM_SUM) {

            double branchP = 0;

            // Probability of k hidden events P(k) under a Poisson(mu)
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

    // Multi-type version
    public double[] getCumulativeProbs(int nodeNr, int type) {

        Node node = bsp.treeInput.get().getNode(nodeNr);

        if (node.isRoot() || node.isDirectAncestor()) {
            return new double[0];
        }

        List<Double> probs = new ArrayList<>();

        double spikeShape;
        if (bsp.spikeShapeInput.get().getDimension() == 1) spikeShape = bsp.spikeShapeInput.get().getValue();
         else spikeShape = bsp.spikeShapeInput.get().getValue(type);

        // Check spikeShape is positive
        if (spikeShape <= 0) {
            throw new IllegalArgumentException("Cannot sample spikes because spikeShape is non-positive " + spikeShape);
        }

        // Spike size
        double branchSpike = bsp.spikesInput.get().getValue(nodeNr * nTypes + type);

        double expNrHiddenEvents = bsp.getExpectedHiddenEvents(nodeNr, type);

        double pi = bsp.getPiVals(nodeNr, type);

        int k = 0;
        double poissonCumSum = 0;

        double branchPSum = 0;
        while (poissonCumSum < MAX_CUM_SUM) {

            double branchP = 0;

            // Probability of k hidden events P(k) under a Poisson(mu)
            double logpk = -expNrHiddenEvents + k*Math.log(expNrHiddenEvents) - Gamma.logGamma(k+1);
            double pk = Math.exp(logpk);

            if (indicatorInput.get() != null && indicatorInput.get().getValue()) {

                // Integrate over type of the observed speciation event
                for (int obsEvent = 0; obsEvent <= 1; obsEvent++) {
                    double pObs = obsEvent == 1 ? pi : (1 - pi);

                    // Integrate across all possible values in poisson distribution
                    int nSpikes = node.getParent().isFake() ? k : k + 1;

                    if (nSpikes == 0) {

                        // Valid zero spike
                        if (branchSpike == 0) branchP += Math.exp(logpk + Math.log(pObs));

                    } else {
                        // Compute log-probability of observed spike under Gamma distribution
                        GammaDistribution gamma = new GammaDistributionImpl(
                                spikeShape * nSpikes, 1 / spikeShape);
                        double gammaLogP = gamma.logDensity(branchSpike);

                        if (branchSpike != 0 && Double.isFinite(gammaLogP)) {
                            branchP += Math.exp(logpk + gammaLogP + Math.log(pObs));
                        }
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
