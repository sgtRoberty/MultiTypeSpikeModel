package multitypespike.distribution;

import bdmmprime.parameterization.Parameterization;
import bdmmprime.util.Utils;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;
import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince54Integrator;

import java.io.PrintStream;



public class MultiTypeHiddenEventsIntegrator implements FirstOrderDifferentialEquations, Loggable {

    private final ContinuousOutputModel[] p0geComArray;

    private final double[][] b;
    private final double[][][] M, b_ij;

    private final double[] state;
    private final double[][] piAtNodes;
    private final double[][] expNrHiddenEvents;

    public int nTypes, nIntervals;
    private int currentNodeNr;
    public double[] intervalEndTimes;
    private final Tree tree;
    protected int interval;

    protected FirstOrderIntegrator integrator;

    protected double integrationMinStep, integrationMaxStep;
    private static final double EPS = 1e-8; // Small value to avoid division by zero

    // Optional storage of π continuous output trajectories for testing purposes
    private final boolean storePiTrajectories;
    private final ContinuousOutputModel[] piIntegrationResults;

    public MultiTypeHiddenEventsIntegrator(Parameterization parameterization, Tree tree, ContinuousOutputModel[] p0geComArray, double absoluteTolerance, double relativeTolerance, boolean storePiTrajectories) {

        this.tree = tree;
        this.p0geComArray = p0geComArray;
        this.b = parameterization.getBirthRates();
        this.M = parameterization.getMigRates();
        this.b_ij = parameterization.getCrossBirthRates();
        this.nTypes = parameterization.getNTypes();
        this.nIntervals = parameterization.getTotalIntervalCount();
        this.intervalEndTimes = parameterization.getIntervalEndTimes();

        this.state = new double[2 * nTypes];
        this.piAtNodes = new double[tree.getNodeCount()][nTypes];
        this.expNrHiddenEvents = new double[tree.getNodeCount()][nTypes];

        integrationMinStep = parameterization.getTotalProcessLength() * 1e-100;
        integrationMaxStep= parameterization.getTotalProcessLength() / 10;

        this.storePiTrajectories = storePiTrajectories;
        this.piIntegrationResults = storePiTrajectories ? new ContinuousOutputModel[tree.getNodeCount()] : null;

        this.integrator = new DormandPrince54Integrator(
                integrationMinStep, integrationMaxStep,
                absoluteTolerance, relativeTolerance
        );
    }


    public void setInterval(int interval) {this.interval = interval;}

    @Override
    public int getDimension() {
        return 2*this.nTypes;
    }

    public void setCurrentNodeNr(int nodeNr) {
        this.currentNodeNr = nodeNr;
    }

    private ContinuousOutputModel getP0GeIntegrationResults(int nodeNr) {
        return p0geComArray[nodeNr];
    }


    /**
     * Returns the interpolated p0 and ge values at a given time for an edge.
     *
     * @param nodeNr
     * @param time
     * @return
     */
    public double[] getP0Ge(int nodeNr, double time) {
        //  p0:  (0 .. dim-1)
        //  ge: (dim .. 2*dim-1)
        ContinuousOutputModel p0geCom = getP0GeIntegrationResults(nodeNr);
        p0geCom.setInterpolatedTime(time);
        double[] p0geRaw = p0geCom.getInterpolatedState();
        double[] p0ge = p0geRaw.clone();

        // Trim away small negative values due to numerical integration errors
        for (int i=0; i<p0ge.length; i++) {
            if (p0ge[i] < 0)
                p0ge[i] = 0.0;
        }
        return p0ge;
    }


    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        double[] p0ge = getP0Ge(currentNodeNr, t);


        for (int i = 0; i< nTypes; i++) {

            /*  π equations (0 .. dim-1)  */
            yDot[i] = 0;

            for (int j = 0; j < nTypes; j++) {
                if (j == i) continue;


                yDot[i] += ((b_ij[interval][j][i] * p0ge[j] + M[interval][j][i]) * (p0ge[nTypes + i] / Math.max(p0ge[nTypes + j], EPS))) * y[j];       // Floor denominator at EPS to avoid division by zero
                yDot[i] -= ((b_ij[interval][i][j] * p0ge[i] + M[interval][i][j]) * (p0ge[nTypes + j] / Math.max(p0ge[nTypes + i], EPS))) * y[i];
            }

            /*  Hidden speciation events ((dim .. 2*dim-1) */
            yDot[nTypes + i] = 2.0 * y[i] * b[interval][i] * p0ge[i];
        }
    }


    public void setInitialConditionsAtRoot(double[] startTypePriorProbs, double rootTime) {
        double total = 0.0;
        double[] p0geInit = getP0Ge(currentNodeNr, rootTime);

        for (int i = 0; i < nTypes; i++) {
            this.state[i] = p0geInit[i + nTypes] * startTypePriorProbs[i];
            total += this.state[i];
        }

        for (int i = 0; i < nTypes; i++) {
            this.state[i] /= total;
            this.state[nTypes + i] = 0.0;
        }
    }


    public void setInitialConditionsAtNode(double[] initialPi) {
        for (int i = 0; i < nTypes; i++) {
            // π initialisation
            state[i] = initialPi[i];

            // Reset hidden events to zero at start of each branch
            state[nTypes + i] = 0.0;
        }
    }

    public void storeResultsAtNode(int nodeNr) {

        for (int i = 0; i < nTypes; i++) {
            // Store π at node
            piAtNodes[nodeNr][i] = state[i];

            // Store hidden events at node
            expNrHiddenEvents[nodeNr][i] = state[nTypes + i];
        }
    }


    public void integrate(double tStart, double tEnd) {
        integrator.integrate(this, tStart, this.state, tEnd, this.state);
    }


    public void integrateAlongEdge(Node node, double tStart, Parameterization parameterization,
                                   double finalSampleOffset) {

        setCurrentNodeNr(node.getNr());
        double thisTime = tStart;
        double tEnd = parameterization.getNodeTime(node, finalSampleOffset);

        int thisInterval = parameterization.getIntervalIndex(thisTime);
        int endInterval = parameterization.getNodeIntervalIndex(node, finalSampleOffset);

        // If storePiTrajectories = false, perform integration without storing π trajectories
        if (!storePiTrajectories) {
            while (thisInterval < endInterval) {
                double nextTime = intervalEndTimes[thisInterval];
                if (Utils.lessThanWithPrecision(thisTime, nextTime)) {
                    setInterval(thisInterval);
                    integrate(thisTime, nextTime);
                }
                thisTime = nextTime;
                thisInterval++;
            }
            if (Utils.lessThanWithPrecision(thisTime, tEnd)) {
                setInterval(thisInterval);
                integrate(thisTime, tEnd);
            }
            storeResultsAtNode(currentNodeNr);
            return;
        }

        // If storePiTrajectories = true, sotre ContinuousOutputModel for π
        ContinuousOutputModel fullModel = new ContinuousOutputModel();

        while (thisInterval < endInterval) {
            double nextTime = intervalEndTimes[thisInterval];

            if (Utils.lessThanWithPrecision(thisTime, nextTime)) {
                // Make a new continuous output model for this interval
                ContinuousOutputModel segment = new ContinuousOutputModel();
                integrator.clearStepHandlers();
                integrator.addStepHandler(segment);

                setInterval(thisInterval);
                integrate(thisTime, nextTime);

                // Append segment to fullModel
                fullModel.append(segment);
            }

            thisTime = nextTime;
            thisInterval++;
        }

        if (Utils.lessThanWithPrecision(thisTime, tEnd)) {
            ContinuousOutputModel segment = new ContinuousOutputModel();
            integrator.clearStepHandlers();
            integrator.addStepHandler(segment);

            setInterval(thisInterval);
            integrate(thisTime, tEnd);

            fullModel.append(segment);
        }

        // store both numeric results and the continuous output model if requested
        storeResultsAtNode(currentNodeNr);
        if (storePiTrajectories) {
            piIntegrationResults[currentNodeNr] = fullModel;
        }
    }


    public void integrateAtNode(Node node, double parentTime, Parameterization parameterization, double finalSampleOffset) {

        // Get the time of this node
        double nodeTime = parameterization.getNodeTime(node, finalSampleOffset);

        // Skip integration on origin and zero length branches
        if (!node.isRoot() && !node.isDirectAncestor()) {
            // Determine intervals and integrate from parent to this node
            integrateAlongEdge(node, parentTime, parameterization, finalSampleOffset);
        }

        // Recurse to children
        for (Node child : node.getChildren()) {

            // Set initial conditions
            setInitialConditionsAtNode(piAtNodes[node.getNr()]);

            integrateAtNode(child, nodeTime, parameterization, finalSampleOffset);
        }
    }

    /**
     * Performs recursive integration of π and expected hidden events along the tree.
     *
     * @param startTypePriorProbs
     * @param parameterization
     * @param finalSampleOffset
     */
    public void integrateHiddenEvents(double[] startTypePriorProbs, Parameterization parameterization, double finalSampleOffset) {

        Node root = this.tree.getRoot();
        setCurrentNodeNr(root.getNr());
        double rootTime = parameterization.getNodeTime(root, finalSampleOffset);

        // Set initial conditions at root
        setInitialConditionsAtRoot(startTypePriorProbs, rootTime);

        // Store initial conditions at root
        storeResultsAtNode(currentNodeNr);

        // Start pre-order traversal integration from root
        integrateAtNode(root, rootTime, parameterization, finalSampleOffset);
    }

    public ContinuousOutputModel getPiIntegrationResultsForNode(int nodeNr) {
        if (!storePiTrajectories) {
            throw new IllegalStateException("π trajectories not stored; set storePiTrajectories to true.");
        }
        return piIntegrationResults[nodeNr];
    }

    /**
     *  Retrieves the expected number of hidden events per type for a branch
     *
     * @param nodeNr
     * @return
     */
    public double[] getExpNrHiddenEventsForNode(int nodeNr) {
        return expNrHiddenEvents[nodeNr];
    }

    /**
     *  Retrieves array of π values at a node
     *
     * @param nodeNr
     * @return
     */
    public double[] getPiAtNode(int nodeNr) {
        return piAtNodes[nodeNr];
    }

    /**
     *  Integrates π and hidden events along a single lineage tree (for testing purposes only).
     *
     * @param startTypePriorProbs
     * @param parameterization
     * @param startTime
     * @param endTime
     */
    public void integrateSingleLineage(double[] startTypePriorProbs, Parameterization parameterization, double startTime, double endTime) {

        Node root = this.tree.getRoot();

        if (!root.isLeaf()) {
            throw new RuntimeException("Tree has more than one lineage! Only single-lineage trees supported.");
        }

        setCurrentNodeNr(root.getNr());

        // Set initial conditions at root
        setInitialConditionsAtRoot(startTypePriorProbs, startTime);

        ContinuousOutputModel com = new ContinuousOutputModel();
        if (storePiTrajectories) {
            // Optional: store continuous output
            integrator.clearStepHandlers();
            integrator.addStepHandler(com);
        }

        // Integrate from startTime to endTime
        integrate(startTime, endTime);

        // Store results
        storeResultsAtNode(root.getNr());

        // If you want to keep the full trajectory for testing
        if (storePiTrajectories) {
            piIntegrationResults[root.getNr()] = com;
        }
    }


    @Override
    public void init(PrintStream out) {
    }

    @Override
    public void log(long sample, PrintStream out) {
    }

    @Override
    public void close(PrintStream out) {
    }

}