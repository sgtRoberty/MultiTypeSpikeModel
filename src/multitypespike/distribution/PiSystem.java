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

public class PiSystem implements FirstOrderDifferentialEquations, Loggable {

    private final ContinuousOutputModel[] p0geComArray;

    private final double[][] b;
    private final double[][][] M, b_ij;

    public int nTypes, nIntervals;
    private int currentNodeNr;
    public double[] intervalEndTimes;
    private final Tree tree;
    protected int interval;

    protected FirstOrderIntegrator piIntegrator;

    protected double integrationMinStep, integrationMaxStep;

    public ContinuousOutputModel[] integrationResults;

    public PiSystem(Parameterization parameterization, Tree tree, ContinuousOutputModel[] p0geComArray, double absoluteTolerance, double relativeTolerance) {

        this.tree = tree;
        this.p0geComArray = p0geComArray;
        this.b = parameterization.getBirthRates();
        this.M = parameterization.getMigRates();
        this.b_ij = parameterization.getCrossBirthRates();

        this.nTypes = parameterization.getNTypes();
        this.nIntervals = parameterization.getTotalIntervalCount();

        this.intervalEndTimes = parameterization.getIntervalEndTimes();

        integrationMinStep = parameterization.getTotalProcessLength() * 1e-100;
        integrationMaxStep= parameterization.getTotalProcessLength() / 10;

        integrationResults = new ContinuousOutputModel[tree.getNodeCount()];

        this.piIntegrator = new DormandPrince54Integrator(
                integrationMinStep, integrationMaxStep,
                absoluteTolerance, relativeTolerance);

    }


    public void setInterval(int interval) {this.interval = interval;}

    public int getDimension() {
        return this.nTypes;
    }

    public void setCurrentNodeNr(int nodeNr) {
        this.currentNodeNr = nodeNr;
    }

    private ContinuousOutputModel getP0GeIntegrationResults(int nodeNr) {
        return p0geComArray[nodeNr];
    }

    public double[] getP0Ge(int nodeNr, double time) {
        //  p0:  (0 .. dim-1)
        //  ge: (dim .. 2*dim-1)
        ContinuousOutputModel p0geCom = getP0GeIntegrationResults(nodeNr);
        p0geCom.setInterpolatedTime(time);
        double[] p0ge = p0geCom.getInterpolatedState();

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

//        System.out.format("*** t=%g pi[0]=%g pi[1]=%g\n", t, y[0], y[1]);

        for (int i = 0; i< nTypes; i++) {
            yDot[i] = 0;

            for (int j = 0; j < nTypes; j++) {
                if (j == i) continue;


//                yDot[i] += ((b_ij[interval][j][i] * p0ge[j] + M[interval][j][i]) * (p0ge[nTypes + i] / p0ge[nTypes + j])) * y[j];
//                yDot[i] -= ((b_ij[interval][i][j] * p0ge[i] + M[interval][i][j]) * (p0ge[nTypes + j] / p0ge[nTypes + i])) * y[i];

                yDot[i] += ((b_ij[interval][j][i] * p0ge[j] + M[interval][j][i]) * (p0ge[nTypes + i] / Math.max(p0ge[nTypes + j], 1e-6))) * y[j];
                yDot[i] -= ((b_ij[interval][i][j] * p0ge[i] + M[interval][i][j]) * (p0ge[nTypes + j] / Math.max(p0ge[nTypes + i], 1e-6))) * y[i];
            }
        }
    }


    public void setInitialConditionsForPi(PiState state, double[] startTypePriorProbs, double rootTime) {
        double total = 0.0;
        double[] p0geInit = getP0Ge(currentNodeNr, rootTime);

        for (int type = 0; type < nTypes; type++) {
            state.pi[type] = p0geInit[type + nTypes] * startTypePriorProbs[type];
            total += state.pi[type];
        }

        for (int type = 0; type < nTypes; type++) {
            state.pi[type] /= total;
        }
    }


    public ContinuousOutputModel getIntegrationResultsForNode(int nodeNr) {
        return integrationResults[nodeNr];
    }

    public void integrate(PiState state, double tStart, double tEnd) {
        piIntegrator.integrate(this, tStart, state.pi, tEnd, state.pi);
    }

//    public void integratePiAlongEdge(Node node, double tStart, PiState state,
//                                  Parameterization parameterization, double finalSampleOffset) {
//        ContinuousOutputModel continuousOutputModel = new ContinuousOutputModel();
////        ContinuousOutputModel com = new ContinuousOutputModel();
////        piIntegrator.clearStepHandlers();
////        piIntegrator.addStepHandler(com);
//
//        setCurrentNodeNr(node.getNr());
//
//        double thisTime = tStart;
//        double tEnd = parameterization.getNodeTime(node, finalSampleOffset);
//
////        if (node.isLeaf()) tEnd = tEnd - 1e-4;
//
//        int thisInterval = parameterization.getIntervalIndex(thisTime);
//        int endInterval = parameterization.getNodeIntervalIndex(node, finalSampleOffset);
//
//        while (thisInterval < endInterval) {
//            // Forward integration across intervals
//            double nextTime = intervalEndTimes[thisInterval];
//
//            if (Utils.lessThanWithPrecision(thisTime, nextTime)) {
////                ContinuousOutputModel com = new ContinuousOutputModel();
////                piIntegrator.clearStepHandlers();
////                piIntegrator.addStepHandler(com);
//
//                setInterval(thisInterval);
//                integrate(state, thisTime, nextTime);
////                com.append(com);
//
//                continuousOutputModel.append(com);
//
////                for (int idx=0; idx<101; idx+=1) {
////                    double thist = tStart + idx*(tEnd-tStart)/101;
////                    com.setInterpolatedTime(thist);
////                    double[] thispi = com.getInterpolatedState();
////                    System.out.format("*** t=%g pi[0]=%g pi[1]=%g\n", thist, thispi[0], thispi[1]);
////                }
////                System.exit(0);
//            }
//
//            thisTime = nextTime;
//            thisInterval += 1;
//        }
//
//        if (Utils.lessThanWithPrecision(thisTime, tEnd)) {
//            setInterval(thisInterval);
//            integrate(state, thisTime, tEnd);
//            continuousOutputModel.append(com);
//
//        }
//
//        // Store integration results for each edge
////        integrationResults[node.getNr()] = com;
//        integrationResults[node.getNr()] = continuousOutputModel;
//    }



    public void integratePiAlongEdge(Node node, double tStart, PiState state,
                                     Parameterization parameterization, double finalSampleOffset) {

        ContinuousOutputModel fullModel = new ContinuousOutputModel();

        setCurrentNodeNr(node.getNr());

        double thisTime = tStart;
        double tEnd = parameterization.getNodeTime(node, finalSampleOffset);

        int thisInterval = parameterization.getIntervalIndex(thisTime);
        int endInterval = parameterization.getNodeIntervalIndex(node, finalSampleOffset);

        while (thisInterval < endInterval) {

            double nextTime = intervalEndTimes[thisInterval];

            if (Utils.lessThanWithPrecision(thisTime, nextTime)) {

                ContinuousOutputModel segment = new ContinuousOutputModel();

                piIntegrator.clearStepHandlers();
                piIntegrator.addStepHandler(segment);

                setInterval(thisInterval);
                integrate(state, thisTime, nextTime);

                fullModel.append(segment);
            }

            thisTime = nextTime;
            thisInterval++;
        }

        if (Utils.lessThanWithPrecision(thisTime, tEnd)) {

            ContinuousOutputModel segment = new ContinuousOutputModel();

            piIntegrator.clearStepHandlers();
            piIntegrator.addStepHandler(segment);

            setInterval(thisInterval);
            integrate(state, thisTime, tEnd);

            fullModel.append(segment);
        }

        integrationResults[node.getNr()] = fullModel;
    }



    public void integratePiAtNode(Node node, double parentTime, PiState state,
                                  Parameterization parameterization, double finalSampleOffset) {

        // Get the time of this node
        double nodeTime = parameterization.getNodeTime(node, finalSampleOffset);

        // Skip integration on origin
        if (!node.isRoot()) {
            // Determine intervals and integrate from parent to this node
            integratePiAlongEdge(node, parentTime, state, parameterization, finalSampleOffset);
        }
        // Recurse to children
        for (Node child : node.getChildren()) {
            PiState childState = new PiState(nTypes);
            System.arraycopy(state.pi, 0, childState.pi, 0, nTypes);
            // Copy current state as starting condition for child
            integratePiAtNode(child, nodeTime, childState, parameterization, finalSampleOffset);
        }
    }


    public void integratePi(double[] startTypePriorProbs,
                            Parameterization parameterization, double finalSampleOffset) {

        Node root = this.tree.getRoot();
        setCurrentNodeNr(root.getNr());
        double rootTime = parameterization.getNodeTime(root, finalSampleOffset);
        PiState state = new PiState(parameterization.getNTypes());

        // Set initial conditions at root
        setInitialConditionsForPi(state, startTypePriorProbs, rootTime);
        // Start pre-order traversal integration from root
        integratePiAtNode(root, rootTime,
                state, parameterization, finalSampleOffset);
    }

    // FOR TESTING PURPOSES ONLY
    public void integratePiSingleLineage(double[] startTypePriorProbs,
                            Parameterization parameterization, double startTime, double endTime) {

        Node root = this.tree.getRoot();
        if(!root.isLeaf()) throw new RuntimeException("Not a single-lineage tree!");

        setCurrentNodeNr(root.getNr());
        PiState state = new PiState(parameterization.getNTypes());

        // Set initial conditions at root
        setInitialConditionsForPi(state, startTypePriorProbs, startTime);

        ContinuousOutputModel com = new ContinuousOutputModel();
        piIntegrator.addStepHandler(com);

        // Integrate from root to tip
        integrate(state, startTime, endTime);

        // Store result
        integrationResults[root.getNr()] = com;
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