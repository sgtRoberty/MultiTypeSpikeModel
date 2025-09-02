package multitypespike.distribution;

import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;


public class TypedHiddenEventsODE implements FirstOrderDifferentialEquations {

    private final int i;
    private final double lambda_i;
    private final double[][] lambda_ij;
    private final int nTypes;
    private final PiState piState;
    private final int nodeNr;
    private final ContinuousOutputModel[] p0geComArray;


    public TypedHiddenEventsODE(int i, int nodeNr, double lambda_i, double[][] lambda_ij,
                                int nTypes, PiState piState, ContinuousOutputModel[] p0geComArray) {
        this.i = i;
        this.lambda_i = lambda_i;
        this.lambda_ij = lambda_ij;
        this.nTypes = nTypes;
        this.piState = piState;
        this.nodeNr = nodeNr;
        this.p0geComArray = p0geComArray;
    }


    @Override
    public int getDimension() {
        return 1;
    }


    private ContinuousOutputModel getP0GeModel(int nodeNr) {
        return p0geComArray[nodeNr];
    }


    public double getP0_i(int nodeNr, double time, int i) {
        ContinuousOutputModel p0geCom = getP0GeModel(nodeNr);
        p0geCom.setInterpolatedTime(time);
        double[] p0ge = p0geCom.getInterpolatedState();

        return ( p0ge[i] < 0) ? 0.0 :  p0ge[i];
    }


    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {

        double p0_i = getP0_i(nodeNr, t, i);
        double pi_i = piState.getPiAtTime(nodeNr, t)[i];

        double sum = 2 * lambda_i * p0_i;

        for (int j = 0; j < nTypes; j++) {
            if (j != i) {
                double p0_j = getP0_i(nodeNr, t, j);
                sum += lambda_ij[i][j] * p0_j;
            }
        }

        yDot[0] = pi_i * sum;
    }
}
