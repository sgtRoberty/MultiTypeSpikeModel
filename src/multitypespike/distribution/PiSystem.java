package multitypespike.distribution;

import bdmmprime.distribution.P0System;
import bdmmprime.parameterization.Parameterization;

public class PiSystem extends P0System {

    public double[][]  p0, ge;

    public PiSystem(Parameterization parameterization, double absoluteTolerance, double relativeTolerance) {
        super(parameterization, absoluteTolerance, relativeTolerance);
    }

    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {

        for (int i = 0; i< nTypes; i++) {

            yDot[i] = 2 * y[i] * (b[interval][i] * p0[i] + M[interval][i][i]);

            for (int j = 0; j < nTypes; j++) {
                if (j == i) continue;

                yDot[i] += ((b_ij[interval][i][j] * p0[i] + M[interval][i][j]) * (ge[j] / ge[i])) * y[i];
                yDot[i] -= ((b_ij[interval][j][i] * p0[j] + M[interval][j][i]) * (ge[i] / ge[j])) * y[j];

            }
        }
    }




}