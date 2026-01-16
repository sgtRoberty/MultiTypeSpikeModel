package multitypespike.operator;

import bdmmprime.parameterization.Parameterization;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.kernel.BactrianScaleOperator;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;


@Description("Moves spikes up and corresponding rates down or vice versa")
public class SpikeUpDownOperator extends BactrianScaleOperator {

    final public Input<RealParameter> spikesInput = new Input<>("spikes", "spikes associated with each branch", Input.Validate.REQUIRED);
    final public Input<Parameterization> parameterizationInput = new Input<>("parameterization",
            "BDMM-prime parameterization object (see BDMM-prime package for available parameterizations)",
            Input.Validate.REQUIRED);
    final public Input<Boolean> reverseInput = new Input<>("reverseDirection", "If true, moves spikes up and corresponding rates down", false);

    private int nTypes;
    private Boolean reverse;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        reverse = reverseInput.get();
        nTypes = parameterizationInput.get().getNTypes();
    }

    @Override
    public double proposal() {


        RealParameter rates = parameterInput.get();
        RealParameter spikes = spikesInput.get();


        // Sample an index
        final int index = Randomizer.nextInt(rates.getDimension());
        final double scale = getScaler(index, 0);


        double r = rates.getValue(index);
        double r_;
        if (!reverse) {
            r_ = r * scale;
        } else r_ = r / scale;

        rates.setValue(index, r_);


        for (int i = 0; i < nTypes; i++) {
            double s = spikes.getValue(index * nTypes + i);

            double s_;
            if (!reverse) {
                s_ = s / scale;
            } else s_ = s * scale;

            spikes.setValue(index, s_);
        }

        return (1 - nTypes) * Math.log(scale);
    }

}
