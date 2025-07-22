package multitypespike.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.kernel.BactrianScaleOperator;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;


/**
 * Bactrian scale operator for spike parameters.
 * Proposes  moves only on valid spikes (excluding zero-valued root and sampled ancestors).
 * The 'parameter' input must reference the spikes parameter.
 * The 'refTree' input is used to ensure moves are only made on valid spikes.
 */

@Description("Bactrian scale operator that proposes multiplicative moves only on valid branches (excluding root, sampled ancestors).")
public class BactrianSpikeScaleOperator extends BactrianScaleOperator {

    public Input<Tree> refTreeInput = new Input<>("refTree", "Reference tree used to determine valid branch spikes",
            Input.Validate.REQUIRED);

    public Input<Integer> maxTriesInput = new Input<>("maxTries", "maximum attempts to find a valid branch", 100);
    int maxTries = maxTriesInput.get();


    @Override
    public double proposal() {
        try {

            double hastingsRatio;

            // not a tree scaler, so scale a parameter
            final boolean scaleAll = scaleAllInput.get();
            final int specifiedDoF = degreesOfFreedomInput.get();
            final boolean scaleAllIndependently = scaleAllIndependentlyInput.get();

            final RealParameter param = (RealParameter)InputUtil.get(parameterInput, this);

            assert param.getLower() != null && param.getUpper() != null;

            final int dim = param.getDimension();

            if (scaleAllIndependently) {
                // update all dimensions independently.
                hastingsRatio = 0;
                final BooleanParameter indicators = indicatorInput.get();
                if (indicators != null) {
                    final int dimCount = indicators.getDimension();
                    final Boolean[] indicator = indicators.getValues();
                    final boolean impliedOne = dimCount == (dim - 1);
                    for (int i = 0; i < dim; i++) {
                        if( (impliedOne && (i == 0 || indicator[i-1])) || (!impliedOne && indicator[i]) )  {
                            final double scaleOne = getScaler(i, param.getValue(i));
                            final double newValue = scaleOne * param.getValue(i);

                            hastingsRatio += Math.log(scaleOne);

                            if (outsideBounds(newValue, param)) {
                                return Double.NEGATIVE_INFINITY;
                            }

                            param.setValue(i, newValue);
                        }
                    }
                }  else {

                    for (int i = 0; i < dim; i++) {

                        final double scaleOne = getScaler(i, param.getValue(i));
                        final double newValue = scaleOne * param.getValue(i);

                        hastingsRatio += Math.log(scaleOne);

                        if( outsideBounds(newValue, param) ) {
                            return Double.NEGATIVE_INFINITY;
                        }

                        param.setValue(i, newValue);
                    }
                }
            } else if (scaleAll) {
                // update all dimensions
                // hasting ratio is dim-2 times of 1dim case. would be nice to have a reference here
                // for the proof. It is supposed to be somewhere in an Alexei/Nicholes article.

                // all Values assumed independent!
                final double scale = getScaler(0, param.getValue(0));
                final int computedDoF = param.scale(scale);
                final int usedDoF = (specifiedDoF > 0) ? specifiedDoF : computedDoF ;
                hastingsRatio = usedDoF * Math.log(scale);
            } else {

                // which position to scale
                final int index;
                final BooleanParameter indicators = indicatorInput.get();
                if (indicators != null) {
                    final int dimCount = indicators.getDimension();
                    final Boolean[] indicator = indicators.getValues();
                    final boolean impliedOne = dimCount == (dim - 1);

                    // available bit locations. there can be hundreds of them. scan list only once.
                    final int[] loc = new int[dimCount + 1];
                    int locIndex = 0;

                    if (impliedOne) {
                        loc[locIndex] = 0;
                        ++locIndex;
                    }
                    for (int i = 0; i < dimCount; i++) {
                        if (indicator[i]) {
                            loc[locIndex] = i + (impliedOne ? 1 : 0);
                            ++locIndex;
                        }
                    }

                    if (locIndex > 0) {
                        final int rand = Randomizer.nextInt(locIndex);
                        index = loc[rand];
                    } else {
                        return Double.NEGATIVE_INFINITY; // no active indicators
                    }

                } else {
                    /**
                     * Changes: Ensures that random index corresponds to valid spike parameters.
                     */
                    // valid spike selection
                    final Tree refTree = refTreeInput.get();
                    index = SpikeOperatorBase.getRandomValidBranch(refTree, parameterInput.get(), maxTries);
                    if (index == -1) {
                        return Double.NEGATIVE_INFINITY;
                }
            }

                final double oldValue = param.getValue(index);

                if (oldValue == 0) {
                    // Error: parameter has value 0 and cannot be scaled
                    return Double.NEGATIVE_INFINITY;
                }

                final double scale = getScaler(index, oldValue);
                hastingsRatio = Math.log(scale);

                final double newValue = scale * oldValue;

                if (outsideBounds(newValue, param)) {
                    // reject out of bounds scales
                    return Double.NEGATIVE_INFINITY;
                }

                param.setValue(index, newValue);
                // provides a hook for subclasses
                //cleanupOperation(newValue, oldValue);
            }

            return hastingsRatio;

        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }
    }


}
