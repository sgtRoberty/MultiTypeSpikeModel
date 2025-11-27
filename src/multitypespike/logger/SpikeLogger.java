package multitypespike.logger;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import multitypespike.clockmodel.PunctuatedClockModel;


@Description("Logs spike values per branch scaled by the spike mean parameter")
public class SpikeLogger extends CalculationNode implements Function, Loggable {
    final public Input<PunctuatedClockModel> clockModelInput =
            new Input<>("clock", "Punctuated clock model input", Validate.REQUIRED);



    @Override
    public void initAndValidate() {
    }

    @Override
    public int getDimension() {
        return clockModelInput.get().getSpikeDimension();
    }

    @Override
    public double getArrayValue() {
        return getArrayValue(0);
    }


    @Override
    public double getArrayValue(int dim) {
        return clockModelInput.get().getSpikeSize(dim);
    }


    /**
     * Loggable interface implementation follows
     */
    @Override
    public void init(PrintStream out) {
        for (int i = 0; i < this.getDimension(); i ++) {
            String id = this.getID();
            if (id == null || id.equals("")) id = "scaledSpike";
            out.print(id + "." + i + "\t");
        }
    }

    @Override
    public void log(long sampleNr, PrintStream out) {
        for (int i = 0; i < this.getDimension(); i ++) {
            out.print(this.getArrayValue(i) + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

}
