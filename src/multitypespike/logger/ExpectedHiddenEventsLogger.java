package multitypespike.logger;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import multitypespike.distribution.BranchSpikePrior;

import java.io.PrintStream;

@Description("Logs expected number of hidden speciation events per branch")
public class ExpectedHiddenEventsLogger extends CalculationNode implements Function, Loggable {
    final public Input<BranchSpikePrior> branchSpikePriorInput =
            new Input<>("branchSpikePrior", "Branch spike prior input", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public void log(long sample, PrintStream out) {
        for (int i = 0; i < this.getDimension(); i ++) {
            out.print(this.getArrayValue(i) + "\t");
        }
    }


    @Override
    public int getDimension() {
        return branchSpikePriorInput.get().nTypes * branchSpikePriorInput.get().nodeCount;
    }

    @Override
    public double getArrayValue(int dim) {
        return branchSpikePriorInput.get().getExpectedHiddenEvents(dim);
    }

    @Override
    public void init(PrintStream out) {
        for (int i = 0; i < this.getDimension(); i ++) {
            String id = this.getID();
            if (id == null || id.equals("")) id = "expectedHiddenEvents";
            out.print(id + "." + i + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

}