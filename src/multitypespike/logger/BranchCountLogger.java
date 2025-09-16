package multitypespike.logger;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Tree;

public class BranchCountLogger extends BEASTObject implements Loggable {

    final public Input<Tree> treeInput = new Input<>("tree", "Tree to be logged", Validate.REQUIRED);

    @Override
    public void init(PrintStream out) {
        out.print(treeInput.get().getID() + ".nbranches\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.print(treeInput.get().getNodeCount() + "\t");
    }

    @Override
    public void close(PrintStream out) {

    }

    @Override
    public void initAndValidate() {

    }

}
