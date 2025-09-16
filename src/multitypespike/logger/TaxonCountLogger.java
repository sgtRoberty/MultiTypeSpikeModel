package multitypespike.logger;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Tree;

// Based on <GammaSpikeModel>  Copyright (C) <2025>  <Jordan Douglas>

public class TaxonCountLogger extends BEASTObject implements Loggable {

    final public Input<Tree> treeInput = new Input<>("tree", "Tree to be logged", Validate.REQUIRED);

    @Override
    public void init(PrintStream out) {
        out.print(treeInput.get().getID() + ".ntaxa\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        out.print(treeInput.get().getLeafNodeCount() + "\t");
    }

    @Override
    public void close(PrintStream out) {

    }

    @Override
    public void initAndValidate() {

    }

}
