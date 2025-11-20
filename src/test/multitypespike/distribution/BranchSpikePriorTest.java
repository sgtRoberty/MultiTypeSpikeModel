package test.multitypespike.distribution;

import bdmmprime.distribution.BirthDeathMigrationDistribution;
import bdmmprime.parameterization.*;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.tree.*;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.evolution.tree.coalescent.RandomTree;
import beast.base.inference.parameter.RealParameter;
import gammaspike.distribution.StumpedTreePrior;
import gammaspike.tree.Stubs;
import multitypespike.distribution.BranchSpikePrior;
import multitypespike.distribution.MultiTypeHiddenEvents;
import multitypespike.distribution.PiSystem;
import org.apache.commons.math3.ode.ContinuousOutputModel;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.assertEquals;

public class BranchSpikePriorTest {


    public static void main(String[] args) {

        String newick = "((2:22.625392002719074,1:1.575392002719063)3:12.930687094526313,0:6.03107909724541)4:0.0";
        Tree tree = new TreeParser(newick, false,false,true,0);
        RealParameter origin = new RealParameter("100.0");
        RealParameter startTypePriorProbs = new RealParameter("0.3333333333333333 0.3333333333333333 0.3333333333333333");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(3),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.1 0.1 0.1"), 3),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.1 0.1 0.1"), 3),
                "birthRateAmongDemes", new SkylineMatrixParameter(
                        null,
                        new RealParameter("0.0 0.0 0.0 0.0 0.0 0.0"), 3),
                "migrationRate", new SkylineMatrixParameter(
                        null,
                        new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0"), 3),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.1 0.1 0.1"), 3),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.0 0.0 0.0"), 3),
                "rhoSampling", new TimedParameter(
                        origin,
                        new RealParameter("0.0 0.0 0.0"), 3));
//              <typeTraitSet id="typeTraitSet.t:amniotes" spec="bdmmprime.util.InitializedTraitSet" traitname="type" value="Gephyrostegus_bohemicus_Eu_309=Eu,Seymouria_spp_NAm_288.51=NAm,Tseajaia_campi_NAm_296.9=NAm,Limnoscelis_paludis_NAm_298.61=NAm,Oedaleops_campi_NAm_296.21=NAm,Eothyris_parkeyi_NAm_285.05=NAm,Vaughnictis_smithae_NAm_294.5=NAm,Eocasea_martini_NAm_300.45=NAm,Euromycter_rutena_Eu_288.51=Eu,Casea_broilii_NAm_283=NAm,Ennatosaurus_tecton_Eu_269.025=Eu,Echinerpeton_intermedium_NAm_309=NAm,Archaeothyris_florensis_NAm_309=NAm,Varanosaurus_acutirostris_NAm_286.8=NAm,Ophiacodon_spp_NAm_294.5=NAm,Cutleria_wilmarthi_NAm_296.21=NAm,Pantelosaurus_saxonicus_Eu_296.21=Eu,Haptodus_garnettensis_NAm_305.35=NAm,Secodontosaurus_obtusidens_NAm_286.8=NAm,Dimetrodon_spp_NAm_286.8=NAm,Ianthasaurus_hardestiorum_NAm_305.35=NAm,Edaphosaurus_boanerges_NAm_286.8=NAm,Hylonomus_lyelli_NAm_316.65=NAm,Anthracodromeus_longipes_NAm_309=NAm,Paleothyris_acadiana_NAm_309=NAm,Protorothyris_archeri_NAm_294.5=NAm,Petrolacosaurus_kansensis_NAm_305.35=NAm,Araeoscelis_spp_NAm_286.8=NAm,Orovenator_mayorum_NAm_289=NAm,Archaeovenator_hamiltonensis_NAm_300.45=NAm,Ascendonanus_nestleri_Eu_290=Eu,Aerosaurus_wellesi_NAm_296.21=NAm,Apsisaurus_witteri_NAm_294.5=NAm,Heleosaurus_scholtzi_SA_260.55=SA,Mesenosaurus_romeri_Eu_267.025=Eu,Mycterosaurus_longiceps_NAm_286.8=NAm,Elliotsmithia_longiceps_SA_260.55=SA,Varanops_brevirostris_NAm_285.05=NAm,Varanodon_agilis_NAm_272.5=NAm,Watongia_meieri_NAm_272.5=NAm,Youngina_capensis_SA_256.62=SA,Acerosodontosaurus_piveteaui_SA_260.35=SA,Lanthanolania_ivakhnenkoi_Eu_269.025=Eu,Claudiosaurus_germaini_SA_263.95=SA,Prolacerta_broomi_SA_249.55=SA,Proterosuchus_spp_SA_249.55=SA,Thuringothyris_mahlendorffae_Eu_286.8=Eu,Captorhinus_aguti_NAm_286.8=NAm,Captorhinus_laticeps_NAm_286.8=NAm,Labidosaurikos_meachami_NAm_278.225=NAm,Labidosaurus_hamatus_NAm_278.225=NAm,Protocaptorhinus_pricei_NAm_286.8=NAm,Romeria_spp_NAm_294.5=NAm,Euconcordia_cunninghami_NAm_300.45=NAm,Reiszorhinus_olsoni_NAm_286.8=NAm,Eudibamus_cursoris_Eu_286.8=Eu,Mesosaurus_tenuidens_SA_280.25=SA,Stereosternum_tumidum_SA_280.25=SA,Erpetonyx_arsenaultorum_NAm_301.15=NAm,Belebey_vegrandis_Eu_270.875=Eu,Macroleter_poezicus_Eu_269.025=Eu,Milleretta_rubidgei_SA_256.62=SA,Acleistorhinus_pteroticus_NAm_284.75=NAm,Colobomycter_pholeter_NAm_289=NAm,Nyctiphruretus_acudens_Eu_258.5=Eu,Emeroleter_levis_Eu_258.5=Eu,Deltavjatia_rossicus_Eu_258.5=Eu,Procolophon_trigoniceps_SA_249.55=SA,Candelaria_barbouri_NAm_239.5=NAm,Owenetta_kitchingorum_SA_253.15=SA">

//        typeTraitSet typeTraitSet =

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();

        density.initByName(
                "parameterization", parameterization,
                "startTypePriorProbs", startTypePriorProbs,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false,
                "useAnalyticalSingleTypeSolution", false,
                "storeIntegrationResults", true
        );

        density.calculateLogP();  // Calculate LogP to call integration method
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();

        int nodeNr = 132;

        PiSystem piSystem = new PiSystem(parameterization, tree, p0geResults, 1e-100, 1e-7);
        piSystem.integratePiSingleLineage(startTypePriorProbs.getDoubleValues(), parameterization, 0.0, 1.0);


        int steps = 100;
        double branchLength = 1;
        double stepSize = branchLength / steps;
        double parentTime = 0;
        ContinuousOutputModel model = piSystem.integrationResults[nodeNr];

        for (int i = 0; i < steps; i++) {
            double t = parentTime + (i + 0.5) * stepSize;  // Midpoint of interval
            model.setInterpolatedTime(t);
            double[] state = model.getInterpolatedState();  // [π0, π1, ge0, ge1]
            System.out.printf(" π₀=%.4f π₁=%.4f %n",
                    state[0], state[1]);
        }


        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "bdmDistr", density);

        MultiTypeHiddenEvents multitypeHiddenEvents = new MultiTypeHiddenEvents(nodeNr, parameterization,
                p0geResults, piSystem.getIntegrationResultsForNode(nodeNr) ,1e-100,1e-7);

        double[] hiddenEvents = multitypeHiddenEvents.integrateSingleLineage(new double[2], 0.0, 1.0);

        System.out.println("expected number of type A hidden events = " + hiddenEvents[0]);
        System.out.println("expected number of type B hidden events = " + hiddenEvents[1]);
        System.out.println("sum of expected number of hidden events = " +  (hiddenEvents[0] + hiddenEvents[1]) );
    }




//    public static void main(String[] args) {
//
//        String newick = "(t1:1.0, t2:1.0);";
//        Tree tree = new TreeParser(newick, false,false,true,0);
//        RealParameter origin = new RealParameter("2.5");
//
//        Parameterization parameterization = new CanonicalParameterization();
//        parameterization.initByName(
//                "processLength", origin,
//                "birthRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("3.0"), 1),
//                "deathRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.5"), 1),
//                "samplingRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.0"), 1),
//                "removalProb", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.0"), 1),
//                "rhoSampling", new TimedParameter(
//                        origin,
//                        new RealParameter("0.2"))
//        );
//
//        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
//
//        density.initByName(
//                "parameterization", parameterization,
//                "conditionOnSurvival", false,
//                "tree", tree,
//                "typeLabel", "state",
//                "parallelize", false,
//                "useAnalyticalSingleTypeSolution", false,
//                "storeIntegrationResults", true);
//
//        // nodeNr 0 = t1
//        // nodeNr 1 = t2
//        // nodeNr 2 = root
//
//        int nodeNr = 0;
//        Node node = tree.getNode(nodeNr);
//
//        BranchSpikePrior bsp = new BranchSpikePrior();
//        bsp.initByName("parameterization", parameterization,
//                "tree", tree,
//                "spikeShape", "1.0",
//                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
//                "bdmDistr", density);
//
//        System.out.println("single-type method = " + bsp.getExpNrHiddenEventsForBranch(node));
//
//        Stubs stub = new Stubs();
//        gammaspike.distribution.BranchSpikePrior gamma_bsp = new gammaspike.distribution.BranchSpikePrior();
//        StumpedTreePrior stp = new StumpedTreePrior();
//        stp.initByName("lambda", "3",
//                "samplingProportion", "0.5",
//                "turnover", "0.16666666666667",
//                "tree", tree);
//        stub.initByName("tree", tree, "prior", stp);
//        gamma_bsp.initByName("spikes", "1.0 0.5 0.1 0.2 0.7 0.1", "shape", "1.0", "stubs", stub, "tree", tree);
////        System.out.println(node.getHeight());
////        System.out.println(node.getParent().getHeight());
////        System.out.println("GS Model = " + stp.getMeanStubNumber(node.getHeight(), node.getParent().getHeight()));
//    }



    @Test
    public void BSPMultiTypeSingleTypeComparisonTest() {

        String newick = "(t1[&state=1]:1.5, t2[&state=1]:0.5);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealParameter origin = new RealParameter("2.5");
        RealParameter startTypePriorProbs = new RealParameter("1.0");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(null, new RealParameter("2.0"), 1),
                "deathRate", new SkylineVectorParameter(null, new RealParameter("1.0"), 1),
                "samplingRate", new SkylineVectorParameter(null, new RealParameter("0.5"), 1),
                "removalProb", new SkylineVectorParameter(null, new RealParameter("0.0"), 1)
        );

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "startTypePriorProbs", startTypePriorProbs,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false,
                "useAnalyticalSingleTypeSolution", false,
                "storeIntegrationResults", true
        );
        density.calculateLogP();
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();

        PiSystem piSystem = new PiSystem(parameterization, tree, p0geResults, 1e-100, 1e-7);
        piSystem.integratePi(startTypePriorProbs.getDoubleValues(), parameterization, 0.0);

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "useAnalyticalSingleTypeSolution", false,
                "bdmDistr", density
        );

        double multiTypeResult = bsp.calculateLogP();

        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "useAnalyticalSingleTypeSolution", true,
                "bdmDistr", density
        );
        double singleTypeResult = bsp.calculateLogP();

        double tolerance = 1e-3;

        assertEquals("Mismatch:"  + "single-type result = " + singleTypeResult
                        + " does not match multi-type result = " +  multiTypeResult,
                singleTypeResult, multiTypeResult, tolerance);

    }

    @Test
    public void BSPMultiTypeSingleTypeComparisonSkylineTest() {

        String newick = "(t1[&state=1]:1.5, t2[&state=1]:0.5);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealParameter origin = new RealParameter("2.5");
        RealParameter startTypePriorProbs = new RealParameter("1.0");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(new RealParameter("1.0"), new RealParameter("2.0 0.5"), 1),
                "deathRate", new SkylineVectorParameter(new RealParameter("1.5"), new RealParameter("1.0 0.2"), 1),
                "samplingRate", new SkylineVectorParameter(new RealParameter("2.0"), new RealParameter("0.5 1.8"), 1),
                "removalProb", new SkylineVectorParameter(null, new RealParameter("0.0"), 1)
        );

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "startTypePriorProbs", startTypePriorProbs,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false,
                "useAnalyticalSingleTypeSolution", false,
                "storeIntegrationResults", true
        );
        density.calculateLogP();
        ContinuousOutputModel[] p0geResults = density.getIntegrationResults();

        PiSystem piSystem = new PiSystem(parameterization, tree, p0geResults, 1e-100, 1e-7);
        piSystem.integratePi(startTypePriorProbs.getDoubleValues(), parameterization, 0.0);

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "useAnalyticalSingleTypeSolution", false,
                "bdmDistr", density
        );

        double multiTypeResult = bsp.calculateLogP();

        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "startTypePriorProbs", startTypePriorProbs,
                "useAnalyticalSingleTypeSolution", true,
                "bdmDistr", density
        );
        double singleTypeResult = bsp.calculateLogP();

        double tolerance = 1e-3;

        assertEquals("Mismatch:"  + "single-type result = " + singleTypeResult
                        + " does not match multi-type result = " +  multiTypeResult,
                singleTypeResult, multiTypeResult, tolerance);

    }

    /**
     * Single-type test for expected number of hidden events on a branch with no rate shifts
     * Compares with simulations of birth-death trajectories performed in R
     */
    @Test
    public void singleTypeExpectedHiddenEventsTest() {
        String newick = "(t1:1.0, t2:1.0);";
        Tree tree = new TreeParser(newick, false, false, true, 0);
        RealParameter origin = new RealParameter("2.5");

        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "processLength", origin,
                "birthRate", new SkylineVectorParameter(null, new RealParameter("3.0"), 1),
                "deathRate", new SkylineVectorParameter(null, new RealParameter("0.5"), 1),
                "samplingRate", new SkylineVectorParameter(null, new RealParameter("0.0"), 1),
                "removalProb", new SkylineVectorParameter(null, new RealParameter("0.0"), 1),
                "rhoSampling", new TimedParameter(origin, new RealParameter("0.2"))
        );

        BirthDeathMigrationDistribution density = new BirthDeathMigrationDistribution();
        density.initByName(
                "parameterization", parameterization,
                "conditionOnSurvival", false,
                "tree", tree,
                "typeLabel", "state",
                "parallelize", false,
                "useAnalyticalSingleTypeSolution", false,
                "storeIntegrationResults", true
        );

        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName(
                "parameterization", parameterization,
                "tree", tree,
                "spikeShape", "1.0",
                "spikes", "1.0 0.5 0.1 0.2 0.7 0.1",
                "bdmDistr", density
        );

        // Get expected hidden events for node 0 (t1)
        Node node = tree.getNode(0);
        double expectedHiddenEvents = bsp.getExpNrHiddenEventsForBranch(node);

        // Compare to R simulation result (hiddenEventsSim.R)
        double expectedFromSimulation = 3.39054;
        double tolerance = 5e-3;

        assertEquals("Expected hidden events should match simulation", expectedFromSimulation, expectedHiddenEvents, tolerance);
    }

    /**
     * Single-type test for expected number of hidden events on a branch with no rate shifts
     * Compares with the output of GammaSpike Model - 01/03/2025
     */
    @Test
    public void noRateShiftsSingleTypeCaseTest() {

        String newick = "((0:1.0,1:1.0)4:1.0,(2:1.0,3:1.0)5:0.5)6:0.0;";
        TreeParser treeParser = new TreeParser(newick, false, false, false, 0);
        Tree tree = treeParser;

        RealParameter originParam = new RealParameter("2.0");
        Parameterization parameterization = new CanonicalParameterization();
        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.75"), 1),
                "deathRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.3"), 1),
                "samplingRate", new SkylineVectorParameter(
                        null,
                        new RealParameter("0.1"), 1),
                "removalProb", new SkylineVectorParameter(
                        null,
                        new RealParameter("0"), 1),
                "rhoSampling", new TimedParameter(
                        originParam,
                        new RealParameter("1.0"))
        );
        BranchSpikePrior bsp = new BranchSpikePrior();
        bsp.initByName("parameterization", parameterization, "tree", tree, "spikeShape", "1.0", "spikes", "1.0 0.5 0.1 0.2 0.7 0.1");
        Node node = tree.getNode(5);

        /*
        Stubs stub = new Stubs();
        gammaspike.distribution.BranchSpikePrior gamma_bsp = new gammaspike.distribution.BranchSpikePrior();
        StumpedTreePrior stp = new StumpedTreePrior();
        stp.initByName("lambda", "0.75", "r0", "2.5", "samplingProportion", "0.25", "tree", tree);
        stub.initByName("tree", tree, "prior", stp);
        gamma_bsp.initByName("spikes", "1.0 0.5 0.1 0.2 0.7 0.1", "shape", "1.0", "stubs", stub, "tree", tree);

        stp.getMeanStubNumber(node.getHeight(), node.getParent().getHeight())  = 0.186082405828208.
        gamma_bsp.calculateLogP() = -3.4385926164063445.
        */

        assertEquals(0.186082405828208, bsp.getExpNrHiddenEventsForBranch(node), 1e-5);
//        assertEquals(-3.4385926164063445, bsp.calculateLogP(), 1e-5); // MTS and GS logPs do not match due to differences in handling of the root and sampled ancestors

    }


    /**
     * Test to check the logic of the branch traversal process used for calculating the
     * expected number of hidden events along a branch (by computing the lengths of branches instead of
     * the expected number of hidden events)
     */
    @Test
    public void branchTraversalTest() {
        RealParameter originParam = new RealParameter("2.0");
        Parameterization parameterization = new CanonicalParameterization();

        parameterization.initByName(
                "typeSet", new TypeSet(1),
                "processLength", originParam,
                "birthRate", new SkylineVectorParameter(
                        new RealParameter("0.1 0.2 0.3 0.4 0.5 0.6"),
                        new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0"), 1),
                "deathRate", new SkylineVectorParameter(
                        new RealParameter("0.7 0.8 0.9 1.0 1.1 1.2"),
                        new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0"), 1),
                "samplingRate", new SkylineVectorParameter(
                        new RealParameter("1.3 1.4 1.5 1.6 1.7 1.8"),
                        new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0"), 1),
                "removalProb", new SkylineVectorParameter(
                        new RealParameter("1.9 2.0 2.1 2.2 2.5 3.0"),
                        new RealParameter("1.0 1.0 1.0 1.0 1.0 1.0 1.0"), 1),
                "rhoSampling", new TimedParameter(
                        originParam,
                        new RealParameter("0.0")));

        // Code to generate random trees for testing

        List<Sequence> seqList = new ArrayList<Sequence>();
        int Nleaves = 10;
        int Niter = 10000;

        for (int i = 0; i < Nleaves; i++) {
            String taxonID = "t " + i;
            seqList.add(new Sequence(taxonID, "?"));
        }
        Alignment alignment = new Alignment(seqList, "nucleotide");

        // Population model
        ConstantPopulation populationModel = new ConstantPopulation();
        populationModel.initByName("popSize", new RealParameter("1.0"));

        // Create RandomTree
        RandomTree randomTree = new RandomTree();
        for (int k = 0; k < Niter; k++) {
            randomTree.initByName("taxa", alignment, "populationModel", populationModel);

            for (int n = 1; n < randomTree.getNodeCount() - 1; n++) {
                Node node = randomTree.getNode(n);

                double trueBranchTime = parameterization.getNodeTime(node, 0)
                        - parameterization.getNodeTime(node.getParent(), 0);

                double[] intervalEndTimes = parameterization.getIntervalEndTimes();
                double branchTime = 0;
                int nodeIndex = parameterization.getNodeIntervalIndex(node, 0);
                int parentIndex = parameterization.getNodeIntervalIndex(node.getParent(), 0);
                double t0 = parameterization.getNodeTime(node.getParent(), 0);
                double T = parameterization.getNodeTime(node, 0);

                if (nodeIndex == parentIndex) branchTime = (T - t0);
                else {
                    for (int i = parentIndex; i <= nodeIndex - 1; i++) {

                        double t1 = intervalEndTimes[i];
                        branchTime += (t1 - t0);
                        t0 = t1;
                    }

                    branchTime += (T - t0);
                }

                assertEquals(branchTime, trueBranchTime, 1e-9);
            }
        }
    }




//    public static void main(String[] args) {
//
//        String newick = "((0:1.0,1:1.0)4:1.0,(2:1.0,3:1.0)5:0.5)6:0.0;";
//        TreeParser treeParser = new TreeParser(newick, false, false, false, 0);
//        Tree tree = treeParser;
//
//        RealParameter originParam = new RealParameter("2.0");
//
//        Parameterization parameterization = new CanonicalParameterization();
//        parameterization.initByName(
//                "typeSet", new TypeSet(1),
//                "processLength", originParam,
//                "birthRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.75"), 1),
//                "deathRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.3"), 1),
//                "samplingRate", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0.1"), 1),
//                "removalProb", new SkylineVectorParameter(
//                        null,
//                        new RealParameter("0"), 1),
//                "rhoSampling", new TimedParameter(
//                        originParam,
//                        new RealParameter("1.0"))
//        );
//
//        Node node = tree.getNode(5);
//        BranchSpikePrior bsp = new BranchSpikePrior();
//        bsp.initByName("parameterization", parameterization, "tree", tree, "spikeShape", "1.0", "spikes", "1.0 0.5 0.1 0.2 0.7 0.1");
//        System.out.println(bsp.getExpNrHiddenEventsForBranch(node));
//        System.out.println(bsp.calculateLogP());

//        System.out.println(parameterization.getNodeTime(node,0));
//        System.out.println(parameterization.getNodeTime(node.getParent(),0));
//
//        System.out.println(node.getParent().getHeight());
//        System.out.println(node.getHeight());



//        Stubs stub = new Stubs();
//        gammaspike.distribution.BranchSpikePrior gamma_bsp = new gammaspike.distribution.BranchSpikePrior();
//        StumpedTreePrior stp = new StumpedTreePrior();
//        stp.initByName("lambda", "0.75", "r0", "2.5", "samplingProportion", "0.25", "tree", tree);
//        stub.initByName("tree", tree, "prior", stp);
//        gamma_bsp.initByName("spikes", "1.0 0.5 0.1 0.2 0.7 0.1", "shape", "1.0", "stubs", stub, "tree", tree);
//        System.out.println(stp.getMeanStubNumber(node.getHeight(), node.getParent().getHeight()));
//        System.out.println(gamma_bsp.calculateLogP());


//}


}

