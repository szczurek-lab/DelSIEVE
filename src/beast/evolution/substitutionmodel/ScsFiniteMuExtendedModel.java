package beast.evolution.substitutionmodel;

import beast.evolution.datatype.CovSup;
import beast.evolution.datatype.FullSupsCov;
import beast.evolution.datatype.ReadCounts;

import java.lang.reflect.InvocationTargetException;

public class ScsFiniteMuExtendedModel extends ScsSubstitutionModelBase {

    /**
     * constructor for testing purpose
     */
    public ScsFiniteMuExtendedModel() {
        super();
        initAndValidate();
    } // constructor


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // only 2 alleles modeled
        modeledAlleles = new int[]{2};

        // fix the number of genotypes
        nrOfStates = 4;

        // fix the genotype of the tree root
        rootGenotype = 0;

        // fix the genotype of the constant site
        constGenotype = 0;

        // initialize the nrOfExistingAlleles
        nrOfExistingAlleles = new int[1];
        nrOfExistingAlleles[0] = 2;

        nrOfAltAlleles = new int[]{0, 1, 2, 2};

        /*
         * 0: 0/0
         * 1: 0/1
         * 2: 1/1
         * 3: 1/2 (namely 1/1')
         */
        ternaryCodes = new int[]{0, 1, 2, 3};

        try {
            eigenSystem = createEigenSystem();
        } catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException
                 | IllegalArgumentException | InvocationTargetException e) {
            throw new IllegalArgumentException(e.getMessage());
        }
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[nrOfStates * (nrOfStates - 1)];
        storedRelativeRates = new double[nrOfStates * (nrOfStates - 1)]; // maybe not of any use
    } // initAndValidate

    @Override
    protected void setupRelativeRates() {
        /*
         * 0/0 0/1 1/1 1/2
         */

        /* 0/0 */
        relativeRates[0] = 1.0; // -> 0/1
        relativeRates[1] = 0.0; // -> 1/1
        relativeRates[2] = 0.0; // -> 1/2

        /* 0/1 */
        relativeRates[3] = 1.0 / 6; // -> 0/0
        relativeRates[4] = 1.0 / 6; // -> 1/1
        relativeRates[5] = 1.0 / 3; // -> 1/2

        /* 1/1 */
        relativeRates[6] = 0.0; // -> 0/0
        relativeRates[7] = 1.0 / 3; // -> 0/1
        relativeRates[8] = 2.0 / 3; // -> 1/2

        /* 1/2 */
        relativeRates[9] = 0.0; // -> 0/0
        relativeRates[10] = 1.0 / 3; // -> 0/1
        relativeRates[11] = 1.0 / 3; // -> 1/1
    } // setupRelativeRates

    @Override
    public String getAlphabetGenotype(int index) {
        if (index < 0 || index > nrOfStates - 1) {
            throw new IllegalArgumentException("Index exceeds the boundary (0 - " + (nrOfStates - 1) + ")");
        }

        if (index == 0) {
            return "0/0";
        } else if (index == 1) {
            return "0/1";
        } else if (index == 2) {
            return "1/1";
        } else {
            return "1/2";
        }
    } // getAlphabetGenotype

    @Override
    public boolean canHandleDataType(ReadCounts readCountsDataType) {
        return readCountsDataType instanceof CovSup || readCountsDataType instanceof FullSupsCov;
    } // canHandleDataType

    @Override
    public EvolutionaryEventType[] getEvolutionaryEvents(int parentGenotype, int childGenotype) {
        switch (parentGenotype) {
            case 0:
                switch (childGenotype) {
                    case 0:
                        return null;
                    case 1:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_MUTATION};
                    case 2:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.COIN_HOMO_DOUBLE_MUTATION};
                    case 3:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.COIN_HETERO_DOUBLE_MUTATION};
                    default:
                        throw new IllegalArgumentException("Error! Unsupported genotype for the child: " + childGenotype);
                }
            case 1:
                switch (childGenotype) {
                    case 0:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_BACK_MUTATION};
                    case 1:
                        return null;
                    case 2:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.HOMO_SINGLE_MUTATION_ADDITION};
                    case 3:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.HETERO_SINGLE_MUTATION_ADDITION};
                    default:
                        throw new IllegalArgumentException("Error! Unsupported genotype for the child: " + childGenotype);
                }
            case 2:
                switch (childGenotype) {
                    case 0:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.COIN_DOUBLE_BACK_MUTATION};
                    case 1:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_BACK_MUTATION};
                    case 2:
                        return null;
                    case 3:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.HETERO_SUBST_SINGLE_MUTATION};
                    default:
                        throw new IllegalArgumentException("Error! Unsupported genotype for the child: " + childGenotype);
                }
            case 3:
                switch (childGenotype) {
                    case 0:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.COIN_DOUBLE_BACK_MUTATION};
                    case 1:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_BACK_MUTATION};
                    case 2:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.HOMO_SUBST_SINGLE_MUTATION};
                    case 3:
                        return null;
                    default:
                        throw new IllegalArgumentException("Error! Unsupported genotype for the child: " + childGenotype);
                }
            default:
                throw new IllegalArgumentException("Error! Unsupported genotype for the parent: " + parentGenotype);
        }
    } // getEvolutionaryEvents


    //***********************************************
    //*              Getter and Setter              *
    //***********************************************

    @Override
    public int getNrOfAlleles(int genotypeIndex) {
        return nrOfExistingAlleles[0];
    } // getNrOfAlleles

}
