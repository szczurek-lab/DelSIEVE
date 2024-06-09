package beast.evolution.substitutionmodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.CovSup;
import beast.evolution.datatype.FullSupsCov;
import beast.evolution.datatype.ReadCounts;

import java.lang.reflect.InvocationTargetException;

@Description(value = "A finite mutation and deletion substitution model with 7 genotypes")
public class ScsFiniteMuDelModel extends ScsSubstitutionModelBase {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    public final Input<RealParameter> deletionRateInput = new Input<>("deletionRate", "deletion rate " +
            "(non-negative; default 0); measured relatively to the mutation rate", Input.Validate.REQUIRED);

    protected RealParameter deletionRate;


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    /**
     * constructor for testing purpose
     */
    public ScsFiniteMuDelModel() {
        super();
        initAndValidate();
    } // constructor

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        // model 0, 1, 2 alleles
        modeledAlleles = new int[]{0, 1, 2};

        // fix the number of genotypes
        nrOfStates = 7;

        // fix the genotype of the tree root
        rootGenotype = 0;

        // fix the genotype of the constant site
        constGenotype = 0;

        // initialize the nrOfExistingAlleles
        nrOfExistingAlleles = new int[nrOfStates];
        for (int i = 0; i < nrOfStates; i++) {
            if (i <= 3) {
                nrOfExistingAlleles[i] = 2;
            } else if (i <= 5) {
                nrOfExistingAlleles[i] = 1;
            } else {
                nrOfExistingAlleles[i] = 0;
            }
        }

        nrOfAltAlleles = new int[]{0, 1, 2, 2, 0, 1, 0};

        /*
         * -3: -
         * -2: 1/-
         * -1: 0/-
         * 0: 0/0
         * 1: 0/1
         * 2: 1/1
         * 3: 1/2 (namely 1/1')
         */
        ternaryCodes = new int[]{0, 1, 2, 3, -1, -2, -3};

        try {
            eigenSystem = createEigenSystem();
        } catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException
                 | IllegalArgumentException | InvocationTargetException e) {
            throw new IllegalArgumentException(e.getMessage());
        }
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[nrOfStates * (nrOfStates - 1)];
        storedRelativeRates = new double[nrOfStates * (nrOfStates - 1)]; // maybe not of any use

        // parameter sanity check
        if (getParameter(deletionRateInput).getArrayValue() < 0.0) {
            throw new IllegalArgumentException("deletionRate is out of bound, which should not be smaller than 0 (" +
                    this.getClass().getName() + ")");
        }
        deletionRate = getParameter(deletionRateInput);
    } // initAndValidate

    @Override
    protected void setupRelativeRates() {
        /*
         * 0/0 0/1 1/1 1/2 0/- 1/- -
         */

        /* 0/0 */
        relativeRates[0] = 1.0; // -> 0/1
        relativeRates[1] = 0.0; // -> 1/1
        relativeRates[2] = 0.0; // -> 1/2
        relativeRates[3] = deletionRate.getArrayValue(); // -> 0/-
        relativeRates[4] = 0.0; // -> 1/-
        relativeRates[5] = 0.0; // -> -

        /* 0/1 */
        relativeRates[6] = 1.0 / 6; // -> 0/0
        relativeRates[7] = 1.0 / 6; // -> 1/1
        relativeRates[8] = 1.0 / 3; // -> 1/2
        relativeRates[9] = deletionRate.getArrayValue() / 2; // -> 0/-
        relativeRates[10] = deletionRate.getArrayValue() / 2; // -> 1/-
        relativeRates[11] = 0.0; // -> -

        /* 1/1 */
        relativeRates[12] = 0.0; // -> 0/0
        relativeRates[13] = 1.0 / 3; // -> 0/1
        relativeRates[14] = 2.0 / 3; // -> 1/2
        relativeRates[15] = 0.0; // -> 0/-
        relativeRates[16] = deletionRate.getArrayValue(); // -> 1/-
        relativeRates[17] = 0.0; // -> -

        /* 1/2 */
        relativeRates[18] = 0.0; // -> 0/0
        relativeRates[19] = 1.0 / 3; // -> 0/1
        relativeRates[20] = 1.0 / 3; // -> 1/1
        relativeRates[21] = 0.0; // -> 0/-
        relativeRates[22] = deletionRate.getArrayValue(); // -> 1/-
        relativeRates[23] = 0.0; // -> -

        /* 0/- */
        relativeRates[24] = 0.0; // -> 0/0
        relativeRates[25] = 0.0; // -> 0/1
        relativeRates[26] = 0.0; // -> 1/1
        relativeRates[27] = 0.0; // -> 1/2
        relativeRates[28] = 1.0 / 2; // -> 1/-
        relativeRates[29] = deletionRate.getArrayValue() / 2; // -> -

        /* 1/- */
        relativeRates[30] = 0.0; // -> 0/0
        relativeRates[31] = 0.0; // -> 0/1
        relativeRates[32] = 0.0; // -> 1/1
        relativeRates[33] = 0.0; // -> 1/2
        relativeRates[34] = 1.0 / 6; // -> 0/-
        relativeRates[35] = deletionRate.getArrayValue() / 2; // -> -

        /* - */
        relativeRates[36] = 0.0; // -> 0/0
        relativeRates[37] = 0.0; // -> 0/1
        relativeRates[38] = 0.0; // -> 1/1
        relativeRates[39] = 0.0; // -> 1/2
        relativeRates[40] = 0.0; // -> 0/-
        relativeRates[41] = 0.0; // -> 1/-
    } // setupRelativeRates

    @Override
    public String getAlphabetGenotype(final int index) {
        if (index < 0 || index > nrOfStates - 1) {
            throw new IllegalArgumentException("Index exceeds the boundary (0 - " + (nrOfStates - 1) + ")");
        }

        if (index == 0) {
            return "0/0";
        } else if (index == 1) {
            return "0/1";
        } else if (index == 2) {
            return "1/1";
        } else if (index == 3) {
            return "1/2";
        } else if (index == 4) {
            return "0/.";
        } else if (index == 5) {
            return "1/.";
        } else {
            return "./.";
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
                    case 4:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_DELETION_NOT_LOH};
                    case 5:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.COIN_DELETION_AND_MUTATION};
                    case 6:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.COIN_DOUBLE_DELETION};
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
                    case 4:
                    case 5:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_DELETION_LOH};
                    case 6:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.COIN_DOUBLE_DELETION};
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
                    case 4:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.COIN_DELETION_AND_BACK_MUTATION};
                    case 5:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_DELETION_NOT_LOH};
                    case 6:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.COIN_DOUBLE_DELETION};
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
                    case 4:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.COIN_DELETION_AND_BACK_MUTATION};
                    case 5:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_DELETION_LOH};
                    case 6:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.COIN_DOUBLE_DELETION};
                    default:
                        throw new IllegalArgumentException("Error! Unsupported genotype for the child: " + childGenotype);
                }
            case 4:
                switch (childGenotype) {
                    case 4:
                        return null;
                    case 5:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_DELETION_MUTATION_ADDITION};
                    case 6:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_DELETION_ADDITION};
                    default:
                        throw new IllegalArgumentException("Error! Unsupported genotype for the child: " + childGenotype);
                }
            case 5:
                switch (childGenotype) {
                    case 4:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_DELETION_BACK_MUTATION_ADDITION};
                    case 5:
                        return null;
                    case 6:
                        return new EvolutionaryEventType[]{EvolutionaryEventType.SINGLE_DELETION_ADDITION};
                    default:
                        throw new IllegalArgumentException("Error! Unsupported genotype for the child: " + childGenotype);
                }
            case 6:
                if (childGenotype == 6)
                    return null;

                throw new IllegalArgumentException("Error! Unsupported genotype for the child: " + childGenotype);
            default:
                throw new IllegalArgumentException("Error! Unsupported genotype for the parent: " + parentGenotype);
        }
    } // getEvolutionaryEvents


    //***********************************************
    //*              Getter and Setter              *
    //***********************************************

    public RealParameter getDeletionRate() {
        return deletionRate;
    } // getDeletionRate

    public int[] getNrOfExistingAlleles() {
        return nrOfExistingAlleles;
    } // getNrOfAlleles

    @Override
    public int getNrOfAlleles(final int genotypeIndex) {
        return nrOfExistingAlleles[genotypeIndex];
    } // getNrOfAlleles

} // ScsFiniteMuDelModel
