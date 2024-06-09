package beast.evolution.substitutionmodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.ReadCounts;
import beast.evolution.variantsinfo.vcfentry.CandidateAltNuc;
import com.google.common.primitives.Chars;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Each extended class of ScsSubstitutionModelBase should implement:
 * canHandleDataType(ReadCounts readCountsDataType)
 * setupRelativeRates()
 * setupRateMatrix()
 */
@Description(value = "Base implementation of a substitution model for single cell DNA sequencing read counts data")
public abstract class ScsSubstitutionModelBase extends GeneralSubstitutionModel {


    //***********************************************
    //*                  Variables                  *
    //***********************************************

    protected final static String INVALID_TRANSITION = "Transition of genotypes from the first to the second is illegal: ";

    public enum EvolutionaryEventType {

        // single mutation (0/0 -> 0/1, 000 -> 001)
        SINGLE_MUTATION("SM"),

        // coincident homozygous double mutation (0/0 -> 1/1, 000 -> 011)
        COIN_HOMO_DOUBLE_MUTATION("CHoDM"),

        // coincident heterozygous double mutation (0/0 -> 1/1', 000 -> 011')
        COIN_HETERO_DOUBLE_MUTATION("CHeDM"),

        // coincident homozygous triple mutation (000 -> 111)
        COIN_HOMO_TRIPLE_MUTATION("CHoTM"),

        // hybrid homo- and heterozygous coincident triple mutation (000 -> 111')
        HYBRID_HOMO_HETERO_TRIPLE_MUTATION("HoHeCTM"),

        // heterozygous coincident triple mutation (000 -> 11'1'')
        HETERO_COIN_TRIPLE_MUTATION("HeCTM"),

        // single back mutation (0/1 -> 0/0, 1/1 -> 0/1, 1/1' -> 0/1, 001 -> 000, 011 -> 001, 011' -> 001, 111 -> 011, 111' -> 011, 111' -> 011', 11'1'' -> 011, 11'1'' -> 011')
        SINGLE_BACK_MUTATION("SB"),

        // coincident double back mutation (1/1 -> 0/0, 1/1' -> 0/0, 011 -> 000, 011' -> 000, 111 -> 001, 111' -> 001, 11'1'' -> 001)
        COIN_DOUBLE_BACK_MUTATION("CDB"),

        // double back mutation (111 -> 000, 111' -> 000, 11'1'' -> 000)
        TRIPLE_BACK_MUTATION("TB"),

        // hybrid homo- and heterozygous single substitute and back mutation (111 -> 011')
        HYBRID_HOMO_HETERO_SINGLE_SUBST_BACK_MUTATION("HoHeSSB"),

        // homozygous single mutation addition (0/1 -> 1/1, 001 -> 011, 011 -> 111)
        HOMO_SINGLE_MUTATION_ADDITION("HoSMA"),

        // homozygous double mutation addition (001 -> 111)
        HOMO_DOUBLE_MUTATION_ADDITION("HoDMA"),

        // heterozygous single mutation addition (0/1 -> 1/1', 001 -> 011', 011' -> 11'1'')
        HETERO_SINGLE_MUTATION_ADDITION("HeSMA"),

        // heterozygous double mutation addition (001 -> 11'1'')
        HETERO_DOUBLE_MUTATION_ADDITION("HeDMA"),

        // hybrid homo- and heterozygous single mutation addition (011 -> 111', 011' -> 111')
        HYBRID_HOMO_HETERO_SINGLE_MUTATION_ADDITION("HoHeSMA"),

        // hybrid homo- and heterozygous double mutation addition (001 -> 111')
        HYBRID_HOMO_HETERO_DOUBLE_MUTATION_ADDITION("HoHeDMA"),

        // homozygous substitute single mutation (1/1' -> 1/1, 111' -> 111)
        HOMO_SUBST_SINGLE_MUTATION("HoSubSM"),

        // homozygous substitute single mutation (11'1'' -> 111)
        HOMO_SUBST_DOUBLE_MUTATION("HoSubDM"),

        // heterozygous substitute single mutation (1/1 -> 1/1', 011 -> 011', 111' -> 11'1'')
        HETERO_SUBST_SINGLE_MUTATION("HeSubSM"),

        // heterozygous substitute double mutation (111 -> 11'1'')
        HETERO_SUBST_DOUBLE_MUTATION("HeSubDM"),

        // hybrid homo- and heterozygous substitute single mutation (011' -> 011, 111 -> 111', 11'1'' -> 111')
        HYBRID_HOMO_HETERO_SUBST_SINGLE_MUTATION("HoHeSubSM"),

        // homozygous single substitute and mutation (011' -> 111)
        HOMO_SINGLE_SUBST_AND_MUTATION("HoSSubM"),

        // heterozygous single substitute and mutation (011 -> 11'1'')
        HETERO_SINGLE_SUBST_AND_MUTATION("HeSSubM"),

        // single deletion loss of heterozygosity (0/1 -> 0/- or 1/-, 1/1' -> 1/-, 001 -> 0/0, 111' -> 1/1)
        SINGLE_DELETION_LOH("SDLOH"),

        // loss of heterozygosity by double deletion (001 -> 0/-, 001 -> 1/-, 011 -> 0/-, 011 -> 1/-, 011' -> 0/-, 011' -> 1/-)
        LOSS_OF_HETEROZYGOSITY_DOUBLE_DELETION("LOHDD"),

        // single deletion not LOH (0/0 -> 0/-, 1/1 -> 1/-, 000 -> 0/0, 001 -> 0/1, 011 -> 0/1, 011 -> 1/1, 011' -> 0/1, 011' -> 1/1', 111 -> 1/1, 111' -> 1/1', 11'1'' -> 1/1')
        SINGLE_DELETION_NOT_LOH("SDNLOH"),

        // coincident deletion and mutation (0/0 -> 1/-, 000 -> 0/1)
        COIN_DELETION_AND_MUTATION("CDM"),

        // single deletion and homozygous single mutation addition (001 -> 1/1)
        SINGLE_DELETION_AND_HOMO_SINGLE_MUTATION_ADDITION("SDHoSMA"),

        // single deletion and heterozygous single mutation addition (001 -> 1/1')
        SINGLE_DELETION_AND_HETERO_SINGLE_MUTATION_ADDITION("SDHeSMA"),

        // single deletion and homozygous coincident double mutation (000 -> 1/1)
        SINGLE_DELETION_AND_HOMO_COIN_DOUBLE_MUTATION("SDHoCDM"),

        // single deletion and heterozygous coincident double mutation (000 -> 1/1')
        SINGLE_DELETION_AND_HETERO_COIN_DOUBLE_MUTATION("SDHeCDM"),

        // coincident deletion and back mutation (1/1 -> 0/-, 1/1' -> 0/-, 011 -> 0/0, 011' -> 0/0, 111 -> 0/1, 111' -> 0/1, 11'1'' -> 0/1)
        COIN_DELETION_AND_BACK_MUTATION("CDBM"),

        // single deletion and double back mutation (111 -> 0/0, 111' -> 0/0, 11'1'' -> 0/0)
        SINGLE_DELETION_AND_DOUBLE_BACK_MUTATION("SDDB"),

        // single deletion mutation addition (0/- -> 1/-)
        SINGLE_DELETION_MUTATION_ADDITION("SDMA"),

        // single deletion and homozygous substitute single mutation (011' -> 1/1, 11'1'' -> 1/1)
        SINGLE_DELETION_AND_HOMO_SUBST_SINGLE_MUTATION("SDHoSubSA"),

        // single deletion and heterozygous substitute single mutation (011 -> 1/1', 111 -> 1/1')
        SINGLE_DELETION_AND_HETERO_SUBST_SINGLE_MUTATION("SDHeSubSA"),

        // single deletion back mutation addition (1/- -> 0/-)
        SINGLE_DELETION_BACK_MUTATION_ADDITION("SDBA"),

        // single deletion addition (0/- -> -, 1/- -> -)
        SINGLE_DELETION_ADDITION("SDA"),

        // coincident double deletion (0/0 -> -, 0/1 -> -, 1/1 -> -, 1/1' -> -, 000 -> 0/-, 111 -> 1/-, 111' -> 1/-, 11'1'' -> 1/-)
        COIN_DOUBLE_DELETION("CDD"),

        // double deletion and single mutation (000 -> 1/-)
        DOUBLE_DELETION_AND_SINGLE_MUTATION("DDSM"),

        // double deletion and single back mutation (111 -> 0/-, 111' -> 0/-, 11'1'' -> 0/-)
        DOUBLE_DELETION_AND_SINGLE_BACK_MUTATION("DDSB"),

        // triple deletion (000 -> -, 001 -> -, 011 -> -, 011' -> -, 111 -> -, 111' -> -, 11'1'' -> -)
        TRIPLE_DELETION("TD"),

        // single insertion (0/0 -> 000, 0/1 -> 001, 0/1 -> 011, 1/1 -> 111, 1/1' -> 111', 0/- -> 0/0, 1/- -> 1/1)
        SINGLE_INSERTION("SI"),

        // double insertion (0/- -> 000, 1/- -> 111)
        DOUBLE_INSERTION("DI"),

        // single insertion and mutation (0/0 -> 001, 0/1 -> 011', 0/0 -> 011, 0/0 -> 011', 0/0 -> 111, 0/0 -> 11'1'', 0/0 -> 111', 0/- -> 0/1, 0/- -> 1/1, 0/- -> 1/1')
        SINGLE_INSERTION_AND_MUTATION("SIM"),

        // double insertion and mutation (0/- -> 001, 0/- -> 011, 0/- -> 012, 0/- -> 111, 0/- -> 112, 0/- -> 123)
        DOUBLE_INSERTION_AND_MUTATION("DIM"),

        // single insertion and back mutation (0/1 -> 000, 1/1 -> 000, 1/1 -> 001, 1/1 -> 011, 1/1' -> 000, 1/1' -> 001, 1/1' -> 011, 1/1' -> 011', 1/- -> 0/0, 1/- -> 0/1)
        SINGLE_INSERTION_AND_BACK_MUTATION("SIB"),

        // double insertion and back mutation (1/- -> 000, 1/- -> 001, 1/- -> 011)
        DOUBLE_INSERTION_AND_BACK_MUTATION("DIB"),

        // single insertion and mutation addition (0/1 -> 111, 0/1 -> 111', 0/1 -> 11'1'')
        SINGLE_INSERTION_AND_MUTATION_ADDITION("SIMA"),

        // single insertion, back and substitute mutation (1/1 -> 011', 1/- -> 1/1')
        SINGLE_INSERTION_AND_BACK_AND_SUBST_MUTATION("SIBSubM"),

        // double insertion, back and substitute mutation (1/- -> 011')
        DOUBLE_INSERTION_AND_BACK_AND_SUBST_MUTATION("DIBSubM"),

        // single insertion and substitute mutation (1/1 -> 111', 1/1 -> 11'1'', 1/1' -> 111, 1/1' -> 11'1'')
        SINGLE_INSERTION_AND_SUBST_MUTATION("SISubM"),

        // double insertion and substitute mutation (1/- -> 111', 1/- -> 11'1'')
        DOUBLE_INSERTION_AND_SUBST_MUTATION("DISubM");

        final String desc;

        EvolutionaryEventType(String s) {
            desc = s;
        }

        @Override
        public String toString() {
            return desc;
        }

        public static String arrayToString(final EvolutionaryEventType[] arr) {
            return arrayToString(arr, "");
        } // arrayToString

        public static String arrayToString(final EvolutionaryEventType[] arr, final String separator) {
            StringBuilder stringBuilder = new StringBuilder();

            for (int index = 0; index < arr.length; index++) {
                stringBuilder.append(arr[index].toString());

                if (index < arr.length - 1)
                    stringBuilder.append(separator);
            }

            return stringBuilder.toString();
        } // arrayToString

        /**
         * Whether the evolutionary event violates infinite-sites assumption or not?
         *
         * @return yes or no
         */
        public boolean violateISA() {
            return this != EvolutionaryEventType.SINGLE_MUTATION;
        }

        public int compare(@NotNull EvolutionaryEventType e) {
            return this.desc.compareTo(e.desc);
        }

        public static class EventsComparator implements Comparator<EvolutionaryEventType> {

            @Override
            public int compare(EvolutionaryEventType e1, EvolutionaryEventType e2) {
                return e1.compare(e2);
            }

        }

    }

    /**
     * genotype of the root
     */
    protected int rootGenotype;

    /**
     * genotype of the constant site
     */
    protected int constGenotype;

    /**
     * the number of alternative alleles each genotype contains
     */
    protected int[] nrOfAltAlleles;

    /**
     * corresponding number of existing alleles for each genotype
     */
    protected int[] nrOfExistingAlleles;

    /**
     * the number of existing alleles this model studies
     */
    protected int[] modeledAlleles;

    /**
     * ternary codes corresponding to each genotype
     * 0: homogeneous reference
     * 1: heterogeneous alternative
     * 2: homogeneous alternative
     * 3: missing data
     */
    protected int[] ternaryCodes;


    //***********************************************
    //*                   Methods                   *
    //***********************************************

    /**
     * constructor for testing purpose
     */
    public ScsSubstitutionModelBase() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
        try {
            ratesInput.setValue(null, this);
        } catch (Exception e) {
            throw new IllegalArgumentException(e.getMessage());
        }

        frequenciesInput.setRule(Input.Validate.OPTIONAL);
        try {
            frequenciesInput.setValue(null, this);
        } catch (Exception e) {
            throw new IllegalArgumentException(e.getMessage());
        }
    } // constructor

    /**
     * ratesInput and frequenciesInput should not be specified in the configuration xml document as
     * they are not applied to the derived classes of ScsSubstitutionModelBase
     */
    @Override
    public void initAndValidate() {
        if (ratesInput.get() != null) {
            throw new IllegalArgumentException("the rates attribute should not be used for the selected substitution model (" + this.getClass().getName() + ")");
        }
        if (frequenciesInput.get() != null) {
            throw new IllegalArgumentException("the frequencies attribute should not be used for the selected substitution model (" + this.getClass().getName() + ")");
        }
    } // initAndValidate

    /**
     * @param parameterInput a parameter of type Input<Function>
     * @return either the input value or default 0.0
     */
    protected RealParameter getParameter(Input<RealParameter> parameterInput) {
        if (parameterInput.get() != null) {
            return parameterInput.get();
        }
        return new RealParameter("0.0");
    } // getParameter

    @Override
    protected void setupRateMatrix() {
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = 0.0;
            for (int j = 0; j < i; j++) {
                rateMatrix[i][j] = relativeRates[i * (nrOfStates - 1) + j];
            }
            for (int j = i + 1; j < nrOfStates; j++) {
                rateMatrix[i][j] = relativeRates[i * (nrOfStates - 1) + j - 1];
            }
        }

        // set up diagonal
        for (int i = 0; i < nrOfStates; i++) {
            double sum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    sum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -sum;
        }

        // normalise rate matrix to one expected substitution per unit time
        double subst = 0.0;
        for (int i = 0; i < nrOfStates; i++)
            subst += -rateMatrix[i][i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / subst;
            }
        }
    } // setupRateMatrix

    @Override
    protected abstract void setupRelativeRates();

    public abstract String getAlphabetGenotype(final int index);

    public abstract boolean canHandleDataType(ReadCounts readCountsDataType);

    /**
     * Get the evolutionary events during the process of a parent evolving to a child.
     *
     * @param parentGenotype genotype of parent
     * @param childGenotype  genotype of a child
     * @return an array of {@link EvolutionaryEventType}
     */
    public abstract EvolutionaryEventType[] getEvolutionaryEvents(final int parentGenotype, final int childGenotype);

    /**
     * Get genotype compatible with VCF standard.
     *
     * @param genotype      genotype code
     * @param locusAltNucs  locus-wise alternative nucleotides
     * @param cellAltNucs   cell-wise alternative nucleotides in descending order according to read counts
     * @param missingAllele character to represent a missing allele
     * @return adjusted genotype
     */
    public String getGenotypeForVCF(
            int genotype,
            List<CandidateAltNuc> locusAltNucs,
            char[] cellAltNucs,
            char missingAllele
    ) {
        List<Character> existingAltNucs = new ArrayList<>();

        final String[] allelesOri = getAlphabetGenotype(genotype).split("/");

        String chrom1Str = String.copyValueOf(
                adaptAllelesChrom(
                        allelesOri[0],
                        locusAltNucs,
                        cellAltNucs,
                        missingAllele,
                        existingAltNucs
                )
        );
        String chrom2Str = String.copyValueOf(
                adaptAllelesChrom(
                        allelesOri[1],
                        locusAltNucs,
                        cellAltNucs,
                        missingAllele,
                        existingAltNucs
                )
        );

        return String.join("/", new String[]{chrom1Str, chrom2Str});
    } // getGenotype

    /**
     * Get alleles compatible with VCF standard per chromosome.
     *
     * @param alleles         alleles on a strand of chromosome
     * @param locusAltNucs    locus-wise alternative nucleotides
     * @param cellAltNucs     cell-wise alternative nucleotides in descending order according to read counts
     * @param missingAllele   character to represent a missing allele
     * @param existingAltNucs existing alt nucs in a cell
     * @return adjusted alleles
     */
    protected char[] adaptAllelesChrom(
            String alleles,
            List<CandidateAltNuc> locusAltNucs,
            char[] cellAltNucs,
            char missingAllele,
            List<Character> existingAltNucs
    ) {
        List<Character> results = new ArrayList<>();

        for (int i = 0; i < alleles.length(); i++) {
            results.add(
                    adaptAllele(
                            alleles.charAt(i),
                            locusAltNucs,
                            cellAltNucs,
                            missingAllele,
                            existingAltNucs
                    )
            );
        }

        return Chars.toArray(results);
    } // adaptAllelesChrom

    /**
     * Adjust an allele.
     * <p/>
     * This function only handles these characters: '0', '1', '2', and '-'.
     * <p/>
     * Any substitution model not matching the above rule should override this function.
     *
     * @param allele          allele to be adjusted
     * @param locusAltNucs    locus-wise alternative nucleotides
     * @param cellAltNucs     cell-wise alternative nucleotides in descending order according to read counts
     * @param missingAllele   character to represent a missing allele
     * @param existingAltNucs existing alt nucs in a cell
     * @return an adjusted allele
     */
    protected char adaptAllele(
            char allele,
            List<CandidateAltNuc> locusAltNucs,
            char[] cellAltNucs,
            char missingAllele,
            List<Character> existingAltNucs
    ) {
        switch (allele) {
            case '0':
                return '0';
            case '-':
            case '.':
                return missingAllele;
            case '1':
                return (cellAltNucs == null || cellAltNucs.length == 0) ?
                        adaptAltAllele(locusAltNucs, 'N', existingAltNucs) :
                        adaptAltAllele(locusAltNucs, cellAltNucs[0], existingAltNucs);
            case '2':
                return (cellAltNucs == null || cellAltNucs.length < 2) ?
                        adaptAltAllele(locusAltNucs, 'N', existingAltNucs) :
                        adaptAltAllele(locusAltNucs, cellAltNucs[1], existingAltNucs);
            default:
                throw new IllegalArgumentException("Error! Unsupported character: " + allele +
                        ". Only '0', '1', '2', '-', and '.' are allowed.");
        }
    } // adaptAllele

    /**
     * Adjust the alternative allele: {@param cellAltNuc}.
     * <p/>
     * e.g. {@param allele} is '1', {@param cellAltNuc} is 'C', {@param locusAltNucs} is {'T', 'C', 'G'},
     * then the index of 'C' in {@param locusAltNucs} plus 1 is returned, i.e., 2.
     *
     * @param locusAltNucs    locus-wise alternative nucleotides
     * @param cellAltNuc      cell-wise alternative nucleotides in descending order according to read counts
     * @param existingAltNucs existing alt nucs in a cell
     * @return an adjusted allele
     */
    protected char adaptAltAllele(
            List<CandidateAltNuc> locusAltNucs,
            char cellAltNuc,
            List<Character> existingAltNucs
    ) {
        if (Character.toUpperCase(cellAltNuc) == 'N')
            return '1';

        for (int i = 0; i < locusAltNucs.size(); i++) {
            if (Character.toUpperCase(locusAltNucs.get(i).getNuc()) == Character.toUpperCase(cellAltNuc)) {
                locusAltNucs.get(i).addCount(1);
                if (!existingAltNucs.contains(cellAltNuc)) {
                    existingAltNucs.add(cellAltNuc);
                    locusAltNucs.get(i).addNumOfCells(1);
                }

                return Character.forDigit(i + 1, 10);
            }
        }

        locusAltNucs.add(
                new CandidateAltNuc(cellAltNuc)
        );
        locusAltNucs.get(locusAltNucs.size() - 1).addCount(1);
        existingAltNucs.add(cellAltNuc);
        locusAltNucs.get(locusAltNucs.size() - 1).addNumOfCells(1);
        return Character.forDigit(locusAltNucs.size(), 10);
    } // adaptAltAllele

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return false;
    }

    public int getRootGenotype() {
        return rootGenotype;
    } // getRootGenotype

    public int getConstGenotype() {
        return constGenotype;
    } // getConstGenotype

    public abstract int getNrOfAlleles(final int genotypeIndex); // getNrOfAlleles

    public int getModeledAllelesSize() {
        return modeledAlleles.length;
    } // getModeledAllelesSize

    public int[] getModeledAlleles() {
        return modeledAlleles;
    } // getModeledAlleles

    public void getModeledAlleles(int[] out) {
        System.arraycopy(modeledAlleles, 0, out, 0, modeledAlleles.length);
    } // getModeledAlleles

    public String getAllGenotypes(String delimiter) {
        String[] genotypes = new String[this.nrOfStates];

        for (int i = 0; i < nrOfStates; i++) {
            genotypes[i] = getAlphabetGenotype(i);
        }

        return String.join(delimiter, genotypes);
    } // getAllGenotypes

    public int getNrOfAltAlleles(int genotypeIndex) {
        return this.nrOfAltAlleles[genotypeIndex];
    } // getNrOfAltAlleles

    public int getTernaryCode(int genotypeIndex) {
        return this.ternaryCodes[genotypeIndex];
    } // getTernaryCode

    public boolean isVariant(int genotypeIndex) {
        return genotypeIndex != this.rootGenotype;
    } // isVariant


} // ScsSubstitutionModelBase
