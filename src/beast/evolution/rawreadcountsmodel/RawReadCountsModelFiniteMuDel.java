package beast.evolution.rawreadcountsmodel;

import beast.core.Description;

import static beast.math.util.MathFunctions.logSumExp;

@Description("Model of nucleotide read counts described by Dirichlet-multinomial distribution compatible with the substitution model of mutations and deletions.")
public class RawReadCountsModelFiniteMuDel extends RawReadCountsModelFiniteMu {


    //**********************************************
    //*             Overridden methods             *
    //**********************************************


    @Override
    protected double computeMixedLikelihoodCore(
            int genotypeIndex,
            int taxonIndex,
            int indexToSeqCovLikelihood,
            int indexToNucReadCountsLikelihood,
            double[] comp
    ) {
        final double theta = adoRate.getValue();

        double[] tmp = comp == null ? this.getAllocatedArr() : comp;

        switch (genotypeIndex) {
            case 0:
                // 0/0
                if (singleADO) {
                    // single ado
                    if (useLogPartials) {
                        tmp[0] = Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = Math.log(theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);

                        return logSumExp(new double[]{tmp[0], tmp[1]});
                    } else {
                        tmp[0] = (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = theta * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);

                        return tmp[0] + tmp[1];
                    }
                } else {
                    // locus ado
                    if (useLogPartials) {
                        tmp[0] = 2 * Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = Math.log(2) + Math.log(theta) + Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = 2 * Math.log(theta) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1], tmp[2]});
                    } else {
                        tmp[0] = Math.pow(1 - theta, 2) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = 2 * theta * (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = Math.pow(theta, 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1] + tmp[2];
                    }
                }
            case 1:
                // 0/1
                if (singleADO) {
                    // single ado
                    if (useLogPartials) {
                        tmp[0] = Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 3) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = Math.log(theta) - Math.log(2) + logSumExp(new double[]{nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood), nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1)}) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);

                        return logSumExp(new double[]{tmp[0], tmp[1]});
                    } else {
                        tmp[0] = (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 3) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = (theta / 2) * (nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1)) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);

                        return tmp[0] + tmp[1];
                    }
                } else {
                    // locus ado
                    if (useLogPartials) {
                        tmp[0] = 2 * Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 3) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = Math.log(theta) + Math.log(1 - theta) + logSumExp(new double[]{nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood), nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1)}) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = 2 * Math.log(theta) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1], tmp[2]});
                    } else {
                        tmp[0] = Math.pow(1 - theta, 2) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 3) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = theta * (1 - theta) * (nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1)) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = Math.pow(theta, 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1] + tmp[2];
                    }
                }
            case 2:
                // 1/1
                if (singleADO) {
                    // single ado
                    if (useLogPartials) {
                        tmp[0] = Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = Math.log(theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);

                        return logSumExp(new double[]{tmp[0], tmp[1]});
                    } else {
                        tmp[0] = (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = theta * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);

                        return tmp[0] + tmp[1];
                    }
                } else {
                    // locus ado
                    if (useLogPartials) {
                        tmp[0] = 2 * Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = Math.log(2) + Math.log(theta) + Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = 2 * Math.log(theta) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1], tmp[2]});
                    } else {
                        tmp[0] = Math.pow(1 - theta, 2) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = 2 * theta * (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = Math.pow(theta, 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1] + tmp[2];
                    }
                }
            case 3:
                // 1/1'
                if (singleADO) {
                    // single ado
                    if (useLogPartials) {
                        tmp[0] = Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 2) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = Math.log(theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);

                        return logSumExp(new double[]{tmp[0], tmp[1]});
                    } else {
                        tmp[0] = (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = theta * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);

                        return tmp[0] + tmp[1];
                    }
                } else {
                    // locus ado
                    if (useLogPartials) {
                        tmp[0] = 2 * Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 2) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = Math.log(2) + Math.log(theta) + Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = 2 * Math.log(theta) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1], tmp[2]});
                    } else {
                        tmp[0] = Math.pow(1 - theta, 2) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 2);
                        tmp[1] = 2 * theta * (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[2] = Math.pow(theta, 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1] + tmp[2];
                    }
                }
            case 4:
                // 0/-
                if (singleADO) {
                    // single ado
                    if (useLogPartials) {
                        tmp[0] = Math.log(1 - theta / 2) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = Math.log(theta) - Math.log(2) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1]});
                    } else {
                        tmp[0] = (1 - theta / 2) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = (theta / 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1];
                    }
                } else {
                    // locus ado
                    if (useLogPartials) {
                        tmp[0] = Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = Math.log(theta) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1]});
                    } else {
                        tmp[0] = (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = theta * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1];
                    }
                }
            case 5:
                // 1/-
                if (singleADO) {
                    // single ado
                    if (useLogPartials) {
                        tmp[0] = Math.log(1 - theta / 2) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = Math.log(theta) - Math.log(2) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1]});
                    } else {
                        tmp[0] = (1 - theta / 2) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = (theta / 2) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1];
                    }
                } else {
                    // locus ado
                    if (useLogPartials) {
                        tmp[0] = Math.log(1 - theta) + nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = Math.log(theta) + seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return logSumExp(new double[]{tmp[0], tmp[1]});
                    } else {
                        tmp[0] = (1 - theta) * nucReadCountsModel.getNucReadCountsLikelihood(taxonIndex, indexToNucReadCountsLikelihood + 1) * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood + 1);
                        tmp[1] = theta * seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                        return tmp[0] + tmp[1];
                    }
                }
            case 6:
                // -
                tmp[0] = seqCovModel.getSeqCovLikelihood(taxonIndex, indexToSeqCovLikelihood);

                return tmp[0];
            default:
                throw new IllegalStateException("Unexpected genotype index: " + genotypeIndex + " (" + this.getClass().getName() + ")");
        }
    } // computeMixedLikelihoodCore

}
