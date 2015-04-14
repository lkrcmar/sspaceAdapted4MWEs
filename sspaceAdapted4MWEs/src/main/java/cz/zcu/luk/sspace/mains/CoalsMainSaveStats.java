/*
 * Copyright 2009 Keith Stevens 
 *
 * This file is part of the S-Space package and is covered under the terms and
 * conditions therein.
 *
 * The S-Space package is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as published
 * by the Free Software Foundation and distributed hereunder to you.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS" AND NO REPRESENTATIONS OR WARRANTIES,
 * EXPRESS OR IMPLIED ARE MADE.  BY WAY OF EXAMPLE, BUT NOT LIMITATION, WE MAKE
 * NO REPRESENTATIONS OR WARRANTIES OF MERCHANT- ABILITY OR FITNESS FOR ANY
 * PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION
 * WILL NOT INFRINGE ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER
 * RIGHTS.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

package cz.zcu.luk.sspace.mains;

import cz.zcu.luk.sspace.coals.CoalsSaveStats;
import cz.zcu.luk.sspace.common.SemanticSpaceSaveStats;
import cz.zcu.luk.sspace.matrix.*;
import edu.ucla.sspace.common.ArgOptions;
import edu.ucla.sspace.common.SemanticSpaceIO.SSpaceFormat;
import edu.ucla.sspace.matrix.MatrixFactorization;
import edu.ucla.sspace.matrix.SVD;


/**
 * An executable class for running {@link edu.ucla.sspace.coals.Coals} from the
 * command line.  This class takes in several command line arguments.
 *
 * <ul>
 * <li> {@code --dimensions=<int>} how many dimensions to use for the regular
 *      word vectors. See {@link edu.ucla.sspace.coals.Coals} for a default value
 *
 * <li> {@code --reduce} If present, the word-word matrix will be reduced using
 *      Singular Valued Decomposition otherwise no reduction will be performed.
 *
 * <li> {code --reduceDimensions=<int>} size of the reduced svd vectors. See
 *      {@link edu.ucla.sspace.coals.Coals} for a default value.
 * </ul>
 *
 * <p>
 *
 * An invocation will produce one file as output where the file name will mark
 * the number of words kept in the word-word matrix, and wether or not SVD was
 * used.
 *
 * @see edu.ucla.sspace.coals.Coals
 */
public class CoalsMainSaveStats extends GenericMainSaveStats {

    /**
     * Uninstantiable.
     */
    private CoalsMainSaveStats() {
    }

    /**
     * {@inheritDoc}
     */
    public void addExtraOptions(ArgOptions options) {
          options.addOption('n', "dimensions", 
                            "Set the number of columns to keep in the raw " +
                            "co-occurance matrix.",
                            true, "INT", "Optional"); 
          options.addOption('m', "maxWords",
                            "Set the maximum number of words to keep in the " +
                            "space, ordered by frequency",
                            true, "INT", "Optional");
          options.addOption('s', "reducedDimension", 
                            "Set the number of dimension to reduce to " +
                            "using the Singular Value Decompositon.  This is " +
                            "used if --reduce is set.",
                            true, "INT", "Optional");
          options.addOption('r', "reduce", 
                            "Set to true if the co-occurrance matrix should " +
                            "be reduced using the Singluar Value Decomposition",
                            false, null, "Optional");
    }

    public static void main(String[] args) throws Exception {
        CoalsMainSaveStats coals = new CoalsMainSaveStats();
        coals.run(args);
    }
    
    /**
     * {@inheritDoc}
     */
    public SemanticSpaceSaveStats getSpace() {
        TransformExtended transform = new CorrelationTransformExtendedSerialiazable();
        MatrixFactorization reducer = (argOptions.hasOption("reduce"))
            ? SVD.getFastestAvailableFactorization()
            : null;

        return new CoalsSaveStats(transform, reducer,
                         argOptions.getIntOption("reducedDimension", 0),
                         argOptions.getIntOption("maxWords", 0),
                         argOptions.getIntOption("dimensions", 0));
    }

    /**
     * {@inheritDoc}
     */
    protected SSpaceFormat getSpaceFormat() {
        return SSpaceFormat.SPARSE_BINARY;
    }
}
