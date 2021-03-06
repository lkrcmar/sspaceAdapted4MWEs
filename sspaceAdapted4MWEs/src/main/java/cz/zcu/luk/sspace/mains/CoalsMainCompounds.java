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

import edu.ucla.sspace.common.ArgOptions;
import edu.ucla.sspace.common.SemanticSpace;
import edu.ucla.sspace.common.SemanticSpaceIO.SSpaceFormat;
import edu.ucla.sspace.matrix.MatrixFactorization;
import edu.ucla.sspace.matrix.SVD;
import edu.ucla.sspace.util.FileResourceFinder;
import edu.ucla.sspace.util.ResourceFinder;

import java.io.BufferedReader;
import java.io.IOError;
import java.io.IOException;
import java.util.LinkedHashSet;
import java.util.Set;

import cz.zcu.luk.sspace.coals.CoalsCompounds;
import cz.zcu.luk.sspace.matrix.*;


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
public class CoalsMainCompounds extends GenericMainModified {

    /**
     * Uninstantiable.
     */
    private CoalsMainCompounds() {
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
          options.addOption('c', "compoundsInvestigated", "a file where each line is a " +
                "recognized compound for which a statistic is being done." +
                " No words' vectors are influenced", true, "FILE",
                "Program Options");
    }

    public static void main(String[] args) throws Exception {
        CoalsMainCompounds coals = new CoalsMainCompounds();
        coals.run(args);
    }
    
    /**
     * {@inheritDoc}
     */
    public SemanticSpace getSpace() {
        TransformExtended transform = new CorrelationTransformExtended();
        MatrixFactorization reducer = (argOptions.hasOption("reduce"))
            ? SVD.getFastestAvailableFactorization()
            : null;

        // LK changed - added processing of compounds..
        Set<String> compounds = null;
        if (argOptions.getStringOption('c', null) != null) {
            ResourceFinder resourceFinder = new FileResourceFinder();
            String compoundTokensProp =
                    argOptions.getStringOption("compoundsInvestigated");
            if (compoundTokensProp != null) {
                compounds = new LinkedHashSet<String>();
                // Load the tokens from file
                try {
                    BufferedReader br = resourceFinder.open(compoundTokensProp);
                    for (String line = null; (line = br.readLine()) != null; ) {
                        compounds.add(line);
                    }
                } catch (IOException ioe) {
                    // rethrow
                    throw new IOError(ioe);
                }
            }
        }

        return new CoalsCompounds(transform, reducer,
                         argOptions.getIntOption("reducedDimension", 0),
                         argOptions.getIntOption("maxWords", 0),
                         argOptions.getIntOption("dimensions", 0),
                         compounds);
    }

    /**
     * {@inheritDoc}
     */
    protected SSpaceFormat getSpaceFormat() {
        return SSpaceFormat.SPARSE_BINARY;
    }
}
