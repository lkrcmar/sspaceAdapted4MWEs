/*
 * Copyright 2010 David Jurgens
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

import edu.ucla.sspace.basis.StringBasisMapping;
import edu.ucla.sspace.common.ArgOptions;
import edu.ucla.sspace.common.SemanticSpace;
import edu.ucla.sspace.common.SemanticSpaceIO.SSpaceFormat;
import edu.ucla.sspace.matrix.SvdlibcSparseBinaryMatrixBuilder;
import edu.ucla.sspace.util.FileResourceFinder;
import edu.ucla.sspace.util.ReflectionUtil;
import edu.ucla.sspace.util.ResourceFinder;

import java.io.BufferedReader;
import java.io.IOError;
import java.io.IOException;
import java.util.LinkedHashSet;
import java.util.Set;

import cz.zcu.luk.sspace.matrix.LogEntropyTransformExtended;
import cz.zcu.luk.sspace.matrix.TransformExtended;
import cz.zcu.luk.sspace.vsm.VectorSpaceModelCompoundsStopwords;

/**
 * An executable class for running {@link edu.ucla.sspace.vsm.VectorSpaceModel} (VSM) from the
 * command line.  See the <a
 * href="http://code.google.com/p/airhead-research/wiki/VectorSpaceModel"> wiki
 * page</a> for details on running this class from the command line. <p>
 *
 * This class is desgined to run multi-threaded and performs well with one
 * thread per core, which is the default setting.
 *
 * @see VectorSpaceMdeol
 * @see edu.ucla.sspace.matrix.Transform Transform
 *
 * @author David Jurgens
 */
public class VsmMainCompoundsStopwords extends GenericMainModifiedStopwords {

    private VsmMainCompoundsStopwords() { }

    /**
     * Adds all of the options to the {@link edu.ucla.sspace.common.ArgOptions}.
     */
    protected void addExtraOptions(ArgOptions options) {
        options.addOption('p', "preprocess", "a MatrixTransform class to "
                          + "use for preprocessing", true, "CLASSNAME",
                          "Algorithm Options");
        options.addOption('c', "compoundsInvestigated", "a file where each line is a " +
                "recognized compound for which a statistic is being done." +
                " No words' vectors are influenced", true, "FILE",
                "Program Options");
    }

    public static void main(String[] args) {
        VsmMainCompoundsStopwords vsm = new VsmMainCompoundsStopwords();
        try {
            vsm.run(args);
        }
        catch (Throwable t) {
            t.printStackTrace();
        }
    }
    
    /**
     * {@inheritDoc}
     */
    protected SemanticSpace getSpace() {
        // LK changed
        try {
            TransformExtended transform = new LogEntropyTransformExtended();
            if (argOptions.hasOption("preprocess"))
                transform = ReflectionUtil.getObjectInstance(
                        argOptions.getStringOption("preprocess"));

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

            return new VectorSpaceModelCompoundsStopwords(false, new StringBasisMapping(), transform, new SvdlibcSparseBinaryMatrixBuilder(),
                    new SvdlibcSparseBinaryMatrixBuilder(), compounds);
        } catch (IOException ioe) {
            throw new IOError(ioe);
        }
    }

    /**
     * Returns the {@likn SSpaceFormat.SPARSE_BINARY binary} format as the
     * default format of a {@code VectorSpaceModel} space.
     */
    protected SSpaceFormat getSpaceFormat() {
        return SSpaceFormat.SPARSE_BINARY;
    }

    /**
     * {@inheritDoc}
     */
//    protected Properties setupProperties() {
//        // use the System properties in case the user specified them as
//        // -Dprop=<val> to the JVM directly.
//        Properties props = System.getProperties();
//
//        if (argOptions.hasOption("transform")) {
//            props.setProperty(VectorSpaceModel.MATRIX_TRANSFORM_PROPERTY,
//                              argOptions.getStringOption("transform"));
//        }
//
//        return props;
//    }
}
