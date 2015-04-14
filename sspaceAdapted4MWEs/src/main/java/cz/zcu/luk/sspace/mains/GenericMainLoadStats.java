/*
 * Copyright 2009 David Jurgens
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

import edu.ucla.sspace.mains.GenericMain;
import edu.ucla.sspace.text.Document;
import edu.ucla.sspace.text.IteratorFactory;
import edu.ucla.sspace.util.LoggerUtil;
import edu.ucla.sspace.util.WorkQueue;
import edu.ucla.sspace.common.SemanticSpaceIO;

import java.io.File;
import java.util.Iterator;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import cz.zcu.luk.sspace.common.SemanticSpaceLoadStats;
import cz.zcu.luk.sspace.config.Config;

public abstract class GenericMainLoadStats extends GenericMain {

    private static final Logger LOGGER =
        Logger.getLogger(GenericMainLoadStats.class.getName());

    /**
     * Processes the arguments and begins processing the documents using the
     * {@link edu.ucla.sspace.common.SemanticSpace} returned by {@link #getSpace() getSpace}.
     *
     * @param args arguments used to configure this program and the {@code
     *        SemanticSpace}
     */
    public void run(String[] args) throws Exception {
        // LK added
        LOGGER.info("RUN started");
        Set<Thread> threadSet = Thread.getAllStackTraces().keySet();
        System.out.println(threadSet.toString());

        if (args.length == 0) {
            usage();
            System.exit(1);
        }
        argOptions.parseOptions(args);

        if (argOptions.numPositionalArgs() == 0) {
            throw new IllegalArgumentException("must specify additional information");
        }

        verbose = argOptions.hasOption('v') || argOptions.hasOption("verbose");
        // If verbose output is enabled, update all the loggers in the S-Space
        // package logging tree to output at Level.FINE (normally, it is
        // Level.INFO).  This provides a more detailed view of how the execution
        // flow is proceeding.
        if (verbose)
            LoggerUtil.setLevel(Level.FINE);

        // Check whether this class supports mutlithreading when deciding how
        // many threads to use by default
        int numThreads = (isMultiThreaded)
                ? (System.getenv().get("PBS_NUM_PPN") == null ? Runtime.getRuntime().availableProcessors() : (Integer.parseInt(System.getenv().get("PBS_NUM_PPN")) - 1))
                : 1;
        if (argOptions.hasOption("threads")) {
            numThreads = argOptions.getIntOption("threads");
        }
        // Initialize the work queue so that any system that uses one is
        // limitedto the number of processes specified by the command line.
        WorkQueue.getWorkQueue(numThreads);

        boolean overwrite = true;
        if (argOptions.hasOption("overwrite")) {
            overwrite = argOptions.getBooleanOption("overwrite");
        }

        handleExtraOptions();

        Properties props = setupProperties();

        // Initialize the IteratorFactory to tokenize the documents according to
        // the specified configuration (e.g. filtering, compound words)
        if (argOptions.hasOption("tokenFilter")) {
            props.setProperty(IteratorFactory.TOKEN_FILTER_PROPERTY,
                    argOptions.getStringOption("tokenFilter"));
        }

        // Set any tokenizing options.
        if (argOptions.hasOption("stemmingAlgorithm"))
            props.setProperty(IteratorFactory.STEMMER_PROPERTY,
                    argOptions.getStringOption("stemmingAlgorithm"));

        if (argOptions.hasOption("compoundWords")) {
            props.setProperty(IteratorFactory.COMPOUND_TOKENS_FILE_PROPERTY,
                    argOptions.getStringOption("compoundWords"));
        }
        if (argOptions.hasOption("wordLimit"))
            props.setProperty(IteratorFactory.TOKEN_COUNT_LIMIT_PROPERTY,
                    argOptions.getStringOption("wordLimit"));

        // Get the format to be used when writing the semantic space.
        SemanticSpaceIO.SSpaceFormat format = (argOptions.hasOption("outputFormat"))
                ? SemanticSpaceIO.SSpaceFormat.valueOf(
                argOptions.getStringOption("outputFormat").toUpperCase())
                : getSpaceFormat();

        IteratorFactory.setProperties(props);

        // use the System properties in case the user specified them as
        // -Dprop=<val> to the JVM directly.

        SemanticSpaceLoadStats space = getSpace();

        File outputFile = null;
        // LK change
        String lastPar = (argOptions.getPositionalArg(1) == null) ? "" : argOptions.getPositionalArg(1);
        Map<String, String> configuration = Config.getInstance().configuration;
        String dirPlusSpaceNameNoExtension = configuration.get("dataDir") + File.separator +
                configuration.get("wordSpaceDir") + File.separator + space.getSpaceName() + "_" + argOptions.getPositionalArg(0) + lastPar;
        outputFile = new File(dirPlusSpaceNameNoExtension + EXT);
        System.out.println("output File (generated from source code -" +
                " does not follow the given output file name): " + outputFile);

        String spaceNameLoaded = new String("W" + space.getSpaceName().substring(1));
        String dirPlusSpaceNameNoExtensionLoaded = configuration.get("dataDir") + File.separator +
                configuration.get("wordSpaceDir") + File.separator + spaceNameLoaded + "_" + argOptions.getPositionalArg(0);
        loadSStats(space, dirPlusSpaceNameNoExtensionLoaded);

        // all the documents are listed in one file, with one document per line
        Iterator<Document> docIter = getDocumentIterator();

        processDocumentsAndSpace(space, docIter, numThreads, props);

        long startTime = System.currentTimeMillis();
        saveSSpace(space, outputFile, format);
        long endTime = System.currentTimeMillis();
        verbose("printed space in %.3f seconds",
                ((endTime - startTime) / 1000d));

        Set<Thread> threadSet3 = Thread.getAllStackTraces().keySet();
        System.out.println(threadSet3.toString());

        postProcessing();
    }

    protected void loadSStats(SemanticSpaceLoadStats space, String dirPlusSpaceNameNoExtensionLoaded) {
        space.loadStatistics(dirPlusSpaceNameNoExtensionLoaded);
    }

    protected abstract SemanticSpaceLoadStats getSpace();
}
