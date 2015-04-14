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

import edu.ucla.sspace.basis.BasisMapping;
import edu.ucla.sspace.basis.StringBasisMapping;
import edu.ucla.sspace.common.ArgOptions;
import edu.ucla.sspace.common.SemanticSpace;
import edu.ucla.sspace.common.SemanticSpaceIO;
import edu.ucla.sspace.common.SemanticSpaceIO.SSpaceFormat;
import edu.ucla.sspace.matrix.LogEntropyTransform;
import edu.ucla.sspace.matrix.SVD;
import edu.ucla.sspace.matrix.SVD.Algorithm;
import edu.ucla.sspace.matrix.Transform;
import edu.ucla.sspace.matrix.factorization.SingularValueDecomposition;
import edu.ucla.sspace.util.ReflectionUtil;
import edu.ucla.sspace.util.SerializableUtil;

import java.io.File;
import java.io.IOError;
import java.io.IOException;
import java.util.Map;

import cz.zcu.luk.sspace.common.DocumentSemanticSpace;
import cz.zcu.luk.sspace.config.Config;
import cz.zcu.luk.sspace.lsa.LatentSemanticAnalysisModified;


/**
 * An executable class for running {@link edu.ucla.sspace.lsa.LatentSemanticAnalysis} (LSA) from the
 * command line.  This class takes in several command line arguments.
 *
 * <ul>
 *
 * <li><u>Required (at least one of)</u>:
 *   <ul>
 *
 *   <li> {@code -d}, {@code --docFile=FILE[,FILE...]} a file where each line is
 *        a document.  This is the preferred input format for large corpora
 *
 *   <li> {@code -f}, {@code --fileList=FILE[,FILE...]} a list of document files
 *        where each file is specified on its own line.
 *
 *   </ul>
 *
 * <li><u>Algorithm Options</u>:
 *   <ul>
 *
 *   <li> {@code --dimensions=<int>} how many dimensions to use for the LSA
 *        vectors.  See {@link edu.ucla.sspace.lsa.LatentSemanticAnalysis} for default value
 *
 *   <li> {@code --preprocess=<class name>} specifies an instance of {@link
 *        edu.ucla.sspace.lsa.MatrixTransformer} to use in preprocessing the
 *        word-document matrix compiled by LSA prior to computing the SVD.  See
 *        {@link edu.ucla.sspace.lsa.LatentSemanticAnalysis} for default value
 *
 *   <li> {@code -F}, {@code --tokenFilter=FILE[include|exclude][,FILE...]}
 *        specifies a list of one or more files to use for {@link
 *        edu.ucla.sspace.text.TokenFilter filtering} the documents.  An option
 *        flag may be added to each file to specify how the words in the filter
 *        filter should be used: {@code include} if only the words in the filter
 *        file should be retained in the document; {@code exclude} if only the
 *        words <i>not</i> in the filter file should be retained in the
 *        document.
 *
 *   <li> {@code -S}, {@code --svdAlgorithm}={@link
 *        edu.ucla.sspace.matrix.SVD.Algorithm} species a specific {@code
 *        SVD.Algorithm} method to use when reducing the dimensionality in LSA.
 *        In general, users should not need to specify this option, as the
 *        default setting will choose the fastest algorithm available on the
 *        system.  This is only provided as an advanced option for users who
 *        want to compare the algorithms' performance or any variations between
 *        the SVD results.
 *
 *   </ul>
 *
 * <li><u>Program Options</u>:
 *   <ul>
 *
 *   <li> {@code -o}, {@code --outputFormat=}<tt>text|binary}</tt> Specifies the
 *        output formatting to use when generating the semantic space ({@code
 *        .sspace}) file.  See {@link edu.ucla.sspace.common.SemanticSpaceUtils
 *        SemanticSpaceUtils} for format details.
 *
 *   <li> {@code -t}, {@code --threads=INT} how many threads to use when
 *        processing the documents.  The default is one per core.
 *
 *   <li> {@code -w}, {@code --overwrite=BOOL} specifies whether to overwrite
 *        the existing output files.  The default is {@code true}.  If set to
 *        {@code false}, a unique integer is inserted into the file name.
 *
 *   <li> {@code -v}, {@code --verbose}  specifies whether to print runtime
 *        information to standard out
 *
 *   </ul>
 *
 * </ul>
 *
 * <p>
 *
 * An invocation will produce one file as output {@code
 * lsa-semantic-space.sspace}.  If {@code overwrite} was set to {@code true},
 * this file will be replaced for each new semantic space.  Otherwise, a new
 * output file of the format {@code lsa-semantic-space<number>.sspace} will be
 * created, where {@code <number>} is a unique identifier for that program's
 * invocation.  The output file will be placed in the directory specified on the
 * command line.
 *
 * <p>
 *
 * This class is desgined to run multi-threaded and performs well with one
 * thread per core, which is the default setting.
 *
 * @see edu.ucla.sspace.lsa.LatentSemanticAnalysis
 * @see edu.ucla.sspace.matrix.Transform Transform
 *
 * @author David Jurgens
 */
public class LSAMainDocumentsSim extends GenericMainModified {

    private BasisMapping<String, String> basis;

    private LatentSemanticAnalysisModified lsaInst;

    private LSAMainDocumentsSim() {
    }

    /**
     * Adds all of the options to the {@link edu.ucla.sspace.common.ArgOptions}.
     */
    protected void addExtraOptions(ArgOptions options) {
        options.addOption('n', "dimensions", 
                          "the number of dimensions in the semantic space",
                          true, "INT", "Algorithm Options"); 
        options.addOption('p', "preprocess", "a MatrixTransform class to "
                          + "use for preprocessing", true, "CLASSNAME",
                          "Algorithm Options");
        options.addOption('S', "svdAlgorithm", "a specific SVD algorithm to use"
                          , true, "SVD.Algorithm", 
                          "Advanced Algorithm Options");
        options.addOption('B', "saveTermBasis",
                          "If true, the term basis mapping will be stored " +
                          "to the given file name",
                          true, "FILE", "Optional");
    }

    public static void main(String[] args) throws Exception {
        LSAMainDocumentsSim lsa = new LSAMainDocumentsSim();
        lsa.run(args);
    }
    
    protected SemanticSpace getSpace() {
        try {
            int dimensions = argOptions.getIntOption("dimensions", 300);
            Transform transform = new LogEntropyTransform();
            if (argOptions.hasOption("preprocess"))
                transform = ReflectionUtil.getObjectInstance(
                        argOptions.getStringOption("preprocess"));
            String algName = argOptions.getStringOption("svdAlgorithm", "ANY");
            SingularValueDecomposition factorization = SVD.getFactorization(
                    Algorithm.valueOf(algName.toUpperCase()));
            basis = new StringBasisMapping();

            // LK change - first param set to true - preserve doc. space
            // and store in object..
            lsaInst = new LatentSemanticAnalysisModified(
                    true, dimensions, transform, factorization, true, basis);
            return lsaInst;
        } catch (IOException ioe) {
            throw new IOError(ioe);
        }
    }

    /**
     * Returns the {@likn SSpaceFormat.BINARY binary} format as the default
     * format of a {@code LatentSemanticAnalysis} space.
     */
    protected SSpaceFormat getSpaceFormat() {
        return SSpaceFormat.BINARY;
    }

    protected void postProcessing() {
        System.out.println("Post processing runs!!!");
        //System.out.println("Doc space size is: " + lsaInst.documentSpaceSize());
        if (argOptions.hasOption('B'))
            SerializableUtil.save(basis, argOptions.getStringOption('B'));

        // LK change - perform document similarity task..
      //  System.out.println("size of headerToIndex: " + lsaInst.getHeaderToIndex().size());
        //System.out.println("-----------------------------------");
        // load course pairs from file
//        List<CoursePairScored> coursePairsScored = CoursePairScoredIO
//                .loadCoursePairsScored("/storage/home/lkrcmar/baseline.txt");
//
//        for (CoursePairScored onePair : coursePairsScored) {
//            Integer indexFirstDoc = lsaInst.getHeaderToIndex().get(Integer.parseInt(onePair.getIdFirstCourse()));
//            Integer indexSecondDoc = lsaInst.getHeaderToIndex().get(Integer.parseInt(onePair.getIdSecondCourse()));
//            double sim = Similarity.getSimilarity(Similarity.SimType.COSINE,
//                    lsaInst.getDocumentVector(indexFirstDoc), lsaInst.getDocumentVector(indexSecondDoc));
//            onePair.setSimilarity((new Double(sim)).toString());
//            //System.out.println(onePair.toString() + "\t" + lsaInst.getDocumentVector(indexFirstDoc).length());
//        }
//
//        // store course pairs to file
//        CoursePairScoredIO.saveCoursePairsScored("/storage/home/lkrcmar/lsaCoursesSim.txt", coursePairsScored);

//        System.out.println("-----------------------------------");
//        System.out.println(Similarity.getSimilarity(Similarity.SimType.COSINE,
//                lsaInst.getDocumentVector((0)),
//                lsaInst.getDocumentVector((299400))));
//        System.out.println(Similarity.getSimilarity(Similarity.SimType.COSINE,
//                lsaInst.getDocumentVector((1000)),
//                lsaInst.getDocumentVector((299438))));
//        System.out.println(Similarity.getSimilarity(Similarity.SimType.COSINE,
//                lsaInst.getDocumentVector((1000)),
//                lsaInst.getDocumentVector((299439))));
//        System.out.println(Similarity.getSimilarity(Similarity.SimType.COSINE,
//                lsaInst.getDocumentVector((1000)),
//                lsaInst.getDocumentVector((299440))));

        DocumentSemanticSpace docSemSpace = new DocumentSemanticSpace(lsaInst);
        // cannot get number of processed documents in the following way
        // the reason is that empty documents are not processed,
        // however, counted -> they increase headerToIndex size..
        //System.out.println("Number of processed documents: " + lsaInst.getHeaderToIndex().size());
        System.out.println("Number of processed documents: " + lsaInst.getDocumentCounter());
        System.out.println("Document space size: " + docSemSpace.getWords().size());
        try {

            Map<String, String> configuration = Config.getInstance().configuration;
            String lsaDocSpaces = configuration.get("dataDir") + "/" + configuration.get("lsaDocSpacesDN");

            String outputFN = lsaDocSpaces + "/" + docSemSpace.getSpaceName() + EXT;
            System.out.println("Docs space output file name: " + outputFN);
            SemanticSpaceIO.save(docSemSpace, new File(outputFN), SSpaceFormat.BINARY);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    /**
     * {@inheritDoc}
     */
    protected String getAlgorithmSpecifics() {
        return 
            "The --svdAlgorithm provides a way to manually specify which " + 
            "algorithm should\nbe used internally.  This option should not be" +
            " used normally, as LSA will\nselect the fastest algorithm " +
            "available.  However, in the event that it\nis needed, valid" +
            " options are: SVDLIBC, SVDLIBJ, MATLAB, OCTAVE, JAMA and COLT\n";
    }
}
