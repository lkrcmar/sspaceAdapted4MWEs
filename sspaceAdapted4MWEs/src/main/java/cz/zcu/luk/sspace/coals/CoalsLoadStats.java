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

package cz.zcu.luk.sspace.coals;


import edu.ucla.sspace.matrix.Matrices;
import edu.ucla.sspace.matrix.Matrix;
import edu.ucla.sspace.matrix.MatrixFactorization;
import edu.ucla.sspace.matrix.SparseMatrix;
import edu.ucla.sspace.text.IteratorFactory;
import edu.ucla.sspace.vector.*;
import edu.ucla.sspace.vector.Vector;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Logger;

import cz.zcu.luk.sspace.common.SemanticSpaceLoadStats;
import cz.zcu.luk.sspace.common.Serializer;
import cz.zcu.luk.sspace.matrix.*;


/**
 * An implementation of the COALS Semantic Space model.  This implementation is
 * based on:
 *
 * <p style="font-family:Garamond, Georgia, serif"> Rohde, D. L. T.,
 * Gonnerman, L. M., Plaut, D. C. (2005).  An Improved Model of Semantic
 * Similarity Based on Lexical Co-Occurrence. <i>Cognitive Science</i>
 * <b>(submitted)</b>.  Available <a
 * href="http://www.cnbc.cmu.edu/~plaut/papers/pdf/RohdeGonnermanPlautSUB-CogSci.COALS.pdf">here</a></p>
 *
 * COALS first computes a term by term co-occurrance using a ramped 4-word
 * window. Once all documents have been processed, the co-occurrence matrix will
 * be re ordered such that only the {@code N} most frequent terms have their
 * semantic vectors retained and only the {@code M} most frequent terms are used
 * as co-occurrence features.  These values can be set by the {@value
 *#MAX_WORDS_PROPERTY} and {@value MAX_DIMENSIONS_PROPERTY} properties,
 * resepctively.  After re ordering the semantic vectors and features, {@link
 * edu.ucla.sspace.matrix.CorrelationTransform} is used to rerank all co-occurrence scores.  As part of
 * this transform, all negative correlations are dropped and replaced with a 0.
 * Finally, and optionally, the {@link edu.ucla.sspace.matrix.SVD} is used to reduce the semantic
 * space.  To set the number of retained dimensions via {@link edu.ucla.sspace.matrix.SVD}, set the
 * {@value REDUCE_DIMENSION_PROPERTY} property.
 *
 * @author Keith Stevens
 */
public class CoalsLoadStats implements SemanticSpaceLoadStats {

    /**
     * The property prefix for other settings.
     */
    public static final String PROPERTY_PREFIX =
        "edu.ucla.sspace.coals.Coals";

    /**
     * Specifies whether or not the co-occurance matrix should be reduced.
     */
    public static final String REDUCE_MATRIX_PROPERTY =
        PROPERTY_PREFIX + ".reduce";

    /**
     * Specifies the number of dimensions the co-occurance matrix should be
     * reduced to.
     */
    public static final String REDUCE_DIMENSION_PROPERTY =
        PROPERTY_PREFIX + ".dimension";

    /**
     * Specifies the number of dimensions in the raw co-occurrance matrix to
     * maintain.
     */
    public static final String MAX_DIMENSIONS_PROPERTY =
        PROPERTY_PREFIX + ".maxDimensions";

    /**
     * Specifies the number of words to build semantics for.
     */
    public static final String MAX_WORDS_PROPERTY =
        PROPERTY_PREFIX + ".maxWords";

    /**
     * Specifies if Coals should not normalize the co-occurance matrix.
     */
    public static final String DO_NOT_NORMALIZE_PROPERTY =
        PROPERTY_PREFIX + ".doNotNormalize";

    /**
     * The default number of dimensions to reduce to.
     */
    private static final int DEFAULT_REDUCE_DIMENSIONS = 800;

    /**
     * The default number of dimensions to save in the co-occurrance matrix.
     */
    private static final int DEFAULT_MAX_DIMENSIONS = 14000;

    /**
     * The default number of rows to save in the co-occurrance matrix.
     */
    private static final int DEFAULT_MAX_WORDS = 15000;

    /**
     * The name of this {@code SemanticSpace}
     */
    public static final String COALS_SSPACE_NAME =
        "O_COALS";

    /**
     * The logger used to record all output
     */
    private static final Logger COALS_LOGGER =
        Logger.getLogger(CoalsLoadStats.class.getName());

    /**
     * A mapping from each word to the vector the represents its semantics
     */
    //private Map<String, SparseDoubleVector> wordToSemantics;
    private Map<String, SparseDoubleVector> otherToSemantics;

    /**
     * A mapping from word to index number.
     */
    private Map<String, Integer> termToIndex;
    private Map<String, Integer> otherToIndex;

    /**
     * A map containg the total frequency counts of each word.
     */
    //private ConcurrentMap<String, AtomicInteger> totalWordFreq;

    /**
     * The final reduced matrix.
     */
    //private Matrix finalCorrelation;
    private Matrix otherFinalCorrelation;

    /**
     * Specifies the number of reduced dimensions if the matrix is reduced by
     * SVD.
     */
    private final int reducedDimensions;

    /**
     * The maximum number of words that will be retained by {@link CoalsLoadStats}.
     */
    private final int maxWords;

    /**
     * The maximum number of co-occurring words that will be retained by {@link
     * CoalsLoadStats}.
     */
    private final int maxDimensions;

    /**
     * A counter for keeping track of the index values of words.
     */
    //private int wordIndexCounter;
    private int otherIndexCounter;

    /**
     * The {@link edu.ucla.sspace.matrix.MatrixFactorization} algorithm that will decompose the word by
     * document feature space into two smaller feature spaces: a word by class
     * feature space and a class by feature space.
     */
    private final MatrixFactorization reducer;

    /**
     * The {@link edu.ucla.sspace.matrix.Transform} applied to the co-ocucrrence counts, if not {@code
     * null}.
     */
    private TransformExtended transform; // final keyword removed to allow load statistics..

    public CoalsLoadStats(TransformExtended transform, MatrixFactorization reducer) {
        this(transform, reducer, DEFAULT_REDUCE_DIMENSIONS,
             DEFAULT_MAX_WORDS, DEFAULT_MAX_DIMENSIONS);
    }

    /**
     * Creats a {@link CoalsLoadStats} instance.
     */
    public CoalsLoadStats(TransformExtended transform,
                          MatrixFactorization reducer,
                          int reducedDimensions,
                          int maxWords,
                          int maxDimensions) {
        //termToIndex = new HashMap<String, Integer>();
        //totalWordFreq = new ConcurrentHashMap<String, AtomicInteger>();
        //wordToSemantics = new HashMap<String, SparseDoubleVector>(1024, 4f);
        otherToIndex = new HashMap<String, Integer>();
        otherToSemantics = new HashMap<String, SparseDoubleVector>(1024, 4f);
        otherFinalCorrelation = null;
        //this.transform = transform;
        this.reducer = reducer;
        this.reducedDimensions = (reducedDimensions == 0)
            ? DEFAULT_REDUCE_DIMENSIONS
            : reducedDimensions;
        this.maxWords = (maxWords == 0)
            ? DEFAULT_MAX_WORDS
            : maxWords;
        this.maxDimensions = (maxDimensions == 0)
            ? DEFAULT_MAX_DIMENSIONS
            : maxDimensions;
    }

    /**
     * {@inheritDoc}
     */
    public Set<String> getWords() {
        //return termToIndex.keySet();
        return otherToIndex.keySet();
    }

    /**
     * {@inheritDoc}
     */
    public Vector getVector(String term) {
        Integer index = otherToIndex.get(term);
        if (index == null)
            return null;
        return Vectors.immutable(
                otherFinalCorrelation.getRowVector(index.intValue()));
    }

    public String getSpaceName() {
        String ret = COALS_SSPACE_NAME;
        ret += "_M" + maxWords;
        ret += "_N" + maxDimensions;
        if (reducer != null)
            //ret += "_svd_" + reducedDimensions;
            ret += "_D" + reducedDimensions;
        return ret;
    }

    public int getVectorLength() {
        return otherFinalCorrelation.columns();
    }

    /**
     * {@inheritDoc}
     */
    public void processDocument(BufferedReader document) throws IOException {
        Map<String, Integer> wordFreq = new HashMap<String, Integer>();
        Map<String, SparseDoubleVector> otherDocSemantics =
            new HashMap<String, SparseDoubleVector>();
        // System.out.println(termToIndex.toString());
        // Setup queues to track the set of previous and next words in a
        // context.
        Queue<String> prevWords = new ArrayDeque<String>();
        Queue<String> nextWords = new ArrayDeque<String>();

        Iterator<String> it = IteratorFactory.tokenizeOrdered(document);

        for (int i = 0; i < 4 && it.hasNext(); ++i)
            nextWords.offer(it.next());

        // Compute the co-occurrance statistics of each focus word in the
        // document.
        while (!nextWords.isEmpty()) {

            // Slide over the context by one word.
            if (it.hasNext())
                nextWords.offer(it.next());

            // Get the focus word
            String focusWord = nextWords.remove();
            if (!focusWord.equals(IteratorFactory.EMPTY_TOKEN)) {
                if (getIndexFor(focusWord) == null) { // means processing other word.. continue..
                    getIndexForOther(focusWord); // links words through indexes to otherFinalCorr.. matrix..

                    // Get the temprorary semantics for the focus word, create a new
                    // vector for them if needed.
                    SparseDoubleVector focusSemantics = otherDocSemantics.get(
                            focusWord);
                    if (focusSemantics == null) {
                        focusSemantics = new SparseHashDoubleVector(
                                termToIndex.size());
                                //Integer.MAX_VALUE);

                        otherDocSemantics.put(focusWord, focusSemantics);
                    }

                    // Process the previous words.
                    int offset = 4 - prevWords.size();
                    for (String word : prevWords) {
                        offset++;
                        if (word.equals(IteratorFactory.EMPTY_TOKEN))
                            continue;
                        Integer index = getIndexFor(word);
                        if (index != null) { // null means - no main word with no index processed..
                            focusSemantics.add(index, offset);
                        }
                    }

                    // Process the next words.
                    offset = 5;
                    for (String word : nextWords) {
                        offset--;
                        if (word.equals(IteratorFactory.EMPTY_TOKEN))
                            continue;
                        Integer index = getIndexFor(word); // null means - no main word with no index processed..
                        if (index != null) {
                            focusSemantics.add(index, offset);
                        }
                    }
                }
            }

            prevWords.offer(focusWord);
            if (prevWords.size() > 4)
                prevWords.remove();
        }

        // Add the temporary vectors for each word in this document to the
        // actual semantic fectors.
        for (Map.Entry<String, SparseDoubleVector> e :
                otherDocSemantics.entrySet()) {
            SparseDoubleVector focusSemantics = getSemanticVector(
                    e.getKey());
            // Get the non zero indices before hand so that they are cached
            // during the synchronized section.
            focusSemantics.getNonZeroIndices();
            synchronized (focusSemantics) {
                VectorMath.add(focusSemantics, e.getValue());
            }
        }

        // Store the total frequency counts of the words seen in this document
        // so far.
//        for (Map.Entry<String, Integer> entry : wordFreq.entrySet()) {
//            int count = entry.getValue().intValue();
//            AtomicInteger freq = totalWordFreq.putIfAbsent(
//                    entry.getKey(), new AtomicInteger(count));
//            if (freq != null)
//                freq.addAndGet(count);
//        }
    }

    /**
     * Returns the current semantic vector for the provided word, or if the word
     * is not currently in the semantic space, a vector is added for it and
     * returned.
     *
     * @param word a word
     *
     * @return the {@code SemanticVector} for the provide word.
     */
    private SparseDoubleVector getSemanticVector(String other) {
        SparseDoubleVector v = otherToSemantics.get(other);
        if (v == null) {
            // lock on the word in case multiple threads attempt to add it at
            // once
            synchronized(this) {
                // recheck in case another thread added it while we were waiting
                // for the lock
                v = otherToSemantics.get(other);
                if (v == null) {
                    //v = new CompactSparseVector();
                    v = new CompactSparseVector(termToIndex.size()); // ?????????????????????????????
                    otherToSemantics.put(other, v);
                }
            }
        }
        return v;
    }

    /**
     * Returns the index in the co-occurence matrix for this word.  If the word
     * was not previously assigned an index, this method adds one for it and
     * returns that index.
     */
    private Integer getIndexFor(String word) {
        Integer index = termToIndex.get(word);
        return index;
    }

    /**
     * Returns the index in the co-occurence matrix for this word.  If the word
     * was not previously assigned an index, this method adds one for it and
     * returns that index.
     */
    private int getIndexForOther(String other) {
        Integer index = otherToIndex.get(other);
        if (index == null) {
            synchronized(this) {
                // recheck to see if the term was added while blocking
                index = otherToIndex.get(other);
                // if another thread has not already added this word while the
                // current thread was blocking waiting on the lock, then add it.
                if (index == null) {
                    int i = otherIndexCounter++;
                    otherToIndex.put(other, i);
                    return i; // avoid the auto-boxing to assign i to index
                }
            }
        }
        return index;
    }

    /**
     * {@inheritDoc}
     */
    public void processSpace(Properties props) {
        COALS_LOGGER.info("Droppring dimensions from co-occurrance matrix.");
        // Read in the matrix from a file with dimensions dropped.
        buildMatrix(maxWords, maxDimensions);
//        System.out.println();
//        for (int i = 0; i < otherFinalCorrelation.rows(); i++) {
//            for (int j = 0; j < otherFinalCorrelation.columns(); j++) {
//                System.out.print(otherFinalCorrelation.get(i, j) + " ");
//            }
//            System.out.println();
//        }
//        System.out.println();
        COALS_LOGGER.info("Done dropping dimensions.");

        if (transform != null) {
            COALS_LOGGER.info("Normalizing co-occurrance matrix.");

            // Normalize the matrix using correlation.
            SparseDoubleVector[] compoundVectors = new SparseDoubleVector[otherFinalCorrelation.rows()];
            for (int i = 0; i < otherFinalCorrelation.rows(); i++) {
                DoubleVector dv = transform.transformRow(otherFinalCorrelation.getRowVector(i));
                if (dv instanceof SparseDoubleVector) {
                    compoundVectors[i] = (SparseDoubleVector)dv;
                }
                else {
                    throw new IllegalStateException("Expected sparse format of vector, but it is dense!");
                }
            }
            otherFinalCorrelation = Matrices.asSparseMatrix(Arrays.asList(compoundVectors));
            COALS_LOGGER.info("Done normalizing co-occurrance matrix.");
        }
//        System.out.println();
//        for (int i = 0; i < otherFinalCorrelation.rows(); i++) {
//            for (int j = 0; j < otherFinalCorrelation.columns(); j++) {
//                System.out.print(otherFinalCorrelation.get(i, j) + " ");
//            }
//            System.out.println();
//        }
//        System.out.println();
        if (reducer != null) {
            throw new NotImplementedException();
//            if (reducedDimensions > otherFinalCorrelation.columns())
//                throw new IllegalArgumentException(
//                        "Cannot reduce to more dimensions than exist");
//
//            COALS_LOGGER.info("Reducing using SVD.");
//            try {
//                File coalsMatrixFile =
//                    File.createTempFile("coals-term-doc-matrix", "dat");
//                coalsMatrixFile.deleteOnExit();
//                MatrixIO.writeMatrix(otherFinalCorrelation,
//                                     coalsMatrixFile,
//                                     Format.SVDLIBC_SPARSE_BINARY);
//
//                MatrixFile processedSpace = new MatrixFile(
//                        coalsMatrixFile, Format.SVDLIBC_SPARSE_BINARY);
//                reducer.factorize(processedSpace, reducedDimensions);
//
//                otherFinalCorrelation = reducer.dataClasses();
//            } catch (IOException ioe) {
//                throw new IOError(ioe);
//            }
//            COALS_LOGGER.info("Done reducing using SVD.");
        }
    }

    /**
     * Returns a {@link edu.ucla.sspace.matrix.Matrix} that contains {@code maxWords} and {@code
     * maxDimensions} columns.  If {@code maxWords} is 0, then all words will be
     * returned in the semantic {@link edu.ucla.sspace.matrix.Matrix}.  If {@code maxDimensions} is
     * larger than the number of observed features, then all observed features
     * will be maintained.  The resulting rows and columns are both ordred based
     * on the frequency of each term, in descending order, {@code termToIndex}
     * is modified to account for these changed.
     */
    private void buildMatrix(int maxWords, int maxDimensions) {
        // Convert the vectors in the semantic map to a matrix.
        SparseDoubleVector[] vectorListOther =
                new SparseDoubleVector[otherToSemantics.size()];
        for (Map.Entry<String, SparseDoubleVector> e :
                otherToSemantics.entrySet())
            vectorListOther[getIndexForOther(e.getKey())] = e.getValue();
        SparseMatrix matrixOther = Matrices.asSparseMatrix(
                Arrays.asList(vectorListOther));
//        int[] allRows = new int[otherToSemantics.size()];
//        for (int i = 0; i < allRows.length; i++) {
//            allRows[i] = i;
//        }
//        int[] allCols = new int[otherToSemantics.size()];
//        for (int i = 0; i < allRows.length; i++) {
//            allRows[i] = i;
//        }
        //otherFinalCorrelation = new CellMaskedSparseMatrix(matrixCompound, allRows, colMask);
        otherFinalCorrelation = matrixOther;
        otherToSemantics = null;
        // Return a masked version of the original matrix.
        //return new CellMaskedSparseMatrix(matrix, rowMask, colMask);
    }

    public void loadStatistics(String dirPlusSpaceNameNoExtensionLoaded) {
        this.transform = (CorrelationTransformExtendedSerialiazable)Serializer.deserialiazeData(dirPlusSpaceNameNoExtensionLoaded + ".trans");
        this.termToIndex = (Map<String, Integer>)Serializer.deserialiazeData(dirPlusSpaceNameNoExtensionLoaded + ".indmap");
    }

    private class EntryComp
            implements Comparator<Map.Entry<String,AtomicInteger>> {
        public int compare(Map.Entry<String, AtomicInteger> o1,
                           Map.Entry<String, AtomicInteger> o2) {
            int diff = o2.getValue().get() - o1.getValue().get();
            return (diff != 0) ? diff : o2.getKey().compareTo(o1.getKey());
        }
    }
}
