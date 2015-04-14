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


import edu.ucla.sspace.common.SemanticSpace;
import edu.ucla.sspace.matrix.ArrayMatrix;
import edu.ucla.sspace.matrix.CellMaskedSparseMatrix;
import edu.ucla.sspace.matrix.DiagonalMatrix;
import edu.ucla.sspace.matrix.Matrices;
import edu.ucla.sspace.matrix.Matrix;
import edu.ucla.sspace.matrix.MatrixFile;
import edu.ucla.sspace.matrix.MatrixIO;
import edu.ucla.sspace.matrix.SparseMatrix;
import edu.ucla.sspace.matrix.MatrixIO.Format;
import edu.ucla.sspace.matrix.factorization.SingularValueDecomposition;
import edu.ucla.sspace.vector.*;
import edu.ucla.sspace.vector.Vector;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOError;
import java.io.IOException;
import java.lang.ref.WeakReference;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Logger;

import cz.zcu.luk.sspace.matrix.*;
import cz.zcu.luk.sspace.text.IteratorFactoryStopwords;


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
public class CoalsCompoundsStopwords implements SemanticSpace {

    // LK added
    private transient WeakReference<Matrix> SigmaInvTimesVtRef;
    private Matrix sigma;

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
    public static final String COALS_SSPACE_NAME = "xCOALSC";
        //"coals-semantic-space";

    /**
     * The logger used to record all output
     */
    private static final Logger COALS_LOGGER =
        Logger.getLogger(CoalsCompoundsStopwords.class.getName());

    /**
     * A mapping from each word to the vector the represents its semantics
     */
    private Map<String, SparseDoubleVector> wordToSemantics;

    /**
     * A mapping from word to index number.
     */
    private Map<String, Integer> termToIndex;

    /**
     * LK added..
     * A mapping from compound to index number.
     */
    private Map<String, Integer> compoundToIndex;
    /**
     * LK added..
     * A mapping from each compound to the vector the represents its semantics
     */
    private Map<String, SparseDoubleVector> compoundToSemantics;
    /**
     * LK added..
     * A counter for keeping track of the index values of compounds.
     */
    private int compoundIndexCounter;
    /**
     * LK added..
     * The final reduced matrix for compounds
     */
    private Matrix finalCorrelationCompounds;

    /**
     * A map containg the total frequency counts of each word.
     */
    private ConcurrentMap<String, AtomicInteger> totalWordFreq;

    /**
     * The final reduced matrix.
     */
    private Matrix finalCorrelation;

    /**
     * Specifies the number of reduced dimensions if the matrix is reduced by
     * SVD.
     */
    private final int reducedDimensions;

    /**
     * The maximum number of words that will be retained by {@link cz.zcu.luk.sspace.coals.CoalsCompoundsStopwords}.
     */
    private final int maxWords;

    /**
     * The maximum number of co-occurring words that will be retained by {@link
     * cz.zcu.luk.sspace.coals.CoalsCompoundsStopwords}.
     */
    private final int maxDimensions;

    /**
     * A counter for keeping track of the index values of words.
     */
    private int wordIndexCounter;

    /**
     * The {@link edu.ucla.sspace.matrix.MatrixFactorization} algorithm that will decompose the word by
     * document feature space into two smaller feature spaces: a word by class
     * feature space and a class by feature space.
     */
    private final SingularValueDecomposition reducer;

    /**
     * The {@link edu.ucla.sspace.matrix.Transform} applied to the co-ocucrrence counts, if not {@code
     * null}.
     */
    private final TransformExtended transform;

    // LK added
    private Set<String> compounds;

    public CoalsCompoundsStopwords(TransformExtended transform, SingularValueDecomposition reducer) {
        this(transform, reducer, DEFAULT_REDUCE_DIMENSIONS,
             DEFAULT_MAX_WORDS, DEFAULT_MAX_DIMENSIONS);
    }

    /**
     * Creats a {@link cz.zcu.luk.sspace.coals.CoalsCompoundsStopwords} instance.
     */
    public CoalsCompoundsStopwords(TransformExtended transform,
                                   SingularValueDecomposition reducer,
                                   int reducedDimensions,
                                   int maxWords,
                                   int maxDimensions) {
        // LK changed
        this(transform, reducer, reducedDimensions, maxWords, maxDimensions, null);
    }

    public CoalsCompoundsStopwords(TransformExtended transform, SingularValueDecomposition reducer,
                                   int reducedDimensions, int maxWords, int maxDimensions, Set<String> compounds) {
        // LK changed
        termToIndex = new HashMap<String, Integer>();
        compoundToIndex = new HashMap<String, Integer>(); // LK added
        totalWordFreq = new ConcurrentHashMap<String, AtomicInteger>();
        wordToSemantics = new HashMap<String, SparseDoubleVector>(1024, 4f);
        compoundToSemantics = new HashMap<String, SparseDoubleVector>(1024, 4f); // LK added
        finalCorrelation = null;
        this.transform = transform;
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
        this.compounds = compounds;
    }
//
//    /**
//     * {@inheritDoc}
//     */
//    public Set<String> getWords() {
//        return termToIndex.keySet();
//    }

    /**
     * LK changed.. see the original version above..
     * {@inheritDoc}
     */
    public Set<String> getWords() {
        Set<String> wordsAndCompounds = new LinkedHashSet<String>();
        wordsAndCompounds.addAll(termToIndex.keySet());
        if (compounds != null) {
            wordsAndCompounds.addAll(compoundToIndex.keySet());
        }
        return wordsAndCompounds;
    }

//    /**
//     * {@inheritDoc}
//     */
//    public Vector getVector(String term) {
//        Integer index = termToIndex.get(term);
//        if (index == null)
//            return null;
//        return Vectors.immutable(
//                finalCorrelation.getRowVector(index.intValue()));
//    }

    /**
     * LK changed.. see the original version above..
     * {@inheritDoc}
     */
    public Vector getVector(String term) {
        Integer index = termToIndex.get(term);
        if (index != null) {
            return Vectors.immutable(
                    finalCorrelation.getRowVector(index.intValue()));
        }
        if (compounds != null) {
            index = compoundToIndex.get(term);
            if (index == null) {
                return null;
            }
            else {
                return Vectors.immutable(
                    finalCorrelationCompounds.getRowVector(index.intValue()));
            }
        }
        else {
            return null;
        }
    }

    public String getSpaceName() {
        String ret = COALS_SSPACE_NAME;
        ret += "_M" + maxWords;
        ret += "_N" + maxDimensions;
        if (reducer != null)
            //ret += "-svd-" + reducedDimensions;
            ret += "_D" + reducedDimensions;
        return ret;
    }

    public int getVectorLength() {
        return finalCorrelation.columns();
    }


    /*
     * LK added
     */
    private boolean isCompound(String possibleCompound) {
        if (compounds.contains(possibleCompound)) return true;
        else return false;
    }

    /**
     * {@inheritDoc}
     */
    public void processDocument(BufferedReader document) throws IOException {
        Map<String, Integer> wordFreq = new HashMap<String, Integer>();
        Map<String, SparseDoubleVector> wordDocSemantics =
            new HashMap<String, SparseDoubleVector>();
        Map<String, SparseDoubleVector> compoundDocSemantics =
                new HashMap<String, SparseDoubleVector>();

        // LK added
        String removed = null;
        String removedOld = null;

        // Setup queues to track the set of previous and next words in a
        // context.
        Queue<String> prevWords = new ArrayDeque<String>();
        Queue<String> nextWords = new ArrayDeque<String>();

        Iterator<String> it = IteratorFactoryStopwords.tokenizeOrdered(document);

        for (int i = 0; i < 4 && it.hasNext(); ++i)
            nextWords.offer(it.next());

        // Compute the co-occurrance statistics of each focus word in the
        // document.
        while (!nextWords.isEmpty()) {

            // Slide over the context by one word.
            if (it.hasNext())
                nextWords.offer(it.next());

            // Get the focus word
            String focusWord;
            // String notStopwordDeleted;
            String possibleFocusWord = nextWords.remove();
            // System.out.println("possibleFocusWord: \""+possibleFocusWord+"\"");
            if (possibleFocusWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
                focusWord = IteratorFactoryStopwords.EMPTY_TOKEN;
            }
            else {
                focusWord = possibleFocusWord;
            }
            // System.out.println("0\""+focusWord+"\" ");

            if (!focusWord.equals(IteratorFactoryStopwords.EMPTY_TOKEN)) {
                getIndexFor(focusWord);

                // Update the frequency count of the focus word.
                Integer focusFreq = wordFreq.get(focusWord);
                wordFreq.put(focusWord, (focusFreq == null)
                        ? 1
                        : 1 + focusFreq.intValue());


                // Get the temprorary semantics for the focus word, create a new
                // vector for them if needed.
                SparseDoubleVector focusSemantics = wordDocSemantics.get(
                        focusWord);
                if (focusSemantics == null) {
                    focusSemantics = new SparseHashDoubleVector(
                            Integer.MAX_VALUE);
                    wordDocSemantics.put(focusWord, focusSemantics);
                }

                // Process the previous words.
                int offset = 4 - prevWords.size();
                for (String possibleWord : prevWords) {
                    String word;
                    if (possibleWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
                        word = IteratorFactoryStopwords.EMPTY_TOKEN;
                    }
                    else {
                        word = possibleWord;
                    }
                    //System.out.print("1\""+word+"\" ");

                    offset++;
                    if (word.equals(IteratorFactoryStopwords.EMPTY_TOKEN))
                        continue;
                    int index = getIndexFor(word);
                    focusSemantics.add(index, offset);
                }
                //System.out.println();

                // Process the next words.
                offset = 5;
                for (String possibleWord : nextWords) {
                    String word;
                    if (possibleWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
                        word = IteratorFactoryStopwords.EMPTY_TOKEN;
                    }
                    else {
                        word = possibleWord;
                    }
                    // System.out.print("2\""+word+"\" ");

                    offset--;
                    if (word.equals(IteratorFactoryStopwords.EMPTY_TOKEN))
                        continue;
                    int index = getIndexFor(word);
                    focusSemantics.add(index, offset);
                }
            }
            // System.out.println();

            // LK added
            if (compounds != null) {
                String[] prevWordsArr = prevWords.toArray(new String[0]);
                String prevPrevWord = null;
                String prevWord = null;
                String prevWordReal = null;
                if (prevWordsArr.length > 1) {
                    String possiblePrevPrevWord = prevWordsArr[prevWordsArr.length - 2];
                    if (possiblePrevPrevWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
                        prevPrevWord = IteratorFactoryStopwords.EMPTY_TOKEN;
                    }
                    else {
                        prevPrevWord = possiblePrevPrevWord;
                    }
                }
                if (prevWordsArr.length > 0) {
                    String possiblePrevWord = prevWordsArr[prevWordsArr.length - 1];
                    if (possiblePrevWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
                        prevWord = IteratorFactoryStopwords.EMPTY_TOKEN;
                        prevWordReal = possiblePrevWord.substring(0, possiblePrevWord.length()-IteratorFactoryStopwords.STOPWORD_FLAG.length());  // store the stopword!
                    }
                    else {
                        prevWord = possiblePrevWord;
                        prevWordReal = possiblePrevWord;
                    }
                }

                // prevPrevWord and prevWord can be null.. it is expected that
                // compound does not contain null strings..
                String possibleCompoundTrigram = prevPrevWord + " " + prevWordReal + " " + focusWord;
                //System.out.println("posComTrigram: \"" +possibleCompoundTrigram+"\"");
                String possibleCompoundBigram = prevWord + " " + focusWord;
                //System.out.println("posComBigram: \"" +possibleCompoundBigram+"\"");
                String compound = null;
                int compoundSize = 0;
                if (isCompound(possibleCompoundTrigram)) {
                    compound = mapTrigramCompound(possibleCompoundTrigram);
                    compoundSize = 3;
                }
                else if (isCompound(possibleCompoundBigram)) {
                    compound = possibleCompoundBigram;
                    compoundSize = 2;
                }
                if (compound != null) {
                    getIndexForCompound(compound);
                    // System.out.println("3\""+compound+"\" ");
//                    // Update the frequency count of the focus word.
//                    Integer focusFreq = wordFreq.get(compound);
//                    wordFreq.put(compound, (focusFreq == null)
//                            ? 1
//                            : 1 + focusFreq.intValue());


                    // Get the temprorary semantics for the focus word, create a new
                    // vector for them if needed.
                    SparseDoubleVector focusSemanticsCompound = compoundDocSemantics.get(
                            compound);
                    if (focusSemanticsCompound == null) {
                        focusSemanticsCompound = new SparseHashDoubleVector(
                                Integer.MAX_VALUE);
                        compoundDocSemantics.put(compound, focusSemanticsCompound);
                    }

                    // Process the previous words.
                    ArrayList<String> prevWordsList = new ArrayList<String>();
                    if (compoundSize == 3) {
                        if (removedOld != null) {
                            prevWordsList.add(removedOld);
                        }
                        if (removed != null) {
                            prevWordsList.add(removed);
                        }
                        for (int i = 0; i < prevWordsArr.length-2; i++) {
                            prevWordsList.add(prevWordsArr[i]);
                        }
                    }
                    else if (compoundSize == 2) {
                        if (removed != null) prevWordsList.add(removed);
                        for (int i = 0; i < prevWordsArr.length-1; i++) {
                            prevWordsList.add(prevWordsArr[i]);
                        }
                    }

                    int offset = 4 - prevWordsList.size();
                    for (String possibleWord : prevWordsList) {
                        String word;
                        if (possibleWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
                            word = IteratorFactoryStopwords.EMPTY_TOKEN;
                        }
                        else {
                            word = possibleWord;
                        }
                        //System.out.println("4\""+word+"\" ");

                        offset++;
                        if (word.equals(IteratorFactoryStopwords.EMPTY_TOKEN))
                            continue;
                        int index = getIndexFor(word);
                        focusSemanticsCompound.add(index, offset);
                    }

                    // Process the next words.
                    offset = 5;
                    for (String possibleWord : nextWords) {
                        String word;
                        if (possibleWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
                            word = IteratorFactoryStopwords.EMPTY_TOKEN;
                        }
                        else {
                            word = possibleWord;
                        }
                        // System.out.println("5\""+word+"\" ");

                        offset--;
                        if (word.equals(IteratorFactoryStopwords.EMPTY_TOKEN))
                            continue;
                        int index = getIndexFor(word);
                        focusSemanticsCompound.add(index, offset);
                    }
                }
            }
            //  System.out.println();
            prevWords.offer(possibleFocusWord);
            if (prevWords.size() > 4) {
                removedOld = removed;
                removed = prevWords.remove();
            }
        }

        // Add the temporary vectors for each word in this document to the
        // actual semantic vectors.
        for (Map.Entry<String, SparseDoubleVector> e :
                wordDocSemantics.entrySet()) {
            SparseDoubleVector focusSemantics = getSemanticVector(
                    e.getKey());
            // Get the non zero indices before hand so that they are cached
            // during the synchronized section.
            focusSemantics.getNonZeroIndices();
            synchronized (focusSemantics) {
                VectorMath.add(focusSemantics, e.getValue());
            }
        }

        // LK added .. add the temporary vectors for each compound in this document to the
        // actual semantic vectors.
        for (Map.Entry<String, SparseDoubleVector> e :
                compoundDocSemantics.entrySet()) {
            SparseDoubleVector focusSemanticsCompound = getSemanticVectorCompound(
                    e.getKey());
            // Get the non zero indices before hand so that they are cached
            // during the synchronized section.
            focusSemanticsCompound.getNonZeroIndices();
            synchronized (focusSemanticsCompound) {
                VectorMath.add(focusSemanticsCompound, e.getValue());
            }
        }

        // Store the total frequency counts of the words seen in this document
        // so far.
        for (Map.Entry<String, Integer> entry : wordFreq.entrySet()) {
            int count = entry.getValue().intValue();
            AtomicInteger freq = totalWordFreq.putIfAbsent(
                    entry.getKey(), new AtomicInteger(count));
            if (freq != null)
                freq.addAndGet(count);
        }
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
    private SparseDoubleVector getSemanticVector(String word) {
        SparseDoubleVector v = wordToSemantics.get(word);
        if (v == null) {
            // lock on the word in case multiple threads attempt to add it at
            // once
            synchronized(this) {
                // recheck in case another thread added it while we were waiting
                // for the lock
                v = wordToSemantics.get(word);
                if (v == null) {
                    v = new CompactSparseVector();
                    wordToSemantics.put(word, v);
                }
            }
        }
        return v;
    }

    private SparseDoubleVector getSemanticVectorCompound(String compound) {
        SparseDoubleVector v = compoundToSemantics.get(compound);
        if (v == null) {
            // lock on the word in case multiple threads attempt to add it at
            // once
            synchronized(this) {
                // recheck in case another thread added it while we were waiting
                // for the lock
                v = compoundToSemantics.get(compound);
                if (v == null) {
                    v = new CompactSparseVector();
                    compoundToSemantics.put(compound, v);
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
    private int getIndexFor(String word) {
        Integer index = termToIndex.get(word);
        if (index == null) {
            synchronized(this) {
                // recheck to see if the term was added while blocking
                index = termToIndex.get(word);
                // if another thread has not already added this word while the
                // current thread was blocking waiting on the lock, then add it.
                if (index == null) {
                    int i = wordIndexCounter++;
                    termToIndex.put(word, i);
                    return i; // avoid the auto-boxing to assign i to index
                }
            }
        }
        return index;
    }

    private int getIndexForCompound(String compound) {
        Integer index = compoundToIndex.get(compound);
        if (index == null) {
            synchronized(this) {
                // recheck to see if the term was added while blocking
                index = compoundToIndex.get(compound);
                // if another thread has not already added this word while the
                // current thread was blocking waiting on the lock, then add it.
                if (index == null) {
                    int i = compoundIndexCounter++;
                    compoundToIndex.put(compound, i);
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
        finalCorrelation = buildMatrix(maxWords, maxDimensions);
        COALS_LOGGER.info("Done dropping dimensions.");

        if (transform != null) {
            COALS_LOGGER.info("Normalizing co-occurrance matrix.");

            // Normalize the matrix using correlation.
            int wordCount = finalCorrelation.rows();
            finalCorrelation = transform.transform(finalCorrelation);

            // LK added.. transform compoundMatrix
            if (compounds != null) {
                SparseDoubleVector[] compoundVectors = new SparseDoubleVector[finalCorrelationCompounds.rows()];
                for (int i = 0; i < finalCorrelationCompounds.rows(); i++) {
                    DoubleVector dv = transform.transformRow(finalCorrelationCompounds.getRowVector(i));
                    if (dv instanceof SparseDoubleVector) {
                        compoundVectors[i] = (SparseDoubleVector)dv;
                    }
                    else {
                        throw new IllegalStateException("Expected sparse format of vector, but it is dense!");
                    }
                }
                finalCorrelationCompounds = Matrices.asSparseMatrix(Arrays.asList(compoundVectors));
            }

            COALS_LOGGER.info("Done normalizing co-occurrance matrix.");
        }

        if (reducer != null) {
            if (reducedDimensions > finalCorrelation.columns())
                throw new IllegalArgumentException(
                        "Cannot reduce to more dimensions than exist");

            COALS_LOGGER.info("Reducing using SVD.");
            try {
                File coalsMatrixFile =
                    File.createTempFile("coals-term-doc-matrix", "dat");
                coalsMatrixFile.deleteOnExit();
                MatrixIO.writeMatrix(finalCorrelation,
                                     coalsMatrixFile,
                                     Format.SVDLIBC_SPARSE_BINARY);

                MatrixFile processedSpace = new MatrixFile(
                        coalsMatrixFile, Format.SVDLIBC_SPARSE_BINARY);
                reducer.factorize(processedSpace, reducedDimensions);

                finalCorrelation = reducer.dataClasses();

                // LK added
                sigma = reducer.getSingularValues();
                Matrix Vt = reducer.getRightVectors();
                finalCorrelationCompounds = computeCompoundSpace(finalCorrelationCompounds, Vt);
            } catch (IOException ioe) {
                throw new IOError(ioe);
            }
            COALS_LOGGER.info("Done reducing using SVD.");
        }
    }

    /**
     * because isCompound on trigram is called always before this method
     * it should not happen that possibleCompound does not consist of 3 words
     *
     * LK added
     */
    private String mapTrigramCompound(String possibleCompound) {
        String[] wordsInTrigram = possibleCompound.split(" ");
        // remove the word in the middle
        return wordsInTrigram[0] + " " + wordsInTrigram[2];
    }

    /**
     * LK added
     */
    private Set<Integer> loadCompoundIndices() {
        Set<Integer> compoundsIndices = new HashSet<Integer>();
        for (String oneCompound : compounds) {
            Integer indexC;
            if (oneCompound.split(" ").length == 3) {
                indexC = termToIndex.get(mapTrigramCompound(oneCompound));
            }
            else {
                indexC = termToIndex.get(oneCompound);
            }
            //System.out.println("B: "+indexC);
            if (indexC != null && indexC >= 0) {
                compoundsIndices.add(indexC);
            }
        }

        return compoundsIndices;
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
    private Matrix buildMatrix(int maxWords, int maxDimensions) {
        // Convert the vectors in the semantic map to a matrix.
        SparseDoubleVector[] vectorList =
            new SparseDoubleVector[wordToSemantics.size()];
        for (Map.Entry<String, SparseDoubleVector> e :
                wordToSemantics.entrySet())
            vectorList[getIndexFor(e.getKey())] = e.getValue();
        SparseMatrix matrix = Matrices.asSparseMatrix(
                Arrays.asList(vectorList));

        // If maxwords was set to 0, save all words.
        if (maxWords == 0 || maxWords > wordToSemantics.size())
            maxWords = wordToSemantics.size();

        COALS_LOGGER.info("Forming the inverse mapping from terms to indices.");
        // Calculate an inverse mapping from index to word since the binary file
        // stores things by index number.
        String[] indexToTerm = new String[termToIndex.size()];
        for (Map.Entry<String, Integer> entry : termToIndex.entrySet())
            indexToTerm[entry.getValue()] = entry.getKey();

        COALS_LOGGER.info("Sorting the terms based on frequency.");
        // Calculate the new indices for each word that will be kept based on
        // the frequency count, where the most frequent word will be first.
        ArrayList<Map.Entry<String, AtomicInteger>> wordCountList =
            new ArrayList<Map.Entry<String, AtomicInteger>>(
                    totalWordFreq.entrySet());
        Collections.sort(wordCountList, new EntryComp());

        // Calculate the new term to index mapping based on the order of the
        // word frequencies.
        COALS_LOGGER.info("Generating the index masks.");

        // Compute the number of dimensions to maintain. 
        int wordCount = (wordCountList.size() > maxDimensions)
            ? maxDimensions 
            : wordCountList.size();

        int[] rowMask = new int[maxWords];
        int[] colMask = new int[wordCount];

        // Create a new vector list to store the word semantics that will be
        // retained.  When this method exits, it will throw away all the other
        // vectors.
        SparseDoubleVector[] newVectorList = new SparseDoubleVector[maxWords];

        // For each of the terms that we have a mapping, add row and column
        // maskings for the indices of the first maxWords terms.  For all other
        // terms, remove the term to index mapping.
        int termCount = 0;
        for (Map.Entry<String, AtomicInteger> entry : wordCountList) {
            Integer oldIndex = termToIndex.get(entry.getKey());

            // Skip any non mapped terms.
            if (oldIndex == null)
                continue;

            // Add a row and/or column mask from the index of this word to it's
            // index in the original matrix.
            if (termCount < maxWords) {
                if (termCount <  wordCount)
                    colMask[termCount] = oldIndex;
                // Add the vector for this reserved word to the new vector list.
                newVectorList[termCount] = vectorList[oldIndex];
                
                // Record the new dimension for this term.
                rowMask[termCount] = termCount;
                termToIndex.put(entry.getKey(), termCount);
                termCount++;
            }
            // Drop all other mappings.
            else
                termToIndex.remove(entry.getKey());
        }

        wordToSemantics = null;
        matrix = Matrices.asSparseMatrix(Arrays.asList(newVectorList));

        // LK added..
        if (compounds != null) {
            // Convert the vectors in the semantic map to a matrix.
            SparseDoubleVector[] vectorListCompounds =
                    new SparseDoubleVector[compoundToSemantics.size()];
            for (Map.Entry<String, SparseDoubleVector> e :
                    compoundToSemantics.entrySet())
                vectorListCompounds[getIndexForCompound(e.getKey())] = e.getValue();
            SparseMatrix matrixCompound = Matrices.asSparseMatrix(
                    Arrays.asList(vectorListCompounds));
            int[] allRows = new int[compoundToSemantics.size()];
            for (int i = 0; i < allRows.length; i++) {
                allRows[i] = i;
            }
            finalCorrelationCompounds = new CellMaskedSparseMatrix(matrixCompound, allRows, colMask);
            compoundToSemantics = null;
        }

        // Return a masked version of the original matrix.
        return new CellMaskedSparseMatrix(matrix, rowMask, colMask);
    }

    /**
     * LK added..
     */
    private Matrix computeCompoundSpace(Matrix processedCompoundSpace, Matrix Vt) {
        System.out.println(" -------------- PROJECTING COMPOUNDS! -------------- ");
        Matrix projectedCompounds = new ArrayMatrix(processedCompoundSpace.rows(), reducedDimensions);

        // Check that we can actually project the document
        if (finalCorrelation == null)
            throw new IllegalStateException(
                    "processSpace has not been called, so the latent document " +
                            "space does not yet exist");

        // Ensure that when we are projecting the new document that we do not
        // add any new terms to this space's basis.
        //termToIndex.setReadOnly(true);
        //super.compoundToIndex.setReadOnly(true);

        for (int k = 0; k < processedCompoundSpace.rows(); k++) {
            DoubleVector compoundInRow = processedCompoundSpace.getRowVector(k);

            SparseDoubleVector compoundInRowSparse = new SparseHashDoubleVector(compoundInRow.length());
            for (int j = 0; j < compoundInRow.length(); j++) {
                if (compoundInRow.get(j) != 0) {
                    compoundInRowSparse.set(j, compoundInRow.get(j));
                }
            }

            // Represent the compound as a 1-column matrix
            Matrix queryAsMatrix = new ArrayMatrix(compoundInRow.length(), 1);
            for (int nz : compoundInRowSparse.getNonZeroIndices())
                queryAsMatrix.set(nz, 0, compoundInRowSparse.get(nz));

            // Project the new compound vector, t, by using
            //
            //   Sigma_k^-1 * Vt_k * t
            //
            // where k is the dimensionality of the LSA space

            Matrix SigmaInvTimesVt = null;

            // We cache the results of the U_k * Sigma_k^-1 multiplication since
            // this will be the same for all projections.
            while (SigmaInvTimesVt == null) {
                if (SigmaInvTimesVtRef != null
                        && ((SigmaInvTimesVt = SigmaInvTimesVtRef.get()) != null))
                    break;

                int rows = sigma.rows();
                double[] sigmaInv = new double[rows];
                for (int i = 0; i < rows; ++i)
                    sigmaInv[i] = 1d / sigma.get(i, i);
                DiagonalMatrix sigmaInvMatrix = new DiagonalMatrix(sigmaInv);

                System.out.println("sigmaInv: " + sigmaInvMatrix.rows() + " " + sigmaInvMatrix.columns());
                System.out.println("Vt: " + Vt.rows() + " " + Vt.columns());
                SigmaInvTimesVt =
                        Matrices.multiply(sigmaInvMatrix, Vt);
                System.out.println("SigmaInvTimesVt: " + SigmaInvTimesVt.rows() + " " + SigmaInvTimesVt.columns());
                // Update the field with the new reference to the precomputed matrix
                SigmaInvTimesVtRef = new WeakReference<Matrix>(SigmaInvTimesVt);
            }

            // Compute the resulting projected vector as a matrix
            Matrix result = Matrices.multiply(SigmaInvTimesVt, queryAsMatrix);

            int rows = result.rows();
            DoubleVector projected = new DenseVector(result.rows());
            for (int i = 0; i < rows; ++i)
                projected.set(i, result.get(i, 0));

            projectedCompounds.setRow(k, projected);
        }

        return projectedCompounds;
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
