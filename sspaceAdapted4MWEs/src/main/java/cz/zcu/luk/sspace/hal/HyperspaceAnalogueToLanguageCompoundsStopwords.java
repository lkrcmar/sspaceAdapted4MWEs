/*
 * Copyright 2009 Alex Nau
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

package cz.zcu.luk.sspace.hal;

import edu.ucla.sspace.basis.BasisMapping;
import edu.ucla.sspace.basis.StringBasisMapping;
import edu.ucla.sspace.common.SemanticSpace;
import edu.ucla.sspace.hal.ConcatenatedSparseDoubleVector;
import edu.ucla.sspace.hal.LinearWeighting;
import edu.ucla.sspace.hal.WeightingFunction;
import edu.ucla.sspace.matrix.MatrixEntropy;
import edu.ucla.sspace.matrix.MatrixEntropy.EntropyStats;
import edu.ucla.sspace.matrix.SparseMatrix;
import edu.ucla.sspace.matrix.YaleSparseMatrix;
import edu.ucla.sspace.util.BoundedSortedMultiMap;
import edu.ucla.sspace.util.MultiMap;
import edu.ucla.sspace.util.Pair;
import edu.ucla.sspace.vector.CompactSparseVector;
import edu.ucla.sspace.vector.SparseDoubleVector;
import edu.ucla.sspace.vector.Vector;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;
import java.util.logging.Logger;

import cz.zcu.luk.sspace.matrix.AtomicGrowingSparseHashMatrix;
import cz.zcu.luk.sspace.text.IteratorFactoryStopwords;


/**
 * A {@code SemanticSpace} implementation of the Hyperspace Analogue to Language
 * (HAL) algorithm described by Lund and Burgess.  This implementation is based
 * on the following paper: <ul>
 *
 *   <li style="font-family:Garamond, Georgia, serif"> K. Lund and C. Burgess
 *     Producing high-dimensional semantic spaces from lexical Co-occurrence
 *     <i>Behavior Research Methods, Instrumentation, and Computers</i>, 28,
 *     pages 203-208, 1996.  Available <a
 *     href="http://locutus.ucr.edu/reprintPDFs/lb96brmic.pdf">here</a>
 *
 *  </ul> See <a href="http://locutus.ucr.edu/Reprints.html">here</a> for
 *  additional papers that use HAL.  <p>
 *
 * HAL is based on recording the co-occurrence of words in a sparse matrix.  HAL
 * also incorporates word order information by treating the co-occurrences of
 * two words <i>x</i> <i>y</i> as being different than <i>y</i> <i>x</i>.  Each
 * word is assigned a unique index in the co-occurrence matrix.  For some word
 * <i>x</i>, when another word <i>x</i> co-occurs before, matrix entry
 * <i>x</i>,<i>y</i> is update.  Similarly, when <i>y</i> co-occurs after, the
 * matrix entry <i>y</i>,<i>x</i> is updated.  Therefore the full semantic
 * vector for any words is its row vector concatenated with its column
 * vector.<p>
 *
 * Typically, the full vectors are used (for an N x N matrix, these are 2*N in
 * length).  However, HAL also offers two posibilities for dimensionality
 * reduction.  Not all columns provide equal amount of information that can be
 * used to distinguish the meanings of the words.  Specifically, the information
 * theoretic <a
 * href="http://en.wikipedia.org/wiki/Information_entropy">entropy</a> of each
 * column can be calculated as a way of ordering the columns by their
 * importance.  Using this ranking, either a fixed number of columns may be
 * retained, or a threshold may be set to filter out low-entropy columns.
 *
 * <p>
 *
 * A {@link cz.zcu.luk.sspace.hal.HyperspaceAnalogueToLanguageCompoundsStopwords} model is defined by four parameters.
 * The default constructor uses reasonable parameters that match those mentioned
 * in the original publication.  For alternate models, appropriate values must
 * be passed in through the full constructor.  The four parameters are:
 *
 * <dl style="margin-left: 1em">
 *
 * <dt> <i>Parameter:</i> {@code windowSize} </br>
 *      <i>Default:</i> 5
 *
 * <dd style="padding-top: .5em">This parameter sets size of the sliding
 * co-occurrence window such that the {@code windowSize} words before and the {@code
 * windowSize} words after the focus word will be used to count co-occurances.
 * This model always uses symmetric windows.
 *
 * <dt> <i>Property:</i> {@code weighting} </br>
 *      <i>Default:</i> {@link edu.ucla.sspace.hal.LinearWeighting}
 *
 * <dd style="padding-top: .5em">This parameter sets the {@link
 * edu.ucla.sspace.hal.WeightingFunction} used to weight co-occurrences between two words based on
 * the number of interleaving words, i.e. the distance between the two words in
 * the sliding window.  HAL traditionally uses a ramped, linear weighting where
 * those words occurring closets receive more weight, with a linear decrease
 * based on distance.
 *
 * <dt> <i>Property:</i> {@code retainColumns} </br>
 *      <i>Default:</i> -1
 *
 * <dd style="padding-top: .5em">If set to a positive value, this
 * parameter enables dimensionality reduction by retaining only {@code
 * retainColumns} columns.  Columns will be ordered according to their entropy,
 * and the {@code retainColumns} columns with the highest entropy will be
 * retained.  This parameter cannot be set in conjunction with {@code columnThreshold}
 *
 * <dt> <i>Property:</i> {@code columnThreshold} </br>
 *      <i>Default:</i> -1
 *
 * <dd style="padding-top: .5em">If set to a positive value, this parameter enables dimensionality
 *      reduction by retaining only those columns whose entropy is above the
 *      specified threshold.  This parameter may not be set concurrently with
 *      {@code retainColumns}.
 *
 * </dl>
 *
 * <p>
 *
 * For models that require a non-symmetric window, a special {@link
 * edu.ucla.sspace.hal.WeightingFunction} can be used which assigns a weight of {@code 0} to
 * co-occurrences that match the non-symmetric window size.
 *
 * @author Alex Nau
 * @author David Jurgens
 *
 * @see edu.ucla.sspace.common.SemanticSpace
 * @see edu.ucla.sspace.hal.WeightingFunction
 */
public class HyperspaceAnalogueToLanguageCompoundsStopwords implements SemanticSpace {

    /**
     * Logger for HAL
     */
    private static final Logger LOGGER =
            Logger.getLogger(HyperspaceAnalogueToLanguageCompoundsStopwords.class.getName());

    /**
     * The default number of words before and after the focus word to include
     */
    public static final int DEFAULT_WINDOW_SIZE = 5;

    /**
     * A mapping from terms to initial indices in the co-occurrence matrix.
     */
    private final BasisMapping<String, String> termToIndex;

    /**
     * The number of words to consider in one direction to create the symmetric
     * window
     */
    private final int windowSize;

    /**
     * If set to a positive value, this parameter sets a threshold on the
     * columns to be retained.  This cannot be used in conjunction with {@code
     * retainColumns}.
     */
    private final double columnThreshold;

    /**
     * If set to a positive value, this parameter sets the number of highest
     * entropy columns to retain.  This cannot be used in conjunction with
     * {@code columnThreshold}.
     */
    private final int retainColumns;

    /**
     * The type of weight to apply to a the co-occurrence word based on its
     * relative location
     */
    private final WeightingFunction weighting;

    /**
     * The number that keeps track of the index values of words
     */
    private int wordIndexCounter;

    /**
     * The matrix used for storing weight co-occurrence statistics of those
     * words that occur both before and after.
     */
    private AtomicGrowingSparseHashMatrix cooccurrenceMatrix;

    /**
     * The reduced matrix, if columns are to be dropped.
     */
    private SparseMatrix reduced;

    private Set<String> compounds;

    /**
     * Constructs a new instance using the default parameters used in the
     * original publication.
     */
    public HyperspaceAnalogueToLanguageCompoundsStopwords() {
        this(new StringBasisMapping(),
                DEFAULT_WINDOW_SIZE,
                new LinearWeighting(),
                -1d,
                -1, null);
    }

    /**
     * Constructs a {@link cz.zcu.luk.sspace.hal.HyperspaceAnalogueToLanguageCompoundsStopwords} instance using the
     * provided parameters.
     *
     * @param basis A mapping from tokens to dimensions.  This may be preset to
     *              ensure a known column ordering or to prevent particular
     *              tokens from being analyzed.
     * @param windowSize a positive integer that specifies the size of a
     *                   symmetric sliding co-occurence window.
     * @param weightFunction A {@link edu.ucla.sspace.hal.WeightingFunction} that scores a
     *                       co-occurrence between two words based on the
     *                       distance between them in the co-occurrence window.
     * @param columnThreshold An optional, positive threshold.  If set, only
     *                        columns with an entropy above this value will be
     *                        retained as dimensions.  Note, this cannot be used
     *                        in conjunction with {@code retainColumns}.
     * @param retainColumns An optional, positive number of columns.  If set,
     *                      only the {@code retainColumns} columns with the
     *                      highest entropy will be retained. Note that this
     *                      cannot be used in conjunction with {@code
     *                      columnThreshold}.
     * @throws IllegalArgumentException If: either{@code basis} or {@code
     *                                  weightFunction} are {@code null}.
     *                                  {@code windowSize} is non-positive, or
     *                                  {@code columnThreshold} and {@code
     *                                  retainColumns} are both set to positive
     *                                  values.
     */
    public HyperspaceAnalogueToLanguageCompoundsStopwords(BasisMapping<String, String> basis,
                                                          int windowSize,
                                                          WeightingFunction weightFunction,
                                                          double columnThreshold,
                                                          int retainColumns,
                                                          Set<String> compounds
    ) {
        this.cooccurrenceMatrix = new AtomicGrowingSparseHashMatrix();
        this.termToIndex = basis;
        this.windowSize = windowSize;
        this.weighting = weightFunction;
        this.columnThreshold = columnThreshold;
        this.retainColumns = retainColumns;
        reduced = null;
        wordIndexCounter = 0;
        this.compounds = compounds;

        // Validate all the parameters early to prevent processing with invalid
        // parameters.
        if (basis == null)
            throw new IllegalArgumentException(
                    "basis mapping must not be null");
        if (weightFunction == null)
            throw new IllegalArgumentException(
                    "weightFunction must not be null");
        if (windowSize <= 0)
            throw new IllegalArgumentException(
                    "Window size must be a positive, non-negative integer.\n" +
                            "Given: " + windowSize);
        if (columnThreshold > -1d && retainColumns > 0)
            throw new IllegalArgumentException(
                    "columnThreshold and retainColumns cannot both be active.\n" +
                            "columnThreshold: " + columnThreshold + "\n" +
                            "retainColumns: " + retainColumns+ "\n");
    }

    private boolean isCompound(String possibleCompound) {
        if (compounds.contains(possibleCompound)) return true;
        else return false;
    }

    private String mapTrigramCompound(String possibleCompound) {
        String[] wordsInTrigram = possibleCompound.split(" ");
        // remove the word in the middle
        return wordsInTrigram[0] + " " + wordsInTrigram[2];
    }

    /**
     * {@inheritDoc}
     */
    public void  processDocument(BufferedReader document) throws IOException {
        Queue<String> nextWords = new ArrayDeque<String>();
        Queue<String> prevWords = new ArrayDeque<String>();

        // LK added
        String removed = null;
        String removedOld = null;

        Iterator<String> documentTokens =
                IteratorFactoryStopwords.tokenizeOrdered(document);

        String focus = null;

        // Rather than updating the matrix every time an occurrence is seen,
        // keep a thread-local count of what needs to be modified in the matrix
        // and update after the document has been processed.  This saves
        // potential contention from concurrent writes.
        Map<Pair<Integer>,Double> matrixEntryToCount =
                new HashMap<Pair<Integer>,Double>();

        //Load the first windowSize words into the Queue
        for(int i = 0;  i < windowSize && documentTokens.hasNext(); i++)
            nextWords.offer(documentTokens.next());

        while(!nextWords.isEmpty()) {

            // Load the top of the nextWords Queue into the focus word
            // focus = nextWords.remove();
            String possibleWord = nextWords.remove();
            focus = getWord(possibleWord);

            // Add the next word to nextWords queue (if possible)
            if (documentTokens.hasNext())
                nextWords.offer(documentTokens.next());

            // If the filter does not accept this word, skip the semantic
            // processing, continue with the next word
            if (!focus.equals(IteratorFactoryStopwords.EMPTY_TOKEN)) {
                int focusIndex = termToIndex.getDimension(focus);
                // Only process co-occurrences with words with non-negative
                // dimensions.
                if (focusIndex >= 0) {
                    // in front of the focus word
                    int wordDistance = -windowSize + (windowSize - prevWords.size());

                    // LK changed.. added removed, removedOld
                    addTokens(prevWords, focusIndex, wordDistance, matrixEntryToCount, removed, removedOld);
                }
            }

            // LK change - added
            // if focus word and words following create a compound.. add prev words to the compound..
            if (!nextWords.isEmpty() && compounds != null) { // if some words follow and if compounds should be taken into account
                String[] nextWordsArr = nextWords.toArray(new String[0]);
                if (nextWordsArr.length > 1) { // trigram compound
                    String possibleCompoundTrigram = focus + " " + getRealWord(nextWordsArr[0]) + " " + getWord(nextWordsArr[1]);
                    if (isCompound(possibleCompoundTrigram)) {
                        String compoundMapped = mapTrigramCompound(possibleCompoundTrigram);
                        int focusIndex = termToIndex.getDimension(compoundMapped);

                        if (focusIndex >= 0) {
                            // in front of the focus word
                            int wordDistance = -windowSize + (windowSize - prevWords.size());
                            addTokens(prevWords, focusIndex, wordDistance, matrixEntryToCount, removed, removedOld);
                        }
                    }
                }
                if (isCompound((focus + " " + getWord(nextWordsArr[0])))) { // bigram compound
                    int focusIndex = termToIndex.getDimension(focus + " " + getWord(nextWordsArr[0]));
                    if (focusIndex >= 0) {
                        // in front of the focus word
                        int wordDistance = -windowSize + (windowSize - prevWords.size());
                        addTokens(prevWords, focusIndex, wordDistance, matrixEntryToCount, removed, removedOld);
                    }
                }
            }
//            for (int i = 0; i < termToIndex.keySet().size(); i++) {
//                System.out.print(i + " " + termToIndex.getDimensionDescription(i) + " X ");
//            }
//            System.out.println();

            // last, put this focus word in the prev words and shift off the
            // front if it is larger than the window
            prevWords.offer(possibleWord);
            if (prevWords.size() > windowSize) {
                // LK change
                removedOld = removed;
                removed = prevWords.remove();
            }
        }

        // Once the document has been processed, update the co-occurrence matrix
        // accordingly.
        for (Map.Entry<Pair<Integer>,Double> e : matrixEntryToCount.entrySet()){
            Pair<Integer> p = e.getKey();
            cooccurrenceMatrix.addAndGet(p.x, p.y, e.getValue());
        }
    }

    private String getRealWord(String possibleWord) {
        if (possibleWord == null) {
            return null;
        }

        String wordReal;
        if (possibleWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
            wordReal = possibleWord.substring(0, possibleWord.length()-IteratorFactoryStopwords.STOPWORD_FLAG.length()); // store the stopword!
        }
        else {
            wordReal = possibleWord;
        }
        return wordReal;
    }

    private String getWord(String possibleWord) {
        if (possibleWord == null) {
            return null;
        }

        String word;
        if (possibleWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
            word = IteratorFactoryStopwords.EMPTY_TOKEN;
        }
        else {
            word = possibleWord;
        }
        return word;
    }

    /**
     * LK change - added removed and oldRemoved
     *
     * Adds co-occurrence counts between the list of previous words in {@code
     * words} and the focus word represented by {@code focusIndex} which start
     * at {@code distance} tokens away from the focus word.  All Counts will be
     * added into {@code matrixEntryToCount}.
     */
    private void addTokens(Queue<String> words,
                           int focusIndex,
                           int distance,
                           Map<Pair<Integer>, Double> matrixEntryToCount,
                           String removed,
                           String removedOld) {
        String[] wordsArr = words.toArray(new String[0]);
        for (int i = 0; i < wordsArr.length; i++) {
            String word = getWord(wordsArr[i]);
            //System.out.print(word + " ");
            // skip adding co-occurence values for words that are not
            // accepted by the filter
            if (!word.equals(IteratorFactoryStopwords.EMPTY_TOKEN)) {
                int index = termToIndex.getDimension(word);
                if (index >= 0) {
                    addPairToMatrixEntryToCount(focusIndex, index, distance, windowSize, matrixEntryToCount);
                }
            }

            if (compounds != null) { // compounds should be taken into account
                String possibleCompound = null;
                if (i == 0) {
                    possibleCompound = getWord(removedOld) + " " + getRealWord(removed) + " " + getWord(wordsArr[0]);
                    if (isCompound(possibleCompound)) {
                        int compoundIndex = termToIndex.getDimension(mapTrigramCompound(possibleCompound));
                        if (compoundIndex >= 0) {
                            addPairToMatrixEntryToCount(focusIndex, compoundIndex, distance, windowSize, matrixEntryToCount);
                        }
                    }
                    possibleCompound = getWord(removed) + " " + getWord(wordsArr[0]);
                    if (isCompound(possibleCompound)) {
                        int compoundIndex = termToIndex.getDimension(possibleCompound);
                        if (compoundIndex >= 0) {
                            addPairToMatrixEntryToCount(focusIndex, compoundIndex, distance, windowSize, matrixEntryToCount);
                        }
                    }
                }
                else if (i == 1) {
                    possibleCompound = getWord(removed) + " " + getRealWord(wordsArr[0]) + " " + getWord(wordsArr[1]);
                    if (isCompound(possibleCompound)) {
                        int compoundIndex = termToIndex.getDimension(mapTrigramCompound(possibleCompound));
                        if (compoundIndex >= 0) {
                            addPairToMatrixEntryToCount(focusIndex, compoundIndex, distance, windowSize, matrixEntryToCount);
                        }
                    }
                    possibleCompound = getWord(wordsArr[0]) + " " + getWord(wordsArr[1]);
                    if (isCompound(possibleCompound)) {
                        int compoundIndex = termToIndex.getDimension(possibleCompound);
                        if (compoundIndex >= 0) {
                            addPairToMatrixEntryToCount(focusIndex, compoundIndex, distance, windowSize, matrixEntryToCount);
                        }
                    }
                }
                else {
                    possibleCompound = getWord(wordsArr[i-2]) + " " + getRealWord(wordsArr[i-1]) + " " + getWord(wordsArr[i]);
                    if (isCompound(possibleCompound)) {
                        int compoundIndex = termToIndex.getDimension(mapTrigramCompound(possibleCompound));
                        if (compoundIndex >= 0) {
                            addPairToMatrixEntryToCount(focusIndex, compoundIndex, distance, windowSize, matrixEntryToCount);
                        }
                    }
                    possibleCompound = getWord(wordsArr[i-1]) + " " + getWord(wordsArr[i]);
                    if (isCompound(possibleCompound)) {
                        int compoundIndex = termToIndex.getDimension(possibleCompound);
                        if (compoundIndex >= 0) {
                            addPairToMatrixEntryToCount(focusIndex, compoundIndex, distance, windowSize, matrixEntryToCount);
                        }
                    }
                }
            }
            distance++;
        }
        //System.out.println();
    }

    private void addPairToMatrixEntryToCount(int focusIndex, int neighborIndex, int distance, int windowSize,
                                             Map<Pair<Integer>,Double> matrixEntryToCount) {
        // Get the current number of times that the focus word has
        // co-occurred with this word before after it.  Weight the
        // word appropriately baed on distance
        // LK change - swapped index and focusIndex to store predecessors before followers into matrix..
        Pair<Integer> p = new Pair<Integer>(focusIndex, neighborIndex);
        double value = weighting.weight(distance, windowSize);
        Double curCount = matrixEntryToCount.get(p);
        matrixEntryToCount.put(p, (curCount == null)
                ? value : value + curCount);
    }

    /**
     * {@inheritDoc}
     */
    public Set<String> getWords() {
        // If no documents have been processed, it will be empty
        return Collections.unmodifiableSet(termToIndex.keySet());
    }

    /**
     * {@inheritDoc}
     */
    public Vector getVector(String word) {
        Integer index = termToIndex.getDimension(word);
        if (index == null)
            return null;
            // If the matrix hasn't had columns dropped then the returned vector
            // will be the combination of the word's row and column
        else if (reduced == null) {
            // NOTE: the matrix could be asymmetric if the a word has only
            // appeared on one side of a context (its row or column vector would
            // never have been set).  Therefore, check the index with the matrix
            // size first.
            SparseDoubleVector rowVec = (index < cooccurrenceMatrix.rows())
                    ? cooccurrenceMatrix.getRowVectorUnsafe(index)
                    : new CompactSparseVector(termToIndex.numDimensions());
            SparseDoubleVector colVec = (index < cooccurrenceMatrix.columns())
                    ? cooccurrenceMatrix.getColumnVectorUnsafe(index)
                    : new CompactSparseVector(termToIndex.numDimensions());

            return new ConcatenatedSparseDoubleVector(rowVec, colVec);
        }

        // The co-occurrence matrix has had columns dropped so the vector is
        // just the word's row
        return reduced.getRowVector(index);
    }

    /**
     * {@inheritDoc}
     */
    public int getVectorLength() {
        if (cooccurrenceMatrix != null)
            return cooccurrenceMatrix.columns() + cooccurrenceMatrix.rows();
        return reduced.columns();
    }

    /**
     * {@inheritDoc}
     */
    public void processSpace(Properties properties) {
        // Ensure that the bottom right corner of the matrix has a valid value
        // so that we always create a 2 * n set of dimensions in the default
        // case.
        if (cooccurrenceMatrix.get(termToIndex.numDimensions() - 1,
                termToIndex.numDimensions() - 1) == 0d)
            cooccurrenceMatrix.set(termToIndex.numDimensions() - 1,
                    termToIndex.numDimensions() - 1,
                    0d);

        if (columnThreshold > -1d)
            thresholdColumns(columnThreshold);
        if (retainColumns > 0)
            retainOnly(retainColumns);

//        System.out.println();
//        System.out.println();
//        for (int i = 0; i < cooccurrenceMatrix.rows(); i++) {
//            System.out.print(termToIndex.getDimensionDescription(i) + " | ");
//            for (int j = 0; j < cooccurrenceMatrix.columns(); j++) {
//                System.out.print(cooccurrenceMatrix.get(i, j) + " ");
//            }
//            System.out.println();
//        }
//        System.out.println();
//        System.out.println();

        if (compounds != null && retainColumns <= 0) {
            termToIndex.setReadOnly(true);
            retainOnlyNonCompounds();

//            System.out.println();
//            System.out.println();
//            for (int i = 0; i < reduced.rows(); i++) {
//                System.out.print(termToIndex.getDimensionDescription(i) + " | ");
//                for (int j = 0; j < reduced.columns(); j++) {
//                    System.out.print(reduced.get(i, j) + " ");
//                }
//                System.out.println();
//            }
//            System.out.println();
//            System.out.println();
        }
//        for (int i = 0; i < termToIndex.keySet().size(); i++) {
//            System.out.println(i + " " + termToIndex.getDimensionDescription(i));
//        }
    }

    /**
     * LK added
     */
    private Set<Integer> loadCompoundIndices(int words) {
        termToIndex.setReadOnly(true);
        Set<Integer> compoundsIndices = new HashSet<Integer>();
        for (String oneCompound : compounds) {
            int indexC;
            if (oneCompound.split(" ").length == 3) {
                indexC = termToIndex.getDimension(mapTrigramCompound(oneCompound));
            }
            else {
                indexC = termToIndex.getDimension(oneCompound);
            }
            //System.out.println("B: "+indexC);
            if (indexC >= 0) {
                compoundsIndices.add(indexC);
                compoundsIndices.add(indexC + words);
            }
        }
        return compoundsIndices;
    }

    /**
     * LK added
     */
    private void retainOnlyNonCompounds() {
        int words = termToIndex.numDimensions();
        //System.out.println("WORDS COUNT");
        //System.out.println(words);
        //System.out.println("WORDS COUNT");
        Set<Integer> indicesToKeep =
                new HashSet<Integer>();

        // load compounds indices
        Set<Integer> compoundsIndices = loadCompoundIndices(words);

        for (int i = 0; i < words; i++) {
            if (!compoundsIndices.contains(i)) {
                indicesToKeep.add(i);
                indicesToKeep.add(i+words);
            }
        }
//        System.out.println("KEEPING");
//        System.out.println(indicesToKeep.toString());
//        System.out.println("KEEPING");


        LOGGER.info("Retain only NON-COMPOUNDS");

        reduced = retainColumns(indicesToKeep);
        cooccurrenceMatrix = null;
    }

    /**
     * Drops all but the specified number of columns, retaining those that have
     * the highest information theoretic entropy.
     *
     * @param columns the number of columns to keep
     */
    private void retainOnly(int columns) {
        LOGGER.info("Sorting the columns by entropy and computing the top " +
                columns + " columns to retain");

        int words = termToIndex.numDimensions();
        MultiMap<Double,Integer> entropyToIndex =
                new BoundedSortedMultiMap<Double,Integer>(
                        columns, false, true, true);

        // Compute the entropy of each row and column.
        EntropyStats stats = MatrixEntropy.entropy(cooccurrenceMatrix);

        // LK changed.. added condition for compounds.. - do not retain indices of compounds
        // load compounds indices
        if (compounds != null) {
            Set<Integer> compoundsIndices = loadCompoundIndices(words);
            // Add the entropy values for each column.  Since the rows will be
            // concatenated as columns, they represent currently non-existing
            // columns beyond the number of words.
            for (int col = 0; col < words; ++col) {
                if (!compoundsIndices.contains(col)) {
                    entropyToIndex.put(stats.colEntropy[col], col);
                }
            }
            for (int row = 0; row < words; ++row) {
                if (!compoundsIndices.contains(row)) {
                    entropyToIndex.put(stats.rowEntropy[row], row+words);
                }
            }
        }
        else {
            // Add the entropy values for each column.  Since the rows will be
            // concatenated as columns, they represent currently non-existing
            // columns beyond the number of words.
            for (int col = 0; col < words; ++col)
                entropyToIndex.put(stats.colEntropy[col], col);
            for (int row = 0; row < words; ++row)
                entropyToIndex.put(stats.rowEntropy[row], row+words);
        }

        Set<Integer> indicesToKeep =
                new HashSet<Integer>(entropyToIndex.values());

//        System.out.println("KEEPING");
//        System.out.println(indicesToKeep.toString());
//        System.out.println("KEEPING");

        // LK change.. comment extended..
        LOGGER.info("Reducing to " + columns + " highest entropy columns. Retain only NON-COMPOUNDS");

        reduced = retainColumns(indicesToKeep);
        cooccurrenceMatrix = null;
    }

    /**
     * Calculates the entropy of all the columns in the co-occurrence matrix and
     * removes those columns that are below the threshold, setting {@link
     * #reduced} to the remaining columns.
     *
     * @param threshold
     */
    private void thresholdColumns(double threshold) {
        LOGGER.info("Computing the columns which are equal to or above the " +
                "specified threshold");

        int words = termToIndex.numDimensions();
        Set<Integer> colsToRetain = new HashSet<Integer>();

        // Compute the entropy of each row and column.
        EntropyStats stats = MatrixEntropy.entropy(cooccurrenceMatrix);

        // Compare the entropy of each column to the threshold and save and
        // indices that pass the threshold. Since the rows will be concatenated
        // as columns, they represent currently non-existing columns beyond the
        // number of words.
        for (int col = 0; col < words; ++col)
            if (stats.colEntropy[col] >= threshold)
                colsToRetain.add(col);

        for (int row = 0; row < words; ++row)
            if (stats.rowEntropy[row] >= threshold)
                colsToRetain.add(row+words);

        LOGGER.info("Retaining " + colsToRetain.size() + "/" + (words*2) +
                " columns, which passed the threshold of " + threshold);

        reduced = retainColumns(colsToRetain);
        cooccurrenceMatrix = null;
    }

    /**
     * Returns a reduced and concatenated version of {@code cooccurrenceMatrix}
     * which has only the columns specified in {@code indicesToKeep}.
     */
    private SparseMatrix retainColumns(Set<Integer> indicesToKeepSet) {
        // LK change - added - do not want to shuffle columns' indices..
        ArrayList<Integer> indicesToKeep = new ArrayList<Integer>(indicesToKeepSet);
        Collections.sort(indicesToKeep);

        int words = termToIndex.numDimensions();
        int cols = indicesToKeep.size();

        // Create a mapping from old indices to their new index value.
        Map<Integer, Integer> indexMap = new HashMap<Integer, Integer>();
        int newIndex = 0;
        for (Integer index : indicesToKeep)
            indexMap.put(index, newIndex++);

        // Create a reduced matrix that will have only the selected columns in
        // the final space.
        SparseMatrix reduced = new YaleSparseMatrix(
                words, indicesToKeep.size());

        // Iterate over the sparse values in the matrix for added efficiency.
        for (int row = 0; row < words; ++ row) {
            if ((row % 10000) == 0) {
                LOGGER.info("reducing of " + row + " rows done..");
            }
            SparseDoubleVector sv = cooccurrenceMatrix.getRowVector(row);
            for (int col : sv.getNonZeroIndices()) {
                double v = cooccurrenceMatrix.get(row, col);

                // If the original column was retained, get it's new index
                // value and add it to the reduced matrix.
                Integer newColIndex = indexMap.get(col);
                if (newColIndex != null)
                    reduced.set(row, newColIndex, v);

                // If the transposed row column was retained, get it's new index
                // value and add it to the reduced matrix.  This turns the col
                // value into the row and the new index as the column.
                newColIndex = indexMap.get(row + words);
                if (newColIndex != null)
                    reduced.set(col, newColIndex, v);
            }
        }

        // LK added
        LOGGER.info("reducing columns finished..");

        return reduced;
    }

    /**
     * {@inheritDoc}
     */
    public String getSpaceName() {
        return "xHALC" + "_W" + windowSize + "_R" + retainColumns;
    }
}
