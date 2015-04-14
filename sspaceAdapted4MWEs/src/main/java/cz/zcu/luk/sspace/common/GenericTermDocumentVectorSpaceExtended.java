/*
 *  LK Copy of GenericTermDocumentVectorSpaceCompounds
 *  Only but processDocument changed!
 */

package cz.zcu.luk.sspace.common;

import edu.ucla.sspace.basis.BasisMapping;
import edu.ucla.sspace.basis.StringBasisMapping;
import edu.ucla.sspace.common.SemanticSpace;
import edu.ucla.sspace.matrix.Matrices;
import edu.ucla.sspace.matrix.Matrix;
import edu.ucla.sspace.matrix.MatrixBuilder;
import edu.ucla.sspace.matrix.MatrixFile;
import edu.ucla.sspace.matrix.SparseHashMatrix;
import edu.ucla.sspace.util.*;
import edu.ucla.sspace.vector.DoubleVector;
import edu.ucla.sspace.vector.SparseHashDoubleVector;
import edu.ucla.sspace.vector.Vector;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOError;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Logger;

import cz.zcu.luk.sspace.matrix.*;
import cz.zcu.luk.sspace.text.IteratorFactoryStopwords;


/**
 * This base class centralizes much of the common text processing needed for
 * term-document based {@link edu.ucla.sspace.common.SemanticSpace}s.  It processes a document by
 * tokenizing all of the provided text and counting the term occurrences within
 * the document.  Each column in these spaces represent a document, and the
 * column values initially represent the number of occurrences for each word.
 * After all documents are processed, the word space can be modified with one of
 * the many {@link edu.ucla.sspace.matrix.Matrix} {@link edu.ucla.sspace.matrix.Transform} classes.  The transform, if
 * provided, will be used to rescore each term document occurrence count.
 * Typically, this reweighting is typically done to increase the score for
 * important and distinguishing terms while less salient terms, such as stop
 * words, are given a lower score.  After calling {@link
 * #processSpace(edu.ucla.sspace.matrix.Transform) processSpace}, sub classes should call assign thei
 * final data matrix to {@code wordSpace}.  This final matrix should maintain
 * the same row ordering, but the column ordering and dimensionality can be
 * modified in any way.
 *
 * <p>
 *
 * This class is thread-safe for concurrent calls of {@link
 * #processDocument(java.io.BufferedReader) processDocument}.  Once {@link
 * #processSpace(edu.ucla.sspace.matrix.Transform) processSpace} has been called, no further calls to
 * {@link #processDocument(java.io.BufferedReader) processDocument} should be made.
 * This implementation does not support access to the semantic vectors until
 * after {@link #processSpace(java.util.Properties) processSpace} has been called.
 *
 * @see edu.ucla.sspace.matrix.Transform
 * @see edu.ucla.sspace.matrix.SVD
 *
 * @author Keith Stevens
 */
public abstract class GenericTermDocumentVectorSpaceExtended
        implements SemanticSpace, java.io.Serializable {

    // LK added..
    private Set<String> addedWords;

    private static final long serialVersionUID = 1L;

    protected static final Logger LOG =
        Logger.getLogger(GenericTermDocumentVectorSpaceExtended.class.getName());

    /**
     * A mapping from a word to the row index in the that word-document matrix
     * that contains occurrence counts for that word.
     */
    protected final BasisMapping<String, String> termToIndex;

    /**
     * The counter for recording the current number of documents observed.
     * Subclasses can use this for any reporting.
     */
    protected final AtomicInteger documentCounter;

    /**
     * The builder used to construct the term-document matrix as new documents
     * are processed.
     */
    private transient MatrixBuilder termDocumentMatrixBuilder;

    /**
     * If true, the first token in each document is considered to be a document
     * header.
     */
    private final boolean readHeaderToken;

    /**
     * The word space of the term document based word space model.  If the word
     * space is reduced, it is the left factor matrix of the SVD of the
     * word-document matrix.  This matrix is only available after the {@link
     * #processSpace(edu.ucla.sspace.matrix.Transform) processSpace}
     * method has been called.
     */
    protected Matrix wordSpace;

    /**
     * LK added..
     * A mapping from compound to index number.
     */
    protected BasisMapping<String, String> addedWordToIndex;
    /**
     * LK added..
     * The builder used to construct the compound-document matrix as new documents
     * are processed.
     */
    private transient MatrixBuilder addedWordsDocumentMatrixBuilder;
    /**
     * LK added..
     */
    protected Matrix addedWordsSpace;

    /**
     * Constructs the {@code GenericTermDocumentVectorSpace}.
     *
     * @throws java.io.IOException if this instance encounters any errors when creatng
     *         the backing array files required for processing
     */
    public GenericTermDocumentVectorSpaceExtended() throws IOException {
        this(false,new StringBasisMapping(),Matrices.getMatrixBuilderForSVD(), Matrices.getMatrixBuilderForSVD(), null);
    }

    /**
     * Constructs the {@code GenericTermDocumentVectorSpace} using the provided
     * objects for processing.
     *
     * @param readHeaderToken If true, the first token of each document will be
     *        read and passed to {@link #handleDocumentHeader(int, String)
     *        handleDocumentHeader}, which by default discards the header.
     * @param termToIndex The {@link edu.ucla.sspace.basis.BasisMapping} used to map strings to
     *        indices.
     * @param termDocumentMatrixBuilder The {@link edu.ucla.sspace.matrix.MatrixBuilder} used to write
     *        document vectors to disk which later get processed in {@link
     *        #processSpace(java.util.Properties) processSpace}.
     *
     * @throws java.io.IOException if this instance encounters any errors when creatng
     *         the backing array files required for processing
     */
    public GenericTermDocumentVectorSpaceExtended(
            boolean readHeaderToken,
            BasisMapping<String, String> termToIndex,
            MatrixBuilder termDocumentMatrixBuilder,
            MatrixBuilder addedWordsDocumentMatrixBuilder,
            Set<String> addedWords) throws IOException {
        this.readHeaderToken = readHeaderToken;
        this.termToIndex = termToIndex;
        this.addedWords = addedWords;
        documentCounter = new AtomicInteger(0);

        // LK added..
        this.addedWordToIndex = new StringBasisMapping();
        this.addedWordsDocumentMatrixBuilder = addedWordsDocumentMatrixBuilder;

        this.termDocumentMatrixBuilder = termDocumentMatrixBuilder;

        wordSpace = null;
    }

    private boolean isAddedWord(String possibleAddedWord) {
        if (addedWords.contains(possibleAddedWord)) return true;
        else return false;
    }

    /**
     * Tokenizes the document using the {@link edu.ucla.sspace.text.IteratorFactory} and updates the
     * term-document frequency counts.
     *
     * <p>
     *
     * This method is thread-safe and may be called in parallel with separate
     * documents to speed up overall processing time.
     *
     * @param document {@inheritDoc}
     */
    public void processDocument(BufferedReader document) throws IOException {
        // Create a mapping for each term that is seen in the document to the
        // number of times it has been seen.  This mapping would more elegantly
        // be a SparseArray<Integer> however, the length of the sparse array
        // isn't known ahead of time, which prevents it being used by the
        // MatrixBuilder.  Note that the SparseArray implementation would also
        // incur an additional performance hit since each word would have to be
        // converted to its index form for each occurrence, which results in a
        // double Map look-up.
        Counter<String> termCounts = new ObjectCounter<String>();

        // LK added
        Counter<String> addedWordsCounts = new ObjectCounter<String>();

        Iterator<String> documentTokens = IteratorFactoryStopwords.tokenize(document);

        // Increaes the count of documents observed so far.
        int docCount = documentCounter.getAndAdd(1);

        // If the first token is to be interpreted as a document header read it.
        if (readHeaderToken)
            handleDocumentHeader(docCount, documentTokens.next());

        // If the document is empty, skip it
        if (!documentTokens.hasNext())
            return;

        // For each word in the text document, keep a count of how many times it
        // has occurred
        while (documentTokens.hasNext()) {
            // LK changed
            //String word = documentTokens.next();
            String word;
            String possibleWord = documentTokens.next();
            boolean wordInSpace = false;
            if (possibleWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
                word = IteratorFactoryStopwords.EMPTY_TOKEN;
            }
            else {
                word = possibleWord;
                wordInSpace = true;
            }

            // LK added.. for addedWords..
            if (!wordInSpace && addedWords != null) {
                String realWord = possibleWord.substring(0, possibleWord.length()-IteratorFactoryStopwords.STOPWORD_FLAG.length()); // store the stopword!
                if (isAddedWord(realWord)) {
                    addedWordToIndex.getDimension(realWord);

                    addedWordsCounts.count(realWord);
                }
            }

            // Skip added empty tokens for words that have been filtered out
            if (word.equals(IteratorFactoryStopwords.EMPTY_TOKEN))
                continue;

            // Add the term to the total list of terms to ensure it has a proper
            // index.  If the term was already added, this method is a no-op
            termToIndex.getDimension(word);

            termCounts.count(word);
        }

        document.close();

        // Check that we actually loaded in some terms before we increase the
        // documentIndex. This is done after increasing the document count since
        // some configurations may need the document order preserved, for
        // example, if each document corresponds to some cluster assignment.
        // LK changed
        if (termCounts.size() == 0 && addedWordsCounts.size() == 0)
            return;

        // Get the total number of terms encountered so far, including any new
        // unique terms found in the most recent document
        int totalNumberOfUniqueWords = termToIndex.numDimensions();

        // Convert the Map count to a SparseArray
        SparseArray<Integer> documentColumn =
            new SparseIntHashArray(totalNumberOfUniqueWords);
        for (Map.Entry<String,Integer> e : termCounts)
            documentColumn.set(
                    termToIndex.getDimension(e.getKey()), e.getValue());

        // Update the term-document matrix with the results of processing the
        // document.
        termDocumentMatrixBuilder.addColumn(documentColumn);

        // LK added.. for compounds..
        int totalNumberOfUniqueCompounds = addedWordToIndex.numDimensions();
        // Convert the Map count to a SparseArray
        SparseArray<Integer> documentColumnForAddedWords =
                new SparseIntHashArray(totalNumberOfUniqueCompounds);
        for (Map.Entry<String,Integer> e : addedWordsCounts) {
            documentColumnForAddedWords.set(
                    addedWordToIndex.getDimension(e.getKey()), e.getValue());
        }
        addedWordsDocumentMatrixBuilder.addColumn(documentColumnForAddedWords);
    }

//    /**
//     * {@inheritDoc}
//     */
//    public Set<String> getWords() {
//        return Collections.unmodifiableSet(termToIndex.keySet());
//    }

    /**
     * LK changed.. see the original method above
     * {@inheritDoc}
     */
    public Set<String> getWords() {
        Set<String> words = new LinkedHashSet<String>();
        words.addAll(termToIndex.keySet());
        if (addedWords != null) {
            words.addAll(addedWordToIndex.keySet());
        }
        return words;
    }

//    /**
//     * {@inheritDoc}
//     */
//    public Vector getVector(String word) {
//        // determine the index for the word
//        int index = termToIndex.getDimension(word);
//
//        return (index < 0) ? null : wordSpace.getRowVector(index);
//    }

    /**
     * LK changed.. see the original method above
     * {@inheritDoc}
     */
    public Vector getVector(String word) {
        // determine the index for the word
        int index = termToIndex.getDimension(word);

        if (index >= 0) {
            return wordSpace.getRowVector(index);
        }
        if (addedWords != null) {
            index = addedWordToIndex.getDimension(word);
            if (index < 0) {
                return null;
            }
            else {
                return addedWordsSpace.getRowVector(index);
            }
        }
        else {
            return null;
        }
    }

    /**
     * {@inheritDoc}
     */
    public int getVectorLength() {
        return wordSpace.columns();
    }

    /**
     * Processes the {@link cz.zcu.luk.sspace.common.GenericTermDocumentVectorSpaceExtended} with the provided
     * {@link edu.ucla.sspace.matrix.Transform} if it is not {@code null} as a {@link edu.ucla.sspace.matrix.MatrixFile}.
     * Otherwise, the raw term document counts are returned.  Sub classes must
     * call this in order to access the term document counts before doing any
     * other processing.
     *
     * @param transform A matrix transform used to rescale the original raw
     *        document counts.  If {@code null} no transform is done.
     */
    protected MatrixFile processSpace(TransformExtended transform) {
        try {
            // first ensure that we are no longer writing to the matrix
            termDocumentMatrixBuilder.finish();

            // Get the finished matrix file from the builder
            File termDocumentMatrix = termDocumentMatrixBuilder.getFile();

            // If a transform was specified, perform the matrix transform.
            if (transform != null) {
                LoggerUtil.info(LOG, "performing %s transform", transform);

                LoggerUtil.verbose(
                        LOG,"stored term-document matrix in format %s at %s",
                        termDocumentMatrixBuilder.getMatrixFormat(),
                        termDocumentMatrix.getAbsolutePath());

                // Convert the raw term counts using the specified transform
                termDocumentMatrix = transform.transform(
                        termDocumentMatrix, 
                        termDocumentMatrixBuilder.getMatrixFormat());

                LoggerUtil.verbose(
                        LOG, "transformed matrix to %s",
                        termDocumentMatrix.getAbsolutePath());
            }

            return new MatrixFile(
                    termDocumentMatrix, 
                    termDocumentMatrixBuilder.getMatrixFormat());
        } catch (IOException ioe) {
            throw new IOError(ioe);
        }
    }

    /**
     * LK added
     * @param transform
     * @return
     */
    protected Matrix processAddedWordsSpace(TransformExtended transform) {
        // first ensure that we are no longer writing to the matrix
        addedWordsDocumentMatrixBuilder.finish();

        // Get the finished matrix file from the builder
        File compoundDocumentMatrix = addedWordsDocumentMatrixBuilder.getFile();

        Matrix matrixAddedWords = null;
        Matrix matrixAddedWordsTransformed = null;

        // If a transform was specified, perform the matrix transform.
        if (transform != null) {
            LoggerUtil.info(LOG, "performing %s transform", transform);

            LoggerUtil.verbose(
                    LOG,"stored term-document matrix in format %s at %s",
                    addedWordsDocumentMatrixBuilder.getMatrixFormat(),
                    compoundDocumentMatrix.getAbsolutePath());

            // Convert the raw term counts using the specified transform
            MatrixFile mFile = new MatrixFile(
                    compoundDocumentMatrix,
                    addedWordsDocumentMatrixBuilder.getMatrixFormat());
            matrixAddedWords = mFile.load();
            matrixAddedWordsTransformed = new SparseHashMatrix(matrixAddedWords.rows(), matrixAddedWords.columns());
            for (int i = 0; i < matrixAddedWords.rows(); i++) {
                DoubleVector dvs = new SparseHashDoubleVector(matrixAddedWords.columns());
                for (int j = 0; j < matrixAddedWords.columns(); j++) {
                    if (matrixAddedWords.get(i, j) != 0) dvs.set(j, matrixAddedWords.get(i, j));
                }

                transform.transformRow(dvs);
                matrixAddedWordsTransformed.setRow(i, dvs);
            }

            LoggerUtil.verbose(
                    LOG, "transformed matrix to %s",
                    compoundDocumentMatrix.getAbsolutePath());
        }
        return matrixAddedWordsTransformed;
    }

    /**
     * Subclasses should override this method if they need to utilize a header
     * token for each document.  Implementations of this method <b>must</b> be
     * thread safe.  The default action is a no-op.
     *
     * @param docIndex The document id assigned to the current document
     * @param documentName The name of the current document.
     */
    protected void handleDocumentHeader(int docIndex, String header) {
    }
}
