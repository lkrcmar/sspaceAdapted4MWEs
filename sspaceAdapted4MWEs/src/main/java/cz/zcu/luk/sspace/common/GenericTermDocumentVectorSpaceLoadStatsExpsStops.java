/*
 * Copyright 2010 Keith Stevens
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

package cz.zcu.luk.sspace.common;

import edu.ucla.sspace.basis.BasisMapping;
import edu.ucla.sspace.basis.StringBasisMapping;
import edu.ucla.sspace.common.SemanticSpace;
import edu.ucla.sspace.matrix.Matrices;
import edu.ucla.sspace.matrix.Matrix;
import edu.ucla.sspace.matrix.MatrixBuilder;
import edu.ucla.sspace.matrix.MatrixFile;
import edu.ucla.sspace.matrix.SparseHashMatrix;
import edu.ucla.sspace.matrix.Transform;
import edu.ucla.sspace.util.*;
import edu.ucla.sspace.vector.SparseDoubleVector;
import edu.ucla.sspace.vector.Vector;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOError;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Logger;

import cz.zcu.luk.sspace.matrix.*;
import cz.zcu.luk.sspace.text.IteratorFactoryStopwords;


/**
 * This base class centralizes much of the common text processing needed for
 * term-document based {@link SemanticSpace}s.  It processes a document by
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
public abstract class GenericTermDocumentVectorSpaceLoadStatsExpsStops
        implements SemanticSpaceLoadStats, java.io.Serializable {

    private static final long serialVersionUID = 1L;

    protected static final Logger LOG =
        Logger.getLogger(GenericTermDocumentVectorSpaceLoadStatsExpsStops.class.getName());

    protected ArrayList<Integer> processedDocNumbers;

    /**
     * A mapping from a word to the row index in the that word-document matrix
     * that contains occurrence counts for that word.
     */
    //protected BasisMapping<String, String> termToIndex; // final keyword removed to allow load statistics..

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
    protected BasisMapping<String, String> expressionToIndex;
    /**
     * LK added..
     * The builder used to construct the compound-document matrix as new documents
     * are processed.
     */
    private transient MatrixBuilder expressionDocumentMatrixBuilder;
    /**
     * LK added..
     */
    protected Matrix expressionSpace;
    // LK added
    private Set<String> expressions;

    /**
     * Constructs the {@code GenericTermDocumentVectorSpace}.
     *
     * @throws java.io.IOException if this instance encounters any errors when creatng
     *         the backing array files required for processing
     */
    public GenericTermDocumentVectorSpaceLoadStatsExpsStops() throws IOException {
        this(false,new StringBasisMapping(),Matrices.getMatrixBuilderForSVD(), null);
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
    public GenericTermDocumentVectorSpaceLoadStatsExpsStops(
            boolean readHeaderToken,
            BasisMapping<String, String> termToIndex,
            MatrixBuilder expressionDocumentMatrixBuilder, Set<String> expressions) throws IOException {
        this.readHeaderToken = readHeaderToken;
        //this.termToIndex = termToIndex;
        documentCounter = new AtomicInteger(0);

        // LK added and modified
        this.expressionToIndex = new StringBasisMapping();
        this.expressionDocumentMatrixBuilder = expressionDocumentMatrixBuilder;

        wordSpace = null;

        this.expressions = expressions;
    }

    /*
    * LK added
    */
    private boolean isDeterminer(String term) {
        // when term is null (e.g. prevprevword not set.. return false
        // false means no testing for trigram compound.. null cannot be in (only null_XX could be)
        if (term == null) return false;

        return (term.equals("the_XX") || term.equals("a_XX") || term.equals("an_XX"));
    }
    private boolean isExpression(String possibleExpression) {
        if (expressions.contains(possibleExpression)) return true;
        else return false;
    }

    /**
     * Tokenizes the document using the {@link cz.zcu.luk.sspace.text.IteratorFactoryStopwords} and updates the
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
        // LK added
        String removed = null;
        String removedOld = null;

        // Create a mapping for each term that is seen in the document to the
        // number of times it has been seen.  This mapping would more elegantly
        // be a SparseArray<Integer> however, the length of the sparse array
        // isn't known ahead of time, which prevents it being used by the
        // MatrixBuilder.  Note that the SparseArray implementation would also
        // incur an additional performance hit since each word would have to be
        // converted to its index form for each occurrence, which results in a
        // double Map look-up.
        Counter<String> expressionCounts = new ObjectCounter<String>();
        //Iterator<String> documentTokens = IteratorFactoryStopwords.tokenize(document);
        // LK change.. otherwise document is considered as empty even if it contains stopwords!
        Iterator<String> documentTokens = IteratorFactoryStopwords.tokenizeOrdered(document);
        //System.out.println(documentTokens.getClass().toString());

        // Increases the count of documents observed so far.
        int docCount = documentCounter.getAndAdd(1);

        // LK added..
        if (!processedDocNumbers.contains(docCount)) {
            return;
        }

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
            if (possibleWord.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
                word = IteratorFactoryStopwords.EMPTY_TOKEN;
            }
            else {
                word = possibleWord;
            }

            // LK added.. for compounds..
            if (expressions != null) {
                String removedOldWord = null;
                if (removedOld != null) {
                    if (removedOld.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
                        removedOldWord = IteratorFactoryStopwords.EMPTY_TOKEN;
                    }
                    else {
                        removedOldWord = removedOld;
                    }
                }
                String removedReal = null;
                String removedWord = null;
                if (removed != null) {
                    if (removed.endsWith(IteratorFactoryStopwords.STOPWORD_FLAG)) {
                        removedWord = IteratorFactoryStopwords.EMPTY_TOKEN;
                        removedReal = removed.substring(0, removed.length()-IteratorFactoryStopwords.STOPWORD_FLAG.length()); // store the stopword!
                    }
                    else {
                        removedWord = removed;
                        removedReal = removed;
                    }
                }

                String possibleExpressionTrigram = removedOldWord + /*" " + prevWordReal +*/ " " + word;
                //System.out.println("posComTrigram: \"" +possibleExpressionTrigram+"\"");
                String possibleExpressionBigram = removedWord + " " + word;
                //System.out.println("posComBigram: \"" +possibleExpressionBigram+"\"");
                String expression = null;
                if (isDeterminer(removedReal) && isExpression(possibleExpressionTrigram)) {
                    expression = possibleExpressionTrigram;
                    expressionToIndex.getDimension(expression);
                    expressionCounts.count(expression);
                }
                else if (isExpression(possibleExpressionBigram)) {
                    expression = possibleExpressionBigram;
                    expressionToIndex.getDimension(expression);
                    expressionCounts.count(expression);
                }
            }

            // do not skip.. set removedOld and removed..
            // Skip added empty tokens for words that have been filtered out
            //if (word.equals(IteratorFactoryStopwords.EMPTY_TOKEN))
            //    continue;

            // LK added
            removedOld = removed;
            removed = possibleWord;
        }

        document.close();

//        // Check that we actually loaded in some terms before we increase the
//        // documentIndex. This is done after increasing the document count since
//        // some configurations may need the document order preserved, for
//        // example, if each document corresponds to some cluster assignment.
//        if (expressionCounts.size() == 0)
//            return;
        // LK changed.. see up - use processedDocumentSet to decide which docs to skip..
        // empty columns also have to be added!

        // Get the total number of terms encountered so far, including any new
        // unique terms found in the most recent document
        int totalNumberOfUniqueExpressions = expressionToIndex.numDimensions();

        // Convert the Map count to a SparseArray
        SparseArray<Integer> documentColumn =
            new SparseIntHashArray(totalNumberOfUniqueExpressions);
        for (Map.Entry<String,Integer> e : expressionCounts)
            documentColumn.set(
                    expressionToIndex.getDimension(e.getKey()), e.getValue());

        // Update the term-document matrix with the results of processing the
        // document.
        expressionDocumentMatrixBuilder.addColumn(documentColumn);
    }

    /**
     * {@inheritDoc}
     */
    public Set<String> getWords() {
        return Collections.unmodifiableSet(expressionToIndex.keySet());
    }

    /**
     * {@inheritDoc}
     */
    public Vector getVector(String word) {
        // determine the index for the word
        int index = expressionToIndex.getDimension(word);

        return (index < 0) ? null : expressionSpace.getRowVector(index);
    }

    /**
     * {@inheritDoc}
     */
    public int getVectorLength() {
        return expressionSpace.columns();
    }

    /**
     * Processes the {@link cz.zcu.luk.sspace.common.GenericTermDocumentVectorSpaceLoadStatsExpsStops} with the provided
     * {@link edu.ucla.sspace.matrix.Transform} if it is not {@code null} as a {@link edu.ucla.sspace.matrix.MatrixFile}.
     * Otherwise, the raw term document counts are returned.  Sub classes must
     * call this in order to access the term document counts before doing any
     * other processing.
     *
     * @param transform A matrix transform used to rescale the original raw
     *        document counts.  If {@code null} no transform is done.
     */
    protected MatrixFile processSpace(Transform transform) {
        // should not run.. overwritten in child..
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
     * Subclasses should override this method if they need to utilize a header
     * token for each document.  Implementations of this method <b>must</b> be
     * thread safe.  The default action is a no-op.
     *
     * @param docIndex The document id assigned to the current document
     * @param documentName The name of the current document.
     */
    protected void handleDocumentHeader(int docIndex, String header) {
    }

    /**
     * LK added
     * @param transform
     * @return
     */
    protected Matrix processExpressionSpace(TransformExtended transform) {
        // first ensure that we are no longer writing to the matrix
        expressionDocumentMatrixBuilder.finish();

        // Get the finished matrix file from the builder
        File expressionDocumentMatrix = expressionDocumentMatrixBuilder.getFile();

        Matrix matrixExpressions = null;
        Matrix matrixExpressionsTransformed = null;

        // If a transform was specified, perform the matrix transform.
        if (transform != null) {
            LoggerUtil.info(LOG, "performing %s transform - expressions", transform);

            LoggerUtil.verbose(
                    LOG,"stored term-document matrix in format %s at %s",
                    expressionDocumentMatrixBuilder.getMatrixFormat(),
                    expressionDocumentMatrix.getAbsolutePath());

            // Convert the raw term counts using the specified transform
            MatrixFile mFile = new MatrixFile(
                    expressionDocumentMatrix,
                    expressionDocumentMatrixBuilder.getMatrixFormat());
            matrixExpressions = mFile.load();
            LoggerUtil.info(LOG, "Matrix loaded in memory!");

            // new LK hack.. since the matrix size does not have to be stored into file
            // (consider adding last empty column to Matlab sparse matrix builder and its finish method..)
            // add 0 (does not change anything since the value is 0) to ensure the right matrix size..
            if (matrixExpressions.columns() != processedDocNumbers.size()) {
                int oldColsCount = matrixExpressions.columns();
                matrixExpressions.set(0, processedDocNumbers.size()-1, 0);
                System.out.println("Fixing column size!!! From: " + oldColsCount + " to: " + matrixExpressions.columns() + "=" + (processedDocNumbers.size()));
            }
            else {
                System.out.println("Fixing column size not needed!!! " + matrixExpressions.columns() + " = " + (processedDocNumbers.size()));
            }

            matrixExpressionsTransformed = new SparseHashMatrix(matrixExpressions.rows(), matrixExpressions.columns());
            for (int i = 0; i < matrixExpressions.rows(); i++) {
                SparseDoubleVector transformedRowVec = transform.transformRow(matrixExpressions.getRowVector(i));
                for (int nz : transformedRowVec.getNonZeroIndices()) {
                    matrixExpressionsTransformed.set(i, nz, transformedRowVec.get(nz));
                }
                if (i % 10000 == 0) LoggerUtil.info(LOG, "Transformed row number: %s", i);
            }

            LoggerUtil.verbose(
                    LOG, "transformed matrix to %s",
                    expressionDocumentMatrix.getAbsolutePath());
        }
        return matrixExpressionsTransformed;
    }

//    private String findWordWithDim(BasisMapping<String, String> expressionToIndex, int i) {
//        for (String word : expressionToIndex.keySet()) {
//            if (expressionToIndex.getDimension(word) == i) {
//                 return word;
//            }
//        }
//        throw new IllegalStateException();
//    }
}
