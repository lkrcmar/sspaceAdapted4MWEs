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
import edu.ucla.sspace.matrix.Matrices;
import edu.ucla.sspace.matrix.Matrix;
import edu.ucla.sspace.matrix.MatrixBuilder;
import edu.ucla.sspace.matrix.MatrixFile;
import edu.ucla.sspace.matrix.SparseHashMatrix;
import edu.ucla.sspace.matrix.Transform;
import edu.ucla.sspace.text.IteratorFactory;
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
public abstract class GenericTermDocumentVectorSpaceLoadStats
        implements SemanticSpaceLoadStats, java.io.Serializable {

    private static final long serialVersionUID = 1L;

    protected static final Logger LOG =
        Logger.getLogger(GenericTermDocumentVectorSpaceLoadStats.class.getName());

    protected ArrayList<Integer> processedDocNumbers;

    /**
     * A mapping from a word to the row index in the that word-document matrix
     * that contains occurrence counts for that word.
     */
    protected BasisMapping<String, String> termToIndex; // final keyword removed to allow load statistics..

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
    protected BasisMapping<String, String> otherToIndex;
    /**
     * LK added..
     * The builder used to construct the compound-document matrix as new documents
     * are processed.
     */
    private transient MatrixBuilder otherDocumentMatrixBuilder;
    /**
     * LK added..
     */
    protected Matrix otherSpace;

    /**
     * Constructs the {@code GenericTermDocumentVectorSpace}.
     *
     * @throws java.io.IOException if this instance encounters any errors when creatng
     *         the backing array files required for processing
     */
    public GenericTermDocumentVectorSpaceLoadStats() throws IOException {
        this(false,new StringBasisMapping(),Matrices.getMatrixBuilderForSVD());
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
    public GenericTermDocumentVectorSpaceLoadStats(
            boolean readHeaderToken,
            BasisMapping<String, String> termToIndex,
            MatrixBuilder otherDocumentMatrixBuilder) throws IOException {
        this.readHeaderToken = readHeaderToken;
        this.termToIndex = termToIndex;
        documentCounter = new AtomicInteger(0);

        // LK added and modified
        this.otherToIndex = new StringBasisMapping();
        this.otherDocumentMatrixBuilder = otherDocumentMatrixBuilder;

        wordSpace = null;
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
        Counter<String> otherCounts = new ObjectCounter<String>();
        //Iterator<String> documentTokens = IteratorFactory.tokenize(document);
        // LK change.. otherwise document is considered as empty even if it contains stopwords!
        Iterator<String> documentTokens = IteratorFactory.tokenizeOrdered(document);
        //System.out.println(documentTokens.getClass().toString());

        // Increaes the count of documents observed so far.
        int docCount = documentCounter.getAndAdd(1);

        // LK added..
        if (!processedDocNumbers.contains(docCount)) {
            return;
        }
        //System.out.println(processedDocNumbers.toString());
        // If the first token is to be interpreted as a document header read it.
        if (readHeaderToken)
            handleDocumentHeader(docCount, documentTokens.next());

        // If the document is empty, skip it
        if (!documentTokens.hasNext())
            return;

        // For each word in the text document, keep a count of how many times it
        // has occurred
        while (documentTokens.hasNext()) {
            String word = documentTokens.next();

            // Skip added empty tokens for words that have been filtered out
            if (word.equals(IteratorFactory.EMPTY_TOKEN))
                continue;

            // Add the term to the total list of terms to ensure it has a proper
            // index.  If the term was already added, this method is a no-op
            otherToIndex.getDimension(word);

            otherCounts.count(word);
        }

        document.close();

//        // Check that we actually loaded in some terms before we increase the
//        // documentIndex. This is done after increasing the document count since
//        // some configurations may need the document order preserved, for
//        // example, if each document corresponds to some cluster assignment.
//        if (otherCounts.size() == 0)
//            return;
        // LK changed.. see up - use processedDocumentSet to decide which docs to skip..
        // empty columns also have to be added!

        // Get the total number of terms encountered so far, including any new
        // unique terms found in the most recent document
        int totalNumberOfUniqueOthers = otherToIndex.numDimensions();

        // Convert the Map count to a SparseArray
        SparseArray<Integer> documentColumn =
            new SparseIntHashArray(totalNumberOfUniqueOthers);
        for (Map.Entry<String,Integer> e : otherCounts)
            documentColumn.set(
                    otherToIndex.getDimension(e.getKey()), e.getValue());

        // Update the term-document matrix with the results of processing the
        // document.
        otherDocumentMatrixBuilder.addColumn(documentColumn);
    }

    /**
     * {@inheritDoc}
     */
    public Set<String> getWords() {
        return Collections.unmodifiableSet(otherToIndex.keySet());
    }

    /**
     * {@inheritDoc}
     */
    public Vector getVector(String word) {
        // determine the index for the word
        int index = otherToIndex.getDimension(word);

        return (index < 0) ? null : otherSpace.getRowVector(index);
    }

    /**
     * {@inheritDoc}
     */
    public int getVectorLength() {
        return otherSpace.columns();
    }

    /**
     * Processes the {@link cz.zcu.luk.sspace.common.GenericTermDocumentVectorSpaceLoadStats} with the provided
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
    protected Matrix processOtherSpace(TransformExtended transform) {
        // first ensure that we are no longer writing to the matrix
        otherDocumentMatrixBuilder.finish();

        // Get the finished matrix file from the builder
        File otherDocumentMatrix = otherDocumentMatrixBuilder.getFile();

        Matrix matrixOthers = null;
        Matrix matrixOthersTransformed = null;

        // If a transform was specified, perform the matrix transform.
        if (transform != null) {
            LoggerUtil.info(LOG, "performing %s transform - others", transform);

            LoggerUtil.verbose(
                    LOG,"stored term-document matrix in format %s at %s",
                    otherDocumentMatrixBuilder.getMatrixFormat(),
                    otherDocumentMatrix.getAbsolutePath());

            // Convert the raw term counts using the specified transform
            MatrixFile mFile = new MatrixFile(
                    otherDocumentMatrix,
                    otherDocumentMatrixBuilder.getMatrixFormat());
            matrixOthers = mFile.load();
            LoggerUtil.info(LOG, "Matrix loaded in memory!");

            // new LK hack.. since the matrix size does not have to be stored into file
            // (consider adding last empty column to Matlab sparse matrix builder and its finish method..)
            // add 0 (does not change anything since the value is 0) to ensure the right matrix size..
            if (matrixOthers.columns() != processedDocNumbers.size()) {
                int oldColsCount = matrixOthers.columns();
                matrixOthers.set(0, processedDocNumbers.size()-1, 0);
                System.out.println("Fixing column size!!! From: " + oldColsCount + " to: " + matrixOthers.columns() + "=" + (processedDocNumbers.size()));
            }
            else {
                System.out.println("Fixing column size not needed!!! " + matrixOthers.columns() + " = " + (processedDocNumbers.size()));
            }

            matrixOthersTransformed = new SparseHashMatrix(matrixOthers.rows(), matrixOthers.columns());
            for (int i = 0; i < matrixOthers.rows(); i++) {
                SparseDoubleVector transformedRowVec = transform.transformRow(matrixOthers.getRowVector(i));
                for (int nz : transformedRowVec.getNonZeroIndices()) {
                    matrixOthersTransformed.set(i, nz, transformedRowVec.get(nz));
                }
                if (i % 10000 == 0) LoggerUtil.info(LOG, "Transformed row number: %s", i);
            }

            LoggerUtil.verbose(
                    LOG, "transformed matrix to %s",
                    otherDocumentMatrix.getAbsolutePath());
        }
        return matrixOthersTransformed;
    }

//    private String findWordWithDim(BasisMapping<String, String> otherToIndex, int i) {
//        for (String word : otherToIndex.keySet()) {
//            if (otherToIndex.getDimension(word) == i) {
//                 return word;
//            }
//        }
//        throw new IllegalStateException();
//    }
}
