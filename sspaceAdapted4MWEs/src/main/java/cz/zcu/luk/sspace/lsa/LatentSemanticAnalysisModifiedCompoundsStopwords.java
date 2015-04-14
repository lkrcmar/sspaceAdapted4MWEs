/*
 * LK Copy of LatentSemanticAnalysisModified..
 * only but parent changed!
 */

package cz.zcu.luk.sspace.lsa;

import edu.ucla.sspace.basis.BasisMapping;
import edu.ucla.sspace.basis.StringBasisMapping;
import edu.ucla.sspace.matrix.ArrayMatrix;
import edu.ucla.sspace.matrix.DiagonalMatrix;
import edu.ucla.sspace.matrix.Matrices;
import edu.ucla.sspace.matrix.Matrix;
import edu.ucla.sspace.matrix.MatrixFile;
import edu.ucla.sspace.matrix.SVD;
import edu.ucla.sspace.matrix.factorization.SingularValueDecomposition;
import edu.ucla.sspace.text.Document;
import edu.ucla.sspace.util.LoggerUtil;
import edu.ucla.sspace.vector.DenseVector;
import edu.ucla.sspace.vector.DoubleVector;
import edu.ucla.sspace.vector.SparseDoubleVector;
import edu.ucla.sspace.vector.SparseHashDoubleVector;

import java.io.IOException;
import java.lang.ref.WeakReference;
import java.util.Iterator;
import java.util.Map;
import java.util.Properties;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

import cz.zcu.luk.sspace.common.GenericTermDocumentVectorSpaceCompoundsStopwords;
import cz.zcu.luk.sspace.matrix.*;
import cz.zcu.luk.sspace.text.IteratorFactoryStopwords;


/**
 * An implementation of Latent Semantic Analysis (LSA).  This implementation is
 * based on two papers.
 * <ul>
 *
 *   <li style="font-family:Garamond, Georgia, serif"> Landauer, T. K., Foltz,
 *     P. W., and Laham, D. (1998).  Introduction to Latent Semantic
 *     Analysis. <i>Discourse Processes</i>, <b>25</b>, 259-284.  Available <a
 *     href="http://lsa.colorado.edu/papers/dp1.LSAintro.pdf">here</a> </li>
 *
 * <li style="font-family:Garamond, Georgia, serif"> Landauer, T. K., and
 *    Dumais, S. T. (1997). A solution to Plato's problem: The Latent Semantic
 *    Analysis theory of the acquisition, induction, and representation of
 *    knowledge.  <i>Psychological Review</i>, <b>104</b>, 211-240.  Available
 *    <a href="http://lsa.colorado.edu/papers/plato/plato.annote.html">here</a>
 *    </li>
 *
 * </ul> See the Wikipedia page on <a
 * href="http://en.wikipedia.org/wiki/Latent_semantic_analysis"> Latent Semantic
 * Analysis </a> for an execuative summary.
 *
 * <p>
 *
 * LSA first processes documents into a word-document matrix where each unique
 * word is a assigned a row in the matrix, and each column represents a
 * document.  The values of ths matrix correspond to the number of times the
 * row's word occurs in the column's document.  After the matrix has been built,
 * the <a
 * href="http://en.wikipedia.org/wiki/Singular_value_decomposition">Singular
 * Value Decomposition</a> (SVD) is used to reduce the dimensionality of the
 * original word-document matrix, denoted as <span style="font-family:Garamond,
 * Georgia, serif">A</span>. The SVD is a way of factoring any matrix A into
 * three matrices <span style="font-family:Garamond, Georgia, serif">U &Sigma;
 * V<sup>T</sup></span> such that <span style="font-family:Garamond, Georgia,
 * serif"> &Sigma; </span> is a diagonal matrix containing the singular values
 * of <span style="font-family:Garamond, Georgia, serif">A</span>. The singular
 * values of <span style="font-family:Garamond, Georgia, serif"> &Sigma; </span>
 * are ordered according to which causes the most variance in the values of
 * <span style="font-family:Garamond, Georgia, serif">A</span>. The original
 * matrix may be approximated by recomputing the matrix with only <span
 * style="font-family:Garamond, Georgia, serif">k</span> of these singular
 * values and setting the rest to 0. The approximated matrix <span
 * style="font-family:Garamond, Georgia, serif"> &Acirc; = U<sub>k</sub>
 * &Sigma;<sub>k</sub> V<sub>k</sub><sup>T</sup></span> is the least squares
 * best-ﬁt rank-<span style="font-family:Garamond, Georgia, serif">k</span>
 * approximation of <span style="font-family:Garamond, Georgia, serif">A</span>.
 * LSA reduces the dimensions by keeping only the ﬁrst <span
 * style="font-family:Garamond, Georgia, serif">k</span> dimensions from the row
 * vectors of <span style="font-family:Garamond, Georgia, serif">U</span>.
 * These vectors form the <i>semantic space</i> of the words.
 *
 * <p>
 *
 * This class offers configurable preprocessing and dimensionality reduction.
 * through three parameters.  These properties should be specified in the {@code
 * Properties} object passed to the {@link #processSpace(java.util.Properties)
 * processSpace} method.
 *
 * <dl style="margin-left: 1em">
 *
 * <dt> <i>Property:</i> <code><b>{@value #MATRIX_TRANSFORM_PROPERTY}
 *      </b></code> <br>
 *      <i>Default:</i> {@link edu.ucla.sspace.matrix.LogEntropyTransform}
 *
 * <dd style="padding-top: .5em">This variable sets the preprocessing algorithm
 *      to use on the term-document matrix prior to computing the SVD.  The
 *      property value should be the fully qualified named of a class that
 *      implements {@link edu.ucla.sspace.matrix.Transform}.  The class should be public, not abstract,
 *      and should provide a public no-arg constructor.<p>
 *
 * <dt> <i>Property:</i> <code><b>{@value LSA_DIMENSIONS_PROPERTY}
 *      </b></code> <br>
 *      <i>Default:</i> {@code 300}
 *
 * <dd style="padding-top: .5em">The number of dimensions to use for the
 *       semantic space.  This value is used as input to the SVD.<p>
 *
 * <dt> <i>Property:</i> <code><b>{@value LSA_SVD_ALGORITHM_PROPERTY}
 *      </b></code> <br>
 *      <i>Default:</i> {@link edu.ucla.sspace.matrix.SVD.Algorithm#ANY}
 *
 * <dd style="padding-top: .5em">This property sets the specific SVD algorithm
 *       that LSA will use to reduce the dimensionality of the word-document
 *       matrix.  In general, users should not need to set this property, as the
 *       default behavior will choose the fastest available on the system.<p>
 *
 * <dt> <i>Property:</i> <code><b>{@value RETAIN_DOCUMENT_SPACE_PROPERTY}
 *      </b></code> <br>
 *      <i>Default:</i> {@code false}
 *
 * <dd style="padding-top: .5em">This property indicate whether the document
 *       space should be retained after {@code processSpace}.  Setting this
 *       property to {@code true} will enable the {@link #getDocumentVector(int)
 *       getDocumentVector} method. <p>
 *
 * </dl> <p>
 *
 * <p>
 *
 * This class is thread-safe for concurrent calls of {@link
 * #processDocument(BufferedReader) processDocument}.  Once {@link
 * #processSpace(java.util.Properties) processSpace} has been called, no further calls to
 * {@code processDocument} should be made.  This implementation does not support
 * access to the semantic vectors until after {@code processSpace} has been
 * called.
 *
 * @see edu.ucla.sspace.matrix.Transform
 * @see edu.ucla.sspace.matrix.factorization.SingularValueDecomposition
 *
 * @author David Jurgens
 */
public class LatentSemanticAnalysisModifiedCompoundsStopwords extends GenericTermDocumentVectorSpaceCompoundsStopwords
        implements java.io.Serializable {

    private static final long serialVersionUID = 1L;

    // LK added
    private boolean compoundsInvolved = false;
    private transient WeakReference<Matrix> SigmaInvTimesVtRef;

    /**
     * The prefix for naming publically accessible properties
     */
    private static final String PROPERTY_PREFIX =
            "edu.ucla.sspace.lsa.LatentSemanticAnalysis";

    /**
     * The property to define the {@link edu.ucla.sspace.matrix.Transform} class to be used
     * when processing the space after all the documents have been seen.
     */
    public static final String MATRIX_TRANSFORM_PROPERTY =
            PROPERTY_PREFIX + ".transform";

//    // LK change - added
//    public static String getMatrixTransformPropertyShort() {
//        return MATRIX_TRANSFORM_PROPERTY.substring(PROPERTY_PREFIX.length()).replace(".", "-");
//    }

    /**
     * The property to set the number of dimension to which the space should be
     * reduced using the SVD
     */
    public static final String LSA_DIMENSIONS_PROPERTY =
            PROPERTY_PREFIX + ".dimensions";

    /**
     * The property to set the specific SVD algorithm used by an instance during
     * {@code processSpace}.  The value should be the name of a {@link
     * edu.ucla.sspace.matrix.SVD.Algorithm}.  If this property is unset, any
     * available algorithm will be used according to the ordering defined in
     * {@link edu.ucla.sspace.matrix.SVD}.
     */
    public static final String LSA_SVD_ALGORITHM_PROPERTY =
            PROPERTY_PREFIX + ".svd.algorithm";

    /**
     * The property whose boolean value indicate whether the document space
     * should be retained after {@code processSpace}.  Setting this property to
     * {@code true} will enable the {@link #getDocumentVector(int)
     * getDocumentVector} method.
     */
    public static final String RETAIN_DOCUMENT_SPACE_PROPERTY =
            PROPERTY_PREFIX + ".retainDocSpace";

    /**
     * The name prefix used with {@link #getName()}
     */
    private static final String LSA_SSPACE_NAME = "xLSAC";
            //"lsa-semantic-space";

    /**
     * The document space of the term document based word space If the word
     * space is reduced.  After reduction it is the right factor matrix of the
     * SVD of the word-document matrix.  This matrix is only available after the
     * {@link #processSpace(edu.ucla.sspace.matrix.Transform, edu.ucla.sspace.matrix.SVD.Algorithm, int, boolean)
     * processSpace} method has been called.
     */
    private Matrix documentSpace;

    /**
     * The diagonal matrix of the singular values
     */
    private Matrix sigma;


    /**
     * The precomputed result of U * Sigma^-1, which is used to project new
     * document vectors into the latent document space.  Since this matrix is
     * potentially large and its use only depends on the number of calls to
     * {@link #project(edu.ucla.sspace.text.Document)}, we cache the results with a {@code
     * WeakReference}, letting the result be garbage collected if memory
     * pressure gets too high.
     */
    private transient WeakReference<Matrix> UtimesSigmaInvRef;

    /**
     * The left factor matrix of the SVD operation, which is the word space
     * prior to be multiplied by the singular values.
     */
    private Matrix U;

    /**
     * The {@link edu.ucla.sspace.matrix.factorization.SingularValueDecomposition} algorithm that will decompose the word by
     * document feature space into two smaller feature spaces: a word by class
     * feature space and a class by feature space.
     */
    private final SingularValueDecomposition reducer;

    public LatentSemanticAnalysisModifiedCompoundsStopwords(boolean retainDocumentSpace, int dimensions, TransformExtended transform,
                                                            SingularValueDecomposition reducer, boolean readHeaderToken,
                                                            BasisMapping<String, String> termToIndex, Set<String> compounds) throws IOException {

        super(readHeaderToken, termToIndex, reducer.getBuilder(), reducer.getBuilder(), compounds);
        this.reducer = reducer;
        this.transform = transform;
        this.dimensions = dimensions;
        this.retainDocumentSpace = retainDocumentSpace;
        if (compounds != null) compoundsInvolved = true;

        // LK change - added
        this.headerToIndex = new ConcurrentHashMap<Integer, Integer>();
        System.out.println("header to index inicialized!");

    }

    // LK change - added
    public int getDimensions() {
        return dimensions;
    }
    // LK change - added
    public TransformExtended getTransform() {
        return transform;
    }

    /**
     * The {@link edu.ucla.sspace.matrix.Transform} applied to the term document matrix prior to being
     * reduced.
     */
    private final TransformExtended transform;

    /**
     * The final number of latent classes that will be used to represent the
     * word space.
     */
    private final int dimensions;

    /**
     * Set the true if the reduced document space will be made accessible.
     */
    private final boolean retainDocumentSpace;

    // LK change - added
    private Map<Integer, Integer> headerToIndex;
    public Map<Integer, Integer> getHeaderToIndex() {
        return headerToIndex;
    }

    /**
     * Creates a new {@link cz.zcu.luk.sspace.lsa.LatentSemanticAnalysisModifiedCompoundsStopwords} instance.  This intializes
     * {@Link LatentSemanticAnalysis} with the default parameters set in the
     * original paper.  This construct initializes this instance such that the
     * document space is <i>not</i> retained.
     */
    public LatentSemanticAnalysisModifiedCompoundsStopwords() throws IOException {
        this(false, 300, new LogEntropyTransformExtended(),
                SVD.getFastestAvailableFactorization(),
                false, new StringBasisMapping());
    }

    /**
     * Creates a new {@link cz.zcu.luk.sspace.lsa.LatentSemanticAnalysisModifiedCompoundsStopwords} instance with the specified
     * number of dimensions.  This intializes {@Link LatentSemanticAnalysis}
     * with the default parameters set in the original paper for all other
     * parameter values.  This construct initializes this instance such that the
     * document space is <i>not</i> retained.
     *
     * @param dimensions The number of dimensions to retain in the reduced space
     */
    public LatentSemanticAnalysisModifiedCompoundsStopwords(int numDimensions) throws IOException {
        this(false, numDimensions, new LogEntropyTransformExtended(),
                SVD.getFastestAvailableFactorization(),
                false, new StringBasisMapping());
    }

    /**
     * Creates a new {@link cz.zcu.luk.sspace.lsa.LatentSemanticAnalysisModifiedCompoundsStopwords} instance with the specified
     * number of dimensions, which optionally retains both the word and document
     * spaces.  This intializes {@Link LatentSemanticAnalysis} with the default
     * parameters set in the original paper for all other parameter values.
     *
     * @param dimensions The number of dimensions to retain in the reduced space
     * @param retainDocumentSpace If true, the document space will be made
     *        accessible
     */
    public LatentSemanticAnalysisModifiedCompoundsStopwords(int numDimensions,
                                                            boolean retainDocumentSpace)
            throws IOException {
        this(retainDocumentSpace, numDimensions, new LogEntropyTransformExtended(),
                SVD.getFastestAvailableFactorization(),
                false, new StringBasisMapping());
    }

    /**
     * Constructs a new {@code LatentSemanticAnalysisModifiedCompounds} using the provided
     * objects for processing.
     *
     * @param retainDocumentSpace If true, the document space will be made
     *        accessible
     * @param dimensions The number of dimensions to retain in the reduced space
     * @param transform The {@link edu.ucla.sspace.matrix.Transform} to apply before reduction
     * @param reducer The {@link edu.ucla.sspace.matrix.factorization.SingularValueDecomposition} algorithm to
     *        apply to reduce the transformed term document matrix
     * @param readHeaderToken If true, the first token of each document will be
     *        read and passed to {@link #handleDocumentHeader(int, String)
     *        handleDocumentHeader}, which discards the header
     * @param termToIndex The {@link java.util.concurrent.ConcurrentMap} used to map strings to
     *        indices
     *
     * @throws java.io.IOException if this instance encounters any errors when creatng
     *         the backing array files required for processing
     */
    public LatentSemanticAnalysisModifiedCompoundsStopwords(boolean retainDocumentSpace,
                                                            int dimensions,
                                                            TransformExtended transform,
                                                            SingularValueDecomposition reducer,
                                                            boolean readHeaderToken,
                                                            BasisMapping<String, String> termToIndex)
            throws IOException {
        this (retainDocumentSpace,
                dimensions,
                transform,
                reducer,
                readHeaderToken,
                termToIndex, null);
    }

    /**
     * {@inheritDoc}
     */
    public String getSpaceName() {
        // LK changed
        //return LSA_SSPACE_NAME;
        return  LSA_SSPACE_NAME + "_D"+getDimensions()+"_"+getTransform().toString();
    }

    /**
     * Returns the semantics of the document as represented by a numeric vector.
     * Note that document semantics may be represented in an entirely different
     * space, so the corresponding semantic dimensions in the word space will be
     * completely unrelated.  However, document vectors may be compared to find
     * those document with similar content.
     *
     * </p>
     *
     * Similar to {@code getVector}, this method is only to be used after {@code
     * processSpace} has been called.  By default, the document space is not
     * retained unless {@code retainDocumentSpace} is set to true.
     *
     * </p>
     *
     * Implementation note: If a specific document ordering is needed, caution
     * should be used when using this class in a multi-threaded environment.
     * Beacuse the document number is based on what order it was
     * <i>processed</i>, no guarantee is made that this will correspond with the
     * original document ordering as it exists in the corpus files.  However, in
     * a single-threaded environment, the ordering will be preserved.
     *
     * @param documentNumber the number of the document according to when it was
     *        processed
     *
     * @return the semantics of the document in the document space.
     * @throws IllegalArgumentException If the document space was not retained
     *         or the document number is out of range.
     */
    public DoubleVector getDocumentVector(int documentNumber) {
        if (documentSpace == null)
            throw new IllegalArgumentException(
                    "The document space has not been retained or generated.");

        if (documentNumber < 0 || documentNumber >= documentSpace.rows()) {
            throw new IllegalArgumentException(
                    "Document number is not within the bounds of the number of "
                            + "documents: " + documentNumber);
        }
        return documentSpace.getRowVector(documentNumber);
    }

    /**
     * Returns the number of documents processed by {@link
     * LatentSemanticAnalysisModified} if the document space has been retained.
     *
     * @throws IllegalArgumentException If the document space has not been
     *         retained.
     */
    public int documentSpaceSize() {
        if (documentSpace == null)
            throw new IllegalArgumentException(
                    "The document space has not been retained or generated.");

        return documentSpace.rows();
    }

    /**
     * {@inheritDoc}
     *
     * @param properties {@inheritDoc} See this class's {@link
     *        LatentSemanticAnalysisModified javadoc} for the full list of supported
     *        properties.
     */
    public void processSpace(Properties properties) {
        // Perform any optional transformations (e.g., tf-idf) on the
        // term-document matrix
        MatrixFile processedSpace = processSpace(transform);

        LoggerUtil.info(LOG, "reducing to %d dimensions", dimensions);

        // Compute the SVD on the term-document space
        reducer.factorize(processedSpace, dimensions);

        wordSpace = reducer.dataClasses();

        U = reducer.getLeftVectors();
        sigma = reducer.getSingularValues();

        // LK added
        if (compoundsInvolved) {
            Matrix processedCompoundSpace = processCompoundSpace(transform);
            Matrix Vt = reducer.getRightVectors();
            compoundSpace = computeCompoundSpace(processedCompoundSpace, Vt);
        }

        // Save the reduced document space if requested.
        if (retainDocumentSpace) {
            LoggerUtil.verbose(LOG, "loading in document space");

            // We transpose the document space to provide easier access to
            // the document vectors, which in the un-transposed version are
            // the columns.
            documentSpace = Matrices.transpose(reducer.classFeatures());
        }
    }

    /**
     * LK added..
     */
    private Matrix computeCompoundSpace(Matrix processedCompoundSpace, Matrix Vt) {
        System.out.println(" -------------- PROJECTING COMPOUNDS! -------------- ");
        Matrix projectedCompounds = new ArrayMatrix(processedCompoundSpace.rows(), dimensions);

        // Check that we can actually project the document
        if (wordSpace == null)
            throw new IllegalStateException(
                    "processSpace has not been called, so the latent document " +
                            "space does not yet exist");

        // Ensure that when we are projecting the new document that we do not
        // add any new terms to this space's basis.
        termToIndex.setReadOnly(true);
        super.compoundToIndex.setReadOnly(true);

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

    /**
     * Projects this document into the latent document space based on the
     * distirbution of terms contained within it
     *
     * @param doc A document whose contents are to be mapped into the latent
     *        document space of this instance.  The contents of this document
     *        are tokenized using the existing {@link edu.ucla.sspace.text.IteratorFactory} settings.
     *        Note that tokens that are not in the word space of this {@link
     *        LatentSemanticAnalysisModified} instance will be ignored, so documents
     *        that consist mostly of unseen terms will likely not be represented
     *        well.
     *
     * @return the projected version of {@code doc} in this instances latent
     *         document space, using the recognized terms in the document
     *
     * @throws IllegalStateException if {@link #processSpace(java.util.Properties)} has
     *         not yet been called (that is, no latent document space exists
     *         yet).
     */
    public DoubleVector project(Document doc) {
        // Check that we can actually project the document
        if (wordSpace == null)
            throw new IllegalStateException(
                    "processSpace has not been called, so the latent document " +
                            "space does not yet exist");

        // Tokenize the document using the existing tokenization rules
        Iterator<String> docTokens =
                IteratorFactoryStopwords.tokenize(doc.reader());

        // Ensure that when we are projecting the new document that we do not
        // add any new terms to this space's basis.
        termToIndex.setReadOnly(true);
        int numDims = termToIndex.numDimensions();

        // Iterate through the document's tokens and build the document
        // representation for those terms that have an existing basis in the
        // space
        SparseDoubleVector docVec = new SparseHashDoubleVector(numDims);
        while (docTokens.hasNext()) {
            int dim = termToIndex.getDimension(docTokens.next());
            if (dim >= 0)
                docVec.add(dim, 1d);
        }

        // Transform the vector according to this instance's transform's state,
        // which should normalize the vector as the original vectors were.
        DoubleVector transformed = transform.transform(docVec);
        // LK note!!! further not transformed version of query instead of transformed is used!! wrong!!????

        // Represent the document as a 1-column matrix       // LK note - one column?? storing to row!
        Matrix queryAsMatrix = new ArrayMatrix(1, numDims);
        for (int nz : docVec.getNonZeroIndices())
            queryAsMatrix.set(0, nz, docVec.get(nz));

        // Project the new document vector, d, by using
        //
        //   d * U_k * Sigma_k^-1
        //
        // where k is the dimensionality of the LSA space

        Matrix UtimesSigmaInv = null;

        // We cache the reuslts of the U_k * Sigma_k^-1 multiplication since
        // this will be the same for all projections.
        while (UtimesSigmaInv == null) {
            if (UtimesSigmaInvRef != null
                    && ((UtimesSigmaInv = UtimesSigmaInvRef.get()) != null))
                break;

            int rows = sigma.rows();
            double[] sigmaInv = new double[rows];
            for (int i = 0; i < rows; ++i)
                sigmaInv[i] = 1d / sigma.get(i, i);
            DiagonalMatrix sigmaInvMatrix = new DiagonalMatrix(sigmaInv);

            UtimesSigmaInv =
                    Matrices.multiply(U, sigmaInvMatrix);
            // Update the field with the new reference to the precomputed matrix
            UtimesSigmaInvRef = new WeakReference<Matrix>(UtimesSigmaInv);
        }

        // Compute the resulting projected vector as a matrix
        Matrix result = Matrices.multiply(queryAsMatrix, UtimesSigmaInv);

        // Copy out the vector itself so that we don't retain a reference to the
        // matrix as a result of its getRowVector call, which isn't guaranteed
        // to return a copy.
        int cols = result.columns();
        DoubleVector projected = new DenseVector(result.columns());
        for (int i = 0; i < cols; ++i)
            projected.set(i, result.get(0, i));
        return projected;
    }

    public AtomicInteger getDocumentCounter() {
        return documentCounter;
    }

    protected void handleDocumentHeader(int docIndex, String header) {
        if (docIndex == 0 || docIndex == 100 || docIndex == 20000 || docIndex == 299400 || docIndex == 299438 || docIndex == 299439 || docIndex == 299440) {
            System.out.println("Handled doc with docNumber!!: \"" + docIndex + "\" coded as header: \"" + header + "\" ..in int.. \"" + Integer.parseInt(header)+"\"" );
        }
        headerToIndex.put(Integer.parseInt(header), docIndex);
    }
}
