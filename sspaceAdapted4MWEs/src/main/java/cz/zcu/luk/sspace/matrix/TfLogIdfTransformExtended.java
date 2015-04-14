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

package cz.zcu.luk.sspace.matrix;

import edu.ucla.sspace.matrix.Matrix;
import edu.ucla.sspace.matrix.TransformStatistics;
import edu.ucla.sspace.matrix.MatrixIO.Format;
import edu.ucla.sspace.matrix.TransformStatistics.MatrixStatistics;
import edu.ucla.sspace.vector.DoubleVector;
import edu.ucla.sspace.vector.VectorMath;

import java.io.File;


/**
 * Tranforms a matrix according to the <a
 * href="http://en.wikipedia.org/wiki/Tf%E2%80%93idf">Term frequency-Inverse
 * Document Frequency</a> weighting.  The input matrix is assumed to be
 * formatted as rows representing terms and columns representing documents.
 * Each matrix cell indicates the number of times the row's word occurs within
 * the column's document.  For full details see:
 *
 * <ul><li style="font-family:Garamond, Georgia, serif">Spärck Jones, Karen
 *      (1972). "A statistical interpretation of term specificity and its
 *      application in retrieval". <i>Journal of Documentation</i> <b>28</b>
 *      (1): 11–21.</li></ul>
 *
 * @author David Jurgens
 *
 * @see edu.ucla.sspace.matrix.LogEntropyTransform
 */
public class TfLogIdfTransformExtended extends BaseTransformExtended
        implements java.io.Serializable {

    private static final long serialVersionUID = 1L;

    // LK added
    private double docsWithTermCountCashed;

    /**
     * {@inheritDoc}
     */
    protected GlobalTransformExtended getTransform(Matrix matrix) {
        return new TfIdfLogGlobalTransformExtended(matrix);
    }

    /**
     * {@inheritDoc}
     */
    protected GlobalTransformExtended getTransform(File inputMatrixFile,
                                           Format format) {
        return new TfIdfLogGlobalTransformExtended(inputMatrixFile, format);
    }

    /**
     * Returns the name of this transform.
     */
    public String toString() {
        return "TFLOGIDF";
    }

    public class TfIdfLogGlobalTransformExtended
            implements GlobalTransformExtended, java.io.Serializable {

        private static final long serialVersionUID = 1L;

        /**
         * The total number of documents (columns) that each row occurs in.
         */
        private double[] docTermCount;

        /**
         * The total number of documents (columns) that each term occurs in.
         */
        private double[] termDocCount;

        /**
         * The total number of documents (columns) present in the matrix.
         */
        private int totalDocCount;

        /**
         * Creates an instance of {@code TfIdfGlobalTransform} from a {@link
         * edu.ucla.sspace.matrix.Matrix}.
         */
        public TfIdfLogGlobalTransformExtended(Matrix matrix) {
            MatrixStatistics stats =
                TransformStatistics.extractStatistics(matrix, true, false);
            docTermCount = stats.columnSums;
            termDocCount = stats.rowSums;
            totalDocCount = docTermCount.length;
        }

        /**
         * Creates an instance of {@code TfIdfGlobalTransform} from a {@code
         * File} in the format {@link edu.ucla.sspace.matrix.MatrixIO.Format}.
         */
        public TfIdfLogGlobalTransformExtended(File inputMatrixFile, Format format) {
            MatrixStatistics stats = TransformStatistics.extractStatistics(
                    inputMatrixFile, format, true, false);
            docTermCount = stats.columnSums;
            termDocCount = stats.rowSums;
            totalDocCount = docTermCount.length;
        }

        /**
         * Computes the Term Frequency-Inverse Document Frequency for a given
         * value where {@code value} is the observed frequency of term {@code
         * row} in document {@code column}.
         *
         * @param row The index speicifying the term being observed
         * @param column The index specifying the document being observed
         * @param value The number of occurances of the term in the document.
         *
         * @return the TF-IDF of the observed value
         */
        public double transform(int row, int column, double value) {
            // LK change
            //double tf = value / docTermCount[column];
            double tf = Math.log(value + 1);

            double idf =
                Math.log(totalDocCount / (termDocCount[row] + 1));
            return tf * idf;
        }

        /**
         * Computes the Term Frequency-Inverse Document Frequency for a given
         * value where {@code value} is the observed frequency of term {@code
         * row} in document {@code column}.
         *
         * @param row The index speicifying the term being observed
         * @param column The index specifying the document being observed
         * @param value The number of occurances of the term in the document.
         *
         * @return the TF-IDF of the observed value
         */
        public double transform(int row, DoubleVector column) {
            // Calcuate the term frequencies in this new document
            double sum = VectorMath.sum(column);

            // LK change
            //double tf = column.get(row) / sum;
            double tf = Math.log(column.get(row) + 1);

            double idf =
                Math.log(totalDocCount / (termDocCount[row] + 1));
            return tf * idf;
        }

        public double transformColRow(int column, DoubleVector row) {
            double tf = Math.log(row.get(column) + 1);

            if (column == 0) { // hack - the method is called for the whole row from the first column..
                docsWithTermCountCashed = VectorMath.getNonZeroIndices(row);
            }

            if (tf == 0.0) return 0.0; // clear that result is 0..

            double idf = Math.log(totalDocCount / (docsWithTermCountCashed + 1));
            return tf * idf;
        }
    }
}
