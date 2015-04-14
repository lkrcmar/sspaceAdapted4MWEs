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

package cz.zcu.luk.sspace.matrix;

import edu.ucla.sspace.matrix.Matrix;
import edu.ucla.sspace.matrix.TransformStatistics;
import edu.ucla.sspace.matrix.MatrixIO.Format;
import edu.ucla.sspace.matrix.TransformStatistics.MatrixStatistics;
import edu.ucla.sspace.vector.DoubleVector;
import edu.ucla.sspace.vector.SparseVector;

import java.io.File;


/**
 * Transforms a matrix using row correlation weighting.  The input matrix is
 * assumed to be formatted as rows representing terms and columns representing
 * co-occuring terms.  Each matrix cell indicates the number of times the row's
 * word occurs the other term.  See the following paper for details and
 * analysis:
 *
 * <p style="font-family:Garamond, Georgia, serif"> Rohde, D. L. T., Gonnerman,
 * L. M., Plaut, D. C. (2005).  An Improved Model of Semantic Similarity Based
 * on Lexical Co-Occurrence. <i>Cognitive Science</i> <b>(submitted)</b>.
 * Available <a
 * href="http://www.cnbc.cmu.edu/~plaut/papers/pdf/RohdeGonnermanPlautSUB-CogSci.COALS.pdf">here</a></p>

 * @author Keith Stevens
 */
public class CorrelationTransformExtended extends BaseTransformExtended {

    // LK added
    private double rowSumCashed;

    /**
     * {@inheritDoc}
     */
    protected GlobalTransformExtended getTransform(File inputMatrixFile,
                                           Format format) {
        return new CorrelationGlobalTransformExtended(inputMatrixFile, format);
    }

    /**
     * {@inheritDoc}
     */
    protected CorrelationGlobalTransformExtended getTransform(Matrix matrix) {
        return new CorrelationGlobalTransformExtended(matrix);
    }

    /**
     * Returns the name of this transform.
     */
    public String toString() {
        return "CORR";
    }

    public class CorrelationGlobalTransformExtended implements GlobalTransformExtended {

        /**
         * The summation of the values each row
         */
        private double[] rowSums;

        /**
         * The summation of the values each column
         */
        private double[] colSums;

        /**
         * The total sum of all values in the matrix.
         */
        private double totalSum;

        /**
         * Creates an instance of {@code CorrelationTransform} from a {@link
         * edu.ucla.sspace.matrix.Matrix}.
         */
        public CorrelationGlobalTransformExtended(Matrix matrix) {
            MatrixStatistics stats =
                TransformStatistics.extractStatistics(matrix);
            rowSums = stats.rowSums;
            colSums = stats.columnSums;
            totalSum = stats.matrixSum;
        }

        /**
         * Creates an instance of {@code CorrelationTransform} from a {@code
         * File} for format {@link edu.ucla.sspace.matrix.MatrixIO.Format}.
         */
        public CorrelationGlobalTransformExtended(File inputMatrixFile,
                                          Format format) {
            MatrixStatistics stats =
                TransformStatistics.extractStatistics(inputMatrixFile, format);
            rowSums = stats.rowSums;
            colSums = stats.columnSums;
            totalSum = stats.matrixSum;
        }

        /**
         * Computes the correlation, scaled using the square root, between item
         * {@code row} and feature {@code column} where {@code value} specifies
         * the number of occurances.   If {@code value} is zero, the correlation
         * is zero.
         *
         * @param row The index specifying the item being observed
         * @param column The index specifying the feature being observed
         * @param value The number of occurance of the item and feature
         *
         * @return the square root of the correlation between the item aand
         *         feature
         */
        public double transform(int row, int column, double value) {
            if (value == 0d)
                return 0;

            double newValue = 
                (totalSum * value - rowSums[row] * colSums[column]) / 
                Math.sqrt(rowSums[row] * (totalSum - rowSums[row]) *
                        colSums[column] * (totalSum - colSums[column]));
            return (newValue > 0) ? Math.sqrt(newValue) : 0;
        }

        /**
         * Computes the correlation, scaled using the square root, between item
         * {@code row} and feature {@code column} where {@code value} specifies
         * the number of occurances.   If {@code value} is zero, the correlation
         * is zero.
         *
         * @param row The index specifying the item being observed
         * @param column The index specifying the feature being observed
         * @param value The number of occurance of the item and feature
         *
         * @return the square root of the correlation between the item aand
         *         feature
         */
        public double transform(int row, DoubleVector column) {
            double value = column.get(row);
            if (value == 0d)
                return 0;

            // Calcuate the term frequencies in this new document
            double colSum = 0;
            if (column instanceof SparseVector) {
                SparseVector sv = (SparseVector)column;
                for (int nz : sv.getNonZeroIndices())
                    colSum += column.get(nz);
            }
            else {
                int length = column.length();
                for (int i = 0; i < length; ++i)
                    colSum += column.get(i);
            }

            double newValue = 
                (totalSum * value - rowSums[row] * colSum) / 
                Math.sqrt(rowSums[row] * (totalSum - rowSums[row]) *
                        colSum * (totalSum - colSum));
            return (newValue > 0) ? Math.sqrt(newValue) : 0;
        }

        public double transformColRow(int column, DoubleVector row) {
            double rowSum;
            if (column == 0) { // hack - the method is called for the whole row from the first column..
                // Calcuate the term frequencies in this new term (compound..)
                rowSum = 0;
                if (row instanceof SparseVector) {
                    SparseVector sv = (SparseVector)row;
                    for (int nz : sv.getNonZeroIndices())
                        rowSum += row.get(nz);
                }
                else {
                    int length = row.length();
                    for (int i = 0; i < length; ++i)
                        rowSum += row.get(i);
                }
                rowSumCashed = rowSum;
            }
            else {
                rowSum = rowSumCashed;
            }

            double value = row.get(column);
            if (value == 0d) {
                return 0;
            }

            double newValue =
                    (totalSum * value - rowSum * colSums[column]) /
                            Math.sqrt(rowSum * (totalSum - rowSum) *
                                    colSums[column] * (totalSum - colSums[column]));
            return (newValue > 0) ? Math.sqrt(newValue) : 0;
        }
    }
}
