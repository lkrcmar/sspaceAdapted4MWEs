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
import edu.ucla.sspace.matrix.MatrixIO;
import edu.ucla.sspace.vector.DoubleVector;

import java.io.File;


/**
 * Performs no transform on the input matrix.
 */
public class NoTransformExtended extends BaseTransformExtended
        implements /*TransformExtended.. extends..*/ java.io.Serializable {

    private static final long serialVersionUID = 1L;

    /**
     * {@inheritDoc}
     */
    protected GlobalTransformExtended getTransform(File inputMatrixFile,
                                           MatrixIO.Format format) {
        return new NoOpTransformExtended();
    }
    
    /**
     * {@inheritDoc}
     */
    protected GlobalTransformExtended getTransform(Matrix matrix) {
        return new NoOpTransformExtended();
    }

    /**
     * {@inheritDoc}
     */
    public String toString() {
        return "NO";
    }    

    static class NoOpTransformExtended implements GlobalTransformExtended {

        public double transform(int row, int column, double value) {
            return value;
        }

        public double transform(int row, DoubleVector column) {
            return column.get(row);
        }

        public double transformColRow(int column, DoubleVector row) {
            return row.get(column);
        }
    }
}
