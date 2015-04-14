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

import edu.ucla.sspace.matrix.Transform;
import edu.ucla.sspace.vector.DoubleVector;
import edu.ucla.sspace.vector.SparseDoubleVector;


/**
 * An interface for {@link edu.ucla.sspace.matrix.Matrix} transformations.  These tranformations should
 * support both a {@link edu.ucla.sspace.matrix.Matrix} and a serialized {@link edu.ucla.sspace.matrix.Matrix} stored in a
 * {@code File}, in one of the supported matrix formats.  Implementations are
 * strongly encouraged to implement the {@code toString} method, as many {@link
 * edu.ucla.sspace.common.SemanticSpace} implementations will use the name of
 * the applied transform for describing their state.
 *
 * @author David Jurgens
 */
public interface TransformExtended extends Transform {

    SparseDoubleVector transformRow(DoubleVector row);
}
