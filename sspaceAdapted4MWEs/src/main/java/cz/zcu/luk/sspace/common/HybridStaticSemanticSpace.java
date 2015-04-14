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

package cz.zcu.luk.sspace.common;

import edu.ucla.sspace.common.SemanticSpace;
import edu.ucla.sspace.common.SemanticSpaceIO;
import edu.ucla.sspace.vector.Vector;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Properties;
import java.util.Set;

public class HybridStaticSemanticSpace implements SemanticSpace {
    private SemanticSpace mainWordsSpace;
    private List<SemanticSpace> otherWordsSpace;
    private List<SemanticSpace> expressionsWordsSpace;

    /**
     * The name of this semantic space.
     */
    private String spaceName;

    public HybridStaticSemanticSpace(String mainWordsSpaceFNFull) throws IOException {
        String dirName = mainWordsSpaceFNFull.substring(0, mainWordsSpaceFNFull.lastIndexOf(File.separator)+1);
        String wordSpaceFileName = mainWordsSpaceFNFull.substring(mainWordsSpaceFNFull.lastIndexOf(File.separator)+1);

        if (!wordSpaceFileName.startsWith("H_") ) {
            throw new IllegalStateException("Not a hybrid semantic space file name!");
        }

        this.spaceName = "H" + wordSpaceFileName.substring(1);
        this.mainWordsSpace = SemanticSpaceIO.load(dirName + "W" + wordSpaceFileName.substring(1));
        this.otherWordsSpace = new ArrayList<SemanticSpace>();
        this.expressionsWordsSpace = new ArrayList<SemanticSpace>();
        String otherFileNamePrefix = dirName + "O" + wordSpaceFileName.substring(1).replace(".sspace", "");
        int index = 1;
        File lfile;
        while ((lfile = new File(otherFileNamePrefix +index + ".sspace")).exists()) {
            System.out.println(lfile.toString() + " loaded!");
            otherWordsSpace.add(SemanticSpaceIO.load(otherFileNamePrefix +index + ".sspace"));
            index++;
        }

        String expressionFileNamePrefix = dirName + "E" + wordSpaceFileName.substring(1).replace(".sspace", "");
        index = 1;
        while ((lfile = new File(expressionFileNamePrefix +index + ".sspace")).exists()) {
            System.out.println(lfile.toString() + " loaded!");
            expressionsWordsSpace.add(SemanticSpaceIO.load(expressionFileNamePrefix +index + ".sspace"));
            index++;
        }

        //this.expressionsWordsSpace = SemanticSpaceIO.load(dirName + "E" + wordSpaceFileName.substring(1));

        checkHybridSpace();
    }

    private void checkHybridSpace() {
        // check same vector lengths
        if (mainWordsSpace.getVectorLength() != otherWordsSpace.get(0).getVectorLength() ||
            mainWordsSpace.getVectorLength() != expressionsWordsSpace.get(0).getVectorLength()) {
            throw new IllegalStateException("Different sizes of parallel spaces!");
        }
        for (int i = 0; i < otherWordsSpace.size(); i++) {
            if (otherWordsSpace.get(0).getVectorLength() != otherWordsSpace.get(i).getVectorLength()) {
                throw new IllegalStateException("Different sizes of parallel spaces!");
            }
        }
        for (int i = 0; i < expressionsWordsSpace.size(); i++) {
            if (expressionsWordsSpace.get(0).getVectorLength() != expressionsWordsSpace.get(i).getVectorLength()) {
                throw new IllegalStateException("Different sizes of parallel spaces!");
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    public Set<String> getWords() {
        return mainWordsSpace.getWords();
    }
  
    /**
     * {@inheritDoc}
     */
    public Vector getVector(String term) {
        Vector v = mainWordsSpace.getVector(term);
        if (v == null) {
            for (SemanticSpace otherSemSpa : otherWordsSpace) {
                v = otherSemSpa.getVector(term);
                if (v != null) {
                    break;
                }
            }
        }
        if (v == null) {
            for (SemanticSpace expSemSpac : expressionsWordsSpace) {
                v = expSemSpac.getVector(term);
                if (v != null) {
                    break;
                }
            }
        }

        if (v == null) {
            System.out.println("!!!!!!!!!!!!! Null vector for: " + term);
            throw new IllegalArgumentException("Vector for term: \"" + term + "\" not found!!");
        }
        //else System.out.println("OK vector for: " + term);
        return v;
    }

    /**
     * {@inheritDoc}
     */
    public String getSpaceName() {
        return spaceName;
    }

    /**
     * {@inheritDoc}
     */
    public int getVectorLength() {
        return mainWordsSpace.getVectorLength();
    }

    /**
     * Not supported; throws an {@link UnsupportedOperationException} if called.
     *
     * @throws an {@link UnsupportedOperationException} if called
     */
    public void processDocument(BufferedReader document) { 
        throw new UnsupportedOperationException(
            "StaticSemanticSpace instances cannot be updated");
    }

    /**
     * Not supported; throws an {@link UnsupportedOperationException} if called.
     *
     * @throws an {@link UnsupportedOperationException} if called
     */
    public void processSpace(Properties props) { 
        throw new UnsupportedOperationException(
            "StaticSemanticSpace instances cannot be updated");
    }
}
