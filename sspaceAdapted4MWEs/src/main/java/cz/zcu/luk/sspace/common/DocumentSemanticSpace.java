package cz.zcu.luk.sspace.common;

import edu.ucla.sspace.common.SemanticSpace;
import edu.ucla.sspace.vector.Vector;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import cz.zcu.luk.sspace.lsa.LatentSemanticAnalysisModified;

/**
 * Created with IntelliJ IDEA.
 * User: Lukr
 * Date: 2.10.12
 * Time: 13:47
 * To change this template use File | Settings | File Templates.
 */
public class DocumentSemanticSpace implements SemanticSpace {
    private LatentSemanticAnalysisModified lsaInst;

    public DocumentSemanticSpace (LatentSemanticAnalysisModified lsaInst) {
        this.lsaInst = lsaInst;
    }

    public void processDocument(BufferedReader document) throws IOException {
        // processed already..
    	return;
    }

    public Set<String> getWords() {
        Set<String> words = new HashSet<String>();
        for (Map.Entry<Integer, Integer> e : lsaInst.getHeaderToIndex().entrySet()) {
            words.add(e.getKey()+"");
        }
        return words;
    }

    public Vector getVector(String word) {
        return lsaInst.getDocumentVector(lsaInst.getHeaderToIndex().get(Integer.parseInt(word)));
    }

    public void processSpace(Properties properties) {
        // processed already..
    	return;
    }

    public String getSpaceName() {
        return "LSA-DocsSpace-D"+lsaInst.getDimensions()+"-"+lsaInst.getTransform().toString();
    }

    public int getVectorLength() {
        return lsaInst.getDimensions();
    }
}
