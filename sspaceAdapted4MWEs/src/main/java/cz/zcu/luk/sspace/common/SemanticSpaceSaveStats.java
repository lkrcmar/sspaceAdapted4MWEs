package cz.zcu.luk.sspace.common;

import edu.ucla.sspace.common.SemanticSpace;

/**
 * Created with IntelliJ IDEA.
 * User: Lukr
 * Date: 30.7.13
 * Time: 13:15
 * To change this template use File | Settings | File Templates.
 */
public interface SemanticSpaceSaveStats extends SemanticSpace {

    public void saveStatistics(String dirPlusSpaceNameNoExtension);
}
